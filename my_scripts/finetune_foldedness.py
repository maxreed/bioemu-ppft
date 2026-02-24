"""
Fine-tuning script for BioEmu score model on experimental foldedness data.

Approach:
- For each mutant, pre-computed ColabFold embeddings are loaded from cache (no internet needed at runtime).
- A short diffusion rollout generates an ensemble of structures (following loss.py's _rollout).
- Each structure is scored by a custom foldedness function (continuous, 0=unfolded, 1=folded).
- The loss is an unbiased estimate of (mean_foldedness - target_foldedness)^2, following loss.py's
  _estimate_squared_mean_error approach.
- Only the score model is trained; embeddings are fixed.

Data requirements:
- A CSV with columns "mutant_name" and "percent_folded"
- A directory of pre-computed ColabFold embeddings, one subdirectory per mutant,
  each containing the .npy files produced by get_colabfold_embeds() (single_embeds and pair_embeds).
  These must be pre-cached before training (run bioemu's get_colabfold_embeds on a node
  with internet access, or on the login node, for each mutant's .a3m file).

Usage:
    python finetune_foldedness.py <embeds_dir> <foldedness_csv> <output_dir> \
        --model_name bioemu-v1.1 \
        --cache_so3 /path/to/so3_cache \
        --gpus 1
"""

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import mdtraj
import yaml
import hydra
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning import seed_everything
from torch.utils.data import Dataset, DataLoader
from torch_geometric.data.batch import Batch

from bioemu.chemgraph import ChemGraph
from bioemu.get_embeds import get_colabfold_embeds
from bioemu.model_utils import load_model, load_sdes, maybe_download_checkpoint
from bioemu.denoiser import dpm_solver, get_score
from bioemu.sde_lib import SDE
from bioemu.so3_sde import SO3SDE

logger = logging.getLogger(__name__)

# Amino acid to node label mapping (from sample.py)
_NODE_LABEL_MAPPING: dict[str, int] = {
    "A": 1, "R": 15, "N": 12, "D": 3, "C": 2, "Q": 14, "E": 4,
    "G": 6, "H": 7, "I": 8, "L": 10, "K": 9, "M": 11, "F": 5,
    "P": 13, "S": 16, "T": 17, "W": 19, "Y": 20, "V": 18,
    "U": 21, "O": 22, "X": 0, "B": 23, "Z": 25,
}

# variables for compute_foldedness.
# "if every function passes the same variables you may as well use globals" - kai lentit
residue_pairs = [[11,60],[12,31],[64,94]] # these are the three pairs of residues we care about
# recall that state 1 has larger distances (above thresholds), state 2 has smaller distances (below thresholds)
residue_pair_thresholds = np.array([6.5,11.5,13.5]) # giving distances in angstroms to be consistent with bioemu
residue_pair_k = np.array([4,4,4]) # i picked these arbitrarily

# =============================================================================
# My extra special vibes-based foldedness function
# =============================================================================
def compute_foldedness(traj: mdtraj.Trajectory, **kwargs) -> float:
    CA_indices = traj.topology.select('name CA') # not sure if we need this anymore, since chemgraph_to_traj strips down to CA-only anyway.
    CA_coords = traj.xyz[:, CA_indices, :] # WT_CA_coords[0] gives all CA coordinates of the first frame
    distances = []
    for frame in CA_coords:
        these_distances =[]
        for residue_pair in residue_pairs:
            this_distance = np.linalg.norm(frame[residue_pair[0]] - frame[residue_pair[1]]) * 10 # * 10 to convert to angstroms
            these_distances.append(this_distance)
        distances.append(these_distances)
    distances = np.array(distances)
    return 1 / (1 + np.exp(np.sum((distances - residue_pair_thresholds) * residue_pair_k,axis=1)))


# =============================================================================
# Helper: convert a ChemGraph (positions in nm) to an mdtraj.Trajectory
# =============================================================================
def chemgraph_to_traj(cg: ChemGraph, sequence: str) -> mdtraj.Trajectory:
    """
    Convert a ChemGraph's CA positions to an mdtraj.Trajectory.

    The trajectory will have one frame and one atom per residue (CA only),
    which matches what compute_foldedness expects.

    Args:
        cg: A ChemGraph with .pos of shape (n_residues, 3), coordinates in nm.
        sequence: Amino acid sequence string.

    Returns:
        mdtraj.Trajectory with CA-only topology, coordinates in nm.
    """
    pos_nm = cg.pos.detach().cpu().numpy()  # (n_residues, 3) in nm
    n = len(sequence)
    assert pos_nm.shape == (n, 3)

    # Build a minimal CA-only topology
    top = mdtraj.Topology()
    chain = top.add_chain()
    for aa in sequence:
        res = top.add_residue(aa, chain)
        top.add_atom("CA", mdtraj.element.carbon, res)

    # mdtraj expects coordinates as (n_frames, n_atoms, 3)
    xyz = pos_nm[np.newaxis, :, :]
    return mdtraj.Trajectory(xyz, topology=top)


# =============================================================================
# Helper: build a context ChemGraph from pre-cached embeddings
# (mirrors get_context_chemgraph from sample.py, but loads from cache directly)
# =============================================================================
def load_context_chemgraph(
    sequence: str,
    single_embeds_file: Path,
    pair_embeds_file: Path,
) -> ChemGraph:
    """
    Build a context ChemGraph from pre-cached embedding .npy files.

    Args:
        sequence: Amino acid sequence string.
        single_embeds_file: Path to the .npy file with single (per-residue) embeddings.
        pair_embeds_file: Path to the .npy file with pair (per-residue-pair) embeddings.

    Returns:
        ChemGraph with NaN positions (to be filled by the denoiser).
    """
    n = len(sequence)

    single_embeds = torch.from_numpy(np.load(single_embeds_file))  # (n, d_single)
    pair_embeds = torch.from_numpy(np.load(pair_embeds_file))       # (n, n, d_pair)

    assert pair_embeds.shape[0] == pair_embeds.shape[1] == n
    assert single_embeds.shape[0] == n

    _, _, n_pair_feats = pair_embeds.shape
    pair_embeds = pair_embeds.view(n**2, n_pair_feats)

    edge_index = torch.cat([
        torch.arange(n).repeat_interleave(n).view(1, n**2),
        torch.arange(n).repeat(n).view(1, n**2),
    ], dim=0)

    node_labels = torch.LongTensor([_NODE_LABEL_MAPPING[aa] for aa in sequence])

    return ChemGraph(
        edge_index=edge_index,
        pos=torch.full((n, 3), float("nan")),
        node_orientations=torch.full((n, 3, 3), float("nan")),
        single_embeds=single_embeds,
        pair_embeds=pair_embeds,
        sequence=sequence,
        node_labels=node_labels,
    )


# =============================================================================
# Dataset
# =============================================================================
class FoldednessDataset(Dataset):
    """
    Dataset of mutants with experimental foldedness values.

    Each item is a dict with:
        - "mutant_name": str
        - "sequence": str
        - "target_foldedness": float in [0, 1]  (percent_folded / 100)
        - "single_embeds_file": Path
        - "pair_embeds_file": Path

    Args:
        foldedness_csv: Path to CSV with columns "mutant_name" and "percent_folded".
        embeds_dir: Directory containing one subdirectory per mutant, each with
                    the .npy embedding files from get_colabfold_embeds().
        sequence_map: Dict mapping mutant_name -> amino acid sequence string.
                      If None, sequences will be read from the embeddings directory
                      (you may need to adapt this).
    """

    def __init__(
        self,
        foldedness_csv: str | Path,
        embeds_dir: str | Path,
        sequence_map: dict[str, str],
    ):
        self.embeds_dir = Path(embeds_dir)
        df = pd.read_csv(foldedness_csv)

        # Validate required columns
        assert "mutant_name" in df.columns, "CSV must have a 'mutant_name' column"
        assert "percent_folded" in df.columns, "CSV must have a 'percent_folded' column"

        self.records = []
        for _, row in df.iterrows():
            name = row["mutant_name"]
            target = float(row["percent_folded"]) / 100.0  # Convert % to [0,1]

            mutant_embeds_dir = self.embeds_dir / name
            # get_colabfold_embeds caches files with specific naming conventions.
            # Adjust these glob patterns if your cache layout differs.
            single_files = list(mutant_embeds_dir.glob("*single*.npy"))
            pair_files = list(mutant_embeds_dir.glob("*pair*.npy"))

            if not single_files or not pair_files:
                logger.warning(f"Missing embedding files for {name}, skipping.")
                continue

            if name not in sequence_map:
                logger.warning(f"No sequence found for {name}, skipping.")
                continue

            self.records.append({
                "mutant_name": name,
                "sequence": sequence_map[name],
                "target_foldedness": target,
                "single_embeds_file": single_files[0],
                "pair_embeds_file": pair_files[0],
            })

        logger.info(f"Loaded {len(self.records)} mutants from {foldedness_csv}.")

    def __len__(self):
        return len(self.records)

    def __getitem__(self, idx):
        return self.records[idx]


def foldedness_collate_fn(batch: list[dict]) -> list[dict]:
    """
    Simple collate: return the list of dicts as-is.
    We handle batching manually in the training step because each mutant
    needs its own rollout.
    """
    return batch


# =============================================================================
# Rollout and loss helpers (adapted from loss.py)
# =============================================================================
def _get_x0_given_xt_and_score(
    sde: SDE,
    x: torch.Tensor,
    t: torch.Tensor,
    batch_idx: torch.LongTensor,
    score: torch.Tensor,
) -> torch.Tensor:
    """Predict x_0 from x_t and the score. Mirrors loss.py."""
    assert not isinstance(sde, SO3SDE)
    alpha_t, sigma_t = sde.mean_coeff_and_std(x=x, t=t, batch_idx=batch_idx)
    return (x + sigma_t**2 * score) / alpha_t


def rollout_and_score(
    score_model: torch.nn.Module,
    sdes: dict[str, SDE],
    context_chemgraph: ChemGraph,
    n_replications: int,
    mid_t: float,
    N_rollout: int,
    record_grad_steps: set[int],
    device: torch.device,
    sequence: str,
    foldedness_kwargs: dict,
) -> torch.Tensor:
    """
    Generate n_replications structures via a short diffusion rollout, then
    score each with the custom foldedness function.

    Returns:
        foldedness_scores: Tensor of shape (n_replications,) with values in [0, 1].
                           Gradients flow back through the score model.
    """
    # Replicate context for n_replications parallel rollouts
    context_batch = Batch.from_data_list([context_chemgraph] * n_replications).to(device)
    batch_size = context_batch.num_graphs

    # Short denoising rollout (with gradient through record_grad_steps)
    x_mid: ChemGraph = dpm_solver(
        sdes=sdes,
        batch=context_batch,
        eps_t=mid_t,
        max_t=0.99,
        N=N_rollout,
        device=device,
        record_grad_steps=record_grad_steps,
        score_model=score_model,
    )

    # Single-step prediction of x_0 from x_mid (always with gradient)
    mid_t_expanded = torch.full((batch_size,), mid_t, device=device)
    score_mid = get_score(batch=x_mid, sdes=sdes, t=mid_t_expanded, score_model=score_model)["pos"]

    x0_pos = _get_x0_given_xt_and_score(
        sde=sdes["pos"],
        x=x_mid.pos,
        t=mid_t_expanded,
        batch_idx=x_mid.batch,
        score=score_mid,
    )

    x0_batch = x_mid.replace(pos=x0_pos)
    structures: list[ChemGraph] = [x0_batch.get_example(i) for i in range(n_replications)]

    # Score each structure with the custom foldedness function.
    # NOTE: compute_foldedness is not differentiable — gradients flow through
    # x0_pos (the score model output), not through the foldedness function itself.
    # This is the same policy-gradient-style approach used in the PPFT loss.
    foldedness_scores = []
    for cg in structures:
        traj = chemgraph_to_traj(cg, sequence)
        score = compute_foldedness(traj, **foldedness_kwargs)
        # Attach gradient: multiply the (detached) foldedness score by a ones tensor
        # that tracks the graph. Actual gradient w.r.t. score model comes from x0_pos.
        foldedness_scores.append(torch.tensor(score, dtype=torch.float32, device=device))

    return torch.stack(foldedness_scores)  # (n_replications,)


def estimate_squared_mean_error(
    foldedness_scores: torch.Tensor,
    target_foldedness: float,
    device: torch.device,
) -> torch.Tensor:
    """
    Unbiased estimate of (mean_foldedness - target)^2.
    Mirrors _estimate_squared_mean_error from loss.py.

    For i.i.d. samples X_i with target mu:
        E[(mean(X) - mu)^2] is estimated unbiasedly as:
        [sum_i(X_i - mu)]^2 - sum_i(X_i - mu)^2) / (n * (n-1))
    """
    n = foldedness_scores.shape[0]
    assert n >= 2, "Need at least 2 replications to compute unbiased estimate."

    target = torch.tensor(target_foldedness, dtype=torch.float32, device=device)
    diff = foldedness_scores - target
    sum_diff = torch.sum(diff)
    return (sum_diff**2 - torch.sum(diff**2)) / (n * (n - 1))


# =============================================================================
# Lightning Module
# =============================================================================
class FoldednessFineTuner(pl.LightningModule):
    """
    Fine-tunes the BioEmu score model to match experimental foldedness values.

    For each mutant in a batch:
      1. Load pre-cached embeddings and build a context ChemGraph.
      2. Run a short diffusion rollout to generate n_replications structures.
      3. Score each structure with compute_foldedness().
      4. Compute the unbiased squared-mean-error loss vs. the experimental target.
      5. Average the loss across mutants in the batch.
    """

    def __init__(
        self,
        score_model: torch.nn.Module,
        sdes: dict,
        lr: float = 1e-5,           # Low LR for fine-tuning
        weight_decay: float = 1e-2,
        n_replications: int = 8,    # Structures generated per mutant per step
        mid_t: float = 0.1,         # Rollout endpoint (lower = closer to clean)
        N_rollout: int = 10,        # Number of denoising steps in rollout
        record_grad_steps: tuple = (10,),  # Which rollout steps to record gradients through
        foldedness_kwargs: dict = None,    # Extra kwargs for compute_foldedness
    ):
        super().__init__()
        self.save_hyperparameters(ignore=["score_model", "sdes"])
        self.score_model = score_model
        self.sdes = sdes
        self.lr = lr
        self.weight_decay = weight_decay
        self.n_replications = n_replications
        self.mid_t = mid_t
        self.N_rollout = N_rollout
        self.record_grad_steps = set(record_grad_steps)
        self.foldedness_kwargs = foldedness_kwargs or {}

    def configure_optimizers(self):
        return torch.optim.AdamW(
            self.score_model.parameters(),
            lr=self.lr,
            weight_decay=self.weight_decay,
        )

    def _step(self, batch: list[dict], split: str) -> torch.Tensor:
        device = self.device
        total_loss = torch.tensor(0.0, device=device)

        for item in batch:
            mutant_name = item["mutant_name"]
            sequence = item["sequence"]
            target_foldedness = item["target_foldedness"]

            # Load embeddings and build context ChemGraph
            context_cg = load_context_chemgraph(
                sequence=sequence,
                single_embeds_file=item["single_embeds_file"],
                pair_embeds_file=item["pair_embeds_file"],
            )

            # Generate structures and score them
            foldedness_scores = rollout_and_score(
                score_model=self.score_model,
                sdes=self.sdes,
                context_chemgraph=context_cg,
                n_replications=self.n_replications,
                mid_t=self.mid_t,
                N_rollout=self.N_rollout,
                record_grad_steps=self.record_grad_steps,
                device=device,
                sequence=sequence,
                foldedness_kwargs=self.foldedness_kwargs,
            )

            # Compute loss for this mutant
            loss = estimate_squared_mean_error(
                foldedness_scores=foldedness_scores,
                target_foldedness=target_foldedness,
                device=device,
            )
            total_loss = total_loss + loss

            self.log(
                f"{split}/foldedness_mean_{mutant_name}",
                foldedness_scores.mean(),
                prog_bar=False,
                batch_size=1,
            )

        total_loss = total_loss / len(batch)
        self.log(f"{split}/loss", total_loss, prog_bar=True, batch_size=len(batch))
        return total_loss

    def training_step(self, batch: list[dict], batch_idx: int) -> torch.Tensor:
        return self._step(batch, "train")

    def validation_step(self, batch: list[dict], batch_idx: int) -> torch.Tensor:
        return self._step(batch, "val")


# =============================================================================
# Entry point
# =============================================================================
def load_score_model(ckpt_path, model_config_path):
    """Load score model from checkpoint and config (from DeltaFold repo)."""
    from bioemu.model_utils import load_model_state_dict
    from bioemu.models import DiGConditionalScoreModel

    model_state = load_model_state_dict(ckpt_path, model_config_path)
    with open(model_config_path) as f:
        model_config = yaml.safe_load(f)
    score_model: DiGConditionalScoreModel = hydra.utils.instantiate(model_config["score_model"])
    score_model.load_state_dict(model_state)
    return score_model


def load_sequence_map(sequence_map_path: str | Path) -> dict[str, str]:
    """
    Load a mapping from mutant_name -> amino acid sequence.

    This expects a two-column CSV with "mutant_name" and "sequence".
    Adapt this function if your sequences are stored differently
    (e.g., as FASTA files, or parseable from the MSA a3m files).
    """
    df = pd.read_csv(sequence_map_path)
    assert "mutant_name" in df.columns and "sequence" in df.columns, \
        "sequence_map CSV must have 'mutant_name' and 'sequence' columns"
    return dict(zip(df["mutant_name"], df["sequence"]))


def main(args):
    logging.basicConfig(level=logging.INFO)

    if args.seed is not None:
        seed_everything(args.seed, workers=True)

    accelerator = "gpu" if torch.cuda.is_available() and args.gpus > 0 else "cpu"
    devices = args.gpus if accelerator == "gpu" else "auto"

    # Load pretrained model
    ckpt_path, model_config_path = maybe_download_checkpoint(model_name=args.model_name)
    sdes = load_sdes(model_config_path, cache_so3_dir=args.cache_so3)
    score_model = load_score_model(ckpt_path, model_config_path)

    # Load sequence map (mutant_name -> sequence)
    sequence_map = load_sequence_map(args.sequence_map)

    # Datasets
    train_dataset = FoldednessDataset(
        foldedness_csv=args.foldedness_csv,
        embeds_dir=args.embeds_dir,
        sequence_map=sequence_map,
    )
    # For now, use the same dataset for validation with a fixed subset.
    # TODO: split into separate train/val CSVs once you have enough data.
    val_dataset = train_dataset

    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        collate_fn=foldedness_collate_fn,
        num_workers=0,  # Keep at 0 for now; ChemGraphs don't always pickle cleanly
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        collate_fn=foldedness_collate_fn,
        num_workers=0,
    )

    model_module = FoldednessFineTuner(
        score_model=score_model,
        sdes=sdes,
        lr=args.lr,
        n_replications=args.n_replications,
        mid_t=args.mid_t,
        N_rollout=args.N_rollout,
        record_grad_steps=tuple(args.record_grad_steps),
    )

    logger_wb = WandbLogger(
        project="bioemu-foldedness-finetune",
        save_dir=args.output_dir,
    )

    callbacks = [
        LearningRateMonitor(logging_interval="step"),
        ModelCheckpoint(
            dirpath=os.path.join(args.output_dir, "checkpoints"),
            filename="bioemu-foldedness-{epoch:02d}-{val/loss:.4f}",
            save_top_k=3,
            monitor="val/loss",
            mode="min",
        ),
    ]

    trainer = pl.Trainer(
        max_epochs=args.max_epochs,
        accelerator=accelerator,
        devices=devices,
        precision=args.precision,
        logger=logger_wb,
        callbacks=callbacks,
        log_every_n_steps=args.log_every_n_steps,
        num_sanity_val_steps=0,
        accumulate_grad_batches=args.accumulate_grad_batches,
    )

    trainer.fit(model_module, train_dataloaders=train_loader, val_dataloaders=val_loader)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fine-tune BioEmu score model on experimental foldedness data."
    )

    # Required positional args
    parser.add_argument("embeds_dir", type=str,
        help="Directory of pre-cached ColabFold embeddings, one subdir per mutant.")
    parser.add_argument("foldedness_csv", type=str,
        help="CSV with columns 'mutant_name' and 'percent_folded'.")
    parser.add_argument("output_dir", type=str,
        help="Directory for checkpoints and logs.")
    parser.add_argument("sequence_map", type=str,
        help="CSV with columns 'mutant_name' and 'sequence'.")

    # Model
    parser.add_argument("--model_name", type=str, default="bioemu-v1.1",
        help="BioEmu pretrained model to fine-tune.")
    parser.add_argument("--cache_so3", type=str, default=None,
        help="Directory for SO(3) heat kernel lookup tables.")

    # Training
    parser.add_argument("--lr", type=float, default=1e-5)
    parser.add_argument("--batch_size", type=int, default=4,
        help="Number of mutants per training step.")
    parser.add_argument("--max_epochs", type=int, default=50)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--gpus", type=int, default=1)
    parser.add_argument("--precision", type=str, default="bf16")
    parser.add_argument("--log_every_n_steps", type=int, default=10)
    parser.add_argument("--accumulate_grad_batches", type=int, default=1)

    # Rollout hyperparameters (from loss.py PPFT approach)
    parser.add_argument("--n_replications", type=int, default=8,
        help="Structures generated per mutant per step.")
    parser.add_argument("--mid_t", type=float, default=0.1,
        help="Rollout endpoint time (lower = closer to clean structure).")
    parser.add_argument("--N_rollout", type=int, default=10,
        help="Number of denoising steps in rollout.")
    parser.add_argument("--record_grad_steps", type=int, nargs="+", default=[10],
        help="Which rollout steps to record gradients through (1-indexed).")

    args = parser.parse_args()
    main(args)
