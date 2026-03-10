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
- A flat directory of pre-computed ColabFold embeddings produced by get_colabfold_embeds(),
  containing files named {sha256(sequence)}_single.npy and {sha256(sequence)}_pair.npy.
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

import hashlib
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
import hydra
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning import seed_everything
from torch.utils.data import Dataset, DataLoader
from torch_geometric.data.batch import Batch

from bioemu.chemgraph import ChemGraph
from bioemu.model_utils import load_sdes
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


# =============================================================================
# Foldedness parameters
# Residue pairs are 0-indexed CA atoms. Distances are in angstroms.
# State 1 (unfolded) has larger distances (above thresholds).
# State 2 (folded) has smaller distances (below thresholds).
# =============================================================================
_RESIDUE_PAIRS = torch.tensor([[11, 60], [12, 31], [64, 94]], dtype=torch.long)
_PAIR_THRESHOLDS = torch.tensor([6.5, 11.5, 13.5], dtype=torch.float32)  # angstroms
_PAIR_K = torch.tensor([4.0, 4.0, 4.0], dtype=torch.float32)


def compute_foldedness_differentiable(
    x0_pos_stack: torch.Tensor,
) -> torch.Tensor:
    """
    Compute foldedness scores for a batch of structures, fully differentiable.

    Replaces the previous mdtraj-based compute_foldedness. By operating directly
    on the PyTorch position tensor, gradients flow through the distance computation
    and sigmoid back to the score model — no surrogate trick needed.

    The foldedness score for each structure is:
        f = sigmoid( sum_pairs( (d_pair - threshold) * k ) )
    where d_pair is the CA-CA distance in angstroms. This gives f ≈ 0 for
    unfolded (large distances) and f ≈ 1 for folded (small distances).

    Args:
        x0_pos_stack: Tensor of shape (n_replications, n_residues, 3),
                      CA coordinates in nm, grad_fn intact from score model.

    Returns:
        foldedness: Tensor of shape (n_replications,) in [0, 1], grad_fn intact.
    """
    device = x0_pos_stack.device
    pair_indices = _RESIDUE_PAIRS.to(device)       # (n_pairs, 2)
    thresholds = _PAIR_THRESHOLDS.to(device)        # (n_pairs,)
    k_values = _PAIR_K.to(device)                   # (n_pairs,)

    # Extract positions for each residue in each pair
    pos_i = x0_pos_stack[:, pair_indices[:, 0], :]  # (n_rep, n_pairs, 3)
    pos_j = x0_pos_stack[:, pair_indices[:, 1], :]  # (n_rep, n_pairs, 3)

    # Compute CA-CA distances, converting nm -> angstroms
    distances_ang = torch.linalg.norm(pos_i - pos_j, dim=-1) * 10.0  # (n_rep, n_pairs)

    # Sigmoid over weighted sum of (distance - threshold): low distance -> high foldedness
    foldedness = torch.sigmoid(
        torch.sum((thresholds - distances_ang) * k_values, dim=-1)
    )  # (n_rep,), grad_fn intact

    return foldedness


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
        - "target_foldedness": float in [0, 1]
        - "single_embeds_file": Path
        - "pair_embeds_file": Path

    Args:
        foldedness_csv: Path to CSV with columns "mutant_name" and "percent_folded", or a pre-split DataFrame.
        embeds_dir: Flat cache directory produced by get_colabfold_embeds(), containing
                    files named {sha256(sequence)}_single.npy and {sha256(sequence)}_pair.npy.
        sequence_map: Dict mapping mutant_name -> amino acid sequence string.
    """

    def __init__(
        self,
        foldedness_csv: str | Path,
        embeds_dir: str | Path,
        sequence_map: dict[str, str],
    ):
        self.embeds_dir = Path(embeds_dir)
        # Accept either a file path or a pre-split DataFrame
        df = pd.read_csv(foldedness_csv) if not isinstance(foldedness_csv, pd.DataFrame) else foldedness_csv

        assert "mutant_name" in df.columns, "CSV must have a 'mutant_name' column"
        assert "percent_folded" in df.columns, "CSV must have a 'percent_folded' column"

        self.records = []
        for _, row in df.iterrows():
            name = row["mutant_name"]
            target = float(row["percent_folded"])

            if name not in sequence_map:
                logger.warning(f"No sequence found for {name}, skipping.")
                continue

            seq = sequence_map[name]
            seqsha = hashlib.sha256(seq.encode()).hexdigest()

            # get_colabfold_embeds stores all embeddings flat in cache_embeds_dir,
            # named {sha256(sequence)}_single.npy and {sha256(sequence)}_pair.npy.
            single_file = self.embeds_dir / f"{seqsha}_single.npy"
            pair_file = self.embeds_dir / f"{seqsha}_pair.npy"

            if not single_file.exists() or not pair_file.exists():
                logger.warning(
                    f"Missing embedding files for {name} (sha={seqsha[:8]}...), skipping. "
                    f"Run get_colabfold_embeds() for this sequence on a node with internet access."
                )
                continue

            self.records.append({
                "mutant_name": name,
                "sequence": seq,
                "target_foldedness": target,
                "single_embeds_file": single_file,
                "pair_embeds_file": pair_file,
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
    sequence: str,  # still needed for shape assertion
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Generate n_replications structures via a short diffusion rollout, then
    score each with the custom foldedness function.

    Returns:
        x0_pos_stack: Tensor of shape (n_replications, n_residues, 3).
                      This is differentiable w.r.t. the score model parameters
                      and is what carries the gradient signal.
        foldedness_scores: Tensor of shape (n_replications,) in [0, 1],
                           grad_fn intact — fully differentiable through the
                           distance computation and sigmoid back to the score model.
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
    )  # Shape: (n_replications * n_residues, 3), grad_fn intact

    x0_batch = x_mid.replace(pos=x0_pos)
    structures: list[ChemGraph] = [x0_batch.get_example(i) for i in range(n_replications)]

    # Stack x0_pos per structure: (n_replications, n_residues, 3), grad_fn intact
    n_residues = len(sequence)
    x0_pos_stack = torch.stack([s.pos for s in structures])  # keeps grad_fn
    assert x0_pos_stack.shape == (n_replications, n_residues, 3)

    # Compute foldedness scores differentiably directly from x0_pos_stack.
    # grad_fn is preserved through the distance computation and sigmoid,
    # so no surrogate trick is needed in estimate_squared_mean_error.
    foldedness_scores = compute_foldedness_differentiable(x0_pos_stack)
    assert foldedness_scores.shape == (n_replications,)

    return x0_pos_stack, foldedness_scores


def estimate_squared_mean_error(
    x0_pos_stack: torch.Tensor,
    foldedness_scores: torch.Tensor,
    target_foldedness: float,
    device: torch.device,
) -> torch.Tensor:
    """
    Unbiased estimate of (mean_foldedness - target)^2.

    Since foldedness_scores is now fully differentiable (grad_fn intact through
    the distance computation and sigmoid), this is a straightforward implementation
    of the unbiased cross-product estimator from loss.py — no surrogate trick needed.

    For i.i.d. samples with f_i = foldedness of sample i and target mu:
        E[(mean(f) - mu)^2] estimated unbiasedly as:
        [sum_i(f_i - mu)]^2 - sum_i(f_i - mu)^2) / (n * (n-1))

    Mirrors _estimate_squared_mean_error from loss.py.
    """
    n = foldedness_scores.shape[0]
    assert n >= 2, "Need at least 2 replications to compute unbiased estimate."

    target = torch.tensor(target_foldedness, dtype=torch.float32, device=device)
    diff = foldedness_scores - target  # grad_fn intact
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
        n_replications: int = 8,    # Structures per mutant; minimum 2, paper used 2; more = lower variance, higher cost
        mid_t: float = 0.77,        # Rollout endpoint; ~0.77 matches paper (8/35 explicit steps)
        N_rollout: int = 8,         # Number of explicit denoising steps (paper used 8 out of 35 total)
        record_grad_steps: tuple = (6, 7, 8),  # Last 3 of 8 explicit steps; matches paper (steps 3,4,5 of 5 grad steps)
        trainable_components: list[str] = None # I can now name the individual components of the model I want to train
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

        # Components to train (everything else is frozen).
        # Names are matched as prefixes against model_nn.* parameter names.
        # Default: train only the last 2 transformer layers + diffusion head.
        _TRAINABLE_PREFIXES = set()
        for prefix in trainable_components:
            _TRAINABLE_PREFIXES.add(prefix)

        frozen, trainable = 0, 0
        for name, param in self.score_model.named_parameters():
            if any(name.startswith(p) for p in _TRAINABLE_PREFIXES):
                param.requires_grad = True
                trainable += param.numel()
            else:
                param.requires_grad = False
                frozen += param.numel()
        logger.info(
            f"Trainable components: {trainable_components}\n"
            f"  {trainable:,} trainable params, {frozen:,} frozen params."
        )

    def configure_optimizers(self):
        # Only pass parameters that require gradients to the optimizer.
        # Frozen parameters (requires_grad=False) are automatically excluded.
        trainable_params = [p for p in self.score_model.parameters() if p.requires_grad]
        assert len(trainable_params) > 0, (
            "No trainable parameters! Check --trainable_components prefixes match actual parameter names."
        )
        return torch.optim.AdamW(
            trainable_params,
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
            x0_pos_stack, foldedness_scores = rollout_and_score(
                score_model=self.score_model,
                sdes=self.sdes,
                context_chemgraph=context_cg,
                n_replications=self.n_replications,
                mid_t=self.mid_t,
                N_rollout=self.N_rollout,
                record_grad_steps=self.record_grad_steps,
                device=device,
                sequence=sequence,
            )

            # Compute loss for this mutant
            loss = estimate_squared_mean_error(
                x0_pos_stack=x0_pos_stack,
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
        # Log absolute loss for checkpointing: we want the model closest to 0,
        # not the most negative. val/loss_abs is what ModelCheckpoint monitors.
        self.log(f"{split}/loss_abs", torch.abs(total_loss), prog_bar=False, batch_size=len(batch))
        return total_loss

    def training_step(self, batch: list[dict], batch_idx: int) -> torch.Tensor:
        return self._step(batch, "train")

    def validation_step(self, batch: list[dict], batch_idx: int) -> torch.Tensor:
        return self._step(batch, "val")


# =============================================================================
# Entry point
# =============================================================================
def load_score_model(ckpt_path: Path, model_config_path: Path):
    """Load score model from checkpoint and config."""
    from bioemu.models import DiGConditionalScoreModel

    assert ckpt_path.is_file(), f"Checkpoint not found: {ckpt_path}"
    model_state = torch.load(ckpt_path, map_location="cpu", weights_only=True)
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

    # Load pretrained model directly from provided paths (no download needed)
    ckpt_path = Path(args.ckpt_path)
    model_config_path = Path(args.model_config_path)
    assert ckpt_path.exists(), f"Checkpoint not found: {ckpt_path}"
    assert model_config_path.exists(), f"Model config not found: {model_config_path}"
    sdes = load_sdes(model_config_path, cache_so3_dir=args.cache_so3)
    score_model = load_score_model(ckpt_path, model_config_path)

    # Load sequence map (mutant_name -> sequence)
    sequence_map = load_sequence_map(args.sequence_map)

    # 80/20 random train/val split.
    # If results are noisy, consider stratified splitting by percent_folded
    # to ensure both sets cover the full range of foldedness values.
    full_df = pd.read_csv(args.foldedness_csv)
    train_df = full_df.sample(frac=0.8, random_state=args.seed if args.seed is not None else 42)
    val_df = full_df.drop(train_df.index)
    logger.info(f"Split: {len(train_df)} train / {len(val_df)} val mutants.")

    train_dataset = FoldednessDataset(
        foldedness_csv=train_df,
        embeds_dir=args.embeds_dir,
        sequence_map=sequence_map,
    )
    val_dataset = FoldednessDataset(
        foldedness_csv=val_df,
        embeds_dir=args.embeds_dir,
        sequence_map=sequence_map,
    )

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
        trainable_components=args.trainable_components,
    )

    logger_wb = WandbLogger(
        project="bioemu-foldedness-finetune",
        save_dir=args.output_dir,
    )

    callbacks = [
        LearningRateMonitor(logging_interval="step"),
        ModelCheckpoint(
            dirpath=os.path.join(args.output_dir, "checkpoints"),
            filename="bioemu-foldedness-{epoch:02d}-{val/loss_abs:.4f}",
            save_top_k=3,
            monitor="val/loss_abs",  # save checkpoints closest to loss=0, not most negative
            mode="min",
        ),
        ModelCheckpoint(
            dirpath=os.path.join(args.output_dir, "checkpoints"),
            filename="bioemu-foldedness-final",
            save_last=False,      # don't overwrite every epoch
            save_top_k=1,         # just keep the one final checkpoint
            monitor="epoch",
            mode="max",           # highest epoch number = final epoch
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
        help="Flat cache directory produced by get_colabfold_embeds() (contains {sha256}_single/pair.npy files).")
    parser.add_argument("foldedness_csv", type=str,
        help="CSV with columns 'mutant_name' and 'percent_folded'.")
    parser.add_argument("output_dir", type=str,
        help="Directory for checkpoints and logs.")
    parser.add_argument("sequence_map", type=str,
        help="CSV with columns 'mutant_name' and 'sequence'.")

    # Model paths (provide directly; no download needed)
    parser.add_argument("--ckpt_path", type=str,
        default="checkpoints/bioemu-v1.1/checkpoint.ckpt",
        help="Path to BioEmu checkpoint file (.ckpt).")
    parser.add_argument("--model_config_path", type=str,
        default="checkpoints/bioemu-v1.1/config.yaml",
        help="Path to BioEmu model config file (.yaml).")
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
    parser.add_argument(
        "--trainable_components", type=str, nargs="+",
        default=["model_nn.st_module.encoder.layers.6",
                 "model_nn.st_module.encoder.layers.7",
                 "model_nn.st_module.diff_head"],
        help="Prefixes of parameter groups to leave trainable. Everything else is frozen."
    )

    # Rollout hyperparameters (from loss.py PPFT approach)
    parser.add_argument("--n_replications", type=int, default=8,
        help="Structures generated per mutant per step.")
    parser.add_argument("--mid_t", type=float, default=0.77,
        help="Rollout endpoint time. ~0.77 matches paper (8 explicit steps out of 35 total).")
    parser.add_argument("--N_rollout", type=int, default=8,
        help="Number of explicit denoising steps before the skip-to-x0 jump (paper used 8).")
    parser.add_argument("--record_grad_steps", type=int, nargs="+", default=[6, 7, 8],
        help="Which of the N_rollout steps to record gradients through (1-indexed). Paper used last 3 of 8.")

    args = parser.parse_args()
    main(args)
