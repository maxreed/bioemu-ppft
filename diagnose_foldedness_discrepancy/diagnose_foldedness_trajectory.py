"""
diagnose_foldedness_trajectory.py

Diagnostic script for investigating the discrepancy between foldedness estimated
at mid_t (during fine-tuning rollout) vs. after full denoising.

For each sample:
  1. Runs the full DPM denoising trajectory step by step.
  2. At every N_probe steps, computes the skip-to-x0 prediction from the
     current noisy state and evaluates foldedness on that prediction.
  3. Saves the full CA trajectory as PDB + XTC.
  4. Plots foldedness vs. denoising step (and vs. diffusion time t) across
     all samples, with mean ± std shading.

Usage:
    python diagnose_foldedness_trajectory.py \
        --msa_file data/msa_kras_mutants/WT.a3m \
        --ckpt_path /path/to/checkpoint.ckpt \
        --model_config_path /path/to/config.yaml \
        --cache_so3 /path/to/so3_cache \
        --cache_embeds_dir /path/to/embeds_cache \
        --output_dir foldedness_trajectory_diagnostic \
        --num_samples 20 \
        --N_full 35 \
        --N_probe 1 \
        --max_t 0.99 \
        --eps_t 0.001

The key arguments for reproducing the fine-tuning conditions are:
    --mid_t 0.77   (foldedness at this t is what training monitors)
    --N_full 35    (total denoising steps, matching bioemu default)
"""

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import yaml
from torch_geometric.data.batch import Batch
from tqdm import tqdm

from bioemu.chemgraph import ChemGraph
from bioemu.convert_chemgraph import save_pdb_and_xtc
from bioemu.denoiser import EulerMaruyamaPredictor, get_score
from bioemu.model_utils import load_sdes
from bioemu.sde_lib import CosineVPSDE
from bioemu.so3_sde import SO3SDE
from bioemu.get_embeds import get_colabfold_embeds
from bioemu.seq_io import parse_sequence

logger = logging.getLogger(__name__)

# =============================================================================
# Node label mapping (matches sample.py / finetune_foldedness.py)
# =============================================================================
_NODE_LABEL_MAPPING: dict[str, int] = {
    "A": 1, "R": 15, "N": 12, "D": 3, "C": 2, "Q": 14, "E": 4, "G": 6,
    "H": 7, "I": 8, "L": 10, "K": 9, "M": 11, "F": 5, "P": 13, "S": 16,
    "T": 17, "W": 19, "Y": 20, "V": 18, "U": 21, "O": 22, "X": 0,
    "B": 23, "Z": 25,
}

# =============================================================================
# Foldedness parameters — must match finetune_foldedness.py exactly
# =============================================================================
_RESIDUE_PAIRS = torch.tensor([[11, 60], [12, 31], [64, 94]], dtype=torch.long)
_PAIR_THRESHOLDS = torch.tensor([6.5, 11.5, 13.5], dtype=torch.float32)  # angstroms
_PAIR_K = torch.tensor([4.0, 4.0, 4.0], dtype=torch.float32)


def compute_foldedness(x0_pos_stack: torch.Tensor) -> torch.Tensor:
    """
    Differentiable foldedness from CA positions.

    Args:
        x0_pos_stack: (n_samples, n_residues, 3), CA coordinates in nm.

    Returns:
        foldedness: (n_samples,) in [0, 1].
    """
    device = x0_pos_stack.device
    pair_indices = _RESIDUE_PAIRS.to(device)
    thresholds = _PAIR_THRESHOLDS.to(device)
    k_values = _PAIR_K.to(device)

    pos_i = x0_pos_stack[:, pair_indices[:, 0], :]   # (n, n_pairs, 3)
    pos_j = x0_pos_stack[:, pair_indices[:, 1], :]   # (n, n_pairs, 3)
    distances_ang = torch.linalg.norm(pos_i - pos_j, dim=-1) * 10.0  # nm -> Å

    foldedness = torch.sigmoid(
        -torch.sum((distances_ang - thresholds) * k_values, dim=-1)
    )
    return foldedness


# =============================================================================
# Checkpoint loading (no HuggingFace download — local paths only)
# =============================================================================
def load_score_model(ckpt_path: Path, model_config_path: Path):
    """
    Load score model from a local checkpoint.
    Handles both:
      - Raw bioemu state dicts (flat, keys like "model_nn.*")
      - PyTorch Lightning checkpoints (nested under "state_dict", keys like
        "score_model.model_nn.*")
    """
    import hydra
    from bioemu.models import DiGConditionalScoreModel

    assert ckpt_path.is_file(), f"Checkpoint not found: {ckpt_path}"

    raw = torch.load(ckpt_path, map_location="cpu", weights_only=True)

    # Detect Lightning checkpoint by presence of "state_dict" key
    if "state_dict" in raw and isinstance(raw["state_dict"], dict):
        logger.info("Detected PyTorch Lightning checkpoint — extracting state_dict.")
        state_dict = raw["state_dict"]
        prefix = "score_model."
        state_dict = {
            k[len(prefix):]: v
            for k, v in state_dict.items()
            if k.startswith(prefix)
        }
    else:
        state_dict = raw

    with open(model_config_path) as f:
        model_config = yaml.safe_load(f)

    score_model: DiGConditionalScoreModel = hydra.utils.instantiate(
        model_config["score_model"]
    )
    score_model.load_state_dict(state_dict)
    return score_model


# =============================================================================
# Context ChemGraph from cached embeddings
# =============================================================================
def get_context_chemgraph(
    sequence: str,
    cache_embeds_dir: str | Path | None = None,
    msa_file: str | Path | None = None,
) -> ChemGraph:
    n = len(sequence)
    single_embeds_file, pair_embeds_file = get_colabfold_embeds(
        seq=sequence,
        cache_embeds_dir=cache_embeds_dir,
        msa_file=msa_file,
    )
    single_embeds = torch.from_numpy(np.load(single_embeds_file))
    pair_embeds = torch.from_numpy(np.load(pair_embeds_file))

    assert pair_embeds.shape[0] == pair_embeds.shape[1] == n
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
# x0 prediction from current noisy state (mirrors finetune_foldedness.py)
# =============================================================================
def _get_x0_from_batch(
    batch: ChemGraph,
    sdes: dict,
    score_model: torch.nn.Module,
    t_scalar: float,
    device: torch.device,
) -> torch.Tensor:
    """
    From a batch at time t, predict x0 positions via single-step score jump.

    Returns:
        x0_pos_stack: (n_samples, n_residues, 3) in nm.
    """
    batch_size = batch.num_graphs
    t = torch.full((batch_size,), t_scalar, device=device)

    score = get_score(batch=batch, sdes=sdes, t=t, score_model=score_model)

    pos_sde: CosineVPSDE = sdes["pos"]
    alpha_t, sigma_t = pos_sde.mean_coeff_and_std(
        x=batch.pos, t=t, batch_idx=batch.batch
    )
    # x0 = (x_t - sigma_t * score) / alpha_t
    x0_pos_flat = (batch.pos - sigma_t * score["pos"]) / alpha_t
    # Reshape from (batch_size * n_residues, 3) -> (batch_size, n_residues, 3)
    n_residues = batch.pos.shape[0] // batch_size
    x0_pos_stack = x0_pos_flat.view(batch_size, n_residues, 3)
    return x0_pos_stack


# =============================================================================
# Main trajectory diagnostic function
# =============================================================================
@torch.no_grad()
def run_trajectory_diagnostic(
    score_model: torch.nn.Module,
    sdes: dict,
    sequence: str,
    num_samples: int,
    N_full: int,
    max_t: float,
    eps_t: float,
    N_probe: int,
    mid_t: float,
    cache_embeds_dir: str | Path | None,
    msa_file: str | Path | None,
    output_dir: Path,
    seed: int,
    noise: float = 1.0,
) -> dict:
    """
    Run full denoising trajectory for num_samples, recording foldedness at
    every N_probe steps via a skip-to-x0 prediction.

    Returns a dict with:
        "timesteps":   list of t values at which foldedness was measured
        "step_indices": list of step indices at which foldedness was measured
        "foldedness":  np.array of shape (num_probe_points, num_samples)
        "final_pos":   np.array of shape (num_samples, n_residues, 3)
        "final_node_orientations": np.array of shape (num_samples, n_residues, 3, 3)
    """
    torch.manual_seed(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    score_model = score_model.to(device).eval()

    pos_sde: CosineVPSDE = sdes["pos"]
    so3_sde: SO3SDE = sdes["node_orientations"]
    so3_sde.to(device)

    context_cg = get_context_chemgraph(
        sequence=sequence,
        cache_embeds_dir=cache_embeds_dir,
        msa_file=msa_file,
    )
    batch = Batch.from_data_list([context_cg] * num_samples).to(device)

    # Sample from prior
    batch = batch.replace(
        pos=pos_sde.prior_sampling(batch.pos.shape, device=device),
        node_orientations=so3_sde.prior_sampling(
            batch.node_orientations.shape, device=device
        ),
    )

    timesteps = torch.linspace(max_t, eps_t, N_full, device=device)
    dt = -torch.tensor((max_t - eps_t) / (N_full - 1)).to(device)

    noiser_pos = EulerMaruyamaPredictor(
        corruption=pos_sde, noise_weight=noise, marginal_concentration_factor=1.0
    )
    noiser_ori = EulerMaruyamaPredictor(
        corruption=so3_sde, noise_weight=noise, marginal_concentration_factor=1.0
    )

    # Storage for probed foldedness values
    probe_steps = []
    probe_times = []
    probe_foldedness = []  # list of (num_samples,) arrays

    # Find the step index closest to mid_t for annotation
    mid_t_step = int(torch.argmin(torch.abs(timesteps - mid_t)).item())

    logger.info(
        f"Running {N_full}-step denoising for {num_samples} samples, "
        f"probing every {N_probe} steps. mid_t={mid_t} ≈ step {mid_t_step}."
    )

    for i in tqdm(range(N_full - 1), desc="Denoising"):
        t_scalar = timesteps[i].item()
        t = torch.full((batch.num_graphs,), t_scalar, device=device)
        t_hat = t - noise * dt if (i > 0 and t_scalar > 0.0 and t_scalar < 1.0) else t

        # Apply noise perturbation (DPM noise step)
        pos_hat = noiser_pos.forward_sde_step(
            x=batch.pos, t=t, dt=(t_hat - t)[0], batch_idx=batch.batch
        )[0]
        ori_hat = noiser_ori.forward_sde_step(
            x=batch.node_orientations, t=t, dt=(t_hat - t)[0], batch_idx=batch.batch
        )[0]
        batch_hat = batch.replace(pos=pos_hat, node_orientations=ori_hat)

        # Score at t_hat
        score = get_score(batch=batch_hat, t=t_hat, score_model=score_model, sdes=sdes)

        # DPM-Solver-2 position update
        alpha_t, sigma_t = pos_sde.mean_coeff_and_std(
            x=batch.pos, t=t_hat, batch_idx=batch.batch
        )
        lambda_t = torch.log(alpha_t / sigma_t)
        alpha_t_next, sigma_t_next = pos_sde.mean_coeff_and_std(
            x=batch.pos, t=t + dt, batch_idx=batch.batch
        )
        lambda_t_next = torch.log(alpha_t_next / sigma_t_next)
        h_t = lambda_t_next - lambda_t

        # Intermediate step
        lambda_t_middle = (lambda_t + lambda_t_next) / 2
        from bioemu.denoiser import _t_from_lambda
        t_lambda_val = _t_from_lambda(sde=pos_sde, lambda_t=lambda_t_middle)
        t_lambda = torch.full((batch.num_graphs,), t_lambda_val[0][0], device=device)
        alpha_t_lambda, sigma_t_lambda = pos_sde.mean_coeff_and_std(
            x=batch.pos, t=t_lambda, batch_idx=batch.batch
        )
        u = (
            alpha_t_lambda / alpha_t * batch_hat.pos
            + sigma_t_lambda * sigma_t * (torch.exp(h_t / 2) - 1) * score["pos"]
        )
        batch_u = batch.replace(pos=u)

        # SO3 predictor step
        so3_predictor = EulerMaruyamaPredictor(
            corruption=so3_sde, noise_weight=0.0, marginal_concentration_factor=1.0
        )
        drift, _ = so3_predictor.reverse_drift_and_diffusion(
            x=batch_hat.node_orientations, score=score["node_orientations"],
            t=t_hat, batch_idx=batch.batch,
        )
        sample, _ = so3_predictor.update_given_drift_and_diffusion(
            x=batch_hat.node_orientations, drift=drift, diffusion=0.0,
            dt=t_lambda[0] - t_hat[0],
        )
        batch_u = batch_u.replace(node_orientations=sample)

        # Correction step
        score_u = get_score(batch=batch_u, t=t_lambda, sdes=sdes, score_model=score_model)
        pos_next = (
            alpha_t_next / alpha_t * batch_hat.pos
            + sigma_t_next * sigma_t_lambda * (torch.exp(h_t) - 1) * score_u["pos"]
        )
        batch_next = batch.replace(pos=pos_next)

        # SO3 correction
        dt_hat = t + dt - t_hat
        node_score = (
            score_u["node_orientations"]
            + 0.5
            * (score_u["node_orientations"] - score["node_orientations"])
            / (t_lambda[0] - t_hat[0])
            * dt_hat[0]
        )
        drift2, _ = so3_predictor.reverse_drift_and_diffusion(
            x=batch_u.node_orientations, score=node_score,
            t=t_lambda, batch_idx=batch.batch,
        )
        sample2, _ = so3_predictor.update_given_drift_and_diffusion(
            x=batch_hat.node_orientations, drift=drift2, diffusion=0.0,
            dt=dt_hat[0],
        )
        batch = batch_next.replace(node_orientations=sample2)

        # Probe foldedness every N_probe steps and at the final step
        is_last = (i == N_full - 2)
        if i % N_probe == 0 or is_last:
            t_probe = t_scalar
            x0_stack = _get_x0_from_batch(
                batch=batch, sdes=sdes, score_model=score_model,
                t_scalar=t_probe, device=device,
            )
            fold = compute_foldedness(x0_stack)  # (num_samples,)
            probe_steps.append(i)
            probe_times.append(t_probe)
            probe_foldedness.append(fold.cpu().numpy())
            logger.debug(
                f"Step {i:3d} | t={t_probe:.3f} | "
                f"foldedness mean={fold.mean():.3f} ± {fold.std():.3f}"
            )

    # Also probe the final fully-denoised x0 directly (t=eps_t)
    # This is the "true" final foldedness without skip prediction
    n_residues = batch.pos.shape[0] // num_samples
    final_pos = batch.pos.view(num_samples, n_residues, 3).cpu()
    final_fold = compute_foldedness(final_pos.to(device)).cpu().numpy()
    logger.info(
        f"Final foldedness (direct, no skip): "
        f"mean={final_fold.mean():.3f} ± {final_fold.std():.3f}"
    )

    # Collect final structures for saving
    sampled = batch.to_data_list()
    final_pos_np = torch.stack([s.pos for s in sampled]).cpu().numpy()
    final_ori_np = torch.stack([s.node_orientations for s in sampled]).cpu().numpy()

    return {
        "step_indices": probe_steps,
        "timesteps": probe_times,
        "foldedness": np.stack(probe_foldedness, axis=0),  # (n_probes, n_samples)
        "final_foldedness_direct": final_fold,
        "mid_t_step": mid_t_step,
        "mid_t": mid_t,
        "final_pos": final_pos_np,
        "final_node_orientations": final_ori_np,
    }


# =============================================================================
# Plotting
# =============================================================================
def plot_foldedness_trajectory(results: dict, output_dir: Path, target_foldedness: float | None = None):
    step_indices = results["step_indices"]
    timesteps = results["timesteps"]
    foldedness = results["foldedness"]  # (n_probes, n_samples)
    mid_t_step = results["mid_t_step"]
    mid_t = results["mid_t"]
    final_direct = results["final_foldedness_direct"]

    mean_fold = foldedness.mean(axis=1)
    std_fold = foldedness.std(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, x_vals, xlabel in [
        (axes[0], step_indices, "Denoising step"),
        (axes[1], timesteps, "Diffusion time t"),
    ]:
        ax.plot(x_vals, mean_fold, color="steelblue", linewidth=2, label="Mean foldedness (skip-to-x0)")
        ax.fill_between(
            x_vals,
            mean_fold - std_fold,
            mean_fold + std_fold,
            alpha=0.3,
            color="steelblue",
            label="±1 std",
        )

        # Individual sample traces (faint)
        for s in range(foldedness.shape[1]):
            ax.plot(x_vals, foldedness[:, s], color="steelblue", alpha=0.08, linewidth=0.5)

        # Mark mid_t
        if xlabel == "Denoising step":
            ax.axvline(mid_t_step, color="darkorange", linestyle="--", linewidth=1.5,
                       label=f"mid_t step ({mid_t_step}, t≈{mid_t:.2f})")
        else:
            ax.axvline(mid_t, color="darkorange", linestyle="--", linewidth=1.5,
                       label=f"mid_t={mid_t:.2f}")

        # Mark final direct foldedness
        final_x = step_indices[-1] if xlabel == "Denoising step" else timesteps[-1]
        ax.axhline(
            final_direct.mean(), color="firebrick", linestyle=":",
            linewidth=1.5, label=f"Final direct mean ({final_direct.mean():.3f})"
        )

        if target_foldedness is not None:
            ax.axhline(target_foldedness, color="green", linestyle="-.", linewidth=1.5,
                       label=f"Target ({target_foldedness:.2f})")

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Predicted foldedness (skip-to-x0)")
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Foldedness across denoising trajectory", fontsize=13)
    fig.tight_layout()

    plot_path = output_dir / "foldedness_trajectory.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    logger.info(f"Saved plot to {plot_path}")
    plt.close(fig)

    # Also save a CSV of the trajectory data for further analysis
    import csv
    csv_path = output_dir / "foldedness_trajectory.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["step", "t", "mean_foldedness", "std_foldedness"] + \
                 [f"sample_{i}" for i in range(foldedness.shape[1])]
        writer.writerow(header)
        for j, (step, t) in enumerate(zip(step_indices, timesteps)):
            row = [step, f"{t:.4f}", f"{mean_fold[j]:.4f}", f"{std_fold[j]:.4f}"]
            row += [f"{foldedness[j, s]:.4f}" for s in range(foldedness.shape[1])]
            writer.writerow(row)
    logger.info(f"Saved trajectory CSV to {csv_path}")


# =============================================================================
# Entry point
# =============================================================================
def main(args):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ckpt_path = Path(args.ckpt_path)
    model_config_path = Path(args.model_config_path)

    logger.info(f"Loading score model from {ckpt_path}")
    score_model = load_score_model(ckpt_path, model_config_path)
    score_model.eval()

    logger.info(f"Loading SDEs from {model_config_path}")
    sdes = load_sdes(
        model_config_path=model_config_path,
        cache_so3_dir=args.cache_so3,
    )

    msa_file = args.msa_file if str(args.msa_file).endswith(".a3m") else None
    sequence = parse_sequence(args.msa_file)
    logger.info(f"Sequence length: {len(sequence)} residues")

    results = run_trajectory_diagnostic(
        score_model=score_model,
        sdes=sdes,
        sequence=sequence,
        num_samples=args.num_samples,
        N_full=args.N_full,
        max_t=args.max_t,
        eps_t=args.eps_t,
        N_probe=args.N_probe,
        mid_t=args.mid_t,
        cache_embeds_dir=args.cache_embeds_dir,
        msa_file=msa_file,
        output_dir=output_dir,
        seed=args.seed,
        noise=args.noise,
    )

    # Save the full final structures as PDB + XTC
    logger.info("Saving final structures as PDB + XTC...")
    pos_tensor = torch.tensor(results["final_pos"])
    ori_tensor = torch.tensor(results["final_node_orientations"])
    save_pdb_and_xtc(
        pos_nm=pos_tensor,
        node_orientations=ori_tensor,
        topology_path=output_dir / "topology.pdb",
        xtc_path=output_dir / "samples.xtc",
        sequence=sequence,
        filter_samples=False,  # don't filter — we want all samples for the diagnostic
    )

    # Save raw foldedness arrays
    np.save(output_dir / "foldedness_trajectory.npy", results["foldedness"])
    np.save(output_dir / "foldedness_final_direct.npy", results["final_foldedness_direct"])
    np.save(output_dir / "probe_timesteps.npy", np.array(results["timesteps"]))

    # Plot
    plot_foldedness_trajectory(
        results=results,
        output_dir=output_dir,
        target_foldedness=args.target_foldedness,
    )

    # Print summary
    foldedness = results["foldedness"]
    mid_step = results["mid_t_step"]
    if mid_step < len(results["step_indices"]):
        # Find the probed step closest to mid_t_step
        closest_probe = min(results["step_indices"], key=lambda s: abs(s - mid_step))
        probe_idx = results["step_indices"].index(closest_probe)
        mid_fold = foldedness[probe_idx]
        logger.info(
            f"\n--- Summary ---\n"
            f"  mid_t={args.mid_t:.2f} (step ~{mid_step}): "
            f"mean foldedness = {mid_fold.mean():.3f} ± {mid_fold.std():.3f}\n"
            f"  Final (skip-to-x0): "
            f"mean = {foldedness[-1].mean():.3f} ± {foldedness[-1].std():.3f}\n"
            f"  Final (direct, no skip): "
            f"mean = {results['final_foldedness_direct'].mean():.3f} "
            f"± {results['final_foldedness_direct'].std():.3f}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Diagnose foldedness across the denoising trajectory."
    )
    parser.add_argument("--msa_file", required=True,
        help="Path to .a3m MSA file (or .fasta / plain sequence).")
    parser.add_argument("--ckpt_path", required=True,
        help="Path to model checkpoint (.ckpt). Accepts both raw bioemu and "
             "PyTorch Lightning checkpoints.")
    parser.add_argument("--model_config_path", required=True,
        help="Path to model config YAML.")
    parser.add_argument("--output_dir", default="foldedness_trajectory_diagnostic",
        help="Directory for outputs (plot, XTC, CSV).")
    parser.add_argument("--cache_so3", default=None,
        help="Path to SO3 cache directory.")
    parser.add_argument("--cache_embeds_dir", default=None,
        help="Path to colabfold embeddings cache directory.")
    parser.add_argument("--num_samples", type=int, default=20,
        help="Number of independent denoising trajectories to run.")
    parser.add_argument("--N_full", type=int, default=35,
        help="Total number of denoising steps (bioemu default is 35).")
    parser.add_argument("--max_t", type=float, default=0.99,
        help="Starting diffusion time.")
    parser.add_argument("--eps_t", type=float, default=0.001,
        help="Ending diffusion time.")
    parser.add_argument("--mid_t", type=float, default=0.77,
        help="The mid_t used during fine-tuning, annotated on the plot.")
    parser.add_argument("--N_probe", type=int, default=1,
        help="Probe foldedness every N_probe steps. 1 = every step.")
    parser.add_argument("--target_foldedness", type=float, default=None,
        help="Optional target foldedness value to draw on the plot (e.g. 0.9 for WT).")
    parser.add_argument("--noise", type=float, default=1.0,
        help="Noise weight for DPM sampler (1.0 = stochastic, 0.0 = ODE).")
    parser.add_argument("--seed", type=int, default=42,
        help="Random seed.")

    args = parser.parse_args()
    main(args)
