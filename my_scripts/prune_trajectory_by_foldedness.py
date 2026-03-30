#!/usr/bin/env python3
"""
Prune frames from an XTC trajectory until the ensemble foldedness matches
a target value from a CSV file.

Foldedness is defined as state-2 occupancy using a sigmoid over CA-CA distances
(ported from new_loss.py). Frames are removed by randomly sampling from the
"responsible" half of the foldedness distribution (> 0.5 if overshooting,
< 0.5 if undershooting) until the mean foldedness is within --tolerance of
the target.

Single-mutant usage:
    python prune_trajectory_by_foldedness.py \
        --topology outputs/v13_allVariants/KRAS_WT/topology.pdb \
        --xtc      outputs/v13_allVariants/KRAS_WT/samples.xtc \
        --data_csv data/kras_mutants_data.csv

Batch usage (runs on every KRAS_* subdirectory):
    python prune_trajectory_by_foldedness.py \
        --batch_dir outputs/v13_allVariants \
        --data_csv  data/kras_mutants_data.csv
"""

import argparse
import os

import numpy as np
import mdtraj
import pandas as pd


# ── Foldedness definition (ported from new_loss.py) ──────────────────────────

RESIDUE_PAIRS = [[11, 60], [12, 31], [64, 94]]
THRESHOLDS = np.array([6.5, 11.5, 13.5])   # Å
K_VALS = np.array([4.0, 4.0, 4.0])


def compute_per_frame_foldedness(traj: mdtraj.Trajectory) -> np.ndarray:
    """Return per-frame foldedness scores in [0, 1] (state-2 occupancy proxy)."""
    CA_indices = traj.topology.select("name CA")
    CA_coords = traj.xyz[:, CA_indices, :]          # (n_frames, n_CA, 3) in nm
    distances = np.empty((traj.n_frames, len(RESIDUE_PAIRS)))
    for f, frame in enumerate(CA_coords):
        for p, pair in enumerate(RESIDUE_PAIRS):
            distances[f, p] = np.linalg.norm(frame[pair[0]] - frame[pair[1]]) * 10  # nm→Å
    return 1.0 / (1.0 + np.exp(np.sum((distances - THRESHOLDS) * K_VALS, axis=1)))


# ── Pruning logic ─────────────────────────────────────────────────────────────

def prune_to_target_foldedness(
    state_vals: np.ndarray,
    target: float,
    tolerance: float,
    rng: np.random.Generator = None,
) -> np.ndarray:
    """
    Randomly remove frames to bring mean foldedness within tolerance of target.

    Strategy:
      - If mean > target: randomly remove frames with foldedness > 0.5 one at a time.
      - If mean < target: randomly remove frames with foldedness < 0.5 one at a time.
    This avoids systematically depleting one tail of the foldedness distribution
    (the bias that greedy extreme-value removal would introduce).

    Stops when within tolerance or no eligible candidates remain.
    Returns a boolean mask (True = keep).
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(state_vals)
    kept = np.ones(n, dtype=bool)
    current_sum = float(state_vals.sum())
    current_n = n

    def current_mean():
        return current_sum / current_n if current_n > 0 else float("nan")

    while abs(current_mean() - target) > tolerance:
        if current_mean() > target:
            candidates = np.where(kept & (state_vals > 0.5))[0]
        else:
            candidates = np.where(kept & (state_vals < 0.5))[0]

        if len(candidates) == 0:
            break  # no eligible frames left; report whatever mean we reached

        idx = rng.choice(candidates)
        kept[idx] = False
        current_sum -= state_vals[idx]
        current_n -= 1

    return kept


# ── Per-mutant processing ─────────────────────────────────────────────────────

def process_one(
    topology: str,
    xtc: str,
    mutant_name: str,
    df: pd.DataFrame,
    tolerance: float,
    output_name: str,
    rng: np.random.Generator,
) -> dict:
    """
    Process a single mutant trajectory. Returns a result dict for summary reporting.
    """
    print(f"\n{'═' * 68}")
    print(f"  Mutant: {mutant_name}")
    print(f"{'═' * 68}")

    # Look up target foldedness
    row = df[df["mutant_name"] == mutant_name]
    if row.empty:
        print(f"  SKIPPED: '{mutant_name}' not found in CSV.")
        return {"mutant": mutant_name, "status": "skipped — not in CSV"}
    target_foldedness = float(row["percent_folded"].values[0])
    print(f"  Target foldedness: {target_foldedness:.4f}")

    # Load trajectory
    print(f"  Loading trajectory: {xtc} ...")
    traj = mdtraj.load_xtc(xtc, top=topology)
    n_original = traj.n_frames
    print(f"  {n_original} frames loaded.")

    # Compute per-frame foldedness
    print("  Computing per-frame foldedness ...")
    state_vals = compute_per_frame_foldedness(traj)
    initial_foldedness = float(state_vals.mean())
    print(f"  Initial foldedness: {initial_foldedness:.4f}")

    # Prune
    kept_mask = prune_to_target_foldedness(state_vals, target_foldedness, tolerance, rng)
    n_kept = int(kept_mask.sum())
    n_removed = n_original - n_kept
    final_foldedness = float(state_vals[kept_mask].mean())
    within = abs(final_foldedness - target_foldedness) <= tolerance

    print(f"  Frames removed   : {n_removed} / {n_original}  ({100.0 * n_removed / n_original:.1f}%)")
    print(f"  Final foldedness : {final_foldedness:.4f}  "
          f"(target {target_foldedness:.4f}, diff {abs(final_foldedness - target_foldedness):.4f})")
    if not within:
        print(f"  WARNING: could not reach target within tolerance {tolerance}.")

    # Save
    kept_indices = np.where(kept_mask)[0]
    pruned_traj = traj[kept_indices]
    out_dir = os.path.dirname(os.path.abspath(xtc))
    out_path = os.path.join(out_dir, output_name)
    pruned_traj.save_xtc(out_path)
    print(f"  Saved: {out_path}")

    return {
        "mutant": mutant_name,
        "n_original": n_original,
        "n_kept": n_kept,
        "pct_removed": 100.0 * n_removed / n_original,
        "initial_foldedness": initial_foldedness,
        "final_foldedness": final_foldedness,
        "target_foldedness": target_foldedness,
        "diff": abs(final_foldedness - target_foldedness),
        "status": "OK" if within else "WARNING: outside tolerance",
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Prune XTC trajectory frames to match a target foldedness value."
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--batch_dir",
        metavar="DIR",
        help="Directory containing KRAS_<mutant> subdirectories; runs on all of them.",
    )
    mode.add_argument(
        "--xtc",
        metavar="FILE",
        help="Path to a single input XTC trajectory (single-mutant mode).",
    )

    parser.add_argument(
        "--topology",
        metavar="FILE",
        help="Path to topology PDB file (required in single-mutant mode).",
    )
    parser.add_argument(
        "--data_csv",
        default="data/kras_mutants_data.csv",
        help="CSV with mutant_name and percent_folded columns (default: data/kras_mutants_data.csv)",
    )
    parser.add_argument(
        "--mutant_name",
        default=None,
        help=(
            "Mutant name to look up in the CSV (single-mutant mode only). "
            "If omitted, inferred from the XTC parent directory by stripping 'KRAS_'."
        ),
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.02,
        help="Acceptable absolute deviation from target foldedness (default: 0.02)",
    )
    parser.add_argument(
        "--output_name",
        default="samples_pruned.xtc",
        help="Output XTC filename within each trajectory directory (default: samples_pruned.xtc)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible frame removal (default: None)",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    df = pd.read_csv(args.data_csv)

    # ── Batch mode ───────────────────────────────────────────────────────────
    if args.batch_dir:
        batch_dir = os.path.abspath(args.batch_dir)
        subdirs = sorted(
            d for d in os.listdir(batch_dir)
            if os.path.isdir(os.path.join(batch_dir, d)) and d.startswith("KRAS_")
        )
        if not subdirs:
            raise RuntimeError(f"No KRAS_* subdirectories found in {batch_dir}")
        print(f"Batch mode: found {len(subdirs)} KRAS_* subdirectories in {batch_dir}")

        results = []
        for subdir in subdirs:
            mutant_name = subdir.replace("KRAS_", "")
            traj_dir = os.path.join(batch_dir, subdir)
            topology = os.path.join(traj_dir, "topology.pdb")
            xtc = os.path.join(traj_dir, "samples.xtc")

            if not os.path.isfile(topology) or not os.path.isfile(xtc):
                print(f"\nSKIPPED {subdir}: topology.pdb or samples.xtc not found.")
                results.append({"mutant": mutant_name, "status": "skipped — files missing"})
                continue

            result = process_one(topology, xtc, mutant_name, df, args.tolerance, args.output_name, rng)
            results.append(result)

        # Print batch summary table
        print(f"\n\n{'═' * 88}")
        print("  BATCH SUMMARY")
        print(f"{'═' * 88}")
        header = f"  {'Mutant':<12}  {'Original':>8}  {'Kept':>6}  {'Removed%':>9}  {'Initial':>8}  {'Final':>8}  {'Target':>8}  {'Diff':>6}  Status"
        print(header)
        print(f"  {'-'*12}  {'-'*8}  {'-'*6}  {'-'*9}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*6}  ------")
        for r in results:
            if "n_original" not in r:
                print(f"  {r['mutant']:<12}  {'—':>8}  {'—':>6}  {'—':>9}  {'—':>8}  {'—':>8}  {'—':>8}  {'—':>6}  {r['status']}")
            else:
                print(
                    f"  {r['mutant']:<12}  {r['n_original']:>8}  {r['n_kept']:>6}  "
                    f"{r['pct_removed']:>8.1f}%  {r['initial_foldedness']:>8.4f}  "
                    f"{r['final_foldedness']:>8.4f}  {r['target_foldedness']:>8.4f}  "
                    f"{r['diff']:>6.4f}  {r['status']}"
                )
        print(f"{'═' * 88}")

    # ── Single-mutant mode ───────────────────────────────────────────────────
    else:
        if args.topology is None:
            parser.error("--topology is required in single-mutant mode.")

        mutant_name = args.mutant_name
        if mutant_name is None:
            dir_name = os.path.basename(os.path.dirname(os.path.abspath(args.xtc)))
            mutant_name = dir_name.replace("KRAS_", "")
            print(f"Inferred mutant name from directory: '{mutant_name}'")

        process_one(args.topology, args.xtc, mutant_name, df, args.tolerance, args.output_name, rng)


if __name__ == "__main__":
    main()
