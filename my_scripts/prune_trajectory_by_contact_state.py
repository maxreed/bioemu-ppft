#!/usr/bin/env python3
"""
Prune frames from an XTC trajectory until the ensemble binary state-2
fraction matches a target value from a CSV file, using the
ContactStateClassifier from contact_state_loss.py.

The state score uses differential CA-CA contacts between two reference
structures (state 1 = 6MBQ, state 2 = 6XI7). Score 0 = state 1, 1 = state 2.
The percent_folded column in the CSV maps directly to the target state-2
fraction.

A frame is classified as "state 2" if its continuous score exceeds a
threshold (default 0.5). Pruning adjusts the binary state-2 fraction
(number of state-2 frames / total frames) toward the target by randomly
removing state-2 frames when overshooting or non-state-2 frames when
undershooting, until the fraction is within --tolerance of the target.

Single-mutant usage:
    python prune_trajectory_by_contact_state.py \
        --topology outputs/v13_allVariants/KRAS_WT/topology.pdb \
        --xtc      outputs/v13_allVariants/KRAS_WT/samples.xtc \
        --data_csv data/kras_mutants_data.csv

Batch usage (runs on every KRAS_* subdirectory):
    python prune_trajectory_by_contact_state.py \
        --batch_dir outputs/v13_allVariants \
        --data_csv  data/kras_mutants_data.csv

With custom references:
    python prune_trajectory_by_contact_state.py \
        --batch_dir outputs/v13_allVariants \
        --data_csv  data/kras_mutants_data.csv \
        --state1_ref potential_references/6mbq.cif \
        --state2_ref potential_references/6xi7.cif \
        --contact_selection delta3
"""

import argparse
import os
import sys

import numpy as np
import MDAnalysis as mda
import pandas as pd
import torch

# Import ContactStateClassifier from sibling module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from contact_state_loss import ContactStateClassifier


# ── Scoring ──────────────────────────────────────────────────────────────────

def compute_per_frame_state_scores(
    topology: str,
    xtc: str,
    classifier: ContactStateClassifier,
) -> np.ndarray:
    """
    Score every frame using the ContactStateClassifier (soft sigmoid).

    Loads the trajectory with MDAnalysis, extracts CA coordinates per frame,
    and runs the classifier in inference mode (no gradients).

    Returns:
        scores: (n_frames,) array of state scores in [0, 1].
    """
    u = mda.Universe(topology, xtc)
    ca_atoms = u.select_atoms("name CA")
    n_frames = len(u.trajectory)

    scores = np.empty(n_frames, dtype=np.float64)
    with torch.no_grad():
        for frame_idx, ts in enumerate(u.trajectory):
            # MDAnalysis positions are in Angstroms already
            ca_pos = ca_atoms.positions  # (n_ca, 3), float32, Angstroms
            ca_tensor = torch.tensor(ca_pos, dtype=torch.float32).unsqueeze(0)
            score = classifier(ca_tensor, coords_in_nm=False)
            scores[frame_idx] = score.item()

    return scores


# ── Pruning logic ────────────────────────────────────────────────────────────

def prune_to_target(
    state_vals: np.ndarray,
    target: float,
    tolerance: float,
    rng: np.random.Generator = None,
    state2_threshold: float = 0.5,
) -> np.ndarray:
    """
    Randomly remove frames to bring the binary state-2 fraction within
    tolerance of target.

    A frame is "state 2" if its continuous score > state2_threshold.

    Strategy:
      - If fraction > target: randomly remove a state-2 frame.
      - If fraction < target: randomly remove a non-state-2 frame.
    This avoids systematically depleting one tail of the distribution.

    Stops when within tolerance or no eligible candidates remain.
    Returns a boolean mask (True = keep).
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(state_vals)
    is_state2 = state_vals > state2_threshold
    kept = np.ones(n, dtype=bool)
    n_state2 = int(is_state2.sum())
    current_n = n

    def current_fraction():
        return n_state2 / current_n if current_n > 0 else float("nan")

    while abs(current_fraction() - target) > tolerance:
        if current_fraction() > target:
            # Too many state-2 frames -- remove one
            candidates = np.where(kept & is_state2)[0]
        else:
            # Too few state-2 frames -- remove a non-state-2 frame
            candidates = np.where(kept & ~is_state2)[0]

        if len(candidates) == 0:
            break  # no eligible frames left

        idx = rng.choice(candidates)
        kept[idx] = False
        if is_state2[idx]:
            n_state2 -= 1
        current_n -= 1

    return kept


# ── Per-mutant processing ────────────────────────────────────────────────────

def process_one(
    topology: str,
    xtc: str,
    mutant_name: str,
    df: pd.DataFrame,
    classifier: ContactStateClassifier,
    tolerance: float,
    output_name: str,
    rng: np.random.Generator,
) -> dict:
    """
    Process a single mutant trajectory. Returns a result dict for summary reporting.
    """
    print(f"\n{'=' * 68}")
    print(f"  Mutant: {mutant_name}")
    print(f"{'=' * 68}")

    # Look up target state-2 occupancy
    row = df[df["mutant_name"] == mutant_name]
    if row.empty:
        print(f"  SKIPPED: '{mutant_name}' not found in CSV.")
        return {"mutant": mutant_name, "status": "skipped -- not in CSV"}
    target = float(row["percent_folded"].values[0])
    print(f"  Target state-2 occupancy: {target:.4f}")

    # Load trajectory and score
    print(f"  Loading trajectory: {xtc} ...")
    u = mda.Universe(topology, xtc)
    n_original = len(u.trajectory)
    print(f"  {n_original} frames loaded.")

    print("  Computing per-frame contact state scores ...")
    state_vals = compute_per_frame_state_scores(topology, xtc, classifier)
    initial_frac = float((state_vals > 0.5).mean())
    initial_mean = float(state_vals.mean())
    print(f"  Initial state-2 fraction : {initial_frac:.4f}  (continuous mean: {initial_mean:.4f})")

    # Prune
    kept_mask = prune_to_target(state_vals, target, tolerance, rng)
    n_kept = int(kept_mask.sum())
    n_removed = n_original - n_kept
    kept_vals = state_vals[kept_mask]
    final_frac = float((kept_vals > 0.5).mean()) if n_kept > 0 else float("nan")
    final_mean = float(kept_vals.mean()) if n_kept > 0 else float("nan")
    within = abs(final_frac - target) <= tolerance

    print(f"  Frames removed   : {n_removed} / {n_original}  ({100.0 * n_removed / n_original:.1f}%)")
    print(f"  Final state-2 fraction : {final_frac:.4f}  "
          f"(target {target:.4f}, diff {abs(final_frac - target):.4f})")
    print(f"  Final continuous mean  : {final_mean:.4f}")
    if not within:
        print(f"  WARNING: could not reach target within tolerance {tolerance}.")

    # Save pruned trajectory using MDAnalysis
    kept_indices = np.where(kept_mask)[0]
    out_dir = os.path.dirname(os.path.abspath(xtc))
    out_path = os.path.join(out_dir, output_name)

    # Re-load and write only kept frames
    u2 = mda.Universe(topology, xtc)
    with mda.Writer(out_path, n_atoms=u2.atoms.n_atoms) as writer:
        for frame_idx in kept_indices:
            u2.trajectory[frame_idx]
            writer.write(u2.atoms)
    print(f"  Saved: {out_path}")

    return {
        "mutant": mutant_name,
        "n_original": n_original,
        "n_kept": n_kept,
        "pct_removed": 100.0 * n_removed / n_original,
        "initial_frac": initial_frac,
        "final_frac": final_frac,
        "target_frac": target,
        "diff": abs(final_frac - target),
        "final_mean": final_mean,
        "status": "OK" if within else "WARNING: outside tolerance",
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Prune XTC trajectory frames to match a target state-2 occupancy "
                    "using contact-based state classification."
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
        help="Acceptable absolute deviation from target state score (default: 0.02)",
    )
    parser.add_argument(
        "--output_name",
        default="samples_pruned_v2.xtc",
        help="Output XTC filename within each trajectory directory (default: samples_pruned.xtc)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible frame removal (default: None)",
    )
    parser.add_argument(
        "--state1_ref",
        default="potential_references/6mbq.cif",
        help="Path to state 1 reference mmCIF (default: potential_references/6mbq.cif)",
    )
    parser.add_argument(
        "--state2_ref",
        default="potential_references/6xi7.cif",
        help="Path to state 2 reference mmCIF (default: potential_references/6xi7.cif)",
    )
    parser.add_argument(
        "--contact_selection",
        default="delta3",
        help=(
            "Contact selection strategy controlling delta_threshold passed to "
            "ContactStateClassifier. 'delta3' uses delta_threshold=3.0 A (default). "
            "Numeric values are interpreted as the delta_threshold in Angstroms."
        ),
    )
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    df = pd.read_csv(args.data_csv)

    # Parse contact_selection into delta_threshold
    if args.contact_selection == "delta3":
        delta_threshold = 3.0
    else:
        try:
            delta_threshold = float(args.contact_selection)
        except ValueError:
            parser.error(
                f"Unrecognised --contact_selection '{args.contact_selection}'. "
                f"Use 'delta3' or a numeric delta_threshold in Angstroms."
            )

    # Instantiate classifier once (expensive -- loads CIF files)
    print("Instantiating ContactStateClassifier ...")
    print(f"  State 1 ref: {args.state1_ref}")
    print(f"  State 2 ref: {args.state2_ref}")
    print(f"  Delta threshold: {delta_threshold} A")
    classifier = ContactStateClassifier(
        state1_cif=args.state1_ref,
        state2_cif=args.state2_ref,
        delta_threshold=delta_threshold,
    )
    print(f"  {classifier}")

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
                results.append({"mutant": mutant_name, "status": "skipped -- files missing"})
                continue

            result = process_one(topology, xtc, mutant_name, df, classifier,
                                 args.tolerance, args.output_name, rng)
            results.append(result)

        # Print batch summary table
        print(f"\n\n{'=' * 96}")
        print("  BATCH SUMMARY")
        print(f"{'=' * 96}")
        header = (f"  {'Mutant':<12}  {'Original':>8}  {'Kept':>6}  {'Removed%':>9}  "
                  f"{'InitFrac':>8}  {'FinalFrac':>9}  {'Target':>8}  {'Diff':>6}  Status")
        print(header)
        print(f"  {'-'*12}  {'-'*8}  {'-'*6}  {'-'*9}  {'-'*8}  {'-'*9}  {'-'*8}  {'-'*6}  ------")
        for r in results:
            if "n_original" not in r:
                print(f"  {r['mutant']:<12}  {'--':>8}  {'--':>6}  {'--':>9}  "
                      f"{'--':>8}  {'--':>9}  {'--':>8}  {'--':>6}  {r['status']}")
            else:
                print(
                    f"  {r['mutant']:<12}  {r['n_original']:>8}  {r['n_kept']:>6}  "
                    f"{r['pct_removed']:>8.1f}%  {r['initial_frac']:>8.4f}  "
                    f"{r['final_frac']:>9.4f}  {r['target_frac']:>8.4f}  "
                    f"{r['diff']:>6.4f}  {r['status']}"
                )
        print(f"{'=' * 96}")

    # ── Single-mutant mode ───────────────────────────────────────────────────
    else:
        if args.topology is None:
            parser.error("--topology is required in single-mutant mode.")

        mutant_name = args.mutant_name
        if mutant_name is None:
            dir_name = os.path.basename(os.path.dirname(os.path.abspath(args.xtc)))
            mutant_name = dir_name.replace("KRAS_", "")
            print(f"Inferred mutant name from directory: '{mutant_name}'")

        process_one(args.topology, args.xtc, mutant_name, df, classifier,
                    args.tolerance, args.output_name, rng)


if __name__ == "__main__":
    main()
