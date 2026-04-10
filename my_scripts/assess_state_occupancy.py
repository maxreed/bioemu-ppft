#!/usr/bin/env python3
"""
Assess state 2 occupancy of a KRAS mutant ensemble using the contact-based
state classifier (ContactStateClassifier from contact_state_loss.py).

Reports two occupancy measures per mutant:
  binary_occupancy  -- fraction of frames with state_score > 0.5 (hard classification)
  soft_occupancy    -- mean state_score over all frames (continuous)

Output (to stdout, one line per mutant, CSV-compatible):
  name,n_frames,binary_occupancy,soft_occupancy

Usage example:
  python my_scripts/assess_state_occupancy.py \
      --xtc outputs/v21_allVariants/KRAS_WT/samples.xtc \
      --pdb outputs/v21_allVariants/KRAS_WT/topology.pdb \
      --name KRAS_WT \
      --state1_cif potential_references/6mbq.cif \
      --state2_cif potential_references/6xi7.cif
"""

import argparse
import os
import sys

import numpy as np
import torch
import mdtraj

# Allow importing contact_state_loss from the same directory as this script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from contact_state_loss import ContactStateClassifier


def score_trajectory(
    traj: mdtraj.Trajectory,
    classifier: ContactStateClassifier,
    batch_size: int = 128,
) -> np.ndarray:
    """
    Score all frames in a trajectory using the ContactStateClassifier.

    Args:
        traj:       mdtraj.Trajectory (coordinates in nm, as loaded by mdtraj).
        classifier: Instantiated ContactStateClassifier.
        batch_size: Number of frames to process at once (trade-off memory/speed).

    Returns:
        scores: (n_frames,) numpy array of state scores in [0, 1].
    """
    ca_indices = traj.topology.select("name CA")
    ca_coords_nm = traj.xyz[:, ca_indices, :]  # (n_frames, n_CA, 3), nanometres

    n_frames = ca_coords_nm.shape[0]
    scores = np.empty(n_frames, dtype=np.float32)

    classifier.eval()
    with torch.no_grad():
        for start in range(0, n_frames, batch_size):
            end = min(start + batch_size, n_frames)
            batch = torch.tensor(ca_coords_nm[start:end], dtype=torch.float32)
            # coords_in_nm=True: classifier converts nm -> Angstroms internally
            s = classifier(batch, coords_in_nm=True)
            scores[start:end] = s.numpy()

    return scores


def main():
    p = argparse.ArgumentParser(
        description="Report state 2 occupancy (binary and soft) for a KRAS ensemble."
    )
    p.add_argument("--xtc",        required=True,  help="Input trajectory XTC file")
    p.add_argument("--pdb",        required=True,  help="Topology PDB file")
    p.add_argument("--name",       required=True,  help="Mutant name (used in output)")
    p.add_argument("--state1_cif", required=True,  help="State 1 reference mmCIF (e.g. 6mbq.cif)")
    p.add_argument("--state2_cif", required=True,  help="State 2 reference mmCIF (e.g. 6xi7.cif)")
    p.add_argument("--delta_threshold", type=float, default=3.0,
                   help="Min |d_s2 - d_s1| in Å for discriminating contacts (default 3.0)")
    p.add_argument("--batch_size",  type=int,   default=128,
                   help="Frames per scoring batch (default 128)")
    p.add_argument("--header",      action="store_true",
                   help="Print CSV header line before the data row")
    args = p.parse_args()

    # --- Instantiate classifier ---
    classifier = ContactStateClassifier(
        state1_cif=args.state1_cif,
        state2_cif=args.state2_cif,
        delta_threshold=args.delta_threshold,
    )

    # --- Load trajectory ---
    traj = mdtraj.load_xtc(args.xtc, top=args.pdb)
    n_frames = len(traj)

    # --- Score ---
    scores = score_trajectory(traj, classifier, batch_size=args.batch_size)

    # --- Compute occupancies ---
    binary_occ = float(np.mean(scores > 0.5))   # fraction classified as state 2
    soft_occ   = float(np.mean(scores))          # mean continuous score

    # --- Output ---
    if args.header:
        print("name,n_frames,binary_occupancy,soft_occupancy")
    print(f"{args.name},{n_frames},{binary_occ:.6f},{soft_occ:.6f}")


if __name__ == "__main__":
    main()
