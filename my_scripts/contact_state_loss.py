#!/usr/bin/env python3
"""
Differentiable contact-based state classifier and PPFT loss for KRAS fine-tuning.

Replaces the 3-pair sigmoid foldedness score in finetune_foldedness_v4b_hybrid_loss.py
with a richer metric based on differential CA-CA contacts between two reference
structures (state 1 = 6MBQ, state 2 = 6XI7).

State score:
    0 = state 1 (T35-out, minor/excited state)
    1 = state 2 (effector-binding-competent, active)

The metric:
    1. Identify contacts that discriminate the two states (delta >= delta_threshold
       between reference distances).
    2. For each frame, compute FNC_2 (fraction of state-2-specific native contacts
       formed) and FNC_1 (fraction of state-1-specific native contacts formed)
       using a soft sigmoid contact criterion for differentiability.
    3. state_score = (FNC_2 - FNC_1 + 1) / 2  in [0, 1].
"""

import numpy as np
import torch
import torch.nn as nn

from Bio.PDB import MMCIFParser


# =============================================================================
# Reference structure parsing (one-time, at init)
# =============================================================================

def _extract_ca_from_cif(cif_path: str, chain_id: str, max_resid_auth: int = 169):
    """
    Extract CA coordinates and author residue IDs from an mmCIF file.

    Returns (resids, resnames, coords) where resids and coords are numpy arrays
    and resnames is a list of 3-letter strings.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("ref", cif_path)
    model = structure[0]

    resids, resnames, coords = [], [], []
    for chain in model:
        if chain.id != chain_id:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if "CA" not in residue:
                continue
            auth_seq_id = residue.id[1]
            if auth_seq_id > max_resid_auth:
                continue
            resids.append(auth_seq_id)
            resnames.append(residue.resname)
            coords.append(residue["CA"].get_vector().get_array())

    return np.array(resids, dtype=int), resnames, np.array(coords, dtype=np.float64)


def _identify_kras_chain_6xi7(cif_path: str) -> str:
    """Return the author chain ID corresponding to KRAS in 6XI7 (longest chain)."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("6xi7", cif_path)
    model = structure[0]
    chain_sizes = {}
    for chain in model:
        chain_sizes[chain.id] = sum(1 for r in chain if r.id[0] == " " and "CA" in r)
    return max(chain_sizes, key=chain_sizes.get)


def _compute_distance_map(coords, resids, contact_cutoff, seq_sep_min=4):
    """
    Compute CA-CA contacts and full pairwise distance map.

    Returns:
        contacts: set of (auth_i, auth_j) pairs with d < contact_cutoff
        dist_map: dict of (auth_i, auth_j) -> distance for all pairs with |i-j| >= seq_sep_min
    """
    n = len(resids)
    diff = coords[:, None, :] - coords[None, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=-1))

    contacts = set()
    dist_map = {}
    for i in range(n):
        for j in range(i + 1, n):
            if abs(int(resids[i]) - int(resids[j])) < seq_sep_min:
                continue
            d = dists[i, j]
            pair = (int(resids[i]), int(resids[j]))
            dist_map[pair] = d
            if d < contact_cutoff:
                contacts.add(pair)

    return contacts, dist_map


# =============================================================================
# ContactStateClassifier
# =============================================================================

class ContactStateClassifier(nn.Module):
    """
    Computes a differentiable state score (0=state1, 1=state2) for KRAS
    using differential CA-CA contacts between two reference structures.

    Designed to be instantiated once (expensive -- loads CIF files and
    computes contact pairs) then called repeatedly during training.

    The forward pass is pure PyTorch and fully differentiable w.r.t. input
    CA coordinates.

    Contact "formed" criterion (soft, differentiable):
        contact_formed(d, d_ref) = sigmoid(-k * (d - threshold_factor * d_ref))

    where threshold_factor=1.2 follows BioEmu's FNC convention (a contact is
    considered formed if the distance is within 120% of the reference distance).
    """

    def __init__(
        self,
        state1_cif: str,
        state2_cif: str,
        delta_threshold: float = 3.0,
        contact_distance_cutoff: float = 10.0,
        seq_sep_min: int = 4,
        max_resid_auth: int = 169,
        sigmoid_k: float = 5.0,
        threshold_factor: float = 1.2,
        state1_chain: str = "A",
        state2_chain: str | None = None,
    ):
        """
        Args:
            state1_cif: Path to state 1 reference mmCIF (e.g. 6MBQ).
            state2_cif: Path to state 2 reference mmCIF (e.g. 6XI7).
            delta_threshold: Minimum |d_s2 - d_s1| in Angstroms for a contact to be
                             considered state-discriminating (default 3.0 = "delta3").
            contact_distance_cutoff: CA-CA distance cutoff for defining contacts (A).
            seq_sep_min: Minimum sequence separation |i-j| for contact pairs.
            max_resid_auth: Maximum author residue ID to include.
            sigmoid_k: Steepness of the soft contact sigmoid (higher = sharper).
            threshold_factor: Contact is "formed" when d < threshold_factor * d_ref.
            state1_chain: Author chain ID for state 1 reference.
            state2_chain: Author chain ID for state 2 reference. If None, auto-detect
                          longest chain (appropriate for 6XI7 KRAS/RAF complex).
        """
        super().__init__()

        self.sigmoid_k = sigmoid_k
        self.threshold_factor = threshold_factor

        # --- Parse reference structures ---
        s1_resids, s1_resnames, s1_coords = _extract_ca_from_cif(
            state1_cif, state1_chain, max_resid_auth
        )
        if state2_chain is None:
            state2_chain = _identify_kras_chain_6xi7(state2_cif)
        s2_resids, s2_resnames, s2_coords = _extract_ca_from_cif(
            state2_cif, state2_chain, max_resid_auth
        )

        # --- Compute contacts and distance maps ---
        contacts_s1, dists_s1 = _compute_distance_map(
            s1_coords, s1_resids, contact_distance_cutoff, seq_sep_min
        )
        contacts_s2, dists_s2 = _compute_distance_map(
            s2_coords, s2_resids, contact_distance_cutoff, seq_sep_min
        )

        common_pairs = set(dists_s1.keys()) & set(dists_s2.keys())

        # State 2-specific: in contact in s2 but NOT in s1
        state2_specific = sorted(
            p for p in contacts_s2
            if p in common_pairs and p not in contacts_s1
        )
        # State 1-specific: in contact in s1 but NOT in s2
        state1_specific = sorted(
            p for p in contacts_s1
            if p in common_pairs and p not in contacts_s2
        )

        # Apply delta threshold filter
        if delta_threshold > 0.0:
            state2_specific = [
                p for p in state2_specific
                if abs(dists_s2[p] - dists_s1[p]) >= delta_threshold
            ]
            state1_specific = [
                p for p in state1_specific
                if abs(dists_s1[p] - dists_s2[p]) >= delta_threshold
            ]

        self.n_s2_contacts = len(state2_specific)
        self.n_s1_contacts = len(state1_specific)

        if self.n_s2_contacts == 0 and self.n_s1_contacts == 0:
            raise ValueError(
                f"No discriminating contacts found with delta_threshold={delta_threshold}. "
                f"Check reference structures or lower the threshold."
            )

        # --- Build residue index mapping ---
        # BioEmu topology uses 0-indexed residues: topology_resid = auth_seq_id - 1
        # The forward() method receives CA coords indexed 0..N-1 in topology order.
        # We need to map auth_seq_id pairs to 0-indexed positions.
        #
        # NOTE on mapping mismatches:
        #   6MBQ has a GG tag at auth_seq_id 0,1 which doesn't match topology resid 0 (MET).
        #   6XI7 has SER at auth_seq_id 118 vs CYS at topology resid 117.
        #   These mismatches are handled by the residue name check below, same as
        #   classify_states_contact_metric.py's build_residue_mapping().

        # For the loss module we don't have the topology PDB available at init time
        # in the general case. Instead, we store the auth_seq_id pairs and apply
        # the offset (auth_seq_id - 1) to get 0-indexed residue indices.
        # This is valid because BioEmu's KRAS topology is numbered 0-168 = auth 1-169.
        # Contact pairs that involve residues missing from the topology (e.g. 6MBQ
        # auth_seq_id 0) naturally won't appear because those residues are not in
        # the contacts set (they are filtered by the common_pairs intersection).

        def _auth_pairs_to_zero_indexed(pairs, ref_dists_dict):
            """Convert auth_seq_id pairs to 0-indexed and extract reference distances."""
            idx_i, idx_j, d_ref = [], [], []
            for (ri, rj) in pairs:
                # 0-indexed = auth_seq_id - 1
                idx_i.append(ri - 1)
                idx_j.append(rj - 1)
                d_ref.append(ref_dists_dict[(ri, rj)])
            return (
                torch.tensor(idx_i, dtype=torch.long),
                torch.tensor(idx_j, dtype=torch.long),
                torch.tensor(d_ref, dtype=torch.float32),
            )

        # State 2-specific contacts: reference distances from state 2 (6XI7)
        s2_idx_i, s2_idx_j, s2_dref = _auth_pairs_to_zero_indexed(
            state2_specific, dists_s2
        )
        # State 1-specific contacts: reference distances from state 1 (6MBQ)
        s1_idx_i, s1_idx_j, s1_dref = _auth_pairs_to_zero_indexed(
            state1_specific, dists_s1
        )

        # Register as buffers so they move with .to(device) and .cuda()
        self.register_buffer("s2_idx_i", s2_idx_i)
        self.register_buffer("s2_idx_j", s2_idx_j)
        self.register_buffer("s2_dref", s2_dref)  # Angstroms
        self.register_buffer("s1_idx_i", s1_idx_i)
        self.register_buffer("s1_idx_j", s1_idx_j)
        self.register_buffer("s1_dref", s1_dref)  # Angstroms

        # Store contact pair info for inspection
        self._state2_pairs = state2_specific
        self._state1_pairs = state1_specific

    def forward(
        self,
        ca_coords: torch.Tensor,
        residue_mask: torch.Tensor | None = None,
        coords_in_nm: bool = True,
    ) -> torch.Tensor:
        """
        Compute differentiable state scores for a batch of structures.

        Args:
            ca_coords: (batch, N_residues, 3) CA coordinates.
                       Units determined by coords_in_nm flag.
            residue_mask: (batch, N_residues) bool tensor. If provided, masked
                          residues are excluded from contact calculations. Currently
                          unused but reserved for future multi-length batching.
            coords_in_nm: If True (default, matching BioEmu convention), coordinates
                          are in nanometres and will be converted to Angstroms
                          internally. Set False if already in Angstroms.

        Returns:
            state_score: (batch,) tensor, values in [0, 1].
                         0 = state 1, 1 = state 2.
        """
        if coords_in_nm:
            ca_coords = ca_coords * 10.0  # nm -> Angstroms

        # --- Compute pairwise distances for state 2-specific contacts ---
        fnc2 = self._compute_fnc(
            ca_coords, self.s2_idx_i, self.s2_idx_j, self.s2_dref
        )  # (batch,)

        # --- Compute pairwise distances for state 1-specific contacts ---
        fnc1 = self._compute_fnc(
            ca_coords, self.s1_idx_i, self.s1_idx_j, self.s1_dref
        )  # (batch,)

        # state_score = (FNC_2 - FNC_1 + 1) / 2, maps [-1, 1] -> [0, 1]
        state_score = (fnc2 - fnc1 + 1.0) / 2.0
        return state_score

    def _compute_fnc(
        self,
        ca_coords: torch.Tensor,
        idx_i: torch.Tensor,
        idx_j: torch.Tensor,
        d_ref: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute soft fraction of native contacts for a set of contact pairs.

        Soft contact criterion:
            sigma(-k * (d - threshold_factor * d_ref))
        where sigma is the logistic sigmoid. When d < threshold_factor * d_ref,
        the argument is positive and sigmoid -> 1 (contact formed). When
        d >> threshold_factor * d_ref, sigmoid -> 0 (contact broken).

        Args:
            ca_coords: (batch, N_residues, 3), Angstroms.
            idx_i, idx_j: (n_contacts,) 0-indexed residue indices.
            d_ref: (n_contacts,) reference distances in Angstroms.

        Returns:
            fnc: (batch,) mean soft contact formation fraction.
        """
        if len(idx_i) == 0:
            return torch.zeros(ca_coords.shape[0], device=ca_coords.device)

        # Extract positions: (batch, n_contacts, 3)
        pos_i = ca_coords[:, idx_i, :]
        pos_j = ca_coords[:, idx_j, :]

        # Pairwise distances: (batch, n_contacts)
        # Using norm with a small epsilon inside sqrt for numerical stability
        diff = pos_i - pos_j
        dists = torch.sqrt(torch.sum(diff ** 2, dim=-1) + 1e-8)

        # Soft contact formation: sigmoid(-k * (d - threshold * d_ref))
        # d_ref is (n_contacts,), broadcast over batch dimension
        contact_formed = torch.sigmoid(
            -self.sigmoid_k * (dists - self.threshold_factor * d_ref.unsqueeze(0))
        )  # (batch, n_contacts)

        # Fraction of native contacts formed: mean over contacts
        fnc = contact_formed.mean(dim=-1)  # (batch,)
        return fnc

    def __repr__(self) -> str:
        return (
            f"ContactStateClassifier("
            f"n_s2_contacts={self.n_s2_contacts}, "
            f"n_s1_contacts={self.n_s1_contacts}, "
            f"sigmoid_k={self.sigmoid_k}, "
            f"threshold_factor={self.threshold_factor})"
        )


# =============================================================================
# PPFT loss: unbiased squared-mean-error estimator
# =============================================================================

def compute_ppft_loss(
    state_scores: torch.Tensor,
    target_state2_fraction: float,
    n_replications: int,
) -> torch.Tensor:
    """
    Unbiased squared-mean-error PPFT loss.
    Matches the formulation in finetune_foldedness_v4b_hybrid_loss.py.

    For independent samples with s_i = state_score of sample i and target mu:
        E[(mean(s) - mu)^2] is estimated unbiasedly by the cross-product form:
        ( [sum_i (s_i - mu)]^2 - sum_i (s_i - mu)^2 ) / (n * (n - 1))

    This is mathematically equivalent to:
        (1 / (n*(n-1))) * sum_{i != j} (s_i - mu)(s_j - mu)

    which is an unbiased estimator of (E[s] - mu)^2 when samples are i.i.d.
    from the model.

    Args:
        state_scores: (n_replications,) tensor of per-sample state scores.
        target_state2_fraction: float, target fraction in state 2 (in [0, 1]).
        n_replications: int, number of replications (must be >= 2).

    Returns:
        Scalar loss tensor with grad_fn intact.
    """
    n = state_scores.shape[0]
    assert n >= 2, "Need at least 2 replications for the unbiased estimator."
    assert n == n_replications, (
        f"state_scores has {n} elements but n_replications={n_replications}"
    )

    target = torch.tensor(
        target_state2_fraction, dtype=state_scores.dtype, device=state_scores.device
    )
    diff = state_scores - target
    sum_diff = torch.sum(diff)
    # Unbiased cross-product estimator: (sum^2 - sum_of_squares) / (n * (n-1))
    return (sum_diff ** 2 - torch.sum(diff ** 2)) / (n * (n - 1))


# =============================================================================
# Test block
# =============================================================================

if __name__ == "__main__":
    import MDAnalysis as mda

    STATE1_CIF = "potential_references/6mbq.cif"
    STATE2_CIF = "potential_references/6xi7.cif"
    TOPOLOGY = "outputs/bioemu1p1_KRAS/KRAS_WT/topology.pdb"
    TRAJECTORY = "outputs/bioemu1p1_KRAS/KRAS_WT/samples.xtc"

    print("=" * 70)
    print("Instantiating ContactStateClassifier (delta3, 29 contacts)")
    print("=" * 70)

    classifier = ContactStateClassifier(
        state1_cif=STATE1_CIF,
        state2_cif=STATE2_CIF,
        delta_threshold=3.0,
        contact_distance_cutoff=10.0,
    )
    print(classifier)
    print(f"  State 2-specific contacts: {classifier.n_s2_contacts}")
    print(f"  State 1-specific contacts: {classifier.n_s1_contacts}")
    print(f"  State 2 pairs (auth_seq_id): {classifier._state2_pairs}")
    print(f"  State 1 pairs (auth_seq_id): {classifier._state1_pairs}")

    # Load trajectory with MDAnalysis
    print(f"\nLoading trajectory: {TRAJECTORY}")
    u = mda.Universe(TOPOLOGY, TRAJECTORY)
    ca_atoms = u.select_atoms("name CA")
    n_ca = len(ca_atoms)
    n_frames = len(u.trajectory)
    print(f"  {n_ca} CA atoms, {n_frames} frames")

    # Score all frames
    scores_soft = np.zeros(n_frames)
    scores_hard = np.zeros(n_frames)

    # Also compute hard-threshold scores for comparison with classify_states_contact_metric.py
    # (which uses d < 1.2 * d_ref as a hard cutoff)
    s2_idx_i = classifier.s2_idx_i.numpy()
    s2_idx_j = classifier.s2_idx_j.numpy()
    s2_dref = classifier.s2_dref.numpy()
    s1_idx_i = classifier.s1_idx_i.numpy()
    s1_idx_j = classifier.s1_idx_j.numpy()
    s1_dref = classifier.s1_dref.numpy()

    for frame_idx, ts in enumerate(u.trajectory):
        ca_pos = ca_atoms.positions  # (n_ca, 3) in Angstroms

        # --- Soft scores (differentiable sigmoid) ---
        ca_tensor = torch.tensor(ca_pos, dtype=torch.float32).unsqueeze(0)  # (1, n_ca, 3)
        with torch.no_grad():
            score_soft = classifier(ca_tensor, coords_in_nm=False)
        scores_soft[frame_idx] = score_soft.item()

        # --- Hard scores (for comparison with classify_states_contact_metric.py) ---
        if len(s2_idx_i) > 0:
            diff_s2 = ca_pos[s2_idx_i] - ca_pos[s2_idx_j]
            d_s2 = np.sqrt(np.sum(diff_s2 ** 2, axis=1))
            fnc2_hard = np.mean(d_s2 < 1.2 * s2_dref)
        else:
            fnc2_hard = 0.0
        if len(s1_idx_i) > 0:
            diff_s1 = ca_pos[s1_idx_i] - ca_pos[s1_idx_j]
            d_s1 = np.sqrt(np.sum(diff_s1 ** 2, axis=1))
            fnc1_hard = np.mean(d_s1 < 1.2 * s1_dref)
        else:
            fnc1_hard = 0.0
        scores_hard[frame_idx] = (fnc2_hard - fnc1_hard + 1.0) / 2.0

        if (frame_idx + 1) % 200 == 0 or frame_idx == 0:
            print(f"  Frame {frame_idx + 1:5d}/{n_frames}: "
                  f"soft={scores_soft[frame_idx]:.4f}  hard={scores_hard[frame_idx]:.4f}")

    print(f"\n{'=' * 70}")
    print("Results (delta3 selection, 29 contacts)")
    print(f"{'=' * 70}")
    print(f"  Soft sigmoid (k={classifier.sigmoid_k}):")
    print(f"    Mean:  {np.mean(scores_soft):.4f}")
    print(f"    Std:   {np.std(scores_soft):.4f}")
    print(f"    Range: [{np.min(scores_soft):.4f}, {np.max(scores_soft):.4f}]")
    print(f"  Hard threshold (d < 1.2 * d_ref):")
    print(f"    Mean:  {np.mean(scores_hard):.4f}")
    print(f"    Std:   {np.std(scores_hard):.4f}")
    print(f"    Range: [{np.min(scores_hard):.4f}, {np.max(scores_hard):.4f}]")
    print(f"  Correlation (soft vs hard): {np.corrcoef(scores_soft, scores_hard)[0, 1]:.4f}")

    # Test PPFT loss computation
    print(f"\n{'=' * 70}")
    print("PPFT loss test")
    print(f"{'=' * 70}")
    # Simulate a small batch of 8 state scores
    test_scores = torch.tensor(scores_soft[:8], dtype=torch.float32)
    target = 0.5  # example target
    loss = compute_ppft_loss(test_scores, target_state2_fraction=target, n_replications=8)
    print(f"  Test scores (first 8 frames): {test_scores.tolist()}")
    print(f"  Target state 2 fraction: {target}")
    print(f"  PPFT loss: {loss.item():.6f}")

    # Verify gradient flow
    print(f"\n{'=' * 70}")
    print("Gradient flow test")
    print(f"{'=' * 70}")
    ca_test = torch.tensor(
        ca_atoms.positions, dtype=torch.float32
    ).unsqueeze(0).requires_grad_(True)  # (1, n_ca, 3), Angstroms
    score = classifier(ca_test, coords_in_nm=False)
    score.backward()
    grad_norm = ca_test.grad.norm().item()
    print(f"  Score: {score.item():.4f}")
    print(f"  Gradient norm w.r.t. CA coords: {grad_norm:.6f}")
    print(f"  Gradient is nonzero: {grad_norm > 0}")

    print("\nDone.")
