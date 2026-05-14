"""Utility functions extracted from newgraphene.py for band structure projections.

This module provides helper functions for:
- Basis order auto-detection (atom_spin vs spin_first)
- Degenerate eigenvalue identification and gauge fixing
- Atom selection parsing and orbital index mapping
- Projection mode conversion (expectation, weight_up/down, polarization)
"""

import numpy as np
from typing import Tuple, List, Dict, Optional, Union


def guess_basis_order(
    ham_r: Dict[Tuple[int, int, int], np.ndarray],
    n_orbitals: int,
    threshold: float = 1e-3,
) -> str:
    """Auto-detect basis ordering from real-space Hamiltonian H(R).

    Heuristic: For spin_first basis, the first half contains all spin-up states
    and second half contains spin-down. If basis is spin_first, there should be
    weak coupling between halves (no spin-orbit coupling). We detect this by
    comparing cross-half vs within-half coupling strengths.

    Parameters
    ----------
    ham_r : dict
        Real-space Hamiltonian as {(R1, R2, R3): H_matrix}
    n_orbitals : int
        Total number of orbitals (dimension of Hamiltonian matrix)
    threshold : float, default=1e-3
        Threshold for cross/within ratio to classify as spin_first.
        If max(cross_coupling) < threshold * max(within_coupling), classifies as spin_first.

    Returns
    -------
    str
        Either "atom_spin" or "spin_first"

    Notes
    -----
    - "atom_spin": Basis ordered as [atom1_up, atom1_down, atom2_up, atom2_down, ...]
    - "spin_first": Basis ordered as [atom1_up, atom2_up, ..., atom1_down, atom2_down, ...]
    """
    H0 = ham_r.get((0, 0, 0), None)

    # Default to atom_spin if H(0,0,0) not available or odd dimension
    if H0 is None or (n_orbitals % 2 != 0):
        return "atom_spin"

    # Check spin-first hypothesis: split at n_half
    n_half = n_orbitals // 2
    cross = np.max(np.abs(H0[:n_half, n_half:]))
    within = max(
        np.max(np.abs(H0[:n_half, :n_half])), np.max(np.abs(H0[n_half:, n_half:]))
    )

    # If cross coupling is tiny vs within-block coupling -> likely spin_first
    if within > 0 and cross < threshold * within:
        return "spin_first"

    return "atom_spin"


def degeneracy_groups(evals: np.ndarray, tol: float = 1e-8) -> List[Tuple[int, int]]:
    """Identify groups of degenerate eigenvalues.

    Parameters
    ----------
    evals : np.ndarray
        1D array of eigenvalues (assumed sorted)
    tol : float, default=1e-8
        Tolerance for considering eigenvalues equal

    Returns
    -------
    list of tuple
        List of (start_idx, end_idx) for each degeneracy group

    Examples
    --------
    >>> evals = np.array([0.0, 0.0, 0.5, 0.5, 0.5, 2.0])
    >>> degeneracy_groups(evals)
    [(0, 2), (2, 5), (5, 6)]
    """
    groups = []
    start = 0
    for i in range(1, len(evals)):
        if abs(evals[i] - evals[i - 1]) > tol:
            groups.append((start, i))
            start = i
    groups.append((start, len(evals)))
    return groups


def gauge_fix_evecs(
    evals: np.ndarray,
    evecs: np.ndarray,
    sop_sel: np.ndarray,
    sel_idx: np.ndarray,
    deg_tol: float = 1e-8,
) -> np.ndarray:
    """Apply gauge fixing to eigenvectors in degenerate subspaces.

    In each degenerate block, rotates eigenvectors to maximize the spin
    expectation value in the selected atom subspace. This ensures smooth
    band coloring across degeneracies.

    Parameters
    ----------
    evals : np.ndarray
        1D array of eigenvalues at a given k-point
    evecs : np.ndarray
        2D array of eigenvectors (columns are eigenstates)
    sop_sel : np.ndarray
        Spin operator (Sx, Sy, or Sz) restricted to selected atom subspace
    sel_idx : np.ndarray or None
        Indices of selected atoms in orbital basis
    deg_tol : float, default=1e-8
        Tolerance for considering eigenvalues degenerate

    Returns
    -------
    np.ndarray
        Rotated eigenvectors (same shape as input evecs)

    Notes
    -----
    If sel_idx is None, returns evecs unchanged.
    Rotations are performed separately for each degenerate block.
    """
    if sel_idx is None or sop_sel is None:
        return evecs

    # Restrict eigenvectors to selected atom subspace
    Usel_full = evecs[sel_idx, :]  # (n_selected, n_bands)
    groups = degeneracy_groups(evals, tol=deg_tol)

    for start, end in groups:
        m = end - start
        if m <= 1:
            continue

        # Extract degenerate block
        Usel = Usel_full[:, start:end]  # (n_selected, m)

        # Construct projected operator in degenerate subspace
        Osub = Usel.conj().T @ (sop_sel @ Usel)  # (m, m)

        # Diagonalize to find optimal rotation
        w, W = np.linalg.eigh(Osub)  # Hermitian diagonalization

        # Sort by eigenvalue (largest first = highest spin projection)
        order = np.argsort(w.real)[::-1]
        W = W[:, order]

        # Apply rotation to full eigenvectors in this block
        evecs[:, start:end] = evecs[:, start:end] @ W

    return evecs


def atom_orbital_indices(
    atom_indices: List[int],
    orbitals_per_atom: int,
    n_orb_per_spin: int,
    basis_order: str,
    n_atoms_total: int,
) -> np.ndarray:
    """Map selected atoms to their orbital indices in the Hamiltonian basis.

    Handles different basis orderings:
    - "atom_spin": [atom1_up, atom1_down, atom2_up, atom2_down, ...]
    - "spin_first": [atom1_up, atom2_up, ..., atom1_down, atom2_down, ...]

    Parameters
    ----------
    atom_indices : list of int
        0-based atom indices to select
    orbitals_per_atom : int
        Number of orbitals per atom (including spin)
    n_orb_per_spin : int
        Number of orbitals per atom per spin (= orbitals_per_atom / 2)
    basis_order : str
        Either "atom_spin" or "spin_first"
    n_atoms_total : int
        Total number of atoms

    Returns
    -------
    np.ndarray
        1D array of 0-based orbital indices

    Raises
    ------
    ValueError
        If basis_order is unrecognized
    """
    idx = []

    if basis_order == "atom_spin":
        # Each atom's orbitals are contiguous: [up, down]
        for a in atom_indices:
            base = a * orbitals_per_atom
            idx.extend(range(base, base + orbitals_per_atom))
        return np.array(idx, dtype=int)

    elif basis_order == "spin_first":
        # First half: all spin-up, Second half: all spin-down
        Nup = n_atoms_total * n_orb_per_spin
        for a in atom_indices:
            up_base = a * n_orb_per_spin
            dn_base = Nup + a * n_orb_per_spin
            idx.extend(range(up_base, up_base + n_orb_per_spin))
            idx.extend(range(dn_base, dn_base + n_orb_per_spin))
        return np.array(idx, dtype=int)

    else:
        raise ValueError(
            f"Unknown basis_order: {basis_order}. Must be 'atom_spin' or 'spin_first'."
        )


def parse_atom_selection(
    selection: Union[str, int, List[int], np.ndarray],
    n_atoms: int,
    positions: Optional[List[Tuple]] = None,
) -> Tuple[List[int], str]:
    """Parse flexible atom selection specification.

    Parameters
    ----------
    selection : str, int, or list
        Selection specification:
        - "all": Select all atoms
        - "C" (str): Select atoms with given label
        - 1 (int, 1-based): Select single atom
        - [1, 2, 3] (list): Select multiple atoms
    n_atoms : int
        Total number of atoms
    positions : list of tuple, optional
        List of (label, x, y, z) tuples for label-based matching

    Returns
    -------
    tuple
        (atom_indices_0based, description_string)

    Raises
    ------
    ValueError
        If selection is invalid or out of range
    """
    if isinstance(selection, str):
        s = selection.strip()
        if s.lower() == "all":
            return list(range(n_atoms)), "all-atoms"

        # Label-based selection
        if positions is None:
            raise ValueError(f"Cannot match label '{s}' without positions list")

        idx = [a for a, (label, _, _, _) in enumerate(positions) if label == s]
        if not idx:
            available = sorted(set(p[0] for p in positions))
            raise ValueError(f"Label '{s}' matched no atoms. Available: {available}")
        return idx, f"label-{s}"

    if isinstance(selection, (int, np.integer)):
        a = int(selection)
        # Try 1-based indexing first
        if 1 <= a <= n_atoms:
            return [a - 1], f"atom{a}"
        # Fallback to 0-based
        if 0 <= a < n_atoms:
            return [a], f"atom{a + 1}"
        raise ValueError(f"Atom index {a} out of range [0, {n_atoms-1}]")

    if isinstance(selection, (list, tuple, np.ndarray)):
        atoms = [int(x) for x in selection]

        # Check if 1-based
        if all(1 <= a <= n_atoms for a in atoms):
            atoms_0based = [a - 1 for a in atoms]
            return atoms_0based, "atoms" + "-".join(str(a) for a in atoms)

        # Check if 0-based
        if all(0 <= a < n_atoms for a in atoms):
            return atoms, "atoms" + "-".join(str(a + 1) for a in atoms)

        raise ValueError(f"Atom indices {atoms} out of range [0, {n_atoms-1}]")

    raise ValueError(f"Unsupported selection type: {type(selection)}")


def get_projection_vrange(mode: str) -> Tuple[float, float]:
    """Get value range for a projection mode.

    Parameters
    ----------
    mode : str
        Projection mode: "expectation", "weight_up", "weight_down", or "polarization"

    Returns
    -------
    tuple of float
        (vmin, vmax) for colorbar scaling

    Raises
    ------
    ValueError
        If mode is unrecognized
    """
    modes = {
        "expectation": (-0.5, +0.5),
        "weight_up": (0.0, 1.0),
        "weight_down": (0.0, 1.0),
        "polarization": (-1.0, +1.0),
    }

    if mode not in modes:
        raise ValueError(
            f"Unknown projection mode: {mode}. Must be one of {list(modes.keys())}"
        )

    return modes[mode]


def compute_quantity_from_sexp(sexp: np.ndarray, mode: str) -> np.ndarray:
    """Convert spin expectation value to desired projection quantity.

    Parameters
    ----------
    sexp : np.ndarray
        Spin expectation value <S_axis> (typically in [-0.5, +0.5] for 1/2 spin)
    mode : str
        Conversion mode:
        - "expectation": Direct expectation (no conversion)
        - "weight_up": Spin-up weight = 0.5 + <S_axis>
        - "weight_down": Spin-down weight = 0.5 - <S_axis>
        - "polarization": Signed polarization = 2.0 * <S_axis>

    Returns
    -------
    np.ndarray
        Converted quantity (same shape as input)

    Raises
    ------
    ValueError
        If mode is unrecognized

    Notes
    -----
    All conversions preserve NaN values. Mathematically equivalent to:
    - "weight_up" + "weight_down" = 1.0 (conservation of probability)
    - "polarization" = "weight_up" - "weight_down" (definition)
    """
    if mode == "expectation":
        return sexp
    elif mode == "weight_up":
        return 0.5 + sexp
    elif mode == "weight_down":
        return 0.5 - sexp
    elif mode == "polarization":
        return 2.0 * sexp
    else:
        raise ValueError(f"Unknown projection mode: {mode}")


def stabilize_degenerate_projections(
    energies: np.ndarray,
    numerators: np.ndarray,
    denominators: np.ndarray,
    mode: str = "block_average",
    energy_tol: float = 1e-6,
    normalize: bool = True,
    min_abs_value: float = 1e-5,
    sum_tol: float = 5e-4,
    mag_tol: float = 5e-4,
) -> np.ndarray:
    """Stabilize projection values inside degenerate energy blocks.

    For degenerate bands, projection values can oscillate between basis choices.
    This applies stabilization strategies to produce smooth, interpretable results.

    Parameters
    ----------
    energies : np.ndarray
        Band energies at a single k-point (nbands,)
    numerators : np.ndarray
        Numerator values (nbands,), typically :math:`\\langle\\psi|P O P|\\psi\\rangle`
    denominators : np.ndarray
        Denominator values (nbands,), typically :math:`\\langle\\psi|P|\\psi\\rangle`
    mode : str, default="block_average"
        Stabilization mode. One of: "none", "block_average", "pair_cancel".
    energy_tol : float, default=1e-6
        Energy tolerance for grouping degenerate bands
    normalize : bool, default=True
        If True, return num/den; if False, return num
    min_abs_value : float, default=1e-5
        For "pair_cancel": ignore values with :math:`|v| < min\\_abs\\_value`
    sum_tol : float, default=5e-4
        For "pair_cancel": require :math:`|v_i + v_j| < sum\\_tol`
    mag_tol : float, default=5e-4
        For "pair_cancel": require :math:`||v_i| - |v_j|| < mag\\_tol`

    Returns
    -------
    np.ndarray
        Stabilized projections (nbands,)

    Notes
    -----
    - "block_average": average within degenerate manifolds
    - "pair_cancel": zero near-opposite pairs to reduce noise

    Author
    ------
    Ramon Cardias
    """
    mode = mode.lower().strip()
    
    if mode == "none":
        if normalize:
            expv = np.where(
                np.abs(denominators) > 1e-14,
                numerators / denominators,
                np.nan
            )
        else:
            expv = numerators.copy()
        return expv
    
    # Compute base expectation value
    if normalize:
        expv = np.where(
            np.abs(denominators) > 1e-14,
            numerators / denominators,
            np.nan
        )
    else:
        expv = numerators.copy()
    
    # Find degenerate groups
    groups = degeneracy_groups(energies, tol=energy_tol)
    out = expv.copy()
    
    if mode == "block_average":
        for i0, i1 in groups:
            m = i1 - i0
            if m <= 1:
                continue
            
            # Replace entire block with subspace average
            if normalize:
                sden = np.nansum(denominators[i0:i1])
                if sden > 1e-14:
                    block_val = np.nansum(numerators[i0:i1]) / sden
                else:
                    block_val = np.nan
            else:
                block_val = np.nanmean(numerators[i0:i1])
            
            out[i0:i1] = block_val
        
        return out
    
    elif mode == "pair_cancel":
        for i0, i1 in groups:
            m = i1 - i0
            if m <= 1:
                continue
            
            idxs = list(range(i0, i1))
            used = set()
            
            for i in idxs:
                if i in used:
                    continue
                
                qi = out[i]
                if not np.isfinite(qi) or abs(qi) < min_abs_value:
                    continue
                
                # Find best cancellation partner
                best_j = None
                best_score = None
                
                for j in idxs:
                    if j <= i or j in used:
                        continue
                    
                    qj = out[j]
                    if not np.isfinite(qj) or abs(qj) < min_abs_value:
                        continue
                    
                    # Check if qi ≈ -qj (near-opposite)
                    if (abs(qi + qj) < sum_tol and 
                        abs(abs(qi) - abs(qj)) < mag_tol):
                        score = abs(qi + qj)
                        if best_score is None or score < best_score:
                            best_score = score
                            best_j = j
                
                if best_j is not None:
                    out[i] = 0.0
                    out[best_j] = 0.0
                    used.add(i)
                    used.add(best_j)
        
        return out
    
    else:
        raise ValueError(
            f"Unknown stabilization mode: {mode}. "
            "Must be 'none', 'block_average', or 'pair_cancel'."
        )
