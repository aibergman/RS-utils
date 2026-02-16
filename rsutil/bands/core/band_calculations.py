"""
Band calculation utilities.

Provides functions for calculating band energies from Hamiltonians, including
collinear spin-resolved calculations and spin block extraction.

Author: Ramon Cardias
"""

from typing import Callable, Optional, Tuple

import numpy as np
from scipy.linalg import eigh

from ..readers import HamiltonianReader, LatticeReader


def extract_spin_blocks(
    H_full: np.ndarray, n_atoms: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract spin-up and spin-down blocks from full Hamiltonian.

    The full Hamiltonian has structure (for each atom, 18x18 total):
    - Indices 0-8: spin-up orbitals (9 orbitals)
    - Indices 9-17: spin-down orbitals (9 orbitals)

    For collinear magnetism, we extract:
    - H_up: diagonal block [0:9, 0:9] for all atoms
    - H_down: diagonal block [9:18, 9:18] for all atoms

    Off-diagonal blocks represent spin-flip terms (SOC, non-collinear),
    which are ignored in collinear approximation.

    Parameters
    ----------
    H_full : np.ndarray
        Full Hamiltonian matrix (nbas*18 x nbas*18)
    n_atoms : int
        Number of atoms in the basis

    Returns
    -------
    tuple
        (H_up, H_down) spin-resolved Hamiltonians

    Author
    ------
    Ramon Cardias
    """
    n_orb_per_spin = 9
    n_up = n_atoms * n_orb_per_spin
    n_down = n_atoms * n_orb_per_spin

    # Initialize spin-resolved Hamiltonians
    H_up = np.zeros((n_up, n_up), dtype=complex)
    H_down = np.zeros((n_down, n_down), dtype=complex)

    # Extract blocks for each atom
    for iatom in range(n_atoms):
        # Indices in full Hamiltonian
        full_start = iatom * 18
        up_start_full = full_start
        up_end_full = full_start + 9
        down_start_full = full_start + 9
        down_end_full = full_start + 18

        # Indices in spin-resolved Hamiltonians
        spin_start = iatom * 9
        spin_end = spin_start + 9

        # Extract spin-up block
        for jatom in range(n_atoms):
            jfull_start = jatom * 18
            jup_start_full = jfull_start
            jup_end_full = jfull_start + 9
            jdown_start_full = jfull_start + 9
            jdown_end_full = jfull_start + 18

            jspin_start = jatom * 9
            jspin_end = jspin_start + 9

            # Copy spin-up to spin-up block
            H_up[spin_start:spin_end, jspin_start:jspin_end] = H_full[
                up_start_full:up_end_full, jup_start_full:jup_end_full
            ]

            # Copy spin-down to spin-down block
            H_down[spin_start:spin_end, jspin_start:jspin_end] = H_full[
                down_start_full:down_end_full, jdown_start_full:jdown_end_full
            ]

    return H_up, H_down


def calculate_bands(
    kpoints: np.ndarray,
    lattice_reader: LatticeReader,
    hamiltonian_reader: HamiltonianReader,
    hamiltonian_fn: Callable[[np.ndarray], np.ndarray],
    store_eigenvectors: bool = False,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Calculate band structure for given k-points.

    Parameters
    ----------
    kpoints : np.ndarray
        Array of k-points in crystal coordinates (nk x 3)
    lattice_reader : LatticeReader
        Reader with lattice structure
    hamiltonian_reader : HamiltonianReader
        Reader with Hamiltonian data
    hamiltonian_fn : callable
        Function that returns H(k) for a given k-point
    store_eigenvectors : bool, optional
        If True, return eigenvectors for downstream analysis

    Returns
    -------
    tuple
        (band_energies, eigenvectors or None)

    Author
    ------
    Ramon Cardias
    """
    num_kpoints = len(kpoints)
    total_orbitals = lattice_reader.nbas * hamiltonian_reader.orbital_size

    band_energies = np.zeros((num_kpoints, total_orbitals))
    eigenvectors = None

    if store_eigenvectors:
        eigenvectors = np.zeros(
            (num_kpoints, total_orbitals, total_orbitals), dtype=complex
        )

    for i, k in enumerate(kpoints):
        H_k = hamiltonian_fn(k)
        eigenvalues, eigvecs = eigh(H_k)
        band_energies[i, :] = np.real(eigenvalues)

        if store_eigenvectors:
            eigenvectors[i, :, :] = eigvecs

    return band_energies, eigenvectors


def calculate_bands_collinear(
    kpoints: np.ndarray,
    lattice_reader: LatticeReader,
    hamiltonian_fn: Callable[[np.ndarray], np.ndarray],
    store_eigenvectors: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """Calculate spin-polarized band structure for collinear magnetism.

    Extracts spin-up and spin-down Hamiltonians from the full 18x18 matrix
    and diagonalizes them separately.

    Parameters
    ----------
    kpoints : np.ndarray
        Array of k-points in crystal coordinates (nk x 3)
    lattice_reader : LatticeReader
        Reader with lattice structure
    hamiltonian_fn : callable
        Function that returns H(k) for a given k-point
    store_eigenvectors : bool, optional
        If True, return eigenvectors for downstream analysis

    Returns
    -------
    tuple
        (band_energies_up, band_energies_down, eigenvectors_up, eigenvectors_down)

    Author
    ------
    Ramon Cardias
    """
    num_kpoints = len(kpoints)
    n_orb_per_spin = lattice_reader.nbas * 9

    band_energies_up = np.zeros((num_kpoints, n_orb_per_spin))
    band_energies_down = np.zeros((num_kpoints, n_orb_per_spin))

    eigenvectors_up = None
    eigenvectors_down = None

    if store_eigenvectors:
        eigenvectors_up = np.zeros(
            (num_kpoints, n_orb_per_spin, n_orb_per_spin), dtype=complex
        )
        eigenvectors_down = np.zeros(
            (num_kpoints, n_orb_per_spin, n_orb_per_spin), dtype=complex
        )

    for i, k in enumerate(kpoints):
        H_k = hamiltonian_fn(k)
        H_up, H_down = extract_spin_blocks(H_k, lattice_reader.nbas)

        eigs_up, eigvecs_up = eigh(H_up)
        eigs_down, eigvecs_down = eigh(H_down)

        band_energies_up[i, :] = np.real(eigs_up)
        band_energies_down[i, :] = np.real(eigs_down)

        if store_eigenvectors:
            eigenvectors_up[i, :, :] = eigvecs_up
            eigenvectors_down[i, :, :] = eigvecs_down

    return band_energies_up, band_energies_down, eigenvectors_up, eigenvectors_down
