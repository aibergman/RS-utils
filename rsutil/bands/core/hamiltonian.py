"""
Hamiltonian utilities for k-space calculations.

Provides functions for Fourier transforming real-space Hamiltonians to k-space,
and utilities for dumping and analyzing Hamiltonian matrices.

Author: Ramon Cardias (refactoring), Original rslmto_official
"""

import numpy as np
from typing import List, Tuple, Dict

from ..readers import HamiltonianReader, LatticeReader


def fourier_transform_hamiltonian(
    k_point: np.ndarray,
    lattice_reader: LatticeReader,
    hamiltonian_reader: HamiltonianReader,
) -> np.ndarray:
    """Fourier transform the real-space Hamiltonian to k-space.
    
    Computes H(k) = Σ_n H_n * exp(i k·R_n) where H_n are the real-space
    Hamiltonian matrices and R_n are displacements from atomic neighbors.
    
    Parameters
    ----------
    k_point : np.ndarray
        k-point in fractional coordinates (3,)
    lattice_reader : LatticeReader
        Reader with lattice structure information
    hamiltonian_reader : HamiltonianReader
        Reader with Hamiltonian data
        
    Returns
    -------
    np.ndarray
        k-space Hamiltonian matrix (complex)
        
    Author
    ------
    Ramon Cardias (refactored), Original rslmto_official
    """
    # Get reciprocal lattice vectors
    a = lattice_reader.lattice_vectors
    b = 2 * np.pi * np.linalg.inv(a.T)  # Reciprocal lattice vectors as columns

    # Convert k-point to Cartesian coordinates
    k_cart = b @ k_point

    # Initialize k-space Hamiltonian
    total_orbitals = lattice_reader.nbas * hamiltonian_reader.orbital_size
    H_k = np.zeros((total_orbitals, total_orbitals), dtype=complex)

    # Sum over all atoms and their neighbors
    orbital_offset = 0
    for ntype in hamiltonian_reader.hamiltonian_data:
        for neighbor, neighbor_data in hamiltonian_reader.hamiltonian_data[
            ntype
        ].items():
            displacement = neighbor_data["displacement"]
            matrix = neighbor_data["matrix"]

            # Calculate phase factor
            phase = np.exp(1j * np.dot(k_cart, displacement))

            # Add to k-space Hamiltonian
            start_idx = orbital_offset
            end_idx = start_idx + hamiltonian_reader.orbital_size

            H_k[start_idx:end_idx, start_idx:end_idx] += phase * matrix

        orbital_offset += hamiltonian_reader.orbital_size

    return H_k


class HamiltonianDumper:
    """Utilities for exporting and analyzing Hamiltonian matrices.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader with lattice structure
    hamiltonian_reader : HamiltonianReader
        Reader with Hamiltonian data
        
    Author
    ------
    Ramon Cardias
    """
    
    def __init__(
        self,
        lattice_reader: LatticeReader,
        hamiltonian_reader: HamiltonianReader,
    ):
        self.lattice = lattice_reader
        self.hamiltonian = hamiltonian_reader
    
    def dump_at_points(self, output_prefix: str = "H_k") -> Tuple[np.ndarray, np.ndarray]:
        """Dump Hamiltonian matrices at Gamma and H points for BCC structure.
        
        Parameters
        ----------
        output_prefix : str, default="H_k"
            Prefix for output files
            
        Returns
        -------
        tuple
            (H_gamma, H_h) - Hamiltonian matrices at Gamma and H points
            
        Author
        ------
        Ramon Cardias
        """
        # Define high-symmetry points for BCC
        gamma_point = np.array([0.0, 0.0, 0.0])  # Gamma point
        h_point = np.array([0.5, -0.5, 0.5])  # H point for BCC

        print("High-symmetry points for BCC structure:")
        print(
            f"  Gamma: ({gamma_point[0]:.1f}, {gamma_point[1]:.1f}, {gamma_point[2]:.1f}) in units of 2π/a"
        )
        print(
            f"  H: ({h_point[0]:.1f}, {h_point[1]:.1f}, {h_point[2]:.1f}) in units of 2π/a"
        )
        print()

        H_gamma = fourier_transform_hamiltonian(gamma_point, self.lattice, self.hamiltonian)
        H_h = fourier_transform_hamiltonian(h_point, self.lattice, self.hamiltonian)

        # Save to text files
        self._write_matrix_text(H_gamma, gamma_point, f"{output_prefix}_gamma.dat")
        self._write_matrix_text(H_h, h_point, f"{output_prefix}_h.dat")

        # Also save numpy files for programmatic use
        np.save(f"{output_prefix}_gamma.npy", H_gamma)
        np.save(f"{output_prefix}_h.npy", H_h)

        print(
            f"Hamiltonian at Gamma point saved to {output_prefix}_gamma.dat and {output_prefix}_gamma.npy"
        )
        print(
            f"Hamiltonian at H point saved to {output_prefix}_h.dat and {output_prefix}_h.npy"
        )

        # Print some basic info
        print(f"Matrix size: {H_gamma.shape}")
        print(
            f"Gamma point matrix is Hermitian: {np.allclose(H_gamma, H_gamma.conj().T)}"
        )
        print(f"H point matrix is Hermitian: {np.allclose(H_h, H_h.conj().T)}")

        return H_gamma, H_h

    def dump_at_custom_points(
        self, k_points: List[Tuple[str, np.ndarray]], output_prefix: str = "H_k"
    ) -> Dict[str, np.ndarray]:
        """Dump Hamiltonian matrices at custom k-points.
        
        Parameters
        ----------
        k_points : list of tuple
            List of (label, k_point) pairs
        output_prefix : str, default="H_k"
            Prefix for output files
            
        Returns
        -------
        dict
            Dictionary mapping point names to Hamiltonian matrices
            
        Author
        ------
        Ramon Cardias
        """
        results = {}
        for name, k_point in k_points:
            H_k = fourier_transform_hamiltonian(k_point, self.lattice, self.hamiltonian)
            results[name] = H_k

            # Save to text file
            self._write_matrix_text(H_k, k_point, f"{output_prefix}_{name.lower()}.dat")

            # Also save numpy file for programmatic use
            np.save(f"{output_prefix}_{name.lower()}.npy", H_k)

            print(
                f"Hamiltonian at {name} point saved to {output_prefix}_{name.lower()}.dat and {output_prefix}_{name.lower()}.npy"
            )

        return results

    def _write_matrix_text(
        self, H_matrix: np.ndarray, k_point: np.ndarray, filename: str
    ) -> None:
        """Write Hamiltonian matrix to text file.
        
        Parameters
        ----------
        H_matrix : np.ndarray
            Hamiltonian matrix to write
        k_point : np.ndarray
            k-point coordinates
        filename : str
            Output filename
            
        Author
        ------
        Ramon Cardias
        """
        n_rows, n_cols = H_matrix.shape

        with open(filename, "w") as f:
            f.write("# Complex matrix dump for debugging\n")
            f.write(
                f"# k-point: {k_point[0]:12.8f} {k_point[1]:12.8f} {k_point[2]:12.8f}\n"
            )
            f.write(f"# Matrix dimensions: {n_rows} x {n_cols}\n")
            f.write("# Format: row, col, real_part, imag_part\n")
            f.write("#\n")

            for i in range(n_rows):
                for j in range(n_cols):
                    real_part = H_matrix[i, j].real
                    imag_part = H_matrix[i, j].imag
                    f.write(f"{i+1:5d}{j+1:5d}{real_part:23.12E}{imag_part:23.12E}\n")
