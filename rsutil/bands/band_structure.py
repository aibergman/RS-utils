#!/usr/bin/env python3
"""
TB-LMTO Band Structure Calculator

This program reads TB-LMTO Hamiltonian data and lattice information to:
1. Parse the real-space Hamiltonian from rs2txt.dat
2. Read crystal structure from lattice.nml
3. Use spglib to determine crystal symmetry and generate k-path
4. Fourier transform the Hamiltonian to k-space
5. Diagonalize and calculate the band structure
6. Plot the results

Author: Created for rslmto_official
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from pathlib import Path
import argparse
from typing import Tuple, Dict, List
import warnings

from .readers import LatticeReader, HamiltonianReader

try:
    import spglib
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False
    warnings.warn("spglib not found. Install with: pip install spglib")

try:
    from pymatgen.core import Structure
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False
    warnings.warn("pymatgen not found. Install with: pip install pymatgen")


class BandStructureCalculator:
    """Calculate band structure using TB-LMTO Hamiltonian
    
    Supports both non-collinear (default) and collinear spin-polarized calculations.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader for lattice structure
    hamiltonian_reader : HamiltonianReader
        Reader for Hamiltonian data
    collinear : bool, optional
        If True, separate spin-up and spin-down bands for collinear magnetism.
        If False (default), use full non-collinear Hamiltonian.
        
    Attributes
    ----------
    collinear : bool
        Whether to use collinear spin treatment
    """
    
    def __init__(self, lattice_reader: LatticeReader, hamiltonian_reader: HamiltonianReader,
                 collinear: bool = False):
        self.lattice = lattice_reader
        self.hamiltonian = hamiltonian_reader
        self.collinear = collinear
        self.kpath = None
        self.band_energies = None
        self.band_energies_up = None
        self.band_energies_down = None
        
    def setup_kpath(self, num_points: int = 100) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-point path using structure analysis"""
        
        # Import here to avoid circular imports
        from .structure_analysis import StructureAnalyzer
        
        # Analyze structure and generate proper k-path
        analyzer = StructureAnalyzer(self.lattice)
        
        if analyzer.special_points and analyzer.band_path:
            return analyzer.generate_kpath_points(num_points)
        else:
            # Fallback to manual method
            return self._setup_kpath_manual(num_points)
    
    def _setup_kpath_spglib(self, num_points: int) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-path using spglib"""
        cell = (
            self.lattice.lattice_vectors.T,
            self.lattice.atomic_positions,
            self.lattice.atomic_numbers
        )
        
        # Get the symmetry information
        symmetry = spglib.get_symmetry_dataset(cell)
        
        # Get the standardized cell
        mapping, grid = spglib.get_ir_reciprocal_mesh([20, 20, 20], cell)
        
        # Define high symmetry points (this is crystal system dependent)
        # For now, use a simple path
        lattice_type = symmetry['international']
        
        if 'cubic' in lattice_type.lower() or 'P' in lattice_type:
            # Simple cubic path
            special_points = {
                'G': [0, 0, 0],      # Gamma
                'X': [0.5, 0, 0],    # X
                'M': [0.5, 0.5, 0],  # M
                'R': [0.5, 0.5, 0.5] # R
            }
            path = ['G', 'X', 'M', 'G', 'R', 'X']
        else:
            # Generic path
            special_points = {
                'G': [0, 0, 0],
                'X': [0.5, 0, 0],
                'Y': [0, 0.5, 0],
                'Z': [0, 0, 0.5]
            }
            path = ['G', 'X', 'Y', 'G', 'Z']
        
        return self._generate_kpath(special_points, path, num_points)
    
    def _setup_kpath_pymatgen(self, num_points: int) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-path using pymatgen"""
        structure = Structure(
            lattice=self.lattice.lattice_vectors.T,
            species=self.lattice.atomic_numbers,
            coords=self.lattice.atomic_positions
        )
        
        kpath = HighSymmKpath(structure)
        special_points = kpath.kpath['kpoints']
        path = kpath.kpath['path'][0]  # Take first path
        
        return self._generate_kpath(special_points, path, num_points)
    
    def _setup_kpath_manual(self, num_points: int) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Fallback manual k-path setup"""
        warnings.warn("Using manual k-path. Install spglib or pymatgen for automatic detection.")
        
        special_points = {
            'G': [0, 0, 0],
            'X': [0.5, 0, 0],
            'Y': [0, 0.5, 0],
            'Z': [0, 0, 0.5],
            'M': [0.5, 0.5, 0]
        }
        path = ['G', 'X', 'M', 'G', 'Z']
        
        return self._generate_kpath(special_points, path, num_points)
    
    def _generate_kpath(self, special_points: Dict, path: List[str], num_points: int) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Generate k-point path from special points"""
        kpoints = []
        klabels = []
        kdistances = []
        total_distance = 0
        
        for i in range(len(path) - 1):
            start = np.array(special_points[path[i]])
            end = np.array(special_points[path[i + 1]])
            
            # Calculate number of points for this segment
            segment_length = np.linalg.norm(end - start)
            segment_points = max(2, int(num_points * segment_length / 3))  # Rough estimate
            
            # Generate points along segment
            for j in range(segment_points):
                if j == 0 and i == 0:  # First point
                    kpoints.append(start)
                    klabels.append(path[i])
                    kdistances.append(total_distance)
                elif j == 0:  # Skip duplicate points
                    klabels[-1] = f"{klabels[-1]}|{path[i]}"
                else:
                    t = j / (segment_points - 1)
                    k = start + t * (end - start)
                    kpoints.append(k)
                    klabels.append('')
                    total_distance += np.linalg.norm(k - kpoints[-2]) if len(kpoints) > 1 else 0
                    kdistances.append(total_distance)
            
            # Add end point
            if i == len(path) - 2:  # Last segment
                kpoints.append(end)
                klabels.append(path[i + 1])
                total_distance += np.linalg.norm(end - kpoints[-2])
                kdistances.append(total_distance)
        
        return np.array(kpoints), klabels, np.array(kdistances)
    
    def extract_spin_blocks(self, H_full: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Extract spin-up and spin-down blocks from full Hamiltonian
        
        The full Hamiltonian has structure (for each atom, 18×18 total):
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
            Full Hamiltonian matrix (nbas*18 × nbas*18)
            
        Returns
        -------
        H_up : np.ndarray
            Spin-up Hamiltonian (nbas*9 × nbas*9)
        H_down : np.ndarray
            Spin-down Hamiltonian (nbas*9 × nbas*9)
        """
        n_atoms = self.lattice.nbas
        n_orb_per_spin = 9  # 9 orbitals per spin channel
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
                H_up[spin_start:spin_end, jspin_start:jspin_end] = \
                    H_full[up_start_full:up_end_full, jup_start_full:jup_end_full]
                
                # Copy spin-down to spin-down block
                H_down[spin_start:spin_end, jspin_start:jspin_end] = \
                    H_full[down_start_full:down_end_full, jdown_start_full:jdown_end_full]
        
        return H_up, H_down
    
    def fourier_transform_hamiltonian(self, k_point: np.ndarray) -> np.ndarray:
        """Fourier transform the real-space Hamiltonian to k-space"""
        
        # Get reciprocal lattice vectors
        a = self.lattice.lattice_vectors
        b = 2 * np.pi * np.linalg.inv(a.T)  # Reciprocal lattice vectors as columns
        
        # Convert k-point to Cartesian coordinates
        k_cart = b @ k_point
        
        # Initialize k-space Hamiltonian
        total_orbitals = self.lattice.nbas * self.hamiltonian.orbital_size
        H_k = np.zeros((total_orbitals, total_orbitals), dtype=complex)
        
        # Sum over all atoms and their neighbors
        orbital_offset = 0
        for ntype in self.hamiltonian.hamiltonian_data:
            for neighbor, neighbor_data in self.hamiltonian.hamiltonian_data[ntype].items():
                displacement = neighbor_data['displacement']
                matrix = neighbor_data['matrix']
                
                # Calculate phase factor
                phase = np.exp(1j * np.dot(k_cart, displacement))
                
                # Add to k-space Hamiltonian
                start_idx = orbital_offset
                end_idx = start_idx + self.hamiltonian.orbital_size
                
                H_k[start_idx:end_idx, start_idx:end_idx] += phase * matrix
            
            orbital_offset += self.hamiltonian.orbital_size
        
        return H_k
    
    def dump_hamiltonian_at_points(self, output_prefix: str = "H_k"):
        """Dump Hamiltonian matrices at Gamma and H points for BCC structure (temporary function)"""
        
        # Define high-symmetry points for BCC
        gamma_point = np.array([0.0, 0.0, 0.0])  # Gamma point
        h_point = np.array([0.5, -0.5, 0.5])     # H point for BCC
        
        print("High-symmetry points for BCC structure:")
        print(f"  Gamma: ({gamma_point[0]:.1f}, {gamma_point[1]:.1f}, {gamma_point[2]:.1f}) in units of 2π/a")
        print(f"  H: ({h_point[0]:.1f}, {h_point[1]:.1f}, {h_point[2]:.1f}) in units of 2π/a")
        print()
        
        # print("Computing H(k) at Gamma point...")
        H_gamma = self.fourier_transform_hamiltonian(gamma_point)
        
        # print("Computing H(k) at H point...")
        H_h = self.fourier_transform_hamiltonian(h_point)
        
        # Save to text files
        self._write_hamiltonian_text(H_gamma, gamma_point, f"{output_prefix}_gamma.dat")
        self._write_hamiltonian_text(H_h, h_point, f"{output_prefix}_h.dat")
        
        # Also save numpy files for programmatic use
        np.save(f"{output_prefix}_gamma.npy", H_gamma)
        np.save(f"{output_prefix}_h.npy", H_h)
        
        print(f"Hamiltonian at Gamma point saved to {output_prefix}_gamma.dat and {output_prefix}_gamma.npy")
        print(f"Hamiltonian at H point saved to {output_prefix}_h.dat and {output_prefix}_h.npy")
        
        # Print some basic info
        print(f"Matrix size: {H_gamma.shape}")
        print(f"Gamma point matrix is Hermitian: {np.allclose(H_gamma, H_gamma.conj().T)}")
        print(f"H point matrix is Hermitian: {np.allclose(H_h, H_h.conj().T)}")
        
        return H_gamma, H_h
    
    def dump_hamiltonian_at_custom_points(self, k_points: List[Tuple[str, np.ndarray]], output_prefix: str = "H_k"):
        """Dump Hamiltonian matrices at custom k-points"""
        
        # print("Custom k-points:")
        # for name, k_point in k_points:
        #     print(f"  {name}: ({k_point[0]:.6f}, {k_point[1]:.6f}, {k_point[2]:.6f}) in units of 2π/a")
        # print()
        
        results = {}
        for name, k_point in k_points:
            # print(f"Computing H(k) at {name} point...")
            H_k = self.fourier_transform_hamiltonian(k_point)
            results[name] = H_k
            
            # Save to text file
            self._write_hamiltonian_text(H_k, k_point, f"{output_prefix}_{name.lower()}.dat")
            
            # Also save numpy file for programmatic use
            np.save(f"{output_prefix}_{name.lower()}.npy", H_k)
            
            print(f"Hamiltonian at {name} point saved to {output_prefix}_{name.lower()}.dat and {output_prefix}_{name.lower()}.npy")
        
        return results
    
    def _write_hamiltonian_text(self, H_matrix: np.ndarray, k_point: np.ndarray, filename: str):
        """Write Hamiltonian matrix to text file in the specified format"""
        n_rows, n_cols = H_matrix.shape
        
        with open(filename, 'w') as f:
            f.write("# Complex matrix dump for debugging\n")
            f.write(f"# k-point: {k_point[0]:12.8f} {k_point[1]:12.8f} {k_point[2]:12.8f}\n")
            f.write(f"# Matrix dimensions: {n_rows} x {n_cols}\n")
            f.write("# Format: row, col, real_part, imag_part\n")
            f.write("#\n")
            
            for i in range(n_rows):
                for j in range(n_cols):
                    real_part = H_matrix[i, j].real
                    imag_part = H_matrix[i, j].imag
                    f.write(f"{i+1:5d}{j+1:5d}{real_part:23.12E}{imag_part:23.12E}\n")
        
        # print(f"Text dump saved to {filename}")
    
    def calculate_bands(self, kpoints: np.ndarray) -> np.ndarray:
        """Calculate band structure for given k-points
        
        For collinear=False (default): Returns all bands from full 18×18 Hamiltonian
        For collinear=True: Use calculate_bands_collinear() instead
        
        Parameters
        ----------
        kpoints : np.ndarray
            Array of k-points in crystal coordinates (nk × 3)
            
        Returns
        -------
        band_energies : np.ndarray
            Band energies at each k-point (nk × nbands)
        """
        if self.collinear:
            raise ValueError("For collinear calculations, use calculate_bands_collinear() "
                           "which returns separate up and down bands")
        
        num_kpoints = len(kpoints)
        total_orbitals = self.lattice.nbas * self.hamiltonian.orbital_size
        
        band_energies = np.zeros((num_kpoints, total_orbitals))
        
        for i, k in enumerate(kpoints):
            H_k = self.fourier_transform_hamiltonian(k)
            eigenvalues, _ = eigh(H_k)
            band_energies[i, :] = np.real(eigenvalues)  # Take real part
        
        return band_energies
    
    def calculate_bands_collinear(self, kpoints: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate spin-polarized band structure for collinear magnetism
        
        Extracts spin-up and spin-down Hamiltonians from the full 18×18 matrix
        and diagonalizes them separately. This gives two sets of band energies
        that can differ due to exchange splitting in magnetic materials.
        
        Parameters
        ----------
        kpoints : np.ndarray
            Array of k-points in crystal coordinates (nk × 3)
            
        Returns
        -------
        band_energies_up : np.ndarray
            Spin-up band energies (nk × 9*nbas)
        band_energies_down : np.ndarray
            Spin-down band energies (nk × 9*nbas)
            
        Notes
        -----
        - Only valid for collinear magnetic calculations
        - Ignores spin-flip terms (SOC, non-collinear magnetism)
        - Each spin channel has 9*nbas bands (9 orbitals per atom)
        """
        if not self.collinear:
            raise ValueError("calculate_bands_collinear() requires collinear=True. "
                           "Use calculate_bands() for non-collinear calculations.")
        
        num_kpoints = len(kpoints)
        n_orb_per_spin = self.lattice.nbas * 9  # 9 orbitals per atom per spin
        
        band_energies_up = np.zeros((num_kpoints, n_orb_per_spin))
        band_energies_down = np.zeros((num_kpoints, n_orb_per_spin))
        
        for i, k in enumerate(kpoints):
            # Get full Hamiltonian at this k-point
            H_k = self.fourier_transform_hamiltonian(k)
            
            # Extract spin-resolved Hamiltonians
            H_up, H_down = self.extract_spin_blocks(H_k)
            
            # Diagonalize each spin channel
            eigs_up, _ = eigh(H_up)
            eigs_down, _ = eigh(H_down)
            
            band_energies_up[i, :] = np.real(eigs_up)
            band_energies_down[i, :] = np.real(eigs_down)
        
        # Store for later access
        self.band_energies_up = band_energies_up
        self.band_energies_down = band_energies_down
        
        return band_energies_up, band_energies_down
    
    def plot_bands(self, kdistances: np.ndarray, band_energies: np.ndarray, 
                   klabels: List[str], output_file: str = None):
        """Plot the band structure"""
        
        plt.figure(figsize=(10, 8))
        
        # Plot all bands
        for band in range(band_energies.shape[1]):
            plt.plot(kdistances, band_energies[:, band], 'b-', linewidth=1)
        
        # Set up k-point labels
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:  # Non-empty label
                tick_positions.append(kdistances[i])
                tick_labels.append(label)
        
        plt.xticks(tick_positions, tick_labels)
        
        # Add vertical lines at high-symmetry points
        for pos in tick_positions:
            plt.axvline(x=pos, color='k', linestyle='--', alpha=0.3)
        
        plt.xlabel('k-path')
        plt.ylabel('Energy (eV)')
        plt.title('TB-LMTO Band Structure')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Band structure saved to {output_file}")
        
        plt.show()
    
    def plot_bands_collinear(self, kdistances: np.ndarray, 
                            band_energies_up: np.ndarray, 
                            band_energies_down: np.ndarray,
                            klabels: List[str], 
                            output_file: str = None):
        """Plot spin-polarized band structure
        
        Parameters
        ----------
        kdistances : np.ndarray
            Distances along k-path for x-axis
        band_energies_up : np.ndarray
            Spin-up band energies (nk × nbands)
        band_energies_down : np.ndarray
            Spin-down band energies (nk × nbands)
        klabels : List[str]
            Labels for high-symmetry k-points
        output_file : str, optional
            Output filename for saving the plot
            
        Notes
        -----
        Spin-up bands shown in red, spin-down in blue.
        Exchange splitting visible as energy difference between channels.
        """
        plt.figure(figsize=(10, 8))
        
        # Plot spin-up bands in red
        for band in range(band_energies_up.shape[1]):
            plt.plot(kdistances, band_energies_up[:, band], 'r-', 
                    linewidth=1, alpha=0.7, label='↑' if band == 0 else '')
        
        # Plot spin-down bands in blue
        for band in range(band_energies_down.shape[1]):
            plt.plot(kdistances, band_energies_down[:, band], 'b-', 
                    linewidth=1, alpha=0.7, label='↓' if band == 0 else '')
        
        # Set up k-point labels
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:  # Non-empty label
                tick_positions.append(kdistances[i])
                tick_labels.append(label)
        
        plt.xticks(tick_positions, tick_labels)
        
        # Add vertical lines at high-symmetry points
        for pos in tick_positions:
            plt.axvline(x=pos, color='k', linestyle='--', alpha=0.3)
        
        plt.xlabel('k-path')
        plt.ylabel('Energy (eV)')
        plt.title('TB-LMTO Spin-Polarized Band Structure')
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Spin-polarized band structure saved to {output_file}")
        
        plt.show()


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Calculate TB-LMTO band structure')
    parser.add_argument('--lattice', default='lattice.nml', 
                       help='Lattice file (default: lattice.nml)')
    parser.add_argument('--hamiltonian', default='rs2txt.dat',
                       help='Hamiltonian file (default: rs2txt.dat)')
    parser.add_argument('--output', default='band_structure.png',
                       help='Output plot file (default: band_structure.png)')
    parser.add_argument('--nk', type=int, default=100,
                       help='Number of k-points (default: 100)')
    parser.add_argument('--dump-hamiltonian', action='store_true',
                       help='Dump H(k) matrices at Gamma and H points')
    
    args = parser.parse_args()
    
    # Check if files exist
    if not Path(args.lattice).exists():
        print(f"Error: Lattice file {args.lattice} not found")
        return
    
    if not Path(args.hamiltonian).exists():
        print(f"Error: Hamiltonian file {args.hamiltonian} not found")
        return
    
    print("Reading lattice information...")
    lattice_reader = LatticeReader(args.lattice)
    
    print("Reading Hamiltonian data...")
    hamiltonian_reader = HamiltonianReader(args.hamiltonian)
    
    print("Setting up band structure calculator...")
    calculator = BandStructureCalculator(lattice_reader, hamiltonian_reader)
    
    print(f"Generating k-path with {args.nk} points...")
    kpoints, klabels, kdistances = calculator.setup_kpath(args.nk)
    
    print("Calculating band structure...")
    band_energies = calculator.calculate_bands(kpoints)
    
    if args.dump_hamiltonian:
        print("Dumping Hamiltonian matrices...")
        calculator.dump_hamiltonian_at_points()
    
    print("Plotting results...")
    calculator.plot_bands(kdistances, band_energies, klabels, args.output)
    
    print("Done!")


if __name__ == "__main__":
    main()
