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
from typing import Tuple, Dict, List, Optional
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
    """Calculate band structure using TB-LMTO Hamiltonian"""
    
    def __init__(self, lattice_reader: LatticeReader, hamiltonian_reader: HamiltonianReader):
        self.lattice = lattice_reader
        self.hamiltonian = hamiltonian_reader
        self.kpath = None
        self.band_energies = None
        
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
        """Calculate band structure for given k-points"""
        num_kpoints = len(kpoints)
        total_orbitals = self.lattice.nbas * self.hamiltonian.orbital_size
        
        band_energies = np.zeros((num_kpoints, total_orbitals))
        
        for i, k in enumerate(kpoints):
            H_k = self.fourier_transform_hamiltonian(k)
            eigenvalues, _ = eigh(H_k)
            band_energies[i, :] = np.real(eigenvalues)  # Take real part
        
        return band_energies
    
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
