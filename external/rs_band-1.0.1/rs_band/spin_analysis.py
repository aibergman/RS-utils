#!/usr/bin/env python3
"""
Non-collinear spin analysis for TB-LMTO Hamiltonians

This script provides tools to analyze the spin structure of non-collinear
TB-LMTO Hamiltonians with the structure:
[H0+Hz   Hx+iHy]
[Hx-iHy  H0-Hz ]

Where H0 is the charge part and Hx, Hy, Hz are the spin components.
"""

import numpy as np
import matplotlib.pyplot as plt
from .projected_bands import ProjectedBandCalculator
from .band_structure_calculator import LatticeReader, HamiltonianReader
from typing import Dict, Tuple
import argparse


class SpinAnalyzer:
    """Analyze spin components of non-collinear Hamiltonian"""
    
    def __init__(self, calculator: ProjectedBandCalculator):
        self.calculator = calculator
        
    def extract_spin_components(self, k_point: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Extract H0, Hx, Hy, Hz components from the k-space Hamiltonian
        
        Args:
            k_point: k-vector for Fourier transform
            
        Returns:
            Dictionary with 'H0', 'Hx', 'Hy', 'Hz' components
        """
        
        # Get full Hamiltonian at k-point
        H_k = self.calculator.fourier_transform_hamiltonian(k_point)
        
        # Extract blocks for each site
        num_sites = self.calculator.lattice.nbas
        spin_components = {'H0': [], 'Hx': [], 'Hy': [], 'Hz': []}
        
        for site in range(num_sites):
            site_start = site * 18
            site_end = site_start + 18
            
            # Extract site Hamiltonian block
            H_site = H_k[site_start:site_end, site_start:site_end]
            
            # Split into 9x9 blocks
            # Upper left: H0 + Hz
            # Upper right: Hx + i*Hy  
            # Lower left: Hx - i*Hy
            # Lower right: H0 - Hz
            
            H_uu = H_site[:9, :9]    # up-up block
            H_ud = H_site[:9, 9:]    # up-down block
            H_du = H_site[9:, :9]    # down-up block
            H_dd = H_site[9:, 9:]    # down-down block
            
            # Extract components
            H0 = 0.5 * (H_uu + H_dd)
            Hz = 0.5 * (H_uu - H_dd)
            Hx = 0.5 * (H_ud + H_du)
            Hy = 0.5j * (H_du - H_ud)
            
            spin_components['H0'].append(H0)
            spin_components['Hx'].append(Hx)
            spin_components['Hy'].append(Hy)
            spin_components['Hz'].append(Hz)
        
        # Convert to arrays
        for key in spin_components:
            spin_components[key] = np.array(spin_components[key])
        
        return spin_components
    
    def calculate_spin_moments(self, band_energies: np.ndarray, 
                              projections: Dict, kpoints: np.ndarray) -> Dict:
        """
        Calculate spin moments for each band and k-point
        
        Returns:
            Dictionary with spin moment information
        """
        
        num_kpoints, num_bands = band_energies.shape
        
        spin_moments = {
            'total': np.zeros((num_kpoints, num_bands)),
            'x': np.zeros((num_kpoints, num_bands)),
            'y': np.zeros((num_kpoints, num_bands)), 
            'z': np.zeros((num_kpoints, num_bands))
        }
        
        # Calculate spin polarization for each k-point and band
        for k_idx in range(num_kpoints):
            spin_up = projections['spin']['up'][k_idx, :]
            spin_down = projections['spin']['down'][k_idx, :]
            
            # Simple estimate: z-component from up/down difference
            spin_moments['z'][k_idx, :] = 0.5 * (spin_up - spin_down)
            spin_moments['total'][k_idx, :] = np.sqrt(spin_moments['z'][k_idx, :]**2)
            
            # For full x,y components, we'd need the eigenvectors
            # This is a simplified calculation
            
        return spin_moments
    
    def plot_spin_texture(self, kdistances: np.ndarray, band_energies: np.ndarray,
                         spin_moments: Dict, klabels: list, band_indices: list = None,
                         output_file: str = None):
        """Plot spin texture for selected bands"""
        
        if band_indices is None:
            # Find bands near Fermi level (assuming E_F = 0)
            fermi_bands = []
            for band in range(band_energies.shape[1]):
                if np.any(np.abs(band_energies[:, band]) < 0.5):  # Within 0.5 eV of Fermi
                    fermi_bands.append(band)
            band_indices = fermi_bands[:5]  # Take first 5 bands near Fermi
        
        fig, axes = plt.subplots(len(band_indices), 1, figsize=(12, 3*len(band_indices)))
        if len(band_indices) == 1:
            axes = [axes]
        
        for idx, band in enumerate(band_indices):
            ax = axes[idx]
            
            # Plot band structure
            ax.plot(kdistances, band_energies[:, band], 'k-', linewidth=2, label='Band')
            
            # Color by spin polarization
            spin_z = spin_moments['z'][:, band]
            scatter = ax.scatter(kdistances, band_energies[:, band], 
                               c=spin_z, cmap='RdBu', s=50, alpha=0.8)
            
            # Set up k-point labels
            tick_positions = []
            tick_labels = []
            for i, label in enumerate(klabels):
                if label:
                    tick_positions.append(kdistances[i])
                    tick_labels.append(label)
            
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)
            
            for pos in tick_positions:
                ax.axvline(x=pos, color='gray', linestyle='--', alpha=0.3)
            
            ax.set_ylabel('Energy (eV)')
            ax.set_title(f'Band {band+1} - Spin Texture (z-component)')
            ax.grid(True, alpha=0.3)
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Spin z-component')
        
        axes[-1].set_xlabel('k-path')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Spin texture plot saved to {output_file}")
        
        plt.show()
    
    def analyze_magnetic_anisotropy(self, kpoints: np.ndarray, 
                                   sample_points: int = 10) -> Dict:
        """
        Analyze magnetic anisotropy by comparing different spin directions
        
        This is a simplified analysis - full MAE calculation would require
        self-consistent field calculations.
        """
        
        # Sample a few k-points for analysis
        k_indices = np.linspace(0, len(kpoints)-1, sample_points, dtype=int)
        
        anisotropy_data = {
            'k_points': [],
            'H0_trace': [],
            'Hz_trace': [],
            'Hx_trace': [],
            'Hy_trace': []
        }
        
        for k_idx in k_indices:
            k_point = kpoints[k_idx]
            spin_components = self.extract_spin_components(k_point)
            
            anisotropy_data['k_points'].append(k_point)
            anisotropy_data['H0_trace'].append(np.trace(spin_components['H0'][0]).real)
            anisotropy_data['Hz_trace'].append(np.trace(spin_components['Hz'][0]).real)
            anisotropy_data['Hx_trace'].append(np.trace(spin_components['Hx'][0]).real)
            anisotropy_data['Hy_trace'].append(np.trace(spin_components['Hy'][0]).real)
        
        return anisotropy_data
    
    def plot_orbital_magnetism(self, kdistances: np.ndarray, band_energies: np.ndarray,
                              projections: Dict, klabels: list, output_file: str = None):
        """Plot orbital-resolved magnetic moments"""
        
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))
        
        orbital_groups = ['s', 'p', 'd']
        colors = ['red', 'blue', 'green']
        
        for orbital_idx, (orbital, color) in enumerate(zip(orbital_groups, colors)):
            ax = axes[orbital_idx]
            
            # Calculate orbital-resolved spin polarization
            spin_up_proj = projections['orbital'][orbital]
            # For spin down, we need to sum over the second spin channel
            # This is simplified - in practice you'd separate the projections
            spin_polarization = spin_up_proj  # Simplified
            
            # Plot bands colored by orbital magnetism
            for band in range(band_energies.shape[1]):
                weights = spin_polarization[:, band]
                if np.max(weights) > 0.05:  # Only plot if significant weight
                    sizes = weights * 50
                    ax.scatter(kdistances, band_energies[:, band], 
                             s=sizes, c=color, alpha=0.6, label=f'{orbital}' if band == 0 else "")
            
            # Set up k-point labels
            tick_positions = []
            tick_labels = []
            for i, label in enumerate(klabels):
                if label:
                    tick_positions.append(kdistances[i])
                    tick_labels.append(label)
            
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)
            
            for pos in tick_positions:
                ax.axvline(x=pos, color='gray', linestyle='--', alpha=0.3)
            
            ax.set_ylabel('Energy (eV)')
            ax.set_title(f'{orbital.upper()}-orbital Magnetic Contributions')
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        axes[-1].set_xlabel('k-path')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Orbital magnetism plot saved to {output_file}")
        
        plt.show()


def main():
    """Main function for spin analysis"""
    
    parser = argparse.ArgumentParser(description='Analyze non-collinear TB-LMTO spin structure')
    parser.add_argument('--lattice', default='lattice.nml', 
                       help='Lattice file (default: lattice.nml)')
    parser.add_argument('--hamiltonian', default='rs2txt.dat',
                       help='Hamiltonian file (default: rs2txt.dat)')
    parser.add_argument('--output-prefix', default='spin_analysis',
                       help='Output file prefix (default: spin_analysis)')
    parser.add_argument('--nk', type=int, default=50,
                       help='Number of k-points (default: 50)')
    
    args = parser.parse_args()
    
    print("Setting up projected band calculator...")
    lattice_reader = LatticeReader(args.lattice)
    hamiltonian_reader = HamiltonianReader(args.hamiltonian)
    calculator = ProjectedBandCalculator(lattice_reader, hamiltonian_reader)
    
    print("Generating k-path...")
    kpoints, klabels, kdistances = calculator.setup_kpath(args.nk)
    
    print("Calculating projected bands...")
    band_energies, projections = calculator.calculate_projected_bands(kpoints)
    
    print("Analyzing spin structure...")
    spin_analyzer = SpinAnalyzer(calculator)
    
    # Calculate spin moments
    spin_moments = spin_analyzer.calculate_spin_moments(band_energies, projections, kpoints)
    
    # Plot spin texture
    spin_analyzer.plot_spin_texture(kdistances, band_energies, spin_moments, klabels,
                                   output_file=f"{args.output_prefix}_texture.png")
    
    # Plot orbital magnetism
    spin_analyzer.plot_orbital_magnetism(kdistances, band_energies, projections, klabels,
                                        output_file=f"{args.output_prefix}_orbital_mag.png")
    
    # Analyze magnetic anisotropy
    print("Analyzing magnetic anisotropy...")
    anisotropy = spin_analyzer.analyze_magnetic_anisotropy(kpoints)
    
    print("Magnetic anisotropy analysis:")
    print(f"Average H0 trace: {np.mean(anisotropy['H0_trace']):.6f}")
    print(f"Average Hz trace: {np.mean(anisotropy['Hz_trace']):.6f}")
    print(f"Average Hx trace: {np.mean(anisotropy['Hx_trace']):.6f}")
    print(f"Average Hy trace: {np.mean(anisotropy['Hy_trace']):.6f}")
    
    print("Done!")


if __name__ == "__main__":
    main()
