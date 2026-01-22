#!/usr/bin/env python3
"""
Enhanced TB-LMTO Band Structure Calculator with Orbital and Site Projections

This module extends the basic band structure calculator to include:
- Orbital projected band structure (s, p, d contributions)
- Site projected band structure (atom-specific contributions)
- Spin projections for non-collinear systems

The Hamiltonian structure is 18×18 per site:
- 9 orbitals per spin channel (s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2)
- Non-collinear spin structure: [H0+Hz, Hx+iHy; Hx-iHy, H0-Hz]
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from .band_structure_calculator import BandStructureCalculator, LatticeReader, HamiltonianReader
from typing import Dict, List, Tuple
import matplotlib.colors as mcolors


class ProjectedBandCalculator(BandStructureCalculator):
    """Extended band calculator with orbital and site projections"""
    
    def __init__(self, lattice_reader: LatticeReader, hamiltonian_reader: HamiltonianReader):
        super().__init__(lattice_reader, hamiltonian_reader)
        
        # Define orbital names and indices
        self.orbital_names = ['s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dxz', 'dx2-y2', 'dz2']
        self.orbital_indices = {name: i for i, name in enumerate(self.orbital_names)}
        
        # Orbital groups for plotting
        self.orbital_groups = {
            's': [0],           # s orbital
            'p': [1, 2, 3],     # px, py, pz
            'd': [4, 5, 6, 7, 8] # d orbitals
        }
        
        # Spin channels
        self.spin_channels = ['up', 'down']
        
    def calculate_projected_bands(self, kpoints: np.ndarray) -> Tuple[np.ndarray, Dict]:
        """
        Calculate band structure with orbital and site projections
        
        Returns:
            band_energies: Energy eigenvalues
            projections: Dictionary containing all projection data
        """
        
        num_kpoints = len(kpoints)
        total_orbitals = self.lattice.nbas * self.hamiltonian.orbital_size
        
        # Storage for results
        band_energies = np.zeros((num_kpoints, total_orbitals))
        eigenvectors = np.zeros((num_kpoints, total_orbitals, total_orbitals), dtype=complex)
        
        # Calculate eigenvalues and eigenvectors
        for i, k in enumerate(kpoints):
            H_k = self.fourier_transform_hamiltonian(k)
            eigenvals, eigenvecs = eigh(H_k)
            band_energies[i, :] = np.real(eigenvals)
            eigenvectors[i, :, :] = eigenvecs
        
        # Calculate projections
        projections = self._calculate_projections(eigenvectors)
        
        return band_energies, projections
    
    def _calculate_projections(self, eigenvectors: np.ndarray) -> Dict:
        """Calculate orbital and site projections from eigenvectors"""
        
        num_kpoints, total_orbitals, _ = eigenvectors.shape
        num_sites = self.lattice.nbas
        
        projections = {
            'orbital': {},
            'site': {},
            'spin': {},
            'orbital_site': {}
        }
        
        # Initialize projection arrays
        for orbital_group in self.orbital_groups:
            projections['orbital'][orbital_group] = np.zeros((num_kpoints, total_orbitals))
        
        for site in range(num_sites):
            projections['site'][site] = np.zeros((num_kpoints, total_orbitals))
        
        for spin in self.spin_channels:
            projections['spin'][spin] = np.zeros((num_kpoints, total_orbitals))
        
        # Calculate projections for each k-point and band
        for k_idx in range(num_kpoints):
            for band_idx in range(total_orbitals):
                eigenvec = eigenvectors[k_idx, :, band_idx]
                
                # Calculate projections for each site
                for site in range(num_sites):
                    site_start = site * 18
                    site_end = site_start + 18
                    site_coeffs = eigenvec[site_start:site_end]
                    
                    # Site projection (total weight on this site)
                    site_weight = np.sum(np.abs(site_coeffs)**2)
                    projections['site'][site][k_idx, band_idx] = site_weight
                    
                    # Split into spin channels
                    spin_up_coeffs = site_coeffs[:9]    # First 9 orbitals (spin up)
                    spin_down_coeffs = site_coeffs[9:]  # Last 9 orbitals (spin down)
                    
                    spin_up_weight = np.sum(np.abs(spin_up_coeffs)**2)
                    spin_down_weight = np.sum(np.abs(spin_down_coeffs)**2)
                    
                    # Add to total spin projections
                    projections['spin']['up'][k_idx, band_idx] += spin_up_weight
                    projections['spin']['down'][k_idx, band_idx] += spin_down_weight
                    
                    # Orbital projections for each spin channel
                    for spin_name, coeffs in [('up', spin_up_coeffs), ('down', spin_down_coeffs)]:
                        for orbital_group, orbital_indices in self.orbital_groups.items():
                            orbital_weight = np.sum(np.abs(coeffs[orbital_indices])**2)
                            projections['orbital'][orbital_group][k_idx, band_idx] += orbital_weight
                            
                            # Store orbital-site specific projections
                            key = f"{orbital_group}_{site}_{spin_name}"
                            if key not in projections['orbital_site']:
                                projections['orbital_site'][key] = np.zeros((num_kpoints, total_orbitals))
                            projections['orbital_site'][key][k_idx, band_idx] = orbital_weight
        
        return projections
    
    def plot_projected_bands(self, kdistances: np.ndarray, band_energies: np.ndarray,
                           projections: Dict, klabels: List[str], 
                           projection_type: str = 'orbital', output_file: str = None):
        """
        Plot band structure with projections as colors/sizes
        
        Args:
            projection_type: 'orbital', 'site', 'spin', or specific orbital-site combination
        """
        
        plt.figure(figsize=(12, 8))
        
        if projection_type == 'orbital':
            self._plot_orbital_projected(kdistances, band_energies, projections, klabels)
        elif projection_type == 'site':
            self._plot_site_projected(kdistances, band_energies, projections, klabels)
        elif projection_type == 'spin':
            self._plot_spin_projected(kdistances, band_energies, projections, klabels)
        else:
            self._plot_combined_projected(kdistances, band_energies, projections, klabels)
        
        # Set up k-point labels
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:
                tick_positions.append(kdistances[i])
                tick_labels.append(label)
        
        plt.xticks(tick_positions, tick_labels)
        
        # Add vertical lines at high-symmetry points
        for pos in tick_positions:
            plt.axvline(x=pos, color='k', linestyle='--', alpha=0.3)
        
        plt.xlabel('k-path')
        plt.ylabel('Energy (eV)')
        plt.title(f'TB-LMTO Band Structure - {projection_type.title()} Projected')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Projected band structure saved to {output_file}")
        
        plt.show()
    
    def _plot_orbital_projected(self, kdistances: np.ndarray, band_energies: np.ndarray,
                              projections: Dict, klabels: List[str]):
        """Plot orbital projected bands"""
        
        colors = {'s': 'red', 'p': 'blue', 'd': 'green'}
        
        for band in range(band_energies.shape[1]):
            for orbital_group, color in colors.items():
                weights = projections['orbital'][orbital_group][:, band]
                
                # Create scatter plot with point sizes proportional to weights
                sizes = weights * 50  # Scale factor for visibility
                plt.scatter(kdistances, band_energies[:, band], 
                          s=sizes, c=color, alpha=0.6, 
                          label=f'{orbital_group}' if band == 0 else "")
    
    def _plot_site_projected(self, kdistances: np.ndarray, band_energies: np.ndarray,
                           projections: Dict, klabels: List[str]):
        """Plot site projected bands"""
        
        num_sites = len(projections['site'])
        colors = plt.cm.Set1(np.linspace(0, 1, num_sites))
        
        for band in range(band_energies.shape[1]):
            for site, color in enumerate(colors):
                weights = projections['site'][site][:, band]
                
                # Create scatter plot with point sizes proportional to weights
                sizes = weights * 50
                plt.scatter(kdistances, band_energies[:, band], 
                          s=sizes, c=[color], alpha=0.6,
                          label=f'Site {site+1}' if band == 0 else "")
    
    def _plot_spin_projected(self, kdistances: np.ndarray, band_energies: np.ndarray,
                           projections: Dict, klabels: List[str]):
        """Plot spin projected bands"""
        
        colors = {'up': 'red', 'down': 'blue'}
        
        for band in range(band_energies.shape[1]):
            for spin, color in colors.items():
                weights = projections['spin'][spin][:, band]
                
                # Create scatter plot with point sizes proportional to weights
                sizes = weights * 50
                plt.scatter(kdistances, band_energies[:, band], 
                          s=sizes, c=color, alpha=0.6,
                          label=f'Spin {spin}' if band == 0 else "")
    
    def _plot_combined_projected(self, kdistances: np.ndarray, band_energies: np.ndarray,
                               projections: Dict, klabels: List[str]):
        """Plot combined orbital and site projections"""
        
        # Plot background bands in gray
        for band in range(band_energies.shape[1]):
            plt.plot(kdistances, band_energies[:, band], 'gray', alpha=0.3, linewidth=0.5)
        
        # Overlay orbital projections with different colors
        orbital_colors = {'s': 'red', 'p': 'blue', 'd': 'green'}
        
        for band in range(band_energies.shape[1]):
            for orbital_group, color in orbital_colors.items():
                weights = projections['orbital'][orbital_group][:, band]
                
                # Only plot points where projection is significant
                mask = weights > 0.1
                if np.any(mask):
                    sizes = weights[mask] * 10
                    plt.scatter(kdistances[mask], band_energies[mask, band], 
                              s=sizes, c=color, alpha=0.7,
                              label=f'{orbital_group}' if band == 0 else "")
    
    def plot_projection_analysis(self, kdistances: np.ndarray, band_energies: np.ndarray,
                               projections: Dict, klabels: List[str], 
                               energy_range: Tuple[float, float] = None,
                               output_file: str = None):
        """Create comprehensive projection analysis plots"""
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Filter energy range if specified
        if energy_range:
            energy_mask = ((band_energies >= energy_range[0]) & 
                          (band_energies <= energy_range[1]))
        else:
            energy_mask = np.ones_like(band_energies, dtype=bool)
        
        # Plot 1: Orbital projections
        ax1 = axes[0, 0]
        colors = {'s': 'red', 'p': 'blue', 'd': 'green'}
        for orbital_group, color in colors.items():
            for band in range(band_energies.shape[1]):
                weights = projections['orbital'][orbital_group][:, band]
                mask = energy_mask[:, band] & (weights > 0.05)
                if np.any(mask):
                    sizes = weights[mask] * 3
                    ax1.scatter(kdistances[mask], band_energies[mask, band], 
                              s=sizes, c=color, alpha=0.6)
        
        ax1.set_xlabel('k-path')
        ax1.set_ylabel('Energy (eV)')
        ax1.set_title('Orbital Projections')
        ax1.legend(colors.keys())
        
        # Plot 2: Site projections
        ax2 = axes[0, 1]
        num_sites = len(projections['site'])
        site_colors = plt.cm.Set1(np.linspace(0, 1, num_sites))
        for site, color in enumerate(site_colors):
            for band in range(band_energies.shape[1]):
                weights = projections['site'][site][:, band]
                mask = energy_mask[:, band] & (weights > 0.05)
                if np.any(mask):
                    sizes = weights[mask] * 30
                    ax2.scatter(kdistances[mask], band_energies[mask, band], 
                              s=sizes, c=[color], alpha=0.6)
        
        ax2.set_xlabel('k-path')
        ax2.set_ylabel('Energy (eV)')
        ax2.set_title('Site Projections')
        ax2.legend([f'Site {i+1}' for i in range(num_sites)])
        
        # Plot 3: Spin projections
        ax3 = axes[1, 0]
        spin_colors = {'up': 'red', 'down': 'blue'}
        for spin, color in spin_colors.items():
            for band in range(band_energies.shape[1]):
                weights = projections['spin'][spin][:, band]
                mask = energy_mask[:, band] & (weights > 0.05)
                if np.any(mask):
                    sizes = weights[mask] * 30
                    ax3.scatter(kdistances[mask], band_energies[mask, band], 
                              s=sizes, c=color, alpha=0.6)
        
        ax3.set_xlabel('k-path')
        ax3.set_ylabel('Energy (eV)')
        ax3.set_title('Spin Projections')
        ax3.legend(spin_colors.keys())
        
        # Plot 4: Combined view with background bands
        ax4 = axes[1, 1]
        # Plot all bands in gray background
        for band in range(band_energies.shape[1]):
            ax4.plot(kdistances, band_energies[:, band], 'gray', alpha=0.3, linewidth=0.5)
        
        # Overlay most significant projections
        for band in range(band_energies.shape[1]):
            for orbital_group, color in colors.items():
                weights = projections['orbital'][orbital_group][:, band]
                mask = energy_mask[:, band] & (weights > 0.2)  # Higher threshold
                if np.any(mask):
                    sizes = weights[mask] * 40
                    ax4.scatter(kdistances[mask], band_energies[mask, band], 
                              s=sizes, c=color, alpha=0.8)
        
        ax4.set_xlabel('k-path')
        ax4.set_ylabel('Energy (eV)')
        ax4.set_title('Combined View (Major Contributions)')
        
        # Set k-point labels for all subplots
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:
                tick_positions.append(kdistances[i])
                tick_labels.append(label)
        
        for ax in axes.flat:
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)
            for pos in tick_positions:
                ax.axvline(x=pos, color='k', linestyle='--', alpha=0.3)
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Projection analysis saved to {output_file}")
        
        plt.show()


def main():
    """Main function demonstrating projected band calculations"""
    
    import argparse
    from pathlib import Path
    
    parser = argparse.ArgumentParser(description='Calculate projected TB-LMTO band structure')
    parser.add_argument('--lattice', default='lattice.nml', 
                       help='Lattice file (default: lattice.nml)')
    parser.add_argument('--hamiltonian', default='rs2txt.dat',
                       help='Hamiltonian file (default: rs2txt.dat)')
    parser.add_argument('--output-prefix', default='projected_bands',
                       help='Output file prefix (default: projected_bands)')
    parser.add_argument('--nk', type=int, default=100,
                       help='Number of k-points (default: 100)')
    parser.add_argument('--projection', choices=['orbital', 'site', 'spin', 'all'],
                       default='all', help='Type of projection to calculate')
    parser.add_argument('--energy-range', nargs=2, type=float,
                       help='Energy range for analysis (emin emax)')
    
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
    
    print("Setting up projected band calculator...")
    calculator = ProjectedBandCalculator(lattice_reader, hamiltonian_reader)
    
    print(f"Generating k-path with {args.nk} points...")
    kpoints, klabels, kdistances = calculator.setup_kpath(args.nk)
    
    print("Calculating projected band structure...")
    band_energies, projections = calculator.calculate_projected_bands(kpoints)
    
    energy_range = tuple(args.energy_range) if args.energy_range else None
    
    if args.projection == 'all':
        print("Creating comprehensive projection analysis...")
        calculator.plot_projection_analysis(
            kdistances, band_energies, projections, klabels, 
            energy_range=energy_range,
            output_file=f"{args.output_prefix}_analysis.png"
        )
        
        # Also create individual projection plots
        for proj_type in ['orbital', 'site', 'spin']:
            calculator.plot_projected_bands(
                kdistances, band_energies, projections, klabels,
                projection_type=proj_type,
                output_file=f"{args.output_prefix}_{proj_type}.png"
            )
    else:
        print(f"Creating {args.projection} projected bands...")
        calculator.plot_projected_bands(
            kdistances, band_energies, projections, klabels,
            projection_type=args.projection,
            output_file=f"{args.output_prefix}_{args.projection}.png"
        )
    
    print("Done!")


if __name__ == "__main__":
    main()
