#!/usr/bin/env python3
"""
Non-collinear spin analysis for TB-LMTO Hamiltonians

This script provides tools to analyze the spin structure of non-collinear
TB-LMTO Hamiltonians with the structure:
[H0+Hz   Hx+iHy]
[Hx-iHy  H0-Hz ]

Where H0 is the charge part and Hx, Hy, Hz are the spin components.

SpinMomentAnalyzer is the preferred public name; SpinAnalyzer is kept as a
backward-compatible alias.
"""

import numpy as np
from typing import Dict, Tuple

from ...projections import ProjectedBandCalculator
from .plotting import plot_orbital_magnetism, plot_spin_texture


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
        return plot_spin_texture(
            kdistances,
            band_energies,
            spin_moments,
            klabels,
            band_indices=band_indices,
            output_file=output_file,
        )
    
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
        return plot_orbital_magnetism(
            kdistances,
            band_energies,
            projections,
            klabels,
            output_file=output_file,
        )


class SpinMomentAnalyzer(SpinAnalyzer):
    """Alias for SpinAnalyzer (kept for clarity in magnetic moment workflows)."""
    pass


if __name__ == "__main__":
    from .cli import main

    main()
