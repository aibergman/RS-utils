#!/usr/bin/env python3
"""
Enhanced DOS calculator with tetrahedron method and projections.

Main class for calculating density of states over the full Brillouin zone using
either the tetrahedron method (recommended) or Gaussian broadening.

Supports total, orbital-projected, site-projected, and spin-projected DOS.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional

from .integration import TetrahedronDOS
from .mesh import FullBZCalculator
from .dos_moments import DOSMomentsCalculator, analyze_dos_moments


class DOSCalculator:
    """Enhanced DOS calculator with tetrahedron method and projections
    
    Parameters
    ----------
    full_bz_calculator : FullBZCalculator
        Calculator for full Brillouin zone k-point generation and band structure
    """
    
    def __init__(self, full_bz_calculator: FullBZCalculator):
        self.calculator = full_bz_calculator
    
    def calculate_dos(self, mesh_size: Tuple[int, int, int], 
                     energy_range: Optional[Tuple[float, float]] = None,
                     num_energies: int = 1000,
                     use_symmetry: bool = True,
                     use_tetrahedron: bool = True,
                     gaussian_smearing: float = 0.1) -> Dict:
        """
        Calculate total and projected DOS
        
        Parameters
        ----------
        mesh_size : tuple
            (nx, ny, nz) k-point mesh dimensions
        energy_range : tuple, optional
            (emin, emax) energy range for DOS calculation
        num_energies : int
            Number of points on energy grid
        use_symmetry : bool
            Use symmetry reduction for k-mesh
        use_tetrahedron : bool
            Use tetrahedron method (vs Gaussian broadening)
        gaussian_smearing : float
            Gaussian smearing width (eV), used if not tetrahedron method
            
        Returns
        -------
        dos_data : dict
            Dictionary with DOS data and metadata
        """
        
        print("Generating k-point mesh...")
        # For the tetrahedron method we must use the full uniform mesh so
        # the tetrahedra connectivity corresponds to the regular grid.
        # Symmetry reduction breaks the regular connectivity and leads to
        # invalid tetrahedra (and spurious peaks). Force use_symmetry=False
        # when using tetrahedron.
        mesh_symmetry = use_symmetry if not use_tetrahedron else False
        kpoints, weights, tetrahedra = self.calculator.generate_kpoint_mesh(
            mesh_size, use_symmetry=mesh_symmetry)
        
        print("Calculating band structure...")
        band_energies, projections = self.calculator.calculate_full_bz_bands(
            kpoints, include_projections=True)
        
        # Determine energy range
        if energy_range is None:
            emin = np.min(band_energies) - 1.0
            emax = np.max(band_energies) + 1.0
        else:
            emin, emax = energy_range
        
        energy_grid = np.linspace(emin, emax, num_energies)
        
        dos_data = {
            'energies': energy_grid,
            'total_dos': None,
            'orbital_dos': {},
            'site_dos': {},
            'spin_dos': {},
            'kpoints': kpoints,
            'weights': weights,
            'band_energies': band_energies
        }
        
        if use_tetrahedron and len(tetrahedra) > 0:
            print("Calculating DOS using tetrahedron method...")
            tet_dos = TetrahedronDOS(kpoints, tetrahedra, weights)
            
            # Total DOS
            dos_data['total_dos'] = tet_dos.calculate_dos(band_energies, energy_grid)
            
            # Projected DOS
            if projections:
                self._calculate_projected_dos_tetrahedron(
                    tet_dos, band_energies, energy_grid, projections, dos_data)
        
        else:
            print("Calculating DOS using Gaussian broadening...")
            self._calculate_dos_gaussian(band_energies, energy_grid, weights, 
                                       projections, dos_data, gaussian_smearing)
        
        return dos_data
    
    def calculate_dos_with_moments(self, mesh_size: Tuple[int, int, int],
                                  target_electrons: Optional[float] = None,
                                  temperature: float = 300.0,
                                  energy_range: Optional[Tuple[float, float]] = None,
                                  num_energies: int = 1000,
                                  use_symmetry: bool = True,
                                  use_tetrahedron: bool = True,
                                  gaussian_smearing: float = 0.1,
                                  output_prefix: str = "dos_analysis") -> Tuple[Dict, DOSMomentsCalculator]:
        """
        Calculate DOS and analyze moments (occupations, band centers, widths)
        
        Parameters
        ----------
        mesh_size : tuple
            (nx, ny, nz) k-point mesh dimensions
        target_electrons : float, optional
            Number of electrons to determine Fermi level
        temperature : float
            Temperature in K for Fermi-Dirac distribution
        energy_range : tuple, optional
            (emin, emax) energy range for DOS
        num_energies : int
            Number of energy points
        use_symmetry : bool
            Use symmetry reduction
        use_tetrahedron : bool
            Use tetrahedron method
        output_prefix : str
            Prefix for output files
            
        Returns
        -------
        dos_data : dict
            DOS data dictionary
        moments_calc : DOSMomentsCalculator
            Moments calculator with results
        """
        
        # Calculate DOS
        print("Calculating DOS...")
        dos_data = self.calculate_dos(
            mesh_size, energy_range, num_energies,
            use_symmetry, use_tetrahedron, gaussian_smearing)
        
        if dos_data['total_dos'] is None:
            raise ValueError("DOS calculation failed")
        
        # Prepare projected DOS for moments analysis
        projected_dos = {}
        
        # Add orbital projections
        for orbital_name, orbital_dos in dos_data['orbital_dos'].items():
            projected_dos[f'orbital_{orbital_name}'] = orbital_dos
        
        # Add site projections  
        for site_name, site_dos in dos_data['site_dos'].items():
            projected_dos[f'site_{site_name}'] = site_dos
        
        # Add spin projections
        for spin_name, spin_dos in dos_data['spin_dos'].items():
            projected_dos[f'spin_{spin_name}'] = spin_dos
        
        # Calculate moments
        print("Analyzing DOS moments...")
        moments_calc = analyze_dos_moments(
            dos_data['energies'], 
            dos_data['total_dos'],
            projected_dos,
            target_electrons,
            temperature,
            output_prefix
        )
        
        # Add moments data to dos_data
        dos_data['moments'] = moments_calc.moments
        dos_data['fermi_level'] = moments_calc.fermi_level
        dos_data['temperature'] = temperature
        
        return dos_data, moments_calc
    
    def _calculate_projected_dos_tetrahedron(self, tet_dos: TetrahedronDOS,
                                           band_energies: np.ndarray, energy_grid: np.ndarray,
                                           projections: Dict, dos_data: Dict):
        """Calculate projected DOS using tetrahedron method"""
        
        # Orbital projections
        if 'orbital' in projections:
            for orbital in projections['orbital']:
                proj_weights = projections['orbital'][orbital][:, :, np.newaxis]
                dos_data['orbital_dos'][orbital] = tet_dos.calculate_dos(
                    band_energies, energy_grid, proj_weights)[:, 0]
        
        # Site projections
        if 'site' in projections:
            for site in projections['site']:
                proj_weights = projections['site'][site][:, :, np.newaxis]
                dos_data['site_dos'][site] = tet_dos.calculate_dos(
                    band_energies, energy_grid, proj_weights)[:, 0]
        
        # Spin projections
        if 'spin' in projections:
            for spin in projections['spin']:
                proj_weights = projections['spin'][spin][:, :, np.newaxis]
                dos_data['spin_dos'][spin] = tet_dos.calculate_dos(
                    band_energies, energy_grid, proj_weights)[:, 0]
    
    def _calculate_dos_gaussian(self, band_energies: np.ndarray, energy_grid: np.ndarray,
                              weights: np.ndarray, projections: Dict, dos_data: Dict,
                              broadening: float = 0.1):
        """Fallback Gaussian broadening method"""
        
        dos_data['total_dos'] = np.zeros(len(energy_grid))
        
        # Initialize projected DOS
        if projections:
            if 'orbital' in projections:
                for orbital in projections['orbital']:
                    dos_data['orbital_dos'][orbital] = np.zeros(len(energy_grid))
            if 'site' in projections:
                for site in projections['site']:
                    dos_data['site_dos'][site] = np.zeros(len(energy_grid))
            if 'spin' in projections:
                for spin in projections['spin']:
                    dos_data['spin_dos'][spin] = np.zeros(len(energy_grid))
        
        # Calculate DOS
        for k_idx in range(len(band_energies)):
            weight = weights[k_idx]
            
            for band in range(band_energies.shape[1]):
                energy = band_energies[k_idx, band]
                
                # Gaussian contribution
                gaussian = np.exp(-((energy_grid - energy) / broadening)**2)
                gaussian /= (broadening * np.sqrt(np.pi))
                
                dos_data['total_dos'] += weight * gaussian
                
                # Projected contributions
                if projections:
                    if 'orbital' in projections:
                        for orbital in projections['orbital']:
                            proj_weight = projections['orbital'][orbital][k_idx, band]
                            dos_data['orbital_dos'][orbital] += weight * proj_weight * gaussian
                    
                    if 'site' in projections:
                        for site in projections['site']:
                            proj_weight = projections['site'][site][k_idx, band]
                            dos_data['site_dos'][site] += weight * proj_weight * gaussian
                    
                    if 'spin' in projections:
                        for spin in projections['spin']:
                            proj_weight = projections['spin'][spin][k_idx, band]
                            dos_data['spin_dos'][spin] += weight * proj_weight * gaussian
    
    def plot_dos(self, dos_data: Dict, fermi_energy: float = 0.0, 
                output_file: str = None, show_projections: bool = True):
        """Plot comprehensive DOS analysis"""
        
        energies = dos_data['energies']
        
        if show_projections:
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            axes = [ax]
        
        # Plot total DOS
        ax = axes[0]
        ax.plot(energies, dos_data['total_dos'], 'k-', linewidth=2, label='Total DOS')
        ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7, label='E_F')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS (states/eV)')
        ax.set_title('Total Density of States')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if show_projections and len(axes) > 1:
            # Orbital DOS
            ax = axes[1]
            colors = ['red', 'blue', 'green']
            for i, (orbital, color) in enumerate(zip(['s', 'p', 'd'], colors)):
                if orbital in dos_data['orbital_dos']:
                    ax.plot(energies, dos_data['orbital_dos'][orbital], 
                           color=color, linewidth=2, label=f'{orbital}-orbital')
            ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7)
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('DOS (states/eV)')
            ax.set_title('Orbital Projected DOS')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Site DOS
            if len(axes) > 2:
                ax = axes[2]
                num_sites = len(dos_data['site_dos'])
                colors = plt.cm.Set1(np.linspace(0, 1, num_sites))
                for site, color in zip(dos_data['site_dos'], colors):
                    ax.plot(energies, dos_data['site_dos'][site], 
                           color=color, linewidth=2, label=f'Site {site+1}')
                ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7)
                ax.set_xlabel('Energy (eV)')
                ax.set_ylabel('DOS (states/eV)')
                ax.set_title('Site Projected DOS')
                ax.legend()
                ax.grid(True, alpha=0.3)
            
            # Spin DOS
            if len(axes) > 3:
                ax = axes[3]
                if 'up' in dos_data['spin_dos'] and 'down' in dos_data['spin_dos']:
                    ax.plot(energies, dos_data['spin_dos']['up'], 
                           'r-', linewidth=2, label='Spin up')
                    ax.plot(energies, -dos_data['spin_dos']['down'], 
                           'b-', linewidth=2, label='Spin down')
                    ax.axhline(0, color='k', linestyle='-', alpha=0.3)
                ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7)
                ax.set_xlabel('Energy (eV)')
                ax.set_ylabel('DOS (states/eV)')
                ax.set_title('Spin Projected DOS')
                ax.legend()
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"DOS plot saved to {output_file}")
        
        plt.show()
