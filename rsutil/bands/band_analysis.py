#!/usr/bin/env python3
"""
Band structure analysis utilities

This module provides additional analysis tools for the TB-LMTO band structure results.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from .band_structure import BandStructureCalculator


class BandAnalyzer:
    """Additional analysis tools for band structure data"""
    
    def __init__(self, calculator: BandStructureCalculator):
        self.calculator = calculator
    
    def find_band_gap(self, band_energies: np.ndarray, fermi_energy: float = 0.0) -> Tuple[float, str]:
        """
        Find the band gap
        
        Args:
            band_energies: Band energies array (nk, nbands)
            fermi_energy: Fermi energy reference (default: 0.0)
            
        Returns:
            gap_size: Size of the band gap
            gap_type: 'direct' or 'indirect'
        """
        
        # Find valence and conduction bands
        valence_max = -np.inf
        conduction_min = np.inf
        valence_k_idx = -1
        conduction_k_idx = -1
        
        for k_idx in range(band_energies.shape[0]):
            energies = band_energies[k_idx, :]
            
            # Find highest occupied and lowest unoccupied
            occupied = energies[energies <= fermi_energy]
            unoccupied = energies[energies > fermi_energy]
            
            if len(occupied) > 0:
                vb_max = np.max(occupied)
                if vb_max > valence_max:
                    valence_max = vb_max
                    valence_k_idx = k_idx
            
            if len(unoccupied) > 0:
                cb_min = np.min(unoccupied)
                if cb_min < conduction_min:
                    conduction_min = cb_min
                    conduction_k_idx = k_idx
        
        gap_size = conduction_min - valence_max
        gap_type = 'direct' if valence_k_idx == conduction_k_idx else 'indirect'
        
        return gap_size, gap_type
    
    def calculate_dos(self, band_energies: np.ndarray, energy_range: Tuple[float, float] = None,
                     num_points: int = 1000, broadening: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate density of states using Gaussian broadening
        
        Args:
            band_energies: Band energies array (nk, nbands)
            energy_range: (min_energy, max_energy) or None for auto
            num_points: Number of energy points
            broadening: Gaussian broadening parameter
            
        Returns:
            energies: Energy grid
            dos: Density of states
        """
        
        if energy_range is None:
            emin = np.min(band_energies) - 2 * broadening
            emax = np.max(band_energies) + 2 * broadening
        else:
            emin, emax = energy_range
        
        energies = np.linspace(emin, emax, num_points)
        dos = np.zeros(num_points)
        
        # Weight per k-point (assuming uniform sampling)
        weight = 1.0 / band_energies.shape[0]
        
        for k_idx in range(band_energies.shape[0]):
            for band in range(band_energies.shape[1]):
                energy = band_energies[k_idx, band]
                
                # Add Gaussian contribution
                gaussian = np.exp(-((energies - energy) / broadening)**2)
                dos += weight * gaussian / (broadening * np.sqrt(np.pi))
        
        return energies, dos
    
    def plot_dos(self, energies: np.ndarray, dos: np.ndarray, 
                fermi_energy: float = 0.0, output_file: str = None):
        """Plot density of states"""
        
        plt.figure(figsize=(8, 6))
        plt.plot(energies, dos, 'b-', linewidth=2)
        plt.axvline(x=fermi_energy, color='r', linestyle='--', 
                   label=f'Fermi level ({fermi_energy:.2f} eV)')
        
        plt.xlabel('Energy (eV)')
        plt.ylabel('Density of States (states/eV)')
        plt.title('Density of States')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"DOS plot saved to {output_file}")
        
        plt.show()
    
    def plot_bands_and_dos(self, kdistances: np.ndarray, band_energies: np.ndarray,
                          klabels: list, dos_energies: np.ndarray, dos: np.ndarray,
                          fermi_energy: float = 0.0, output_file: str = None):
        """Plot band structure and DOS side by side"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), 
                                      gridspec_kw={'width_ratios': [2, 1]})
        
        # Plot band structure
        for band in range(band_energies.shape[1]):
            ax1.plot(kdistances, band_energies[:, band], 'b-', linewidth=1)
        
        # Set up k-point labels
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:
                tick_positions.append(kdistances[i])
                tick_labels.append(label)
        
        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels(tick_labels)
        
        # Add vertical lines at high-symmetry points
        for pos in tick_positions:
            ax1.axvline(x=pos, color='k', linestyle='--', alpha=0.3)
        
        ax1.axhline(y=fermi_energy, color='r', linestyle='--', alpha=0.7)
        ax1.set_xlabel('k-path')
        ax1.set_ylabel('Energy (eV)')
        ax1.set_title('Band Structure')
        ax1.grid(True, alpha=0.3)
        
        # Plot DOS
        ax2.plot(dos, dos_energies, 'b-', linewidth=2)
        ax2.axhline(y=fermi_energy, color='r', linestyle='--', alpha=0.7)
        ax2.set_xlabel('DOS (states/eV)')
        ax2.set_ylabel('Energy (eV)')
        ax2.set_title('Density of States')
        ax2.grid(True, alpha=0.3)
        
        # Match y-axis ranges
        y_min = min(np.min(band_energies), np.min(dos_energies))
        y_max = max(np.max(band_energies), np.max(dos_energies))
        ax1.set_ylim(y_min, y_max)
        ax2.set_ylim(y_min, y_max)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Combined plot saved to {output_file}")
        
        plt.show()
    
    def analyze_effective_mass(self, kdistances: np.ndarray, band_energies: np.ndarray,
                              band_index: int, k_point_index: int) -> float:
        """
        Calculate effective mass at a specific band and k-point
        
        Args:
            kdistances: K-point distances
            band_energies: Band energies array
            band_index: Which band to analyze
            k_point_index: Which k-point to analyze around
            
        Returns:
            effective_mass: Effective mass in units of electron mass
        """
        
        # Need at least 5 points for numerical differentiation
        if len(kdistances) < 5:
            raise ValueError("Need at least 5 k-points for effective mass calculation")
        
        # Get energies for the specific band
        energies = band_energies[:, band_index]
        
        # Calculate second derivative using finite differences
        if k_point_index < 2:
            k_point_index = 2
        elif k_point_index >= len(kdistances) - 2:
            k_point_index = len(kdistances) - 3
        
        # Use 5-point stencil for second derivative
        h = kdistances[k_point_index + 1] - kdistances[k_point_index]
        
        d2E_dk2 = (-energies[k_point_index + 2] + 16*energies[k_point_index + 1] 
                  - 30*energies[k_point_index] + 16*energies[k_point_index - 1] 
                  - energies[k_point_index - 2]) / (12 * h**2)
        
        # Convert to effective mass (in units of electron mass)
        # m* = ħ²/(d²E/dk²) with appropriate units
        hbar_eV_s = 6.582119569e-16  # eV⋅s
        m_e = 9.1093837015e-31  # kg
        
        if abs(d2E_dk2) < 1e-10:
            return np.inf  # Flat band
        
        effective_mass = (hbar_eV_s**2) / (d2E_dk2 * m_e)
        
        return effective_mass


def main():
    """Demonstration of analysis tools"""
    
    from band_structure_calculator import LatticeReader, HamiltonianReader
    from pathlib import Path
    
    # Check if files exist
    if not Path("lattice.nml").exists() or not Path("rs2txt.dat").exists():
        print("Please ensure lattice.nml and rs2txt.dat exist in current directory")
        return
    
    print("Loading data...")
    lattice_reader = LatticeReader("lattice.nml")
    hamiltonian_reader = HamiltonianReader("rs2txt.dat")
    calculator = BandStructureCalculator(lattice_reader, hamiltonian_reader)
    
    print("Calculating band structure...")
    kpoints, klabels, kdistances = calculator.setup_kpath(100)
    band_energies = calculator.calculate_bands(kpoints)
    
    print("Analyzing results...")
    analyzer = BandAnalyzer(calculator)
    
    # Find band gap
    gap_size, gap_type = analyzer.find_band_gap(band_energies)
    print(f"Band gap: {gap_size:.3f} eV ({gap_type})")
    
    # Calculate DOS
    dos_energies, dos = analyzer.calculate_dos(band_energies)
    
    # Create combined plot
    analyzer.plot_bands_and_dos(kdistances, band_energies, klabels, 
                               dos_energies, dos, output_file="bands_and_dos.png")


if __name__ == "__main__":
    main()
