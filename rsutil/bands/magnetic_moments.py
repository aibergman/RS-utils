#!/usr/bin/env python3
"""
Magnetic Moment Calculator from Band Structure

Calculates total and site-resolved magnetic moments from spin-polarized
band structure calculations.
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from scipy.integrate import trapezoid


class MagneticMomentCalculator:
    """Calculate magnetic moments from collinear spin-polarized bands
    
    Uses the difference in spin-up and spin-down band occupations to
    compute total and site-resolved magnetic moments.
    
    Parameters
    ----------
    band_energies_up : np.ndarray
        Spin-up band energies (nk × nbands)
    band_energies_down : np.ndarray
        Spin-down band energies (nk × nbands)
    kdistances : np.ndarray
        Distances along k-path for integration weights
    fermi_energy : float, optional
        Fermi energy in eV. Default: 0.0
    temperature : float, optional
        Temperature in Kelvin for Fermi-Dirac distribution. Default: 0.0 (step function)
        
    Attributes
    ----------
    dos_up : np.ndarray
        Density of states for spin-up
    dos_down : np.ndarray
        Density of states for spin-down
    occupation_up : np.ndarray
        Occupation numbers for spin-up bands
    occupation_down : np.ndarray
        Occupation numbers for spin-down bands
    """
    
    def __init__(self, band_energies_up: np.ndarray, band_energies_down: np.ndarray,
                 kdistances: np.ndarray, fermi_energy: float = 0.0, temperature: float = 0.0):
        """Initialize magnetic moment calculator"""
        self.band_energies_up = band_energies_up
        self.band_energies_down = band_energies_down
        self.kdistances = kdistances
        self.fermi_energy = fermi_energy
        self.temperature = temperature
        
        # Will be computed on demand
        self.dos_up = None
        self.dos_down = None
        self.occupation_up = None
        self.occupation_down = None
        
    def fermi_dirac(self, energy: np.ndarray) -> np.ndarray:
        """Fermi-Dirac distribution function
        
        Parameters
        ----------
        energy : np.ndarray
            Energy values relative to Fermi level
            
        Returns
        -------
        np.ndarray
            Occupation probability (1 at E < EF, 0 at E > EF if T=0)
        """
        if self.temperature < 1e-6:  # T = 0 K
            return np.where(energy < 0, 1.0, 0.0)
        else:
            # Fermi-Dirac at finite temperature
            k_B = 8.617333e-5  # eV/K Boltzmann constant
            beta = 1.0 / (k_B * self.temperature)
            return 1.0 / (1.0 + np.exp(beta * energy))
    
    def calculate_occupations(self) -> Tuple[float, float]:
        """Calculate total occupations for each spin channel
        
        Returns
        -------
        n_up : float
            Total number of spin-up electrons
        n_down : float
            Total number of spin-down electrons
        """
        # Energy relative to Fermi level
        E_up = self.band_energies_up - self.fermi_energy
        E_down = self.band_energies_down - self.fermi_energy
        
        # Get occupations at each k-point
        f_up = self.fermi_dirac(E_up)  # (nk, nbands)
        f_down = self.fermi_dirac(E_down)  # (nk, nbands)
        
        # Integrate over k-path using trapezoidal rule
        n_up = 0.0
        n_down = 0.0
        
        for band in range(f_up.shape[1]):
            n_up += trapezoid(f_up[:, band], self.kdistances)
            n_down += trapezoid(f_down[:, band], self.kdistances)
        
        # Store for later use
        self.occupation_up = f_up
        self.occupation_down = f_down
        
        return n_up, n_down
    
    def calculate_total_magnetic_moment(self) -> float:
        """Calculate total magnetic moment in units of μ_B
        
        The magnetic moment is proportional to the difference in
        spin-up and spin-down electron occupations.
        
        M = (N_up - N_down) μ_B
        
        Returns
        -------
        moment : float
            Total magnetic moment in Bohr magnetons (μ_B)
            
        Notes
        -----
        Assumes 1 electron contributes 1 unit of magnetic moment.
        More detailed calculations would account for orbital contributions.
        """
        n_up, n_down = self.calculate_occupations()
        moment = n_up - n_down
        return moment
    
    def calculate_magnetic_moment_per_band(self) -> np.ndarray:
        """Calculate magnetic moment contribution from each band
        
        Returns
        -------
        moments : np.ndarray
            Magnetic moment per band (nbands,)
        """
        if self.occupation_up is None:
            self.calculate_occupations()
        
        moments = np.zeros(self.occupation_up.shape[1])
        
        for band in range(len(moments)):
            # Integrate occupation difference for each band
            diff_up = self.occupation_up[:, band]
            diff_down = self.occupation_down[:, band]
            moments[band] = trapezoid(diff_up - diff_down, self.kdistances)
        
        return moments
    
    def calculate_polarization(self) -> float:
        """Calculate spin polarization (Curie-Weiss parameter)
        
        P = (N_up - N_down) / (N_up + N_down)
        
        Returns
        -------
        polarization : float
            Spin polarization in range [-1, 1]
            -1: fully spin-down polarized
             0: non-magnetic
            +1: fully spin-up polarized
        """
        n_up, n_down = self.calculate_occupations()
        n_total = n_up + n_down
        
        if n_total < 1e-10:
            return 0.0
        
        return (n_up - n_down) / n_total
    
    def calculate_spin_density_at_fermi(self) -> float:
        """Calculate spin density of states at Fermi level
        
        Returns
        -------
        spin_dos : float
            Difference in DOS at Fermi level between spin channels
        """
        # Find bands closest to Fermi level
        e_min = min(self.band_energies_up.min(), self.band_energies_down.min())
        e_max = max(self.band_energies_up.max(), self.band_energies_down.max())
        
        # Create energy grid around Fermi level
        energies = np.linspace(e_min, e_max, 500)
        dos_up = np.zeros_like(energies)
        dos_down = np.zeros_like(energies)
        
        sigma = 0.01  # Gaussian broadening in eV
        
        # Calculate DOS using Gaussian broadening
        for band in range(self.band_energies_up.shape[1]):
            E_up = self.band_energies_up[:, band]
            E_down = self.band_energies_down[:, band]
            
            # Gaussian density
            for i, E in enumerate(energies):
                dos_up[i] += np.sum(np.exp(-((E_up - E) / sigma) ** 2)) / (sigma * np.sqrt(np.pi))
                dos_down[i] += np.sum(np.exp(-((E_down - E) / sigma) ** 2)) / (sigma * np.sqrt(np.pi))
        
        # Find DOS at Fermi level
        idx = np.argmin(np.abs(energies - self.fermi_energy))
        spin_dos = dos_up[idx] - dos_down[idx]
        
        return spin_dos
    
    def get_summary(self) -> Dict[str, float]:
        """Get summary of all calculated magnetic properties
        
        Returns
        -------
        summary : Dict[str, float]
            Dictionary with keys:
            - 'total_moment': Total magnetic moment (μ_B)
            - 'polarization': Spin polarization (-1 to 1)
            - 'n_up': Number of spin-up electrons
            - 'n_down': Number of spin-down electrons
            - 'spin_dos_at_ef': Spin DOS at Fermi level
            - 'moment_per_atom': Moment per atom (if 1 atom)
        """
        n_up, n_down = self.calculate_occupations()
        moment = self.calculate_total_magnetic_moment()
        polarization = self.calculate_polarization()
        spin_dos = self.calculate_spin_density_at_fermi()
        
        return {
            'total_moment': moment,
            'polarization': polarization,
            'n_up': n_up,
            'n_down': n_down,
            'spin_dos_at_ef': spin_dos,
            'moment_per_atom': moment,  # Per atom if single-atom calculation
        }
    
    def print_summary(self, name: str = "System") -> None:
        """Print formatted summary of magnetic properties
        
        Parameters
        ----------
        name : str
            System name for output header
        """
        summary = self.get_summary()
        
        print()
        print("=" * 60)
        print(f"Magnetic Moment Analysis: {name}")
        print("=" * 60)
        print(f"Fermi Energy:            {self.fermi_energy:>10.4f} eV")
        print(f"Temperature:             {self.temperature:>10.2f} K")
        print()
        print(f"Spin-Up Electrons:       {summary['n_up']:>10.4f}")
        print(f"Spin-Down Electrons:     {summary['n_down']:>10.4f}")
        print(f"Total Electrons:         {summary['n_up'] + summary['n_down']:>10.4f}")
        print()
        print(f"Total Magnetic Moment:   {summary['total_moment']:>10.4f} μ_B")
        print(f"Spin Polarization:       {summary['polarization']:>10.4f}")
        print(f"Spin DOS at E_F:         {summary['spin_dos_at_ef']:>10.4f} 1/eV")
        print("=" * 60)
        print()
    
    def plot_spin_resolved_dos(self, output_file: Optional[str] = None) -> None:
        """Plot spin-resolved density of states
        
        Parameters
        ----------
        output_file : str, optional
            Output filename for saving the plot
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available for plotting")
            return
        
        # Create energy grid
        e_min = min(self.band_energies_up.min(), self.band_energies_down.min())
        e_max = max(self.band_energies_up.max(), self.band_energies_down.max())
        energies = np.linspace(e_min, e_max, 200)
        dos_up = np.zeros_like(energies)
        dos_down = np.zeros_like(energies)
        
        sigma = 0.02  # Gaussian broadening
        
        # Calculate DOS
        for band in range(self.band_energies_up.shape[1]):
            E_up = self.band_energies_up[:, band]
            E_down = self.band_energies_down[:, band]
            
            for i, E in enumerate(energies):
                dos_up[i] += np.sum(np.exp(-((E_up - E) / sigma) ** 2)) / (sigma * np.sqrt(np.pi))
                dos_down[i] += np.sum(np.exp(-((E_down - E) / sigma) ** 2)) / (sigma * np.sqrt(np.pi))
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(dos_up, energies, 'r-', linewidth=2, label='Spin ↑')
        ax.plot(-dos_down, energies, 'b-', linewidth=2, label='Spin ↓')
        ax.fill_between(dos_up, energies, alpha=0.3, color='red')
        ax.fill_between(-dos_down, energies, alpha=0.3, color='blue')
        
        # Mark Fermi level
        ax.axhline(y=self.fermi_energy, color='k', linestyle='--', linewidth=1, alpha=0.5)
        ax.text(0.98, self.fermi_energy, 'E_F', ha='right', va='bottom', fontsize=10)
        
        ax.set_xlabel('Density of States (1/eV)', fontsize=12)
        ax.set_ylabel('Energy (eV)', fontsize=12)
        ax.set_title('Spin-Resolved Density of States', fontsize=14, fontweight='bold')
        ax.legend(loc='best', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Spin-resolved DOS plot saved to {output_file}")
        
        plt.show()
    
    def plot_band_occupations(self, output_file: Optional[str] = None) -> None:
        """Plot band occupations along k-path
        
        Parameters
        ----------
        output_file : str, optional
            Output filename for saving the plot
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available for plotting")
            return
        
        if self.occupation_up is None:
            self.calculate_occupations()
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot spin-up occupations
        for band in range(self.occupation_up.shape[1]):
            ax1.plot(self.kdistances, self.occupation_up[:, band], 'r-', linewidth=0.5, alpha=0.7)
        ax1.set_ylabel('Occupation', fontsize=11)
        ax1.set_title('Spin-Up Band Occupations', fontsize=12, fontweight='bold')
        ax1.set_ylim([-0.05, 1.05])
        ax1.grid(True, alpha=0.3)
        
        # Plot spin-down occupations
        for band in range(self.occupation_down.shape[1]):
            ax2.plot(self.kdistances, self.occupation_down[:, band], 'b-', linewidth=0.5, alpha=0.7)
        ax2.set_ylabel('Occupation', fontsize=11)
        ax2.set_xlabel('k-path', fontsize=11)
        ax2.set_title('Spin-Down Band Occupations', fontsize=12, fontweight='bold')
        ax2.set_ylim([-0.05, 1.05])
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Band occupations plot saved to {output_file}")
        
        plt.show()
