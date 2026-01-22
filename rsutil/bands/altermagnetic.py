#!/usr/bin/env python3
"""
Altermagnetic Analysis from Band Structure

Detects and analyzes altermagnetic signatures: spin splitting without net magnetization.
Alternets are a new class of magnetic materials with combined properties of 
ferromagnets (spin splitting) and antiferromagnets (zero net magnetization).
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from scipy.interpolate import interp1d


class AltermagneticAnalyzer:
    """Analyze altermagnetic properties from collinear band structure
    
    Alternets exhibit:
    1. Spin-dependent band structures (E_up ≠ E_down)
    2. Zero net magnetization (N_up ≈ N_down)
    3. Spin-dependent transport (AHE without magnetization)
    4. Spin-orbit-like effects without SOC
    
    Parameters
    ----------
    band_energies_up : np.ndarray
        Spin-up band energies (nk × nbands)
    band_energies_down : np.ndarray
        Spin-down band energies (nk × nbands)
    kpoints : np.ndarray
        k-point coordinates (nk × 3)
    kdistances : np.ndarray
        Distances along k-path for integration
    fermi_energy : float, optional
        Fermi energy in eV. Default: 0.0
        
    Attributes
    ----------
    spin_splitting : np.ndarray
        Spin splitting per band: E_up - E_down
    altermagnetic_strength : float
        Measure of altermagnetic character (0-1)
    """
    
    def __init__(self, band_energies_up: np.ndarray, band_energies_down: np.ndarray,
                 kpoints: np.ndarray, kdistances: np.ndarray, fermi_energy: float = 0.0):
        """Initialize altermagnetic analyzer"""
        self.band_energies_up = band_energies_up
        self.band_energies_down = band_energies_down
        self.kpoints = kpoints
        self.kdistances = kdistances
        self.fermi_energy = fermi_energy
        
        # Store dimensions
        self.nk = len(kdistances)
        self.nbands = band_energies_up.shape[1]
        
        # Compute spin splitting
        self.spin_splitting = band_energies_up - band_energies_down
        
    def calculate_net_magnetization(self) -> float:
        """Calculate net magnetization from band occupations
        
        Returns
        -------
        magnetization : float
            Net magnetic moment (should be ~0 for altermagnets)
        """
        # Simple estimate: count occupied states
        ef = self.fermi_energy
        n_up = np.sum(self.band_energies_up < ef)
        n_down = np.sum(self.band_energies_down < ef)
        return float(n_up - n_down)
    
    def calculate_spin_splitting_strength(self) -> float:
        """Calculate overall spin splitting strength
        
        Measures how much bands are split by spin.
        
        Returns
        -------
        strength : float
            RMS spin splitting magnitude (eV)
        """
        # RMS of spin splitting across all bands and k-points
        rms_splitting = np.sqrt(np.mean(self.spin_splitting ** 2))
        return float(rms_splitting)
    
    def identify_altermagnetic_bands(self, magnetization_threshold: float = 0.1) -> Dict:
        """Identify which bands show altermagnetic character
        
        Parameters
        ----------
        magnetization_threshold : float
            Threshold for net magnetization (in units of electrons)
            
        Returns
        -------
        result : Dict
            Dictionary with keys:
            - 'is_altermagnetic': bool, True if system shows altermagnetic signatures
            - 'net_magnetization': Net magnetic moment
            - 'max_splitting': Maximum spin splitting (eV)
            - 'splitting_range': Energy range of largest splittings
            - 'altermagnetic_bands': List of band indices with splitting
            - 'strength': Overall altermagnetic strength (0-1)
        """
        net_mag = self.calculate_net_magnetization()
        splitting_str = self.calculate_spin_splitting_strength()
        max_splitting = np.max(np.abs(self.spin_splitting))
        
        # Find bands with significant splitting
        mean_splitting = np.mean(np.abs(self.spin_splitting), axis=0)
        alt_bands = np.where(mean_splitting > 0.5 * np.max(mean_splitting))[0]
        
        # Altermagnetic signature: small net moment, large splitting
        is_altermagnetic = (
            abs(net_mag) < magnetization_threshold and 
            splitting_str > 0.01  # At least 10 meV splitting
        )
        
        # Strength metric: normalized by typical band energies
        strength = min(1.0, splitting_str / 0.1)  # Saturate at 0.1 eV
        
        return {
            'is_altermagnetic': is_altermagnetic,
            'net_magnetization': net_mag,
            'max_splitting': max_splitting,
            'splitting_range': (np.min(self.spin_splitting), np.max(self.spin_splitting)),
            'altermagnetic_bands': alt_bands.tolist(),
            'strength': strength,
            'splitting_strength': splitting_str,
        }
    
    def calculate_spin_dependent_fermi_surfaces(self) -> Tuple[int, int]:
        """Count states at Fermi level for each spin
        
        Returns
        -------
        n_up_fermi : int
            Number of bands crossing Fermi level (spin-up)
        n_down_fermi : int
            Number of bands crossing Fermi level (spin-down)
        """
        # Bands crossing Fermi level (different in each spin channel)
        up_crosses = np.sum((self.band_energies_up < self.fermi_energy).any(axis=0))
        down_crosses = np.sum((self.band_energies_down < self.fermi_energy).any(axis=0))
        
        return int(up_crosses), int(down_crosses)
    
    def calculate_altermagnetic_gaps(self, prominence_threshold: float = 0.005) -> Dict:
        """Identify altermagnetic gaps in band structure
        
        Altermagnetic gaps occur where one spin channel is gapped while 
        the other is metallic.
        
        Parameters
        ----------
        prominence_threshold : float
            Minimum gap size to consider (eV)
            
        Returns
        -------
        gaps : Dict
            Dictionary with gap information
        """
        from scipy.signal import find_peaks
        
        gaps = {
            'up_gaps': [],
            'down_gaps': [],
            'opposing_gaps': [],  # Up gapped while down metallic
        }
        
        # Simple gap detection per k-point
        for ik in range(len(self.kdistances)):
            # Find minima and maxima in density to identify gaps
            up_band = self.band_energies_up[ik, :]
            down_band = self.band_energies_down[ik, :]
            
            # Find band extrema
            up_peaks, _ = find_peaks(up_band, prominence=prominence_threshold)
            down_peaks, _ = find_peaks(down_band, prominence=prominence_threshold)
            
            if len(up_peaks) > 0 and len(down_peaks) == 0:
                gaps['opposing_gaps'].append(ik)
        
        return gaps
    
    def calculate_band_pair_splitting(self, band_indices: List[int]) -> Dict:
        """Calculate spin splitting for specific band pairs
        
        Useful for analyzing particular crossing points.
        
        Parameters
        ----------
        band_indices : List[int]
            Indices of bands to analyze
            
        Returns
        -------
        result : Dict
            Splitting analysis for selected bands
        """
        result = {}
        
        for band_idx in band_indices:
            if band_idx >= self.band_energies_up.shape[1]:
                continue
            
            splitting = self.band_energies_up[:, band_idx] - self.band_energies_down[:, band_idx]
            
            result[f'band_{band_idx}'] = {
                'max_splitting': float(np.max(np.abs(splitting))),
                'mean_splitting': float(np.mean(np.abs(splitting))),
                'splitting_at_fermi': float(splitting[np.argmin(np.abs(
                    np.mean([self.band_energies_up[:, band_idx], 
                            self.band_energies_down[:, band_idx]], axis=0) - 
                    self.fermi_energy))]),
                'crossing_points': np.where(np.sign(splitting[:-1]) != np.sign(splitting[1:]))[0].tolist(),
            }
        
        return result
    
    def detect_spin_degeneracies(self, tolerance: float = 0.005) -> List[Tuple[int, int, float]]:
        """Find points where spin-up and spin-down bands are degenerate
        
        Parameters
        ----------
        tolerance : float
            Energy tolerance for degeneracy (eV)
            
        Returns
        -------
        degeneracies : List[Tuple[int, int, float]]
            List of (k_index, band_index, splitting) at degeneracies
        """
        degeneracies = []
        
        for ik in range(len(self.kdistances)):
            for band in range(self.band_energies_up.shape[1]):
                splitting = self.spin_splitting[ik, band]
                if abs(splitting) < tolerance:
                    degeneracies.append((ik, band, float(splitting)))
        
        return degeneracies
    
    def plot_spin_splitting_map(self, output_file: Optional[str] = None) -> None:
        """Plot spin splitting as function of k and band
        
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
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Plot spin splitting for all bands
        im = ax.contourf(self.kdistances, range(self.band_energies_up.shape[1]),
                         np.abs(self.spin_splitting).T, levels=20, cmap='RdYlBu')
        
        # Add contour lines
        cs = ax.contour(self.kdistances, range(self.band_energies_up.shape[1]),
                       np.abs(self.spin_splitting).T, levels=10, colors='k', alpha=0.3)
        ax.clabel(cs, inline=True, fontsize=8)
        
        ax.set_xlabel('k-path', fontsize=12)
        ax.set_ylabel('Band Index', fontsize=12)
        ax.set_title('Spin Splitting |E_up - E_down| (eV)', fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('|ΔE| (eV)', fontsize=11)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Spin splitting map saved to {output_file}")
        
        plt.show()
    
    def plot_altermagnetic_character(self, output_file: Optional[str] = None) -> None:
        """Plot altermagnetic character indicators
        
        Shows: net magnetization vs spin splitting strength
        
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
        
        # Calculate per-band properties
        mag_per_band = np.zeros(self.band_energies_up.shape[1])
        split_per_band = np.zeros(self.band_energies_up.shape[1])
        
        for band in range(self.band_energies_up.shape[1]):
            # Magnetization contribution
            up_occupied = np.sum(self.band_energies_up[:, band] < self.fermi_energy)
            down_occupied = np.sum(self.band_energies_down[:, band] < self.fermi_energy)
            mag_per_band[band] = up_occupied - down_occupied
            
            # Splitting strength
            split_per_band[band] = np.mean(np.abs(self.spin_splitting[:, band]))
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        scatter = ax.scatter(mag_per_band, split_per_band, 
                           c=range(len(mag_per_band)), cmap='viridis',
                           s=100, alpha=0.7, edgecolors='black', linewidth=1)
        
        ax.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', linewidth=0.8, alpha=0.3)
        
        # Add regions
        ax.fill_between([-10, 10], 0, 10, alpha=0.1, color='blue', label='Altermagnetic region')
        
        ax.set_xlabel('Magnetization per band (electrons)', fontsize=12)
        ax.set_ylabel('Spin Splitting Strength (eV)', fontsize=12)
        ax.set_title('Altermagnetic Character Analysis', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-10, 10)
        ax.set_ylim(0, max(split_per_band) * 1.2)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Band Index', fontsize=11)
        
        ax.legend(loc='upper right')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Altermagnetic character plot saved to {output_file}")
        
        plt.show()
    
    def print_altermagnetic_summary(self) -> None:
        """Print formatted altermagnetic analysis summary"""
        result = self.identify_altermagnetic_bands()
        up_fs, down_fs = self.calculate_spin_dependent_fermi_surfaces()
        
        print()
        print("=" * 70)
        print("Altermagnetic Analysis Summary")
        print("=" * 70)
        print()
        print(f"Altermagnetic System:        {result['is_altermagnetic']}")
        print(f"Altermagnetic Strength:      {result['strength']:.4f} (0-1 scale)")
        print()
        print(f"Net Magnetization:           {result['net_magnetization']:.4f} μ_B")
        print(f"Spin Splitting Strength:     {result['splitting_strength']:.4f} eV")
        print(f"Maximum Splitting:           {result['max_splitting']:.4f} eV")
        print()
        print(f"Spin-Up Fermi Surface:       {up_fs} bands")
        print(f"Spin-Down Fermi Surface:     {down_fs} bands")
        print(f"Fermi Surface Difference:    {abs(up_fs - down_fs)} bands")
        print()
        if result['altermagnetic_bands']:
            print(f"Bands with Splitting:        {result['altermagnetic_bands']}")
        print()
        print("=" * 70)
        print()
