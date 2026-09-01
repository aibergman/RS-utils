#!/usr/bin/env python3
"""
DOS Moments Calculator for TB-LMTO

This module calculates:
  - Fermi level from total DOS
  - Spin and orbital projected moments:
    - Zeroth moment: occupation numbers
    - First moment: average energies
    - Second moment: energy spread (width)
"""

import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq
from typing import Dict, List, Tuple, Optional
import warnings


class DOSMomentsCalculator:
    """Calculate DOS moments and Fermi level"""

    def __init__(
        self,
        energies: np.ndarray,
        total_dos: np.ndarray,
        projected_dos: Optional[Dict] = None,
        temperature: float = 300.0,
    ):
        """
        Initialize DOS moments calculator

        Parameters
        ----------
        energies : array
            Energy grid (eV)
        total_dos : array
            Total DOS (states/eV)
        projected_dos : dict, optional
            Dictionary with projected DOS data
        temperature : float
            Temperature in K for Fermi-Dirac distribution
        """
        self.energies = energies
        self.total_dos = total_dos
        self.projected_dos = projected_dos or {}
        self.temperature = temperature
        self.kB = 8.617333e-5  # Boltzmann constant in eV/K

        self.fermi_level = None
        self.moments = {}

    def fermi_dirac(self, energy: float, fermi_level: float) -> float:
        """Fermi-Dirac distribution"""
        if self.temperature == 0:
            return 1.0 if energy <= fermi_level else 0.0

        x = (energy - fermi_level) / (self.kB * self.temperature)
        # Avoid overflow
        if x > 50:
            return 0.0
        elif x < -50:
            return 1.0
        else:
            return 1.0 / (1.0 + np.exp(x))

    def find_fermi_level(
        self,
        target_electrons: float,
        energy_range: Optional[Tuple[float, float]] = None,
    ) -> float:
        """
        Find Fermi level for given number of electrons

        Parameters
        ----------
        target_electrons : float
            Target number of electrons
        energy_range : tuple, optional
            Search range for Fermi level (emin, emax)

        Returns
        -------
        float
            Fermi level in eV
        """
        if energy_range is None:
            # Default search range
            emin = self.energies.min() - 1.0
            emax = self.energies.max() + 1.0
        else:
            emin, emax = energy_range

        def electron_count_difference(ef):
            """Calculate difference between actual and target electron count"""
            fd_factors = np.array([self.fermi_dirac(e, ef) for e in self.energies])
            integrand = (self.total_dos / 2) * fd_factors
            electron_count = simpson(integrand, self.energies)
            return electron_count - target_electrons

        try:
            # Check if solution exists in range
            diff_min = electron_count_difference(emin)
            diff_max = electron_count_difference(emax)

            if diff_min * diff_max > 0:
                warnings.warn(
                    f"No solution found in range [{emin:.3f}, {emax:.3f}] eV. "
                    f"Differences: {diff_min:.3f}, {diff_max:.3f}"
                )
                # Return energy where DOS is maximum near the middle
                mid_range = (emin + emax) / 2
                mask = (self.energies >= mid_range - 2) & (
                    self.energies <= mid_range + 2
                )
                if np.any(mask):
                    return self.energies[mask][np.argmax(self.total_dos[mask])]
                else:
                    return mid_range

            fermi_level = brentq(electron_count_difference, emin, emax, xtol=1e-6)
            self.fermi_level = fermi_level
            return fermi_level

        except ValueError as e:
            warnings.warn(f"Failed to find Fermi level: {e}")
            # Fallback: find energy where integrated DOS equals target
            cumulative_dos = np.cumsum(self.total_dos / 2) * (
                self.energies[1] - self.energies[0]
            )
            idx = np.argmin(np.abs(cumulative_dos - target_electrons))
            self.fermi_level = self.energies[idx]
            return self.fermi_level

    def calculate_moments(self, fermi_level: Optional[float] = None) -> Dict:
        """
        Calculate DOS moments for total and projected DOS

        Parameters
        ----------
        fermi_level : float, optional
            Fermi level. If None, use stored value

        Returns
        -------
        dict
            Dictionary containing moments for each DOS component
        """
        if fermi_level is None:
            fermi_level = self.fermi_level

        if fermi_level is None:
            raise ValueError("Fermi level not set. Call find_fermi_level() first.")

        results = {}

        # Calculate moments for total DOS
        results["total"] = self._calculate_component_moments(
            self.total_dos, fermi_level, "Total DOS"
        )

        # Calculate moments for projected DOS components
        if self.projected_dos:
            for component_name, dos_data in self.projected_dos.items():
                if isinstance(dos_data, np.ndarray):
                    results[component_name] = self._calculate_component_moments(
                        dos_data, fermi_level, component_name
                    )
                elif isinstance(dos_data, dict):
                    # Handle nested structure (e.g., spin-resolved)
                    results[component_name] = {}
                    for sub_name, sub_dos in dos_data.items():
                        if isinstance(sub_dos, np.ndarray):
                            results[component_name][sub_name] = (
                                self._calculate_component_moments(
                                    sub_dos, fermi_level, f"{component_name}_{sub_name}"
                                )
                            )

        self.moments = results
        return results

    def _calculate_component_moments(
        self, dos: np.ndarray, fermi_level: float, component_name: str
    ) -> Dict:
        """Calculate moments for a single DOS component"""

        # Fermi-Dirac factors
        fd_factors = np.array([self.fermi_dirac(e, fermi_level) for e in self.energies])

        # Energy differences from Fermi level
        energy_diff = self.energies - fermi_level

        # Weighted DOS (occupied states)
        weighted_dos = (dos / 2) * fd_factors

        # Zeroth moment: occupation (number of electrons)
        n0 = simpson(weighted_dos, self.energies)

        # First moment: average energy relative to Fermi level
        if n0 > 1e-10:  # Avoid division by zero
            n1_integrand = weighted_dos * energy_diff
            n1 = simpson(n1_integrand, self.energies) / n0
        else:
            n1 = 0.0

        # Second moment: energy spread (width squared)
        if n0 > 1e-10:
            n2_integrand = weighted_dos * (energy_diff - n1) ** 2
            n2 = simpson(n2_integrand, self.energies) / n0
        else:
            n2 = 0.0

        # Band center (absolute energy)
        band_center = fermi_level + n1

        # Band width (RMS width)
        band_width = np.sqrt(n2)

        return {
            "occupation": n0,
            "avg_energy_rel_ef": n1,
            "energy_spread_sq": n2,
            "band_center": band_center,
            "band_width": band_width,
            "fermi_level": fermi_level,
        }

    def print_summary(self, components: Optional[List[str]] = None):
        """Print summary of DOS moments"""

        if not self.moments:
            print("No moments calculated. Run calculate_moments() first.")
            return

        print(f"\n{'='*60}")
        print("DOS MOMENTS SUMMARY")
        print(f"{'='*60}")
        print(f"Fermi Level: {self.fermi_level:.6f} eV")
        print(f"Temperature: {self.temperature:.1f} K")
        print(f"{'='*60}")

        if components is None:
            components = list(self.moments.keys())

        for comp_name in components:
            if comp_name not in self.moments:
                continue

            comp_data = self.moments[comp_name]

            if isinstance(comp_data, dict) and "occupation" in comp_data:
                # Single component
                self._print_component_summary(comp_name, comp_data)
            elif isinstance(comp_data, dict):
                # Nested components (e.g., spin-resolved)
                print(f"\n{comp_name}:")
                print("-" * 40)
                for sub_name, sub_data in comp_data.items():
                    if isinstance(sub_data, dict) and "occupation" in sub_data:
                        self._print_component_summary(
                            f"  {sub_name}", sub_data, indent=True
                        )

    def _print_component_summary(self, name: str, data: Dict, indent: bool = False):
        """Print summary for a single component"""
        prefix = "  " if indent else ""

        print(f"\n{prefix}{name}:")
        if not indent:
            print("-" * 40)

        print(f"{prefix}  Occupation:     {data['occupation']:.6f} electrons")
        print(f"{prefix}  Band Center:    {data['band_center']:.6f} eV")
        print(f"{prefix}  Band Width:     {data['band_width']:.6f} eV")
        print(f"{prefix}  Avg E (rel EF): {data['avg_energy_rel_ef']:.6f} eV")

    def save_moments(self, filename: str):
        """Save moments data to numpy file"""
        if not self.moments:
            print("No moments to save. Run calculate_moments() first.")
            return

        save_data = {
            "fermi_level": self.fermi_level,
            "temperature": self.temperature,
            "moments": self.moments,
            "energies": self.energies,
        }

        np.savez(filename, **save_data)
        print(f"DOS moments saved to {filename}")


def analyze_dos_moments(
    energies: np.ndarray,
    total_dos: np.ndarray,
    projected_dos: Optional[Dict] = None,
    target_electrons: Optional[float] = None,
    temperature: float = 300.0,
    output_prefix: str = "dos_moments",
) -> DOSMomentsCalculator:
    """
    Convenience function to analyze DOS moments

    Parameters
    ----------
    energies : array
        Energy grid
    total_dos : array
        Total DOS
    projected_dos : dict, optional
        Projected DOS components
    target_electrons : float, optional
        Number of electrons (if None, integrate total DOS to midpoint)
    temperature : float
        Temperature in K
    output_prefix : str
        Prefix for output files

    Returns
    -------
    DOSMomentsCalculator
        Calculator with computed moments
    """

    calculator = DOSMomentsCalculator(energies, total_dos, projected_dos, temperature)

    # Determine target electrons if not provided
    if target_electrons is None:
        # Integrate total DOS and use that as target
        total_electrons = simpson(total_dos / 2, energies)
        target_electrons = total_electrons
        print(f"Auto-determined target electrons: {target_electrons:.2f}")

    # Find Fermi level
    fermi_level = calculator.find_fermi_level(target_electrons)
    print(f"Fermi level: {fermi_level:.6f} eV")

    # Calculate moments
    moments = calculator.calculate_moments(fermi_level)

    # Print summary
    calculator.print_summary()

    # Save results
    calculator.save_moments(f"{output_prefix}_moments.npz")

    return calculator


if __name__ == "__main__":
    # Example usage
    import matplotlib.pyplot as plt

    # Generate test DOS
    energies = np.linspace(-10, 5, 1000)
    total_dos = np.exp(-((energies + 2) ** 2) / 4) + 0.5 * np.exp(
        -((energies - 1) ** 2) / 2
    )

    # Example projected DOS
    projected_dos = {
        "spin_up": 0.6 * total_dos,
        "spin_down": 0.4 * total_dos,
        "orbital_d": 0.8 * total_dos,
        "orbital_s": 0.2 * total_dos,
    }

    # Analyze
    calculator = analyze_dos_moments(
        energies, total_dos, projected_dos, target_electrons=10.0, temperature=300.0
    )
