#!/usr/bin/env python3
"""
Effective Mass Calculator from Band Structure

Calculates electron/hole effective masses from band curvature.
Supports both scalar and tensor effective masses.
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline


class EffectiveMassCalculator:
    """Calculate effective masses from band structure

    Uses Taylor expansion around band extrema to calculate effective masses:
    E(k) ≈ E_0 + (ħ²/2m*) |k - k_0|²

    Parameters
    ----------
    band_energies : np.ndarray
        Band energies for one band (nk,)
    kpoints : np.ndarray
        k-points in reciprocal space (nk × 3)
    k_lattice : Optional[np.ndarray]
        Reciprocal lattice vectors (3 × 3). Used for metric tensor.
        If None, assumes orthogonal k-space.

    Attributes
    ----------
    band_energies : np.ndarray
        Input band energies
    kpoints : np.ndarray
        Input k-points
    extrema : List[int]
        Indices of band extrema (maxima and minima)
    """

    def __init__(
        self,
        band_energies: np.ndarray,
        kpoints: np.ndarray,
        k_lattice: Optional[np.ndarray] = None,
    ):
        """Initialize effective mass calculator"""
        self.band_energies = band_energies
        self.kpoints = kpoints
        self.k_lattice = k_lattice
        self.extrema = []

        # Find band extrema
        self._find_extrema()

    def _find_extrema(self, prominence: float = 0.01) -> None:
        """Find local extrema in band structure"""
        from scipy.signal import find_peaks

        # Find maxima
        peaks, _ = find_peaks(self.band_energies, prominence=prominence)
        # Find minima (negative of band)
        valleys, _ = find_peaks(-self.band_energies, prominence=prominence)

        self.extrema = sorted(list(peaks) + list(valleys))

    def _quadratic_fit(
        self, k_range: np.ndarray, E_range: np.ndarray
    ) -> Tuple[float, float, float]:
        """Fit quadratic to energy vs |k|²

        Parameters
        ----------
        k_range : np.ndarray
            k-point distances from extremum
        E_range : np.ndarray
            Band energies near extremum

        Returns
        -------
        E0 : float
            Band energy at extremum
        curvature : float
            Curvature (d²E/dk²)
        R2 : float
            R² goodness of fit
        """

        # Fit: E = E0 + curvature * k²
        def quadratic(k, e0, curv):
            return e0 + curv * k**2

        try:
            popt, pcov = curve_fit(quadratic, k_range, E_range, maxfev=1000)
            e0, curv = popt

            # Calculate R²
            residuals = E_range - quadratic(k_range, e0, curv)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((E_range - np.mean(E_range)) ** 2)
            r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

            return e0, curv, r2
        except:
            return np.nan, np.nan, 0.0

    def calculate_effective_mass_at_extremum(
        self, extremum_idx: int, k_range: float = 0.2
    ) -> Dict[str, float]:
        """Calculate effective mass at a specific band extremum

        Parameters
        ----------
        extremum_idx : int
            Index of k-point at band extremum
        k_range : float
            Range in k-space to use for fitting (in units of 2π/a)

        Returns
        -------
        result : Dict[str, float]
            Dictionary with keys:
            - 'index': k-point index
            - 'energy': Band energy at extremum (eV)
            - 'k_point': k-point coordinates
            - 'effective_mass': Scalar effective mass (in units of free electron mass)
            - 'effective_mass_inv': 1/m* (for mobility calculations)
            - 'type': 'electron' (minimum) or 'hole' (maximum)
            - 'goodness_of_fit': R² value (should be > 0.95)
            - 'curvature': d²E/dk² from fit
        """
        if extremum_idx >= len(self.kpoints):
            raise ValueError(f"Invalid extremum index {extremum_idx}")

        # Get k-point at extremum
        k_ext = self.kpoints[extremum_idx]
        e_ext = self.band_energies[extremum_idx]

        # Find k-points within range
        distances = np.linalg.norm(self.kpoints - k_ext, axis=1)
        mask = distances <= k_range

        if np.sum(mask) < 5:
            raise ValueError(
                f"Not enough k-points near extremum (need ≥5, got {np.sum(mask)})"
            )

        # Extract data near extremum
        k_near = distances[mask]
        E_near = self.band_energies[mask]

        # Fit quadratic: E = E0 + (ħ²/2m*) k²
        e0, d2e_dk2, r2 = self._quadratic_fit(k_near, E_near)

        # Calculate effective mass
        # From: E(k) ≈ E0 + (ħ²/2m*) |k|²
        # We have: dE/d|k|² = ħ²/2m*
        # So: m* = ħ²/(2 * d²E/dk²)

        # Constants (in eV units)
        hbar_squared_over_me = 3.8099821e-2  # eV·Å²

        if abs(d2e_dk2) < 1e-10:
            m_eff = np.inf
        else:
            m_eff = hbar_squared_over_me / (2 * abs(d2e_dk2))

        # Determine type (electron or hole)
        if extremum_idx > 0 and extremum_idx < len(self.band_energies) - 1:
            e_prev = self.band_energies[extremum_idx - 1]
            e_next = self.band_energies[extremum_idx + 1]
            if e_prev < e_ext and e_next < e_ext:
                extremum_type = "hole"
            else:
                extremum_type = "electron"
        else:
            extremum_type = "electron" if d2e_dk2 > 0 else "hole"

        return {
            "index": extremum_idx,
            "energy": e_ext,
            "k_point": k_ext,
            "effective_mass": m_eff,
            "effective_mass_inv": 1.0 / m_eff if m_eff != np.inf else 0.0,
            "type": extremum_type,
            "goodness_of_fit": r2,
            "curvature": d2e_dk2,
        }

    def calculate_all_effective_masses(self, k_range: float = 0.2) -> List[Dict]:
        """Calculate effective masses at all band extrema

        Parameters
        ----------
        k_range : float
            Range in k-space for fitting

        Returns
        -------
        results : List[Dict]
            List of effective mass calculations for each extremum
        """
        results = []

        for ext_idx in self.extrema:
            try:
                result = self.calculate_effective_mass_at_extremum(ext_idx, k_range)
                # Only include if fit is good
                if result["goodness_of_fit"] > 0.9:
                    results.append(result)
            except:
                continue

        return results

    def get_band_gap_masses(
        self,
        band_energies_up: np.ndarray,
        band_energies_down: np.ndarray,
        kpoints: np.ndarray,
        k_range: float = 0.2,
    ) -> Dict[str, float]:
        """Calculate electron and hole masses at band gap

        Parameters
        ----------
        band_energies_up : np.ndarray
            Spin-up band energies (nk × nbands)
        band_energies_down : np.ndarray
            Spin-down band energies (nk × nbands)
        kpoints : np.ndarray
            k-point coordinates
        k_range : float
            Range for fitting

        Returns
        -------
        gap_info : Dict[str, float]
            Dictionary with:
            - 'band_gap': Band gap energy
            - 'vbm_energy': Valence band maximum
            - 'cbm_energy': Conduction band minimum
            - 'me_electron': Electron effective mass
            - 'me_hole': Hole effective mass
        """
        # Find VBM and CBM
        vbm_idx = np.argmax(np.max(band_energies_up, axis=0))
        cbm_idx = np.argmin(np.min(band_energies_down, axis=0))

        # Find VBM k-point
        vbm_k_idx = np.argmax(band_energies_up[:, vbm_idx])
        # Find CBM k-point
        cbm_k_idx = np.argmin(band_energies_down[:, cbm_idx])

        e_vbm = band_energies_up[vbm_k_idx, vbm_idx]
        e_cbm = band_energies_down[cbm_k_idx, cbm_idx]

        band_gap = e_cbm - e_vbm

        # Calculate masses
        try:
            vb_calc = EffectiveMassCalculator(
                band_energies_up[:, vbm_idx], kpoints, None
            )
            vb_mass_result = vb_calc.calculate_effective_mass_at_extremum(
                vbm_k_idx, k_range
            )
            me_hole = vb_mass_result["effective_mass"]
        except:
            me_hole = np.nan

        try:
            cb_calc = EffectiveMassCalculator(
                band_energies_down[:, cbm_idx], kpoints, None
            )
            cb_mass_result = cb_calc.calculate_effective_mass_at_extremum(
                cbm_k_idx, k_range
            )
            me_electron = cb_mass_result["effective_mass"]
        except:
            me_electron = np.nan

        return {
            "band_gap": band_gap,
            "vbm_energy": e_vbm,
            "cbm_energy": e_cbm,
            "me_electron": me_electron,
            "me_hole": me_hole,
        }

    def calculate_mobility(
        self, effective_mass: float, scattering_time: float = 1e-14
    ) -> float:
        """Calculate carrier mobility from effective mass

        μ = e * τ / m*

        Parameters
        ----------
        effective_mass : float
            Effective mass in units of electron mass
        scattering_time : float
            Scattering time in seconds. Default: 1e-14 s

        Returns
        -------
        mobility : float
            Mobility in cm²/(V·s)
        """
        if effective_mass <= 0 or effective_mass == np.inf:
            return 0.0

        # Constants
        e = 1.602e-19  # Elementary charge (C)
        m_e = 9.109e-31  # Electron mass (kg)

        # m* in kg
        m_star = effective_mass * m_e

        # Mobility in m²/(V·s)
        mobility = (e * scattering_time) / m_star

        # Convert to cm²/(V·s)
        return mobility * 1e4

    def print_effective_masses(self) -> None:
        """Print formatted effective mass results"""
        results = self.calculate_all_effective_masses()

        if not results:
            print("No band extrema found or fitting failed")
            return

        print()
        print("=" * 70)
        print("Effective Mass Analysis")
        print("=" * 70)
        print()
        print(f"{'Idx':>4} {'Type':>10} {'Energy':>10} {'m*':>10} {'μ':>12} {'R²':>7}")
        print(f"{'':>4} {'':>10} {'(eV)':>10} {'(m_e)':>10} {'(cm²/Vs)':>12} {'':>7}")
        print("-" * 70)

        for result in results:
            energy = result["energy"]
            m_eff = result["effective_mass"]
            e_type = result["type"]
            r2 = result["goodness_of_fit"]

            if np.isfinite(m_eff):
                mu = self.calculate_mobility(m_eff)
                print(
                    f"{result['index']:4d} {e_type:>10} {energy:>10.4f} {m_eff:>10.3f} {mu:>12.2e} {r2:>7.3f}"
                )
            else:
                print(
                    f"{result['index']:4d} {e_type:>10} {energy:>10.4f} {'∞':>10} {'N/A':>12} {r2:>7.3f}"
                )

        print("=" * 70)
        print()
