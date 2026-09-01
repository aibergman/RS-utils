"""
Spin Texture and Spin-Orbit Coupling Analysis

This module provides tools for analyzing spin textures, spin-momentum locking,
and extracting spin-orbit coupling parameters from band structures with stored
eigenvectors.

Key Features:
- Spin expectation value calculation ⟨ψ|σ|ψ⟩
- Spin texture visualization on Fermi surfaces
- Rashba and Dresselhaus parameter extraction
- Spin-momentum locking identification
- Spin Hall angle estimation

Phase 4 Implementation
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

from .texture_plotting import (
    plot_spin_texture_2d as plot_spin_texture_2d_fn,
    plot_fermi_surface_spin_texture as plot_fermi_surface_spin_texture_fn,
)


class SpinTextureCalculator:
    """
    Calculate and analyze spin textures from band structure eigenvectors.
    
    For non-collinear TB-LMTO with 18×18 Hamiltonian, the spinor structure is:
    - First 9 orbitals: spin-up components
    - Last 9 orbitals: spin-down components
    
    The spin expectation values are computed from the Pauli x, y, and z
    operators applied to the spin-up and spin-down components.
    
    Parameters
    ----------
    band_energies : np.ndarray
        Band energies at k-points (nk × nbands)
    eigenvectors : np.ndarray
        Eigenvectors from Hamiltonian diagonalization (nk × nbands × nbands)
        Complex array where eigenvectors[ik, :, n] is the nth eigenvector at k-point ik
    kpoints : np.ndarray
        K-points in crystal coordinates (nk × 3)
    reciprocal_lattice : np.ndarray
        Reciprocal lattice vectors as columns (3 × 3), units of 2π/a
    fermi_energy : float, optional
        Fermi energy for Fermi surface analysis
    
    Attributes
    ----------
    spin_texture : np.ndarray
        Spin expectation values (nk × nbands × 3) for the x, y, and z axes
    spin_polarization : np.ndarray
        Magnitude of spin polarization (nk × nbands)
    """
    
    def __init__(self, band_energies: np.ndarray, eigenvectors: np.ndarray,
                 kpoints: np.ndarray, reciprocal_lattice: np.ndarray,
                 fermi_energy: Optional[float] = None):
        self.energies = band_energies
        self.eigenvectors = eigenvectors
        self.kpoints = kpoints
        self.reciprocal_lattice = reciprocal_lattice
        self.fermi_energy = fermi_energy
        
        self.nk, self.nbands = band_energies.shape
        
        # Determine orbital structure (assumes 18×18 for non-collinear)
        self.n_orb_per_spin = self.nbands // 2
        
        # Storage for computed properties
        self.spin_texture = None
        self.spin_polarization = None
        self.rashba_params = None
        self.dresselhaus_params = None
    
    def calculate_spin_texture(self) -> np.ndarray:
        """
        Calculate the spin expectation vector for all bands at all k-points.
        
        Returns
        -------
        spin_texture : np.ndarray
            Spin vectors (nk × nbands × 3), with the last axis storing x, y, and z
        
        Notes
        -----
        For each band n at k-point ik, we compute:
            The nth spinor at k is formed by concatenating its spin-up and
            spin-down orbital components.
        
        The spin components are extracted and used to calculate Pauli matrix
        expectation values.
        """
        spin_texture = np.zeros((self.nk, self.nbands, 3), dtype=float)
        
        for ik in range(self.nk):
            for n in range(self.nbands):
                # Extract eigenvector for band n at k-point ik
                psi = self.eigenvectors[ik, :, n]
                
                # Split into spin-up and spin-down components
                psi_up = psi[:self.n_orb_per_spin]
                psi_down = psi[self.n_orb_per_spin:]
                
                # Calculate spin expectation values
                # ⟨σx⟩ = ⟨ψ↑|ψ↓⟩ + ⟨ψ↓|ψ↑⟩ = 2 Re(⟨ψ↑|ψ↓⟩)
                spin_x = 2 * np.real(np.vdot(psi_up, psi_down))
                
                # ⟨σy⟩ = -i(⟨ψ↑|ψ↓⟩ - ⟨ψ↓|ψ↑⟩) = 2 Im(⟨ψ↑|ψ↓⟩)
                spin_y = 2 * np.imag(np.vdot(psi_up, psi_down))
                
                # ⟨σz⟩ = ⟨ψ↑|ψ↑⟩ - ⟨ψ↓|ψ↓⟩
                spin_z = np.real(np.vdot(psi_up, psi_up) - np.vdot(psi_down, psi_down))
                
                spin_texture[ik, n, :] = [spin_x, spin_y, spin_z]
        
        self.spin_texture = spin_texture
        
        # Calculate spin polarization magnitude
        self.spin_polarization = np.linalg.norm(spin_texture, axis=2)
        
        return spin_texture
    
    def identify_spin_momentum_locking(self, band_idx: int, 
                                      energy_window: float = 0.1) -> Dict:
        """
        Identify spin-momentum locking for a specific band.
        
        Spin-momentum locking occurs when the spin direction is locked
        perpendicular to momentum (Rashba) or parallel/antiparallel 
        (Dresselhaus) to momentum direction.
        
        Parameters
        ----------
        band_idx : int
            Band index to analyze
        energy_window : float
            Energy window around Fermi level (in eV) for analysis
            
        Returns
        -------
        result : dict
            Dictionary containing:
            - 'locking_type': 'rashba', 'dresselhaus', 'mixed', or 'none'
            - 'locking_strength': measure of locking quality (0 to 1)
            - 'spin_orientation': dominant spin direction
            - 'k_vectors': k-points within energy window
        """
        if self.spin_texture is None:
            self.calculate_spin_texture()
        
        if self.fermi_energy is None:
            raise ValueError("Fermi energy must be set for Fermi surface analysis")
        
        # Find k-points near Fermi level for this band
        band_energies = self.energies[:, band_idx]
        near_ef = np.abs(band_energies - self.fermi_energy) < energy_window
        
        if not np.any(near_ef):
            return {
                'locking_type': 'none',
                'locking_strength': 0.0,
                'spin_orientation': None,
                'k_vectors': []
            }
        
        # Get k-vectors and spin vectors near Fermi level
        k_fermi = self.kpoints[near_ef]
        spin_fermi = self.spin_texture[near_ef, band_idx, :]
        
        # Convert k-points to Cartesian coordinates
        k_cart = k_fermi @ self.reciprocal_lattice.T
        
        # Normalize k-vectors and spin vectors
        k_norm = k_cart / (np.linalg.norm(k_cart, axis=1, keepdims=True) + 1e-10)
        s_norm = spin_fermi / (np.linalg.norm(spin_fermi, axis=1, keepdims=True) + 1e-10)
        
        # Calculate angle between k and spin
        # Rashba: spin ⊥ k (angle ≈ 90°)
        # Dresselhaus: spin ∥ k or spin ∥ -k (angle ≈ 0° or 180°)
        cos_angle = np.abs(np.sum(k_norm * s_norm, axis=1))
        
        # Average alignment
        avg_cos = np.mean(cos_angle)
        std_cos = np.std(cos_angle)
        
        # Classification
        if avg_cos < 0.3:  # Nearly perpendicular
            locking_type = 'rashba'
            locking_strength = 1.0 - avg_cos
        elif avg_cos > 0.7:  # Nearly parallel/antiparallel
            locking_type = 'dresselhaus'
            locking_strength = avg_cos
        else:
            locking_type = 'mixed'
            locking_strength = 0.5
        
        # Determine dominant spin orientation
        avg_spin = np.mean(spin_fermi, axis=0)
        avg_spin /= (np.linalg.norm(avg_spin) + 1e-10)
        
        return {
            'locking_type': locking_type,
            'locking_strength': locking_strength,
            'spin_orientation': avg_spin,
            'k_vectors': k_fermi,
            'alignment_std': std_cos,
            'n_points': len(k_fermi)
        }
    
    def extract_rashba_parameters(self, band_idx: int, 
                                  k_direction: np.ndarray = np.array([1, 0, 0])) -> Dict:
        """
        Extract Rashba parameter αR from dispersion near k=0.
        
        The Rashba Hamiltonian is:
            H_R = αR (σ × k) · ẑ
        
        For 2D systems (kz=0 plane), this gives:
            E±(k) ≈ E₀ ± αR|k|
        
        Parameters
        ----------
        band_idx : int
            Band index showing Rashba splitting
        k_direction : np.ndarray
            Direction in k-space for linear fit (default: [1,0,0])
            
        Returns
        -------
        params : dict
            Dictionary with:
            - 'alpha_R': Rashba parameter (eV·Å)
            - 'E_0': Band center energy
            - 'split_at_k': Splitting at reference k-point
            - 'fit_quality': R² of linear fit
        """
        # Find k-points near Gamma (|k| < 0.2)
        k_cart = self.kpoints @ self.reciprocal_lattice.T
        k_mag = np.linalg.norm(k_cart, axis=1)
        
        near_gamma = k_mag < 0.2
        
        if np.sum(near_gamma) < 5:
            raise ValueError("Insufficient k-points near Gamma for Rashba fit")
        
        # Get energies and k-magnitudes
        k_vals = k_mag[near_gamma]
        E_vals = self.energies[near_gamma, band_idx]
        
        # Linear fit: E = E₀ + αR*k
        # Use both positive and negative slopes for Rashba pair
        coeffs = np.polyfit(k_vals, E_vals, deg=1)
        alpha_R = coeffs[0]
        E_0 = coeffs[1]
        
        # Calculate fit quality
        E_fit = np.polyval(coeffs, k_vals)
        ss_res = np.sum((E_vals - E_fit)**2)
        ss_tot = np.sum((E_vals - np.mean(E_vals))**2)
        r_squared = 1 - ss_res / (ss_tot + 1e-10)
        
        return {
            'alpha_R': abs(alpha_R),
            'E_0': E_0,
            'split_at_k': 2 * abs(alpha_R) * 0.1,  # Splitting at k=0.1 Å⁻¹
            'fit_quality': r_squared,
            'k_range': [k_vals.min(), k_vals.max()],
            'n_points': len(k_vals)
        }
    
    def calculate_spin_hall_conductivity(self, bands: List[int], 
                                        temperature: float = 300) -> Dict:
        """
        Estimate intrinsic spin Hall conductivity from spin texture and Berry curvature.
        
        The spin Hall conductivity is related to the spin Berry curvature:
            σ_SH = (e/ℏ) ∑_n ∫ f(En) Ω_n^s(k) d²k
        
        Parameters
        ----------
        bands : list of int
            Band indices to include in calculation
        temperature : float
            Temperature in Kelvin for Fermi-Dirac distribution
            
        Returns
        -------
        result : dict
            Spin Hall conductivity components and analysis
            
        Notes
        -----
        This is a simplified estimation. Full calculation requires
        Berry curvature with spin projection operators.
        """
        if self.spin_texture is None:
            self.calculate_spin_texture()
        
        if self.fermi_energy is None:
            raise ValueError("Fermi energy required for conductivity calculation")
        
        # Fermi-Dirac distribution
        kB = 8.617e-5  # eV/K
        beta = 1.0 / (kB * temperature)
        
        sigma_sh = np.zeros(3)
        
        for band_idx in bands:
            energies = self.energies[:, band_idx]
            f_fermi = 1.0 / (1.0 + np.exp(beta * (energies - self.fermi_energy)))
            
            # Weighted sum of spin polarization near Fermi level
            spin_contrib = self.spin_texture[:, band_idx, :] * f_fermi[:, np.newaxis]
            sigma_sh += np.sum(spin_contrib, axis=0)
        
        # Normalize by k-space volume (simplified)
        sigma_sh /= self.nk
        
        # Convert to physical units (rough estimate)
        e_over_hbar = 1.519e15  # e/ℏ in SI units (A·m⁻²)
        sigma_sh *= e_over_hbar
        
        return {
            'sigma_SH': sigma_sh,
            'magnitude': np.linalg.norm(sigma_sh),
            'dominant_component': ['x', 'y', 'z'][np.argmax(np.abs(sigma_sh))],
            'temperature': temperature,
            'bands_included': bands
        }
    
    def plot_spin_texture_2d(self, band_idx: int, component: str = 'z',
                            plane: str = 'xy', plane_position: float = 0.0,
                            output_file: Optional[str] = None):
        """
        Plot 2D spin texture with arrows showing spin direction.
        
        Parameters
        ----------
        band_idx : int
            Band index to visualize
        component : str
            Spin component for coloring: 'x', 'y', 'z', or 'magnitude'
        plane : str
            k-space plane: 'xy', 'xz', or 'yz'
        plane_position : float
            Position along perpendicular direction (in crystal coords)
        output_file : str, optional
            Save plot to file if provided
        """
        if self.spin_texture is None:
            self.calculate_spin_texture()
        
        return plot_spin_texture_2d_fn(
            self.kpoints,
            self.spin_texture,
            self.spin_polarization,
            band_idx,
            component=component,
            plane=plane,
            plane_position=plane_position,
            output_file=output_file,
        )
    
    def plot_fermi_surface_spin_texture(self, band_idx: int,
                                       energy_window: float = 0.05,
                                       output_file: Optional[str] = None):
        """
        Plot 3D Fermi surface colored by spin polarization.
        
        Parameters
        ----------
        band_idx : int
            Band index to visualize
        energy_window : float
            Energy window around Fermi level (eV)
        output_file : str, optional
            Save plot to file if provided
        """
        if self.spin_texture is None:
            self.calculate_spin_texture()
        
        if self.fermi_energy is None:
            raise ValueError("Fermi energy must be set for Fermi surface plotting")
        
        return plot_fermi_surface_spin_texture_fn(
            self.kpoints,
            self.energies,
            self.spin_texture,
            self.fermi_energy,
            band_idx,
            energy_window=energy_window,
            output_file=output_file,
        )
