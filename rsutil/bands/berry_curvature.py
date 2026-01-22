#!/usr/bin/env python3
"""
Berry Curvature and Topological Properties from Band Structure

Calculates Berry curvature Ω(k) = ∇ × A(k) and related topological properties.
The Berry connection A_n(k) = i⟨ψ_n(k)|∇_k|ψ_n(k)⟩ is a geometric phase.
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from scipy.interpolate import LinearNDInterpolator, griddata


class BerryCurvatureCalculator:
    """Calculate Berry curvature and topological invariants
    
    Berry curvature appears in:
    1. Anomalous Hall effect: σ_xy = -e²/ℏ ∫ Ω_z(k) f(E_nk) dk
    2. Orbital magnetization: M = -e/2 ∑_nk f(E_nk) Ω(k)
    3. Chern numbers: C = 1/(2π) ∫_BZ Ω_z dk
    
    Parameters
    ----------
    energies : np.ndarray
        Band energies (nk × nbands) in eV
    eigenvectors : np.ndarray
        Eigenvectors (nk × nbands × nbands), eigenvector j in column j
    kpoints : np.ndarray
        k-point coordinates (nk × 3) in reciprocal lattice
    reciprocal_lattice : np.ndarray
        Reciprocal lattice vectors (3 × 3)
    fermi_energy : float
        Fermi energy in eV
    """
    
    def __init__(self, energies: np.ndarray, eigenvectors: np.ndarray,
                 kpoints: np.ndarray, reciprocal_lattice: np.ndarray,
                 fermi_energy: float = 0.0):
        """Initialize Berry curvature calculator"""
        self.energies = energies
        self.eigenvectors = eigenvectors  # (nk, nbands, nbands)
        self.kpoints = kpoints
        self.reciprocal_lattice = reciprocal_lattice
        self.fermi_energy = fermi_energy
        
        self.nk = energies.shape[0]
        self.nbands = energies.shape[1]
        
        # Compute k-space metric
        self._compute_kspace_metric()
        
    def _compute_kspace_metric(self) -> None:
        """Compute metric tensor in k-space for volume element"""
        # Volume element in k-space: dk = |det(b1, b2, b3)| dk1 dk2 dk3
        # For 2D: dk = |b1 × b2| dk1 dk2
        b1, b2, b3 = self.reciprocal_lattice
        volume = np.abs(np.dot(b1, np.cross(b2, b3)))
        self.kspace_volume = volume
        
    def calculate_berry_connection(self, band_idx: int, 
                                  delta_k: float = 0.001) -> np.ndarray:
        """Calculate Berry connection A_n(k) from eigenvectors
        
        Uses finite differences: A_i(k) = i⟨ψ_n(k)|∂ψ_n/∂k_i|ψ_n(k)⟩
        
        Parameters
        ----------
        band_idx : int
            Band index to calculate Berry connection for
        delta_k : float
            k-space step size for finite differences (atomic units)
            
        Returns
        -------
        berry_conn : np.ndarray
            Berry connection (nk × 3) in atomic units
        """
        berry_conn = np.zeros((self.nk, 3), dtype=complex)
        
        # For each k-point, estimate ∂ψ/∂k_i using finite differences
        # This is a simplified calculation that requires k-point ordering
        
        for ik in range(self.nk):
            psi_n = self.eigenvectors[ik, :, band_idx]
            
            # Get neighboring eigenvectors along each k direction
            for direction in range(3):
                # Forward difference (simplified - assumes regular grid)
                if ik + 1 < self.nk:
                    psi_n_plus = self.eigenvectors[min(ik + 1, self.nk - 1), :, band_idx]
                    dpsi_dk = (psi_n_plus - psi_n) / delta_k
                    
                    # Berry connection: A_i = i⟨ψ|∂ψ/∂k_i⟩
                    # Note: Need gauge fixing for phase consistency
                    berry_conn[ik, direction] = 1j * np.vdot(psi_n, dpsi_dk)
        
        return berry_conn
    
    def calculate_berry_curvature_2d(self, band_idx: int,
                                    delta_k: float = 0.001) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate 2D Berry curvature Ω_xy(k) in xy-plane
        
        Ω_xy = ∂A_y/∂k_x - ∂A_x/∂k_y
        
        Parameters
        ----------
        band_idx : int
            Band index
        delta_k : float
            k-space step size (atomic units)
            
        Returns
        -------
        berry_curv_2d : np.ndarray
            2D Berry curvature (nk,)
        berry_conn : np.ndarray
            Berry connection (nk × 2)
        """
        # Calculate Berry connection
        berry_conn = self.calculate_berry_connection(band_idx, delta_k)[:, :2]
        
        # Calculate derivatives using finite differences
        # This is simplified - requires proper k-point grid structure
        berry_curv_2d = np.zeros(self.nk)
        
        # For each k-point, estimate ∂A/∂k
        dk_grid = delta_k * np.ones(3)
        
        for ik in range(self.nk):
            if ik + 1 < self.nk:
                # dA_y/dk_x - dA_x/dk_y (simplified)
                berry_curv_2d[ik] = (
                    np.imag(berry_conn[ik, 1] - berry_conn[ik, 0])
                ) / delta_k
        
        return berry_curv_2d, berry_conn
    
    def calculate_berry_curvature_tensor(self, band_idx: int,
                                        delta_k: float = 0.001) -> np.ndarray:
        """Calculate full 3D Berry curvature tensor Ω_ij(k)
        
        Ω_ij = ∂A_j/∂k_i - ∂A_i/∂k_j (antisymmetric)
        
        Parameters
        ----------
        band_idx : int
            Band index
        delta_k : float
            k-space step size
            
        Returns
        -------
        berry_curv : np.ndarray
            Berry curvature tensor (nk × 3 × 3)
        """
        # Get Berry connection
        berry_conn = self.calculate_berry_connection(band_idx, delta_k)  # (nk × 3)
        
        # Initialize tensor
        berry_curv = np.zeros((self.nk, 3, 3), dtype=complex)
        
        # Approximate derivatives (simplified for non-uniform k-grid)
        for i in range(3):
            for j in range(3):
                if i != j:
                    # ∂A_j/∂k_i - ∂A_i/∂k_j
                    berry_curv[:, i, j] = (
                        np.roll(berry_conn[:, j], -1) - 
                        np.roll(berry_conn[:, j], 1)
                    ) / (2 * delta_k)
        
        return berry_curv
    
    def calculate_chern_number(self, band_idx: int, 
                              integration_method: str = 'trapezoidal') -> float:
        """Calculate Chern number for a given band
        
        C_n = 1/(2π) ∫_BZ Ω_z,n(k) d²k
        
        The Chern number is a topological invariant that counts the number of
        edge states in a 2D system. It must be an integer for gapped bands.
        
        Parameters
        ----------
        band_idx : int
            Band index
        integration_method : str
            Integration method ('trapezoidal', 'simpsons')
            
        Returns
        -------
        chern_number : float
            Chern number (should be integer for gapped bands)
        """
        # Calculate 2D Berry curvature
        omega_z, _ = self.calculate_berry_curvature_2d(band_idx)
        
        # Integrate over Brillouin zone
        if integration_method == 'trapezoidal':
            # Simple trapezoidal rule
            chern = np.trapz(omega_z) / (2 * np.pi)
        else:
            chern = np.trapz(omega_z) / (2 * np.pi)
        
        return float(chern)
    
    def calculate_anomalous_hall_conductivity(self, band_idx: int,
                                             temperature: float = 0.0,
                                             relaxation_time: float = 1e-13) -> float:
        """Calculate anomalous Hall conductivity from Berry curvature
        
        σ_xy = -e²/ℏ ∫ Ω_z(k) f(E_nk) d³k / (2π)³
        
        Parameters
        ----------
        band_idx : int
            Band index
        temperature : float
            Temperature in K (for Fermi-Dirac)
        relaxation_time : float
            Scattering relaxation time in seconds
            
        Returns
        -------
        sigma_xy : float
            Anomalous Hall conductivity (e²/h in units of quantum conductance)
        """
        # Get Berry curvature
        omega_z, _ = self.calculate_berry_curvature_2d(band_idx)
        
        # Fermi-Dirac occupation factor
        if temperature < 0.1:  # T=0
            f_nk = (self.energies[:, band_idx] < self.fermi_energy).astype(float)
        else:
            # Thermal broadening
            beta = 1.0 / (8.617e-5 * temperature)  # 1/k_B T
            f_nk = 1.0 / (1.0 + np.exp(beta * (self.energies[:, band_idx] - self.fermi_energy)))
        
        # Integrate: σ_xy ∝ ∫ Ω_z(k) f(E_nk) dk
        integral = np.sum(omega_z * f_nk) * self.kspace_volume / self.nk
        
        # Convert to e²/h units
        sigma_xy = -integral / (2 * np.pi)
        
        return float(sigma_xy.real)
    
    def calculate_orbital_magnetization(self, band_idx: int,
                                       temperature: float = 0.0) -> float:
        """Calculate orbital magnetization from Berry curvature
        
        M_orb = -e/(2ℏ) ∫ Ω_z(k) f(E_nk) d³k
        
        Parameters
        ----------
        band_idx : int
            Band index
        temperature : float
            Temperature in K
            
        Returns
        -------
        magnetization : float
            Orbital magnetization (μ_B)
        """
        # Get Berry curvature
        omega_z, _ = self.calculate_berry_curvature_2d(band_idx)
        
        # Fermi-Dirac factor
        if temperature < 0.1:
            f_nk = (self.energies[:, band_idx] < self.fermi_energy).astype(float)
        else:
            beta = 1.0 / (8.617e-5 * temperature)
            f_nk = 1.0 / (1.0 + np.exp(beta * (self.energies[:, band_idx] - self.fermi_energy)))
        
        # Integrate
        m_orb = -np.sum(omega_z * f_nk) * self.kspace_volume / (2 * self.nk * np.pi)
        
        # Convert to μ_B (Bohr magneton)
        # μ_B = e ℏ / (2 m_e) - natural unit in DFT
        return float(m_orb)
    
    def calculate_fermi_surface_berry_phase(self, band_idx: int,
                                           fermi_tolerance: float = 0.05) -> float:
        """Calculate Berry phase at Fermi surface for given band
        
        For systems with Fermi surface crossing the band, calculates
        integrated Berry phase around Fermi level.
        
        Parameters
        ----------
        band_idx : int
            Band index
        fermi_tolerance : float
            Energy tolerance around Fermi level (eV)
            
        Returns
        -------
        berry_phase : float
            Berry phase (0 to π)
        """
        # Find k-points near Fermi level
        energy = self.energies[:, band_idx]
        fermi_mask = np.abs(energy - self.fermi_energy) < fermi_tolerance
        
        if np.sum(fermi_mask) == 0:
            return 0.0
        
        # Get Berry connection at Fermi level
        berry_conn = self.calculate_berry_connection(band_idx)
        
        # Integrate along Fermi surface (simplified)
        berry_phase = np.sum(np.abs(berry_conn[fermi_mask])) / np.sum(fermi_mask)
        
        return float(berry_phase % np.pi)
    
    def identify_topological_bands(self, chern_tolerance: float = 0.1) -> Dict:
        """Identify topologically non-trivial bands
        
        Parameters
        ----------
        chern_tolerance : float
            Tolerance for Chern number quantization
            
        Returns
        -------
        topology : Dict
            Dictionary with band indices and their Chern numbers
        """
        topology = {
            'trivial_bands': [],
            'nontrivial_bands': [],
            'chern_numbers': {},
        }
        
        for band in range(self.nbands):
            try:
                chern = self.calculate_chern_number(band)
                topology['chern_numbers'][band] = chern
                
                # Check if Chern number is close to integer
                chern_int = round(chern)
                
                if abs(chern - chern_int) < chern_tolerance and chern_int != 0:
                    topology['nontrivial_bands'].append(band)
                else:
                    topology['trivial_bands'].append(band)
                    
            except Exception as e:
                print(f"Warning: Could not calculate Chern for band {band}: {e}")
        
        return topology
    
    def plot_berry_curvature_map(self, band_idx: int,
                                output_file: Optional[str] = None) -> None:
        """Plot Berry curvature map in k-space
        
        Parameters
        ----------
        band_idx : int
            Band index to plot
        output_file : str, optional
            Output filename for saving
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available for plotting")
            return
        
        # Calculate Berry curvature
        omega_z, _ = self.calculate_berry_curvature_2d(band_idx)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create 2D projection if k-points have structure
        if self.kpoints.shape[1] >= 2:
            kx = self.kpoints[:, 0]
            ky = self.kpoints[:, 1]
            
            # Plot as scatter with color map
            scatter = ax.scatter(kx, ky, c=np.abs(omega_z), s=50,
                               cmap='RdYlBu', alpha=0.8, edgecolors='black')
            
            ax.set_xlabel('k_x (2π/a)', fontsize=12)
            ax.set_ylabel('k_y (2π/a)', fontsize=12)
        else:
            # 1D representation
            ax.plot(range(len(omega_z)), omega_z, 'o-', linewidth=2, markersize=6)
            ax.set_xlabel('k-point index', fontsize=12)
            ax.set_ylabel('Ω_z (k)', fontsize=12)
            scatter = None
        
        ax.set_title(f'Berry Curvature - Band {band_idx}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        if scatter:
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('|Ω_z| (arbitrary units)', fontsize=11)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Berry curvature map saved to {output_file}")
        
        plt.show()
    
    def plot_chern_number_histogram(self, output_file: Optional[str] = None) -> None:
        """Plot Chern numbers for all bands
        
        Parameters
        ----------
        output_file : str, optional
            Output filename for saving
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available for plotting")
            return
        
        chern_numbers = []
        
        for band in range(min(self.nbands, 10)):  # Calculate first 10 bands
            try:
                chern = self.calculate_chern_number(band)
                chern_numbers.append(chern)
            except:
                chern_numbers.append(0.0)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = ['red' if c != 0 else 'blue' for c in chern_numbers]
        ax.bar(range(len(chern_numbers)), chern_numbers, color=colors, alpha=0.7, edgecolor='black')
        
        # Add integer grid
        for i in range(-3, 4):
            ax.axhline(y=i, color='gray', linestyle='--', alpha=0.3, linewidth=0.5)
        
        ax.set_xlabel('Band Index', fontsize=12)
        ax.set_ylabel('Chern Number', fontsize=12)
        ax.set_title('Chern Numbers vs Band Index', fontsize=14, fontweight='bold')
        ax.set_xticks(range(len(chern_numbers)))
        ax.grid(True, axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Chern number histogram saved to {output_file}")
        
        plt.show()
    
    def print_berry_summary(self, num_bands: int = 5) -> None:
        """Print summary of Berry curvature and topological properties
        
        Parameters
        ----------
        num_bands : int
            Number of bands to analyze
        """
        print()
        print("=" * 70)
        print("Berry Curvature and Topological Properties Summary")
        print("=" * 70)
        print()
        print(f"Number of k-points:          {self.nk}")
        print(f"Number of bands:             {self.nbands}")
        print(f"Fermi energy:                {self.fermi_energy:.4f} eV")
        print()
        print("-" * 70)
        print("Band-Resolved Topological Properties:")
        print("-" * 70)
        print(f"{'Band':<8} {'Chern #':<12} {'Anom. Hall':<15} {'Orb. Mag.':<15}")
        print("-" * 70)
        
        for band in range(min(num_bands, self.nbands)):
            try:
                chern = self.calculate_chern_number(band)
                ahc = self.calculate_anomalous_hall_conductivity(band)
                m_orb = self.calculate_orbital_magnetization(band)
                
                print(f"{band:<8} {chern:>10.4f}   {ahc:>13.4e}   {m_orb:>13.4e} μ_B")
            except Exception as e:
                print(f"{band:<8} {'Error':<12} {str(e):<15}")
        
        print()
        print("=" * 70)
        print()
