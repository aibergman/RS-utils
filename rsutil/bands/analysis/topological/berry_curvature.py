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
        """Plot Berry curvature map in k-space with publication-quality styling
        
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
            
            # Plot as scatter with improved color map
            scatter = ax.scatter(kx, ky, c=omega_z, s=60,
                               cmap='coolwarm', alpha=0.8, edgecolors='black', linewidth=0.5)
            
            ax.set_xlabel('k$_x$ (2$\\pi$/a)', fontsize=12)
            ax.set_ylabel('k$_y$ (2$\\pi$/a)', fontsize=12)
            ax.set_aspect('equal', adjustable='datalim')
        else:
            # 1D representation with improved styling
            ax.plot(range(len(omega_z)), omega_z, 'o-', color='#2E86AB',
                   linewidth=2, markersize=7, markerfacecolor='white', 
                   markeredgewidth=1.5, markeredgecolor='#2E86AB')
            ax.set_xlabel('k-point index', fontsize=12)
            ax.set_ylabel('$\\Omega_z$ (k)', fontsize=12)
            scatter = None
            ax.axhline(0, color='gray', linestyle='-', linewidth=0.8, alpha=0.5, zorder=0)
        
        ax.set_title(f'Berry Curvature - Band {band_idx}', fontsize=14)
        ax.grid(False)
        ax.set_axisbelow(True)
        
        if scatter:
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('$\\Omega_z$ (arbitrary units)', fontsize=11)
        
        fig.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Berry curvature map saved to {output_file}")
            plt.close()
        else:
            plt.show()
    
    def plot_chern_number_histogram(self, output_file: Optional[str] = None) -> None:
        """Plot Chern numbers for all bands with publication-quality styling
        
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
        
        # Color code: non-trivial (non-zero) in red, trivial (zero) in blue
        colors = ['#E63946' if abs(c) > 0.1 else '#457B9D' for c in chern_numbers]
        bars = ax.bar(range(len(chern_numbers)), chern_numbers, 
                 color=colors, alpha=0.8, edgecolor='black', linewidth=1.2)
        
        # Add integer grid lines
        for i in range(-3, 4):
            ax.axhline(y=i, color='gray', linestyle='-', alpha=0.3, linewidth=0.5, zorder=0)
        
        ax.set_xlabel('Band Index', fontsize=12)
        ax.set_ylabel('Chern Number', fontsize=12)
        ax.set_title('Chern Numbers vs Band Index', fontsize=14)
        ax.set_xticks(range(len(chern_numbers)))
        ax.grid(False)
        ax.set_axisbelow(True)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E63946', edgecolor='black', label='Non-trivial (C ≠ 0)'),
            Patch(facecolor='#457B9D', edgecolor='black', label='Trivial (C = 0)')
        ]
        ax.legend(handles=legend_elements, loc='best', frameon=False)
        
        fig.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Chern number histogram saved to {output_file}")
            plt.close()
        else:
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
    
    def calculate_orbital_resolved_berry_curvature(self, band_idx: int,
                                                   n_orb_per_spin: int = 9) -> Dict:
        """
        Decompose Berry curvature by orbital character.
        
        For each orbital α, calculate the contribution:
            Ω^α(k) = ∑_n |⟨α|ψ_n(k)⟩|² Ω_n(k)
        
        Parameters
        ----------
        band_idx : int
            Band index to analyze
        n_orb_per_spin : int
            Number of orbitals per spin channel (default: 9 for TB-LMTO)
            
        Returns
        -------
        orbital_berry : dict
            Dictionary containing:
            - 'total_curvature': Total Berry curvature for the band
            - 'orbital_contributions': Berry curvature by orbital (n_orb,)
            - 'spin_contributions': Berry curvature by spin channel (2,)
            - 'dominant_orbital': Index of dominant orbital
            - 'orbital_weights': Orbital character weights
            
        Notes
        -----
        For non-collinear calculations with 18 orbitals:
        - Orbitals 0-8: spin-up (s, px, py, pz, dxy, dyz, dz2, dxz, dx2-y2)
        - Orbitals 9-17: spin-down (same sequence)
        """
        # Calculate Berry curvature for this band
        try:
            berry_curv_z = self.calculate_berry_curvature_kubo(band_idx)
        except:
            # Fallback to 2D calculation
            berry_curv_z, _ = self.calculate_berry_curvature_2d(band_idx)
        
        # Extract orbital contributions from eigenvectors
        total_orbitals = self.eigenvectors.shape[1]
        orbital_berry = np.zeros(total_orbitals)
        
        # For each k-point, weight Berry curvature by orbital character
        for ik in range(self.nk):
            psi = self.eigenvectors[ik, :, band_idx]
            orbital_weights = np.abs(psi)**2
            
            # Distribute Berry curvature according to orbital character
            orbital_berry += orbital_weights * berry_curv_z[ik]
        
        # Normalize by number of k-points
        orbital_berry /= self.nk
        
        # Separate by spin if applicable
        if total_orbitals >= 2 * n_orb_per_spin:
            spin_up_contrib = np.sum(orbital_berry[:n_orb_per_spin])
            spin_down_contrib = np.sum(orbital_berry[n_orb_per_spin:2*n_orb_per_spin])
            spin_contributions = np.array([spin_up_contrib, spin_down_contrib])
        else:
            spin_contributions = np.array([np.sum(orbital_berry), 0.0])
        
        # Identify dominant orbital
        dominant_orbital = np.argmax(np.abs(orbital_berry))
        
        # Calculate total
        total_curvature = np.sum(berry_curv_z) / self.nk
        
        return {
            'total_curvature': total_curvature,
            'orbital_contributions': orbital_berry,
            'spin_contributions': spin_contributions,
            'dominant_orbital': dominant_orbital,
            'dominant_contribution': orbital_berry[dominant_orbital],
            'orbital_labels': self._get_orbital_labels(n_orb_per_spin)
        }
    
    def _get_orbital_labels(self, n_orb_per_spin: int = 9) -> List[str]:
        """Get orbital labels for TB-LMTO basis"""
        if n_orb_per_spin == 9:
            orb_names = ['s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
            labels = [f"{orb}↑" for orb in orb_names] + [f"{orb}↓" for orb in orb_names]
        else:
            labels = [f"orb{i}↑" for i in range(n_orb_per_spin)] + \
                    [f"orb{i}↓" for i in range(n_orb_per_spin)]
        
        return labels
    
    def plot_orbital_resolved_berry_curvature(self, band_idx: int,
                                             output_file: Optional[str] = None):
        """
        Plot orbital-resolved Berry curvature contributions as bar chart.
        
        Parameters
        ----------
        band_idx : int
            Band index to analyze
        output_file : str, optional
            Save plot to file if provided
        """
        import matplotlib.pyplot as plt
        
        # Calculate orbital-resolved Berry curvature
        orb_data = self.calculate_orbital_resolved_berry_curvature(band_idx)
        
        orb_contrib = orb_data['orbital_contributions']
        orb_labels = orb_data['orbital_labels']
        n_orb = len(orb_contrib)
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Plot 1: All orbital contributions
        colors = ['blue'] * (n_orb // 2) + ['red'] * (n_orb // 2)
        ax1.bar(range(n_orb), orb_contrib, color=colors, alpha=0.7, edgecolor='black')
        ax1.set_xlabel('Orbital Index')
        ax1.set_ylabel('Berry Curvature Contribution (a.u.)')
        ax1.set_title(f'Orbital-Resolved Berry Curvature - Band {band_idx}')
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        
        # Add vertical separator between spin channels
        if n_orb >= 18:
            ax1.axvline(x=n_orb//2 - 0.5, color='gray', linestyle='--', linewidth=2)
            ax1.text(n_orb//4, ax1.get_ylim()[1]*0.9, 'Spin ↑', ha='center', fontsize=12, weight='bold')
            ax1.text(3*n_orb//4, ax1.get_ylim()[1]*0.9, 'Spin ↓', ha='center', fontsize=12, weight='bold')
        
        # Plot 2: Spin-resolved contributions
        spin_contrib = orb_data['spin_contributions']
        ax2.bar(['Spin ↑', 'Spin ↓'], spin_contrib, color=['blue', 'red'], alpha=0.7, edgecolor='black')
        ax2.set_ylabel('Total Berry Curvature (a.u.)')
        ax2.set_title('Spin-Resolved Berry Curvature')
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        
        # Add summary text
        total = orb_data['total_curvature']
        dominant = orb_data['dominant_orbital']
        dom_label = orb_labels[dominant] if dominant < len(orb_labels) else f"Orb {dominant}"
        
        summary = f"Total: {total:.4e} a.u.\n"
        summary += f"Dominant: {dom_label} ({orb_data['dominant_contribution']:.4e})"
        ax2.text(0.98, 0.97, summary, transform=ax2.transAxes,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                fontsize=10, family='monospace')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved orbital-resolved Berry curvature plot to {output_file}")
        else:
            plt.show()
        
        plt.close()

