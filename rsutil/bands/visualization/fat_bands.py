"""
Fat Bands Visualization - Orbital/Spin-Projected Band Structure

This module provides tools for visualizing band structures with orbital
and spin projections using "fat bands" where the line thickness or color
intensity represents the projection weight.

Phase 4 Implementation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from typing import List, Optional, Dict, Tuple


class FatBandsPlotter:
    """
    Create fat bands plots showing orbital/spin character of bands.
    
    Fat bands visualize the contribution of specific orbitals or atomic sites
    to each band by varying the line width or color intensity.
    
    Parameters
    ----------
    band_energies : np.ndarray
        Band energies (nk × nbands) in eV
    eigenvectors : np.ndarray
        Eigenvectors (nk × nbands × nbands) for projection
    kdistances : np.ndarray
        K-path distances for x-axis (nk,)
    klabels : list of str
        High-symmetry point labels
    fermi_energy : float, optional
        Fermi energy for reference line
    """
    
    def __init__(self, band_energies: np.ndarray, eigenvectors: np.ndarray,
                 kdistances: np.ndarray, klabels: List[str],
                 fermi_energy: Optional[float] = None):
        self.energies = band_energies
        self.eigenvectors = eigenvectors
        self.kdistances = kdistances
        self.klabels = klabels
        self.fermi_energy = fermi_energy
        
        self.nk, self.nbands = band_energies.shape
        self.n_orb_per_spin = self.nbands // 2  # Assumes non-collinear
    
    def calculate_orbital_projections(self, orbital_indices: List[int]) -> np.ndarray:
        """
        Calculate projection weights for specified orbitals.
        
        For each band n at k-point k:
            P^α_n(k) = ∑_{i∈α} |⟨i|ψ_n(k)⟩|²
        
        Parameters
        ----------
        orbital_indices : list of int
            Indices of orbitals to project onto
            
        Returns
        -------
        projections : np.ndarray
            Projection weights (nk × nbands)
        """
        projections = np.zeros((self.nk, self.nbands))
        
        for ik in range(self.nk):
            for n in range(self.nbands):
                psi = self.eigenvectors[ik, :, n]
                # Sum |coefficient|² for specified orbitals
                projections[ik, n] = np.sum(np.abs(psi[orbital_indices])**2)
        
        return projections
    
    def calculate_spin_projections(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate spin-up and spin-down projections.
        
        Returns
        -------
        proj_up : np.ndarray
            Spin-up projections (nk × nbands)
        proj_down : np.ndarray
            Spin-down projections (nk × nbands)
        """
        proj_up = np.zeros((self.nk, self.nbands))
        proj_down = np.zeros((self.nk, self.nbands))
        
        for ik in range(self.nk):
            for n in range(self.nbands):
                psi = self.eigenvectors[ik, :, n]
                
                # First half: spin-up, second half: spin-down
                proj_up[ik, n] = np.sum(np.abs(psi[:self.n_orb_per_spin])**2)
                proj_down[ik, n] = np.sum(np.abs(psi[self.n_orb_per_spin:])**2)
        
        return proj_up, proj_down
    
    def calculate_atom_projections(self, atom_idx: int, 
                                   n_orb_per_atom: int = 9) -> np.ndarray:
        """
        Calculate projection onto specific atom.
        
        Parameters
        ----------
        atom_idx : int
            Atom index (0-based)
        n_orb_per_atom : int
            Number of orbitals per atom (default: 9 for spd)
            
        Returns
        -------
        projections : np.ndarray
            Atom projection weights (nk × nbands)
        """
        # For non-collinear: atom orbitals are in blocks
        # Atom i has orbitals [i*9:(i+1)*9] for spin-up
        # and [n_atoms*9 + i*9:n_atoms*9 + (i+1)*9] for spin-down
        
        n_atoms = self.n_orb_per_spin // n_orb_per_atom
        
        # Spin-up block
        orb_up = list(range(atom_idx * n_orb_per_atom, 
                           (atom_idx + 1) * n_orb_per_atom))
        # Spin-down block
        orb_down = list(range(self.n_orb_per_spin + atom_idx * n_orb_per_atom,
                             self.n_orb_per_spin + (atom_idx + 1) * n_orb_per_atom))
        
        orbital_indices = orb_up + orb_down
        
        return self.calculate_orbital_projections(orbital_indices)
    
    def plot_fat_bands_linewidth(self, projections: np.ndarray,
                                 scale_factor: float = 5.0,
                                 energy_range: Optional[Tuple[float, float]] = None,
                                 output_file: Optional[str] = None,
                                 title: str = "Fat Bands"):
        """
        Plot fat bands using line width to show projection weight.
        
        Parameters
        ----------
        projections : np.ndarray
            Projection weights (nk × nbands), values in [0, 1]
        scale_factor : float
            Scale factor for line widths (larger = thicker lines)
        energy_range : tuple, optional
            (E_min, E_max) for y-axis
        output_file : str, optional
            Save plot to file
        title : str
            Plot title
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot each band
        for n in range(self.nbands):
            # Get projection weights for this band
            weights = projections[:, n]
            
            # Variable linewidth based on projection
            # Use segments to vary linewidth along the band
            points = np.array([self.kdistances, self.energies[:, n]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Average weights for each segment
            segment_weights = (weights[:-1] + weights[1:]) / 2
            
            # Create LineCollection with variable linewidths
            lwidths = scale_factor * segment_weights
            lc = LineCollection(segments, linewidths=lwidths, colors='blue', alpha=0.7)
            ax.add_collection(lc)
        
        # Add Fermi level
        if self.fermi_energy is not None:
            ax.axhline(self.fermi_energy, color='red', linestyle='--', 
                      linewidth=1.5, label='$E_F$', alpha=0.7)
        
        # Set up k-point labels
        self._setup_kpoint_labels(ax)
        
        # Defer energy range application to after all artists are added
        if energy_range:
            target_ylim = tuple(energy_range)
        else:
            
            e_min = float(np.min(self.energies))
            e_max = float(np.max(self.energies))
            e_range = e_max - e_min
            e_pad = 0.02 * e_range
            target_ylim = (e_min - e_pad, e_max + e_pad)
        
        ax.set_xlabel('k-path', fontsize=12)
        ax.set_ylabel('Energy (eV)', fontsize=12)
        ax.set_title(title, fontsize=14)
        
        # Remove all grid and minor ticks
        ax.grid(False)
        ax.minorticks_off()
        ax.set_axisbelow(True)
        
        if self.fermi_energy is not None:
            ax.legend()

        # Apply y-limits last to avoid unintended autoscale changes
        ax.set_ylim(*target_ylim)
        ax.autoscale(False)
        
        if output_file:
            fig.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved fat bands plot to {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def plot_fat_bands_color(self, projections: np.ndarray,
                            cmap: str = 'hot',
                            energy_range: Optional[Tuple[float, float]] = None,
                            output_file: Optional[str] = None,
                            title: str = "Fat Bands - Color Scale"):
        """
        Plot fat bands using color intensity to show projection weight.
        
        Parameters
        ----------
        projections : np.ndarray
            Projection weights (nk × nbands)
        cmap : str
            Matplotlib colormap name
        energy_range : tuple, optional
            (E_min, E_max) for y-axis
        output_file : str, optional
            Save plot to file
        title : str
            Plot title
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Normalize projection weights for colormap
        norm = Normalize(vmin=0, vmax=1)
        cmap_obj = plt.get_cmap(cmap)
        
        # Plot each band
        for n in range(self.nbands):
            weights = projections[:, n]
            
            # Create segments with color based on projection
            points = np.array([self.kdistances, self.energies[:, n]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Average weights for each segment
            segment_weights = (weights[:-1] + weights[1:]) / 2
            
            # Create LineCollection with variable colors
            lc = LineCollection(segments, cmap=cmap_obj, norm=norm, linewidths=2)
            lc.set_array(segment_weights)
            ax.add_collection(lc)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Projection Weight', fontsize=11)
        
        # Add Fermi level
        if self.fermi_energy is not None:
            ax.axhline(self.fermi_energy, color='green', linestyle='--', 
                      linewidth=1.5, label='$E_F$', alpha=0.8)
        
        # Set up k-point labels
        self._setup_kpoint_labels(ax)
        
        # Defer energy range application to after all artists are added
        if energy_range:
            target_ylim = tuple(energy_range)
        else:
            e_min = float(np.min(self.energies))
            e_max = float(np.max(self.energies))
            e_range = e_max - e_min
            e_pad = 0.02 * e_range
            target_ylim = (e_min - e_pad, e_max + e_pad)
        
        ax.set_xlabel('k-path', fontsize=12)
        ax.set_ylabel('Energy (eV)', fontsize=12)
        ax.set_title(title, fontsize=14)
        
        # Remove all grid and minor ticks
        ax.grid(False)
        ax.minorticks_off()
        ax.set_axisbelow(True)
        
        if self.fermi_energy is not None:
            ax.legend()

        # Apply y-limits last to avoid unintended autoscale changes
        ax.set_ylim(*target_ylim)
        ax.autoscale(False)
        
        if output_file:
            fig.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved fat bands plot to {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def plot_spin_resolved_fat_bands(self, energy_range: Optional[Tuple[float, float]] = None,
                                    output_file: Optional[str] = None):
        """
        Plot fat bands showing spin-up (blue) and spin-down (red) character.
        
        Parameters
        ----------
        energy_range : tuple, optional
            (E_min, E_max) for y-axis
        output_file : str, optional
            Save plot to file
        """
        proj_up, proj_down = self.calculate_spin_projections()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot each band with spin-resolved coloring
        for n in range(self.nbands):
            weights_up = proj_up[:, n]
            weights_down = proj_down[:, n]
            
            # Create segments
            points = np.array([self.kdistances, self.energies[:, n]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Determine dominant spin for each segment
            segment_up = (weights_up[:-1] + weights_up[1:]) / 2
            segment_down = (weights_down[:-1] + weights_down[1:]) / 2
            
            # Color by dominant spin, width by polarization
            spin_diff = segment_up - segment_down
            colors = ['blue' if diff > 0 else 'red' for diff in spin_diff]
            widths = 2 * np.abs(spin_diff)  # Larger width = more polarized
            
            # Plot segments
            for i, seg in enumerate(segments):
                ax.plot(seg[:, 0], seg[:, 1], color=colors[i], 
                       linewidth=widths[i], alpha=0.7)
        
        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='blue', linewidth=2, label='Spin ↑'),
            Line2D([0], [0], color='red', linewidth=2, label='Spin ↓')
        ]
        if self.fermi_energy is not None:
            ax.axhline(self.fermi_energy, color='green', linestyle='--',
                      linewidth=1.5, alpha=0.8)
            legend_elements.append(
                Line2D([0], [0], color='green', linestyle='--', linewidth=1.5, label='$E_F$')
            )
        ax.legend(handles=legend_elements, loc='best')
        
        # Set up k-point labels
        self._setup_kpoint_labels(ax)
        
        # Defer energy range application to after all artists are added
        if energy_range:
            target_ylim = tuple(energy_range)
        else:
            e_min = float(np.min(self.energies))
            e_max = float(np.max(self.energies))
            e_range = e_max - e_min
            e_pad = 0.02 * e_range
            target_ylim = (e_min - e_pad, e_max + e_pad)
        
        ax.set_xlabel('k-path', fontsize=12)
        ax.set_ylabel('Energy (eV)', fontsize=12)
        ax.set_title('Spin-Resolved Fat Bands', fontsize=14)
        
        # Remove all grid and minor ticks
        ax.grid(False)
        ax.minorticks_off()
        ax.set_axisbelow(True)
        
        # Apply y-limits last to avoid unintended autoscale changes
        ax.set_ylim(*target_ylim)
        ax.autoscale(False)

        if output_file:
            fig.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved spin-resolved fat bands to {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def plot_orbital_character(self, orbital_groups: Dict[str, List[int]],
                              energy_range: Optional[Tuple[float, float]] = None,
                              output_file: Optional[str] = None):
        """
        Plot multiple fat bands showing different orbital characters.
        
        Parameters
        ----------
        orbital_groups : dict
            Dictionary mapping labels to orbital index lists
            Example: {'s': [0, 9], 'p': [1,2,3,10,11,12], 'd': [4,5,6,7,8,13,14,15,16,17]}
        energy_range : tuple, optional
            (E_min, E_max) for y-axis
        output_file : str, optional
            Save plot to file
        """
        n_groups = len(orbital_groups)
        fig, axes = plt.subplots(1, n_groups, figsize=(6*n_groups, 8), sharey=True)
        
        if n_groups == 1:
            axes = [axes]
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        for ax, (label, orb_indices), color in zip(axes, orbital_groups.items(), colors):
            projections = self.calculate_orbital_projections(orb_indices)
            
            # Plot bands with projection weights
            for n in range(self.nbands):
                weights = projections[:, n]
                
                points = np.array([self.kdistances, self.energies[:, n]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                
                segment_weights = (weights[:-1] + weights[1:]) / 2
                lwidths = 5 * segment_weights
                
                lc = LineCollection(segments, linewidths=lwidths, colors=color, alpha=0.7)
                ax.add_collection(lc)
            
            # Add Fermi level
            if self.fermi_energy is not None:
                ax.axhline(self.fermi_energy, color='black', linestyle='--',
                          linewidth=1.5, alpha=0.7)
            
            # Set up k-point labels
            self._setup_kpoint_labels(ax)
            
            # Defer per-axis ylim; apply once after all subplots (sharey)
            
            ax.set_xlabel('k-path', fontsize=12)
            if ax == axes[0]:
                ax.set_ylabel('Energy (eV)', fontsize=12)
            ax.set_title(f'{label} character', fontsize=14)
            
            # Remove all grid and minor ticks
            ax.grid(False)
            ax.minorticks_off()
            ax.set_axisbelow(True)
        
        # Compute target ylim once and apply to the shared y-axis
        if energy_range:
            target_ylim = tuple(energy_range)
        else:
            e_min = float(np.min(self.energies))
            e_max = float(np.max(self.energies))
            e_range = e_max - e_min
            e_pad = 0.02 * e_range
            target_ylim = (e_min - e_pad, e_max + e_pad)

        # Apply to the first axis; sharey propagates to others
        axes[0].set_ylim(*target_ylim)
        for ax in axes:
            ax.autoscale(False)

        plt.tight_layout()
        
        if output_file:
            fig.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved orbital character plot to {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def _setup_kpoint_labels(self, ax):
        """Set up high-symmetry labels using actual kdistances indices.

        Uses positions where labels are non-empty in `self.klabels` and aligns
        major ticks and thin vertical lines to those distances. Minor ticks are
        disabled explicitly to avoid style interference.
        """
        tick_positions = []
        tick_labels = []

        if self.klabels and len(self.klabels) == len(self.kdistances):
            for i, label in enumerate(self.klabels):
                if label:
                    tick_positions.append(self.kdistances[i])
                    tick_labels.append(label)

        if tick_positions:
            ax.xaxis.set_major_locator(FixedLocator(tick_positions))
            ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))
            ax.xaxis.set_minor_locator(NullLocator())

            for pos in tick_positions[1:-1]:
                ax.axvline(pos, color='black', linestyle='-', linewidth=0.5, alpha=0.5, zorder=0)

        ax.set_xlim(self.kdistances[0], self.kdistances[-1])
        ax.margins(x=0)
