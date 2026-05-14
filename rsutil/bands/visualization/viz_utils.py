#!/usr/bin/env python3
"""
Brillouin zone visualization utilities.

Provides functions to plot the 3D Brillouin zone polyhedron with
high-symmetry point labels and optional prime convention points.

Author: Ramon Cardias (integrated into rsutil)
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Dict, List, Optional
import warnings


def plot_brillouin_zone(
    reciprocal_lattice: np.ndarray,
    point_coords: Dict[str, np.ndarray],
    label_hs_points: bool = True,
    label_fontsize: int = 10,
    label_primes: bool = True,
    figsize: tuple = (8, 8),
) -> tuple:
    """Plot 3D Brillouin zone polyhedron with high-symmetry point labels.
    
    Parameters
    ----------
    reciprocal_lattice : np.ndarray
        Reciprocal lattice vectors as rows (3×3)
    point_coords : dict
        Dictionary mapping point names to fractional coordinates in 
        reciprocal primitive basis (e.g., from seekpath.get_path_orig_cell)
    label_hs_points : bool, default=True
        If True, label high-symmetry points on the BZ plot
    label_fontsize : int, default=10
        Font size for point labels
    label_primes : bool, default=True
        If True, also plot and label "primed" versions of points with
        in-plane flipping: (kx, ky, kz) → (-kx, -ky, kz)
    figsize : tuple, default=(8, 8)
        Figure size (width, height)
        
    Returns
    -------
    tuple
        (figure, axes) matplotlib objects
        
    Raises
    ------
    ImportError
        If pymatgen is not available
        
    Notes
    -----
    Requires pymatgen for BZ facet generation. Install with:
        pip install pymatgen
    
    Author
    ------
    Ramon Cardias
    """
    try:
        from pymatgen.core.lattice import Lattice
    except ImportError:
        raise ImportError(
            "pymatgen required for BZ visualization. "
            "Install with: pip install pymatgen"
        )
    
    # Helper: fractional to absolute k-point
    def frac_to_abs(kfrac):
        return np.dot(np.array(kfrac, dtype=float), reciprocal_lattice)
    
    # Helper: apply prime convention
    def prime_k(k_abs):
        k = np.array(k_abs, dtype=float).copy()
        k[0] *= -1.0
        k[1] *= -1.0
        return k
    
    # Create figure and 3D axis
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Get BZ from pymatgen
    try:
        bz_facets = Lattice(reciprocal_lattice).get_brillouin_zone()
    except Exception as e:
        warnings.warn(f"Could not compute BZ facets: {e}")
        bz_facets = []
    
    # Plot BZ edges
    allv = []
    for facet in bz_facets:
        facet = np.array(facet, dtype=float)
        allv.append(facet)
        
        # Draw facet edges
        f2 = np.vstack([facet, facet[0]])
        ax.plot(f2[:, 0], f2[:, 1], f2[:, 2], lw=1.0, color='black', alpha=0.3)
    
    # Compute view bounds from all vertices
    if allv:
        allv = np.vstack(allv)
        vmin = allv.min(axis=0)
        vmax = allv.max(axis=0)
        
        # Add padding and enforce cubic view
        center = 0.5 * (vmin + vmax)
        halfspan = 0.5 * (vmax - vmin)
        span = float(np.max(halfspan))
        pad = 0.12 * span
        span = span + pad
        
        ax.set_xlim(center[0] - span, center[0] + span)
        ax.set_ylim(center[1] - span, center[1] + span)
        ax.set_zlim(center[2] - span, center[2] + span)
    
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass  # Older matplotlib versions don't have this
    
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_zlabel('kz')
    ax.set_title('Brillouin Zone', fontsize=14)
    
    # Plot and label high-symmetry points
    if label_hs_points and point_coords:
        for name, kf in point_coords.items():
            ka = frac_to_abs(kf)
            
            # Regular point
            ax.scatter([ka[0]], [ka[1]], [ka[2]], s=20, c='red', zorder=5)
            
            # Format label (Γ instead of GAMMA)
            txt = "Γ" if name.upper() == "GAMMA" else name
            ax.text(ka[0], ka[1], ka[2], txt, fontsize=label_fontsize)
            
            # Primed point if requested
            if label_primes:
                kap = prime_k(ka)
                ax.scatter([kap[0]], [kap[1]], [kap[2]], s=20, c='blue', zorder=5)
                ax.text(kap[0], kap[1], kap[2], txt + "′", fontsize=label_fontsize)
    
    plt.tight_layout()
    
    return fig, ax


def plot_bands_with_projections(
    kpath_x: np.ndarray,
    band_energies: np.ndarray,
    projections: np.ndarray,
    xticks: List[float],
    xticklabels: List[str],
    fermi_energy: float = 0.0,
    projection_mode: str = "expectation",
    marker_sizes: Optional[np.ndarray] = None,
    ylim: Optional[tuple] = None,
    cmap: str = "coolwarm",
    title: str = "",
) -> tuple:
    """Plot band structure with colored scatter by projection.
    
    Parameters
    ----------
    kpath_x : np.ndarray
        Linear coordinate along k-path (nk,)
    band_energies : np.ndarray
        Band energies (nk, nbands)
    projections : np.ndarray
        Projection values for coloring (nk, nbands)
    xticks : list of float
        x-coordinates for major tick marks
    xticklabels : list of str
        TeX-formatted labels for ticks
    fermi_energy : float, default=0.0
        Energy to subtract from bands (typically 0 for Fermi level)
    projection_mode : str, default="expectation"
        Projection mode for value range interpretation
    marker_sizes : np.ndarray, optional
        Marker sizes per point (nk, nbands). If None, uses fixed size.
    ylim : tuple, optional
        Energy axis limits (emin, emax)
    cmap : str, default="coolwarm"
        Matplotlib colormap name
    title : str, default=""
        Plot title
        
    Returns
    -------
    tuple
        (figure, axes)
        
    Author
    ------
    Ramon Cardias
    """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    
    nk, nbands = band_energies.shape
    
    # Plot background lines for all bands
    for b in range(nbands):
        ax.plot(kpath_x, band_energies[:, b] - fermi_energy, 
                lw=0.7, color='0.6', alpha=0.8, zorder=1)
    
    # Prepare data for scatter
    X = np.repeat(kpath_x, nbands)
    Y = (band_energies - fermi_energy).ravel()
    C = projections.ravel()
    
    # Marker sizes
    if marker_sizes is not None:
        S = marker_sizes.ravel()
    else:
        S = 20  # Fixed size
    
    # Determine color range based on projection mode
    from .projection_utils import get_projection_vrange
    vmin, vmax = get_projection_vrange(projection_mode)
    
    # Plot colored scatter
    scatter = ax.scatter(
        X, Y, s=S, c=C, cmap=cmap,
        vmin=vmin, vmax=vmax,
        linewidths=0, alpha=0.85, zorder=2
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(f"Projection ({projection_mode})", fontsize=11)
    
    # Add symmetry point lines and ticks
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    
    for xp in xticks:
        ax.axvline(x=xp, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_xlim(kpath_x[0], kpath_x[-1])
    
    if ylim is not None:
        ax.set_ylim(*ylim)
    
    if title:
        ax.set_title(title, fontsize=12)
    
    plt.tight_layout()
    
    return fig, ax


def export_kpath(
    filename: str,
    kpath_x: np.ndarray,
    kpoints_abs: np.ndarray,
    kpoint_labels: Optional[List[str]] = None,
    header: str = "x_linearcoord   kx   ky   kz",
) -> None:
    """Export k-path to file for reproducibility.
    
    Parameters
    ----------
    filename : str
        Output file path
    kpath_x : np.ndarray
        Linear coordinates along path (nk,)
    kpoints_abs : np.ndarray
        Absolute (Cartesian) k-point coordinates (nk, 3)
    kpoint_labels : list, optional
        K-point labels (nk,). If provided, added as comment column
    header : str, default="x_linearcoord   kx   ky   kz"
        Header line written to file
        
    Notes
    -----
    Output format: Plain text with columns [x, kx, ky, kz] and optional label comments
    
    Author
    ------
    Ramon Cardias
    """
    if len(kpath_x) != len(kpoints_abs):
        raise ValueError("kpath_x and kpoints_abs must have same length")
    
    arr = np.column_stack([kpath_x, kpoints_abs])
    
    # Save with header
    np.savetxt(
        filename,
        arr,
        fmt="%.12e",
        header=header
    )


__all__ = [
    'plot_brillouin_zone',
    'plot_bands_with_projections',
    'export_kpath',
]
