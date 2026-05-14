"""
Seekpath integration for automatic k-path generation.

Provides interface to seekpath library for generating high-symmetry k-paths
with automatic point spacing and interpolation.

Author: Ramon Cardias
"""

import numpy as np
from typing import Dict, List
import warnings

from .coordinates import get_reciprocal_lattice, fractional_to_cartesian
from .labels import tex_format_label
from .point_allocation import allocate_points_per_segment

try:
    import seekpath
    HAS_SEEKPATH = True
except ImportError:
    HAS_SEEKPATH = False
    warnings.warn("seekpath not found. Install with: pip install seekpath")


def seekpath_kpath(
    lattice_vectors: np.ndarray,
    atomic_positions: np.ndarray,
    atomic_numbers: List[int],
    num_points: int = 100,
    reference_distance: float = 0.025,
    symprec: float = 1e-5,
    angle_tol: float = -1.0,
) -> Dict:
    """Generate explicit k-path using seekpath with resampling to requested num_points.
    
    Returns k-points in **fractional coordinates** (same convention as structure_analysis.py).
    This ensures consistency with the working band calculator.
    
    Parameters
    ----------
    lattice_vectors : np.ndarray
        Direct lattice vectors as rows (3×3)
    atomic_positions : np.ndarray
        Atomic positions in fractional coordinates
    atomic_numbers : list of int
        Atomic numbers for each site
    num_points : int, default=100
        Desired total number of k-points (path will be resampled to this)
    reference_distance : float, default=0.025
        Reference distance for seekpath Brillouin zone sampling
    symprec : float, default=1e-5
        Symmetry precision tolerance
    angle_tol : float, default=-1.0
        Angle tolerance (use default if negative)
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'kpoints_frac': Fractional k-point coordinates (nk, 3) [FRACTIONAL, not absolute!]
        - 'kpoints_abs': Absolute (Cartesian) k-point coordinates (nk, 3) [for reference only]
        - 'klabels': K-point labels (list of str)
        - 'xvals': Linear coordinate along k-path in Cartesian (nk,)
        - 'segments': List of (start_idx, end_idx) for each path segment
        - 'xticks': x-coordinates for major labels
        - 'xticklabels': TeX-formatted labels at major points
        - 'point_coords': Dict of all high-symmetry point fractional coordinates
        - 'reciprocal_lattice': Reciprocal lattice vectors
        
    Author
    ------
    Ramon Cardias
    """
    if not HAS_SEEKPATH:
        raise ImportError("seekpath required. Install with: pip install seekpath")
    
    cell_tuple = (lattice_vectors, atomic_positions.tolist(), atomic_numbers)
    
    # Get path information (point_coords are FRACTIONAL)
    path_info = seekpath.get_path_orig_cell(
        cell_tuple,
        with_time_reversal=True,
        recipe="hpkot",
        symprec=symprec,
        angle_tolerance=angle_tol
    )
    
    # Get explicit k-points (with proper spacing)
    explicit_info = seekpath.get_explicit_k_path_orig_cell(
        cell_tuple,
        with_time_reversal=True,
        reference_distance=reference_distance,
        recipe="hpkot",
        symprec=symprec,
        angle_tolerance=angle_tol
    )
    
    point_coords = dict(path_info["point_coords"])  # FRACTIONAL coordinates
    
    # Compute reciprocal lattice using input lattice (same as structure_analysis.py)
    recip_lat = get_reciprocal_lattice(lattice_vectors)
    
    # Build k-path from explicit seekpath points
    kpoints_frac = []
    klabels = []
    
    # Use explicit k-points for interpolation (more precise spacing)
    kpts_abs_raw = np.array(explicit_info["explicit_kpoints_abs"], dtype=float)
    
    # Convert back to fractional to match working convention
    # inverse operation: k_frac = k_cart @ inv(recip_lat)
    recip_inv = np.linalg.inv(recip_lat)
    kpts_frac_raw = kpts_abs_raw @ recip_inv
    
    klabels_raw = list(explicit_info["explicit_kpoints_labels"])
    
    # Include ALL points from seekpath (including interpolated ones with empty labels)
    # Only skip repeated junction points
    prev_kpt = None
    for i, (kf, label) in enumerate(zip(kpts_frac_raw, klabels_raw)):
        # Skip if exactly same k-point as previous (true duplicate junction)
        if prev_kpt is not None and np.allclose(kf, prev_kpt, atol=1e-10):
            continue
        kpoints_frac.append(kf)
        klabels.append(label)
        prev_kpt = kf.copy()
    
    kpoints_frac = np.array(kpoints_frac, dtype=float)
    
    # Now compute kdistances in Cartesian (for plotting x-axis)
    kpoints_abs = kpoints_frac @ recip_lat
    kdistances = np.zeros(len(kpoints_abs))
    for i in range(1, len(kpoints_abs)):
        dk = np.linalg.norm(kpoints_abs[i] - kpoints_abs[i-1])
        kdistances[i] = kdistances[i-1] + dk
    
    # Resample k-path to requested number of points while preserving high-symmetry points
    if len(kpoints_frac) > num_points:
        # Find indices of high-symmetry points (those with non-empty labels)
        special_indices = [i for i, label in enumerate(klabels) if label]
        
        # Resample: create uniform spacing and preserve special points
        total_distance = kdistances[-1]
        new_kdistances = np.linspace(0, total_distance, num_points)
        
        # Use cubic spline interpolation for smooth resampling
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(kdistances, kpoints_frac, axis=0)
        kpoints_frac_resampled = cs(new_kdistances)
        
        # Recalculate k-point labels for resampled path
        # Keep labels only at positions close to original special points
        klabels_resampled = [""] * num_points
        for spec_idx in special_indices:
            original_dist = kdistances[spec_idx]
            # Find closest point in resampled path
            closest_idx = np.argmin(np.abs(new_kdistances - original_dist))
            klabels_resampled[closest_idx] = klabels[spec_idx]
        
        # Update arrays
        kpoints_frac = kpoints_frac_resampled
        klabels = klabels_resampled
        kdistances = new_kdistances
        kpoints_abs = kpoints_frac @ recip_lat
    
    # Generate x-axis ticks and labels using label positions
    ticks_dict = {}
    for i, label in enumerate(klabels):
        if label:
            xval = round(float(kdistances[i]), 10)
            if xval in ticks_dict:
                ticks_dict[xval] += "|" + label
            else:
                ticks_dict[xval] = label
    
    xticks = sorted(list(ticks_dict.keys()))
    xticklabels = [tex_format_label(ticks_dict[x].split("|")[0]) for x in xticks]
    
    return {
        'kpoints_frac': kpoints_frac,  # FRACTIONAL (main return)
        'kpoints_abs': kpoints_abs,    # Cartesian (for reference)
        'klabels': klabels,
        'xvals': kdistances,
        'segments': list(explicit_info["explicit_segments"]),
        'xticks': xticks,
        'xticklabels': xticklabels,
        'point_coords': point_coords,
        'reciprocal_lattice': recip_lat,
    }
