#!/usr/bin/env python3
"""
Advanced k-path generation utilities with seekpath integration.

This module provides:
- Automatic k-path generation via seekpath (explicit path with direct points)
- Interactive k-path selection (standard or custom)
- Prime convention support (K-point flipping for symmetry analysis)
- Proportional k-point allocation by segment length
- K-path export and metadata management

Author: Ramon Cardias (seekpath features integrated into rsutil)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import warnings

try:
    import seekpath
    HAS_SEEKPATH = True
except ImportError:
    HAS_SEEKPATH = False
    warnings.warn("seekpath not found. Install with: pip install seekpath")


def get_reciprocal_lattice(direct_lattice: np.ndarray) -> np.ndarray:
    """Compute reciprocal lattice vectors from direct lattice.
    
    Uses formula: b = 2π * inv(a) where a is direct lattice as rows.
    This matches the working convention in structure_analysis.py line 522
    
    Parameters
    ----------
    direct_lattice : np.ndarray
        Direct lattice vectors as rows (3×3)
        
    Returns
    -------
    np.ndarray
        Reciprocal lattice vectors as rows (3×3)
        
    Author
    ------
    Ramon Cardias
    """
    # Standard formula used in structure_analysis.generate_kpath_points:
    # b = 2π * inv(a) where a is direct lattice as rows
    return 2 * np.pi * np.linalg.inv(direct_lattice)


def cartesian_to_fractional(cartesian: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Convert Cartesian coordinates to fractional.
    
    Parameters
    ----------
    cartesian : np.ndarray
        Cartesian coordinates
    cell : np.ndarray
        Lattice vectors as rows (3×3)
        
    Returns
    -------
    np.ndarray
        Fractional coordinates
        
    Author
    ------
    Ramon Cardias
    """
    inv_cell = np.linalg.inv(cell)
    frac = np.asarray(cartesian, float) @ inv_cell
    return frac - np.floor(frac)


def fractional_to_cartesian(fractional: np.ndarray, reciprocal_lattice: np.ndarray) -> np.ndarray:
    """Convert fractional k-point coordinates to Cartesian (absolute).
    
    Parameters
    ----------
    fractional : np.ndarray
        Fractional coordinates in reciprocal primitive basis
    reciprocal_lattice : np.ndarray
        Reciprocal lattice vectors as rows (3×3)
        
    Returns
    -------
    np.ndarray
        Cartesian (absolute) k-point coordinates
        
    Author
    ------
    Ramon Cardias
    """
    return np.dot(np.array(fractional, dtype=float), reciprocal_lattice)


def apply_prime_convention(k_absolute: np.ndarray) -> np.ndarray:
    """Apply prime convention: flip in-plane components (kx,ky,kz) → (-kx,-ky,kz).
    
    Used for symmetry-related k-points. Returns a new array without modifying input.
    
    Parameters
    ----------
    k_absolute : np.ndarray
        k-point in absolute (Cartesian) coordinates
        
    Returns
    -------
    np.ndarray
        Primed k-point with (kx, ky) flipped
        
    Author
    ------
    Ramon Cardias
    """
    k = np.array(k_absolute, dtype=float).copy()
    k[0] *= -1.0
    k[1] *= -1.0
    return k


def normalize_label(label: str) -> Tuple[str, bool]:
    """Parse k-point label and extract base name and prime flag.
    
    Parameters
    ----------
    label : str
        Label (e.g., "K", "Γ", "K'", "gamma'")
        
    Returns
    -------
    tuple
        (base_label_uppercase, is_prime)
        
    Author
    ------
    Ramon Cardias
    """
    s = label.strip().replace("Γ", "Gamma").replace("γ", "Gamma").replace("'", "'")
    
    if s.lower() in ("g", "gamma"):
        base = "GAMMA"
    elif s.endswith("'"):
        base = s[:-1].upper()
    else:
        base = s.upper()
    
    is_prime = s.endswith("'")
    return base, is_prime


def tex_format_label(label: str) -> str:
    """Format label for LaTeX rendering in plots.
    
    Parameters
    ----------
    label : str
        K-point label (e.g., "K", "Gamma", "K'")
        
    Returns
    -------
    str
        LaTeX-formatted label (e.g., r"$\\Gamma^{\\prime}$")
        
    Author
    ------
    Ramon Cardias
    """
    if label.upper() == "GAMMA":
        if label.endswith("'"):
            return r"$\Gamma^{\prime}$"
        return r"$\Gamma$"
    
    if label.endswith("'"):
        return r"$" + label[:-1] + r"^{\prime}$"
    
    return r"$" + label + r"$"


def allocate_points_per_segment(total_points: int, segment_lengths: np.ndarray) -> List[int]:
    """Allocate k-points proportionally to segment lengths.
    
    Distributes total_points (as unique points) across segments proportional to their lengths,
    with minimum 2 points per segment.
    
    Parameters
    ----------
    total_points : int
        Total number of unique k-points requested
    segment_lengths : np.ndarray
        Length of each k-path segment in |k| units
        
    Returns
    -------
    list of int
        Number of points for each segment
        
    Author
    ------
    Ramon Cardias
    """
    segment_lengths = np.array(segment_lengths, dtype=float)
    nseg = len(segment_lengths)
    
    # Account for repeated points at segment junctions
    total_with_dups = int(total_points) + (nseg - 1)
    if total_with_dups < 2 * nseg:
        total_with_dups = 2 * nseg
    
    Ltot = float(np.sum(segment_lengths))
    if Ltot <= 0:
        return [2] * nseg
    
    w = segment_lengths / Ltot
    raw = w * total_with_dups
    base = np.floor(raw).astype(int)
    base = np.maximum(base, 2)
    
    # Adjust if over/under allocated
    s = int(np.sum(base))
    while s > total_with_dups:
        idx = np.argmax(base)
        if base[idx] > 2:
            base[idx] -= 1
            s -= 1
        else:
            reducible = np.where(base > 2)[0]
            if len(reducible) == 0:
                break
            base[reducible[0]] -= 1
            s -= 1
    
    rema = raw - np.floor(raw)
    while s < total_with_dups:
        idx = np.argmax(rema)
        base[idx] += 1
        rema[idx] = 0.0
        s += 1
    
    return base.tolist()


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


def build_custom_kpath(
    point_coords: Dict,
    reciprocal_lattice: np.ndarray,
    labels_list: List[str],
    total_kpoints: int,
) -> Tuple[np.ndarray, np.ndarray, List, List, Dict]:
    """Build custom k-path from user-specified high-symmetry points.
    
    Parameters
    ----------
    point_coords : dict
        Dictionary mapping point names to fractional coordinates
    reciprocal_lattice : np.ndarray
        Reciprocal lattice vectors as rows (3×3)
    labels_list : list of str
        Ordered list of k-point labels (e.g., ["K'", "GAMMA", "K"])
    total_kpoints : int
        Desired total number of k-points
        
    Returns
    -------
    tuple
        (kpoints_frac, xvals, xticks, xticklabels, metadata)
        Where kpoints_frac are fractional k-points compatible with calculate_bands()
        
    Author
    ------
    Ramon Cardias
    """
    node_names = []
    node_isprime = []
    node_kabs = []
    node_kfrac = []
    
    for lab in labels_list:
        key, is_prime = normalize_label(lab)
        if key not in point_coords:
            available = sorted(point_coords.keys())
            raise ValueError(
                f"Point '{lab}' not found. Available: {available}"
            )
        
        kf = np.array(point_coords[key], dtype=float)
        ka = fractional_to_cartesian(kf, reciprocal_lattice)
        
        # Apply prime convention if requested
        ka_use = apply_prime_convention(ka) if is_prime else np.array(ka, dtype=float)
        
        node_isprime.append(is_prime)
        node_kfrac.append(kf)
        node_kabs.append(ka_use)
        node_names.append(key + "'" if is_prime else key)
    
    # Compute segment lengths and allocate points
    seg_lengths = [np.linalg.norm(b - a) for a, b in zip(node_kabs, node_kabs[1:])]
    n_per_seg = allocate_points_per_segment(total_kpoints, seg_lengths)
    
    # Build k-point path in Cartesian space (for uniform spacing)
    kpts_abs = []
    xvals = []
    xticks = [0.0]
    xticklabels = [tex_format_label(node_names[0])]
    x = 0.0
    
    for iseg, (A, B, nseg_pts, Lseg) in enumerate(
        zip(node_kabs, node_kabs[1:], n_per_seg, seg_lengths)
    ):
        for i in range(nseg_pts):
            t = i / float(nseg_pts - 1) if nseg_pts > 1 else 0.0
            k = A + (B - A) * t
            if iseg > 0 and i == 0:
                continue  # Skip repeated junction point
            kpts_abs.append(k)
            xvals.append(x + Lseg * t)
        
        x += Lseg
        xticks.append(x)
        xticklabels.append(tex_format_label(node_names[iseg + 1]))
    
    # Convert Cartesian k-points back to fractional
    kpts_abs = np.array(kpts_abs, dtype=float)
    recip_inv = np.linalg.inv(reciprocal_lattice)
    kpts_frac = kpts_abs @ recip_inv
    
    metadata = dict(
        node_names=node_names,
        node_isprime=node_isprime,
        node_kfrac=node_kfrac,
        node_kabs=node_kabs,
        seg_lengths=seg_lengths,
        n_per_seg=n_per_seg,
    )
    
    return np.array(kpts_frac, dtype=float), np.array(xvals, dtype=float), xticks, xticklabels, metadata
