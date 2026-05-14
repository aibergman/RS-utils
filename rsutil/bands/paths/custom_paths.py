"""
Custom k-path builder for user-specified high-symmetry points.

Allows construction of k-paths from user-selected high-symmetry points
with automatic point spacing and interpolation.

Author: Ramon Cardias
"""

import numpy as np
from typing import Dict, List, Tuple

from .coordinates import fractional_to_cartesian
from .symmetry import apply_prime_convention
from .labels import normalize_label, tex_format_label
from .point_allocation import allocate_points_per_segment


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
