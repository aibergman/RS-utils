"""
K-point allocation utilities for path segmentation.

Provides functions for distributing k-points proportionally across path segments.

Author: Ramon Cardias
"""

from typing import List
import numpy as np


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
