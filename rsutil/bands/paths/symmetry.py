"""
Symmetry operations for k-point path generation.

Implements prime convention for symmetry-related k-points used in band structure analysis.

Author: Ramon Cardias
"""

import numpy as np


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
