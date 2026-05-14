"""
Coordinate system conversion utilities for k-point calculations.

Provides functions for converting between Cartesian (absolute) and fractional
k-point coordinates, and computing reciprocal lattice vectors.

Author: Ramon Cardias
"""

import numpy as np


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
