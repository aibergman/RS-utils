"""
K-point path generation utilities.

Provides functions for generating high-symmetry k-point paths using various methods
(seekpath, spglib, pymatgen, or manual fallback).

Author: Ramon Cardias
"""

import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..readers import LatticeReader

try:
    import spglib

    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False

try:
    from pymatgen.core import Structure
    from pymatgen.symmetry.bandstructure import HighSymmKpath

    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False

try:
    import seekpath

    HAS_SEEKPATH = True
except ImportError:
    HAS_SEEKPATH = False


def setup_kpath(
    lattice_reader: LatticeReader, num_points: int = 100
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Setup k-point path using structure analysis.
    
    Uses StructureAnalyzer to detect path automatically, with manual fallback.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader with lattice structure
    num_points : int, default=100
        Desired number of k-points
        
    Returns
    -------
    tuple
        (kpoints, klabels, kdistances)
        
    """
    # Import here to avoid circular imports
    from ..analysis.structure_analysis import StructureAnalyzer

    # Analyze structure and generate proper k-path
    analyzer = StructureAnalyzer(lattice_reader)

    if analyzer.special_points and analyzer.band_path:
        return analyzer.generate_kpath_points(num_points)
    else:
        # Fallback to manual method
        return _setup_kpath_manual(lattice_reader, num_points)


def setup_kpath_seekpath(
    lattice_reader: LatticeReader,
    num_points: int = 100,
    interactive: bool = False,
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Setup k-point path using seekpath with optional interactive selection.
    
    Provides automatic high-symmetry k-path via seekpath package, with options
    for interactive selection between standard and custom paths.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader with lattice structure
    num_points : int, default=100
        Desired total number of k-points
    interactive : bool, default=False
        If True, prompt user to choose standard or custom path
        If False, use standard seekpath k-path
        
    Returns
    -------
    tuple
        (kpoints, klabels, kdistances)
        - kpoints: Array of k-points (nk, 3)
        - klabels: Labels for k-points (list of str)
        - kdistances: Linear coordinates along path (nk,)
        
    Raises
    ------
    ImportError
        If seekpath is not installed
        
    Notes
    -----
    Requires: pip install seekpath
    
    """
    if not HAS_SEEKPATH:
        raise ImportError("seekpath required. Install with: pip install seekpath")

    from ..paths import seekpath_kpath, build_custom_kpath

    # Get lattice info
    lat_vec = lattice_reader.lattice_vectors
    pos_cart = np.array(
        [[x, y, z] for x, y, z in lattice_reader.atomic_positions], dtype=float
    )

    # Convert positions to fractional
    inv_cell = np.linalg.inv(lat_vec)
    pos_frac = pos_cart @ inv_cell
    pos_frac = pos_frac - np.floor(pos_frac)

    # Generate paths with requested number of points
    kpath_data = seekpath_kpath(
        lat_vec,
        pos_frac,
        lattice_reader.atomic_numbers,
        num_points=num_points,
        reference_distance=0.025,
        symprec=1e-5,
    )

    # Interactive selection if requested
    if interactive:
        ans = input("Use standard seekpath k-path? [Y/n]: ").strip().lower()
        if ans == "" or ans.startswith("y"):
            # Use fractional k-points (what band calculator expects)
            kpoints = kpath_data["kpoints_frac"]
            klabels = kpath_data["klabels"]
            kdistances = kpath_data["xvals"]
        else:
            # Custom path
            print("\n" + "=" * 60)
            print("AVAILABLE SYMMETRY POINTS")
            print("=" * 60)
            point_coords = kpath_data["point_coords"]
            available_points = sorted(point_coords.keys())
            for point in available_points:
                coords = point_coords[point]
                print(
                    f"  {point:6s}: ({coords[0]:7.4f}, {coords[1]:7.4f}, {coords[2]:7.4f})"
                )
            print("=" * 60)

            nk_s = (
                input(f"\nTotal k-points for custom path? [default {num_points}]: ")
                .strip()
            )
            total_k = num_points if nk_s == "" else int(nk_s)

            p = input("Enter k-path labels separated by '-' (e.g. K'-Γ-K): ").strip()
            if not p:
                raise ValueError("Empty custom path")
            labels_list = [x.strip() for x in p.split("-") if x.strip()]
            if len(labels_list) < 2:
                raise ValueError("Need at least 2 symmetry points")

            kpoints, kdistances, xticks, xticklabels, metadata = build_custom_kpath(
                point_coords, kpath_data["reciprocal_lattice"], labels_list, total_k
            )
            klabels = [""] * len(kpoints)
    else:
        # Use fractional k-points (what band calculator expects)
        kpoints = kpath_data["kpoints_frac"]
        klabels = kpath_data["klabels"]
        kdistances = kpath_data["xvals"]

    # Print summary if available
    print(f"\nK-path: {len(kpoints)} points")

    return kpoints, klabels, kdistances


def _setup_kpath_spglib(
    lattice_reader: LatticeReader, num_points: int
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Setup k-path using spglib.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader with lattice structure
    num_points : int
        Desired number of k-points
        
    Returns
    -------
    tuple
        (kpoints, klabels, kdistances)
        
    Author
    ------
    Ramon Cardias
    """
    cell = (
        lattice_reader.lattice_vectors.T,
        lattice_reader.atomic_positions,
        lattice_reader.atomic_numbers,
    )

    # Get the symmetry information
    symmetry = spglib.get_symmetry_dataset(cell)

    # Get the standardized cell
    mapping, grid = spglib.get_ir_reciprocal_mesh([20, 20, 20], cell)

    # Define high symmetry points (this is crystal system dependent)
    # For now, use a simple path
    lattice_type = symmetry["international"]

    if "cubic" in lattice_type.lower() or "P" in lattice_type:
        # Simple cubic path
        special_points = {
            "G": [0, 0, 0],  # Gamma
            "X": [0.5, 0, 0],  # X
            "M": [0.5, 0.5, 0],  # M
            "R": [0.5, 0.5, 0.5],  # R
        }
        path = ["G", "X", "M", "G", "R", "X"]
    else:
        # Generic path
        special_points = {
            "G": [0, 0, 0],
            "X": [0.5, 0, 0],
            "Y": [0, 0.5, 0],
            "Z": [0, 0, 0.5],
        }
        path = ["G", "X", "Y", "G", "Z"]

    return _generate_kpath(special_points, path, num_points)


def _setup_kpath_pymatgen(
    lattice_reader: LatticeReader, num_points: int
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Setup k-path using pymatgen.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader with lattice structure
    num_points : int
        Desired number of k-points
        
    Returns
    -------
    tuple
        (kpoints, klabels, kdistances)
        
    Author
    ------
    Ramon Cardias
    """
    structure = Structure(
        lattice=lattice_reader.lattice_vectors.T,
        species=lattice_reader.atomic_numbers,
        coords=lattice_reader.atomic_positions,
    )

    kpath = HighSymmKpath(structure)
    special_points = kpath.kpath["kpoints"]
    path = kpath.kpath["path"][0]  # Take first path

    return _generate_kpath(special_points, path, num_points)


def _setup_kpath_manual(
    lattice_reader: LatticeReader, num_points: int
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Fallback manual k-path setup.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader with lattice structure
    num_points : int
        Desired number of k-points
        
    Returns
    -------
    tuple
        (kpoints, klabels, kdistances)
        
    Author
    ------
    Ramon Cardias
    """
    warnings.warn(
        "Using manual k-path. Install spglib or pymatgen for automatic detection."
    )

    special_points = {
        "G": [0, 0, 0],
        "X": [0.5, 0, 0],
        "Y": [0, 0.5, 0],
        "Z": [0, 0, 0.5],
        "M": [0.5, 0.5, 0],
    }
    path = ["G", "X", "M", "G", "Z"]

    return _generate_kpath(special_points, path, num_points)


def _generate_kpath(
    special_points: Dict, path: List[str], num_points: int
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Generate k-point path from special points.
    
    Parameters
    ----------
    special_points : dict
        Mapping of point names to fractional coordinates
    path : list of str
        Ordered list of point names defining the path
    num_points : int
        Total desired number of k-points
        
    Returns
    -------
    tuple
        (kpoints, klabels, kdistances)
        
    Author
    ------
    Ramon Cardias
    """
    kpoints = []
    klabels = []
    kdistances = []
    total_distance = 0

    for i in range(len(path) - 1):
        start = np.array(special_points[path[i]])
        end = np.array(special_points[path[i + 1]])

        # Calculate number of points for this segment
        segment_length = np.linalg.norm(end - start)
        segment_points = max(2, int(num_points * segment_length / 3))  # Rough estimate

        # Generate points along segment
        for j in range(segment_points):
            if j == 0 and i == 0:  # First point
                kpoints.append(start)
                klabels.append(path[i])
                kdistances.append(total_distance)
            elif j == 0:  # Skip duplicate points
                klabels[-1] = f"{klabels[-1]}|{path[i]}"
            else:
                t = j / (segment_points - 1)
                k = start + t * (end - start)
                kpoints.append(k)
                klabels.append("")
                total_distance += (
                    np.linalg.norm(k - kpoints[-2]) if len(kpoints) > 1 else 0
                )
                kdistances.append(total_distance)

        # Add end point
        if i == len(path) - 2:  # Last segment
            kpoints.append(end)
            klabels.append(path[i + 1])
            total_distance += np.linalg.norm(end - kpoints[-2])
            kdistances.append(total_distance)

    return np.array(kpoints), klabels, np.array(kdistances)
