"""Neighbor discovery and the isolated legacy reduction interface."""

from __future__ import annotations

from typing import Any

import numpy as np
from pymatgen.core import Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def discover_default_neighbors(
    cluster: Molecule,
    recursion_sites: list[int] | tuple[int, ...],
    *,
    cutoff: float = 2.0,
    verbose: bool = False,
) -> dict[int, list[Any]]:
    """Discover all non-unique neighbors for each recursion center."""

    neighbors: dict[int, list[Any]] = {}
    for cluster_index, _ in enumerate(recursion_sites):
        entries = sorted(
            cluster.get_neighbors(cluster[cluster_index], cutoff),
            key=lambda neighbor: neighbor.nn_distance,
        )
        neighbors[cluster_index] = entries
    if verbose:
        print("Neighbors:")
        for center, entries in neighbors.items():
            print(f"Neighbors for {center}:")
            for neighbor in entries:
                print(
                    f"  {neighbor.species_string} at {neighbor.coords} "
                    f"with distance {neighbor.nn_distance}"
                )
    return neighbors


def legacy_symmetry_operations(
    analyzer: SpacegroupAnalyzer, *, soc: bool = False
) -> list[Any]:
    """Return the historical SOC-filtered operation list without private mutation."""

    operations = list(analyzer.get_symmetry_operations())
    if not soc:
        return operations
    seen_xy_planes: set[tuple[tuple[float, ...], int]] = set()
    for operation in operations:
        rotation = operation.rotation_matrix
        if np.allclose(rotation[0:2, 2], 0) and np.allclose(rotation[2, 0:2], 0):
            xy_part = tuple(np.round(rotation[0:2, 0:2], 5).flatten())
            seen_xy_planes.add((xy_part, int(np.sign(rotation[2, 2]))))
    filtered: list[Any] = []
    for operation in operations:
        rotation = operation.rotation_matrix
        if np.allclose(rotation[0:2, 2], 0) and np.allclose(rotation[2, 0:2], 0):
            xy_part = tuple(np.round(rotation[0:2, 0:2], 5).flatten())
            z_component = int(np.sign(rotation[2, 2]))
            if z_component > 0 or (xy_part, 1) not in seen_xy_planes:
                filtered.append(operation)
        else:
            filtered.append(operation)
    return filtered


def discover_legacy_unique_neighbors(
    cluster: Molecule,
    recursion_sites: list[int] | tuple[int, ...],
    *,
    cutoff: float = 2.0,
    shell_width: float = 0.1,
    analyzer: SpacegroupAnalyzer | None = None,
    structure: Structure | None = None,
    soc: bool = False,
    verbose: bool = False,
) -> dict[int, list[Any]]:
    """Run the pre-HT-09 unique-neighbor algorithm behind an explicit hook."""

    if shell_width <= 0:
        raise ValueError("shell_width must be positive")
    result: dict[int, list[Any]] = {}
    operations = legacy_symmetry_operations(analyzer, soc=soc) if analyzer else []

    def equivalent(first: np.ndarray, second: np.ndarray) -> bool:
        if analyzer is None:
            return bool(np.allclose(first, second, atol=1e-5))
        if structure is None:
            return any(
                np.allclose(operation.rotation_matrix @ first, second, atol=1e-4)
                for operation in operations
            )
        frac_first = structure.lattice.get_fractional_coords([first])[0]
        frac_second = structure.lattice.get_fractional_coords([second])[0]
        for operation in operations:
            rotated = operation.rotation_matrix @ frac_first
            difference = np.mod(rotated - frac_second + 0.5, 1.0) - 0.5
            if np.allclose(difference, 0.0, atol=1e-4):
                return True
        return False

    for center_index, _ in enumerate(recursion_sites):
        center = cluster[center_index]
        shell_groups: dict[float, list[Any]] = {}
        for neighbor in cluster.get_neighbors(center, cutoff):
            shell = round(neighbor.nn_distance / shell_width) * shell_width
            shell_groups.setdefault(shell, []).append(neighbor)
        selected: list[Any] = []
        for shell in sorted(shell_groups):
            representatives: list[tuple[str, np.ndarray]] = []
            for neighbor in shell_groups[shell]:
                vector = np.asarray(neighbor.coords - center.coords)
                species = neighbor.species_string
                is_new = not any(
                    species == seen_species and equivalent(vector, seen_vector)
                    for seen_species, seen_vector in representatives
                )
                if is_new:
                    representatives.append((species, vector))
                    selected.append(neighbor)
        result[center_index] = sorted(selected, key=lambda neighbor: neighbor.nn_distance)
        if verbose:
            print(f"Unique neighbors for atom {center_index}: {len(result[center_index])}")
    return result


def discover_neighbors(
    cluster: Molecule,
    recursion_sites: list[int] | tuple[int, ...],
    *,
    cutoff: float = 2.0,
    unique: bool = False,
    shell_width: float = 0.1,
    analyzer: SpacegroupAnalyzer | None = None,
    structure: Structure | None = None,
    soc: bool = False,
    verbose: bool = False,
) -> dict[int, list[Any]]:
    """Dispatch through the isolated legacy reduction interface."""

    if unique:
        return discover_legacy_unique_neighbors(
            cluster,
            recursion_sites,
            cutoff=cutoff,
            shell_width=shell_width,
            analyzer=analyzer,
            structure=structure,
            soc=soc,
            verbose=verbose,
        )
    return discover_default_neighbors(
        cluster, recursion_sites, cutoff=cutoff, verbose=verbose
    )


__all__ = [
    "discover_default_neighbors",
    "discover_legacy_unique_neighbors",
    "discover_neighbors",
    "legacy_symmetry_operations",
]
