"""Neighbor discovery and the isolated legacy reduction interface."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from pymatgen.core import Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from rsutil.core import RSStructure, StructureSite

from .pairs import map_cluster_neighbor
from .symmetry import (
    MagneticSymmetryContext,
    allowed_public_symmetry_operations,
    periodic_symmetry_operations,
    reduce_pair_keys,
)


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
    analyzer: SpacegroupAnalyzer,
    *,
    structure: Structure | None = None,
    soc: bool = False,
) -> list[Any]:
    """Return allowed public operations for compatibility callers.

    The old implementation used an approximate vector heuristic.  Exact
    periodic reduction now uses :func:`periodic_symmetry_operations`; this
    wrapper remains for callers that need raw pymatgen ``SymmOp`` objects.
    """

    if structure is None:
        if soc:
            # Without structure coordinates the magnetic validity of an
            # operation cannot be proven.  An empty result is the safe value
            # for callers that use this helper to select reduction operations.
            return []
        return list(analyzer.get_symmetry_operations())
    return allowed_public_symmetry_operations(analyzer, structure, soc=soc)


def _canonical_from_pymatgen(structure: Structure) -> RSStructure:
    """Build the dependency-light structure needed by pair-key mapping."""

    return RSStructure(
        lattice_vectors=tuple(
            tuple(float(value) for value in row) for row in structure.lattice.matrix
        ),
        length_unit="angstrom",
        coordinate_convention="fractional",
        sites=tuple(
            StructureSite(
                index=index,
                position=tuple(float(value) for value in site.frac_coords),
                element=str(site.specie.symbol),
                rs_label=(
                    str(site.properties["rs_label"])
                    if site.properties.get("rs_label") is not None
                    else str(site.species_string)
                ),
                rs_type=(
                    str(site.properties["rs_type"])
                    if site.properties.get("rs_type") is not None
                    else None
                ),
                rs_type_index=(
                    int(site.properties["rs_type_index"])
                    if site.properties.get("rs_type_index") is not None
                    else None
                ),
            )
            for index, site in enumerate(structure.sites)
        ),
    )


def _reference_indices(
    recursion_sites: Sequence[int], site_count: int, *, index_base: int
) -> tuple[int, ...]:
    """Normalize an explicitly declared index base to zero-based sites."""

    if index_base not in (0, 1):
        raise ValueError("index_base must be 0 or 1")

    result: list[int] = []
    for value in recursion_sites:
        integer = int(value)
        canonical = integer - index_base
        if not 0 <= canonical < site_count:
            raise ValueError(
                f"recursion site {integer} is outside the reference cell for index_base={index_base}"
            )
        result.append(canonical)
    return tuple(result)


def discover_exact_unique_neighbors(
    cluster: Molecule,
    recursion_sites: list[int] | tuple[int, ...],
    *,
    cutoff: float = 2.0,
    analyzer: SpacegroupAnalyzer | None,
    structure: Structure,
    canonical_structure: RSStructure | None = None,
    soc: bool = False,
    magnetic_axis: tuple[float, float, float] | None = None,
    magnetic_context: MagneticSymmetryContext | None = None,
    magnetic_reference: str | None = None,
    tolerance: float = 1e-8,
    index_base: int = 1,
    verbose: bool = False,
) -> dict[int, list[Any]]:
    """Reduce neighbors using canonical zero-based sites internally.

    ``recursion_sites`` uses the explicit ``index_base`` boundary declaration
    (the RS-LMTO file boundary defaults to one-based). It is normalized once
    before any pair mapping; the pair-key symmetry machinery never performs
    index-base inference.
    """

    if index_base not in (0, 1):
        raise ValueError("index_base must be 0 or 1")

    if analyzer is None:
        # No crystallographic operation source means no scientifically valid
        # reduction.  Returning all pairs is the conservative behavior.
        return discover_default_neighbors(cluster, recursion_sites, cutoff=cutoff, verbose=verbose)
    canonical = canonical_structure or _canonical_from_pymatgen(structure)
    center_sites = _reference_indices(
        recursion_sites, len(canonical.sites), index_base=index_base
    )
    all_neighbors = discover_default_neighbors(cluster, recursion_sites, cutoff=cutoff, verbose=False)
    operations = periodic_symmetry_operations(
        structure,
        analyzer,
        soc=soc,
        magnetic_axis=magnetic_axis,
        magnetic_context=magnetic_context,
        magnetic_reference=magnetic_reference,
        site_identities=tuple(
            (site.element, site.rs_label, site.rs_type)
            for site in canonical.sites
        ),
        tolerance=tolerance,
    )
    mapped: dict[int, list[tuple[Any, Any]]] = {}
    all_keys = []
    for center_index, entries in all_neighbors.items():
        if center_index >= len(center_sites):
            raise ValueError(f"no reference site for recursion center {center_index}")
        mapped_entries: list[tuple[Any, Any]] = []
        for neighbor in entries:
            pair = map_cluster_neighbor(
                canonical,
                [cluster[center_index]],
                0,
                neighbor,
                center_site_index=center_sites[center_index],
                tolerance=tolerance,
            )
            mapped_entries.append((neighbor, pair))
            all_keys.append(pair.key)
        mapped[center_index] = mapped_entries
    representatives = set(reduce_pair_keys(all_keys, operations))
    result = {
        center_index: sorted(
            [neighbor for neighbor, pair in entries if pair.key in representatives],
            key=lambda neighbor: neighbor.nn_distance,
        )
        for center_index, entries in mapped.items()
    }
    if verbose:
        print(
            f"Exact symmetry reduction retained {len(representatives)} of "
            f"{len(all_keys)} periodic pair keys"
        )
        for center, entries in result.items():
            print(f"Unique neighbors for atom {center}: {len(entries)}")
    return result


def discover_legacy_unique_neighbors(
    cluster: Molecule,
    recursion_sites: list[int] | tuple[int, ...],
    *,
    cutoff: float = 2.0,
    shell_width: float = 0.1,
    analyzer: SpacegroupAnalyzer | None = None,
    structure: Structure | None = None,
    canonical_structure: RSStructure | None = None,
    soc: bool = False,
    magnetic_axis: tuple[float, float, float] | None = None,
    magnetic_context: MagneticSymmetryContext | None = None,
    magnetic_reference: str | None = None,
    index_base: int = 1,
    verbose: bool = False,
) -> dict[int, list[Any]]:
    """Compatibility wrapper for exact periodic symmetry reduction.

    ``shell_width`` is retained and validated for CLI compatibility, but it
    is presentation metadata only and never defines equivalence.
    """

    if shell_width <= 0:
        raise ValueError("shell_width must be positive")
    if structure is None:
        return discover_default_neighbors(cluster, recursion_sites, cutoff=cutoff, verbose=verbose)
    return discover_exact_unique_neighbors(
        cluster,
        recursion_sites,
        cutoff=cutoff,
        analyzer=analyzer,
        structure=structure,
        canonical_structure=canonical_structure,
        soc=soc,
        magnetic_axis=magnetic_axis,
        magnetic_context=magnetic_context,
        magnetic_reference=magnetic_reference,
        index_base=index_base,
        verbose=verbose,
    )


def discover_neighbors(
    cluster: Molecule,
    recursion_sites: list[int] | tuple[int, ...],
    *,
    cutoff: float = 2.0,
    unique: bool = False,
    shell_width: float = 0.1,
    analyzer: SpacegroupAnalyzer | None = None,
    structure: Structure | None = None,
    canonical_structure: RSStructure | None = None,
    soc: bool = False,
    magnetic_axis: tuple[float, float, float] | None = None,
    magnetic_context: MagneticSymmetryContext | None = None,
    magnetic_reference: str | None = None,
    full_pairs: bool = False,
    index_base: int = 1,
    verbose: bool = False,
) -> dict[int, list[Any]]:
    """Dispatch through the isolated legacy reduction interface."""

    if unique and not full_pairs:
        return discover_legacy_unique_neighbors(
            cluster,
            recursion_sites,
            cutoff=cutoff,
            shell_width=shell_width,
            analyzer=analyzer,
            structure=structure,
            canonical_structure=canonical_structure,
            soc=soc,
            magnetic_axis=magnetic_axis,
            magnetic_context=magnetic_context,
            magnetic_reference=magnetic_reference,
            index_base=index_base,
            verbose=verbose,
        )
    return discover_default_neighbors(
        cluster, recursion_sites, cutoff=cutoff, verbose=verbose
    )


__all__ = [
    "discover_default_neighbors",
    "discover_exact_unique_neighbors",
    "discover_legacy_unique_neighbors",
    "discover_neighbors",
    "legacy_symmetry_operations",
]
