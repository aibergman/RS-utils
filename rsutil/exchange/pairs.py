"""Canonical periodic pair identities and cluster-to-cell mapping.

The legacy exchange mapper works with a pymatgen ``Molecule`` and Cartesian
neighbor vectors.  This module keeps that boundary small while making the
scientific identity explicit: a reference-cell site pair and an integer
translation in the reference lattice basis.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from rsutil.core import ExchangePair, PairKey, RSStructure
from rsutil.core.errors import ExchangeError


class PairMappingError(ExchangeError, ValueError):
    """A cluster entry cannot be mapped to exactly one periodic pair."""


def _lattice_matrix(structure: RSStructure) -> np.ndarray:
    matrix = np.asarray(structure.lattice_vectors, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise PairMappingError("reference structure lattice must be a finite 3x3 matrix")
    if abs(float(np.linalg.det(matrix))) <= 1e-12:
        raise PairMappingError("reference structure lattice is singular")
    return matrix


def _cartesian_positions(structure: RSStructure) -> tuple[np.ndarray, ...]:
    matrix = _lattice_matrix(structure)
    if structure.coordinate_convention.value == "fractional":
        return tuple(np.asarray(site.position, dtype=float) @ matrix for site in structure.sites)
    return tuple(np.asarray(site.position, dtype=float) for site in structure.sites)


def _position(entry: Any, description: str) -> np.ndarray:
    value: Any = None
    if isinstance(entry, Mapping):
        for name in ("coords", "position", "cartesian_position"):
            if name in entry:
                value = entry[name]
                break
    else:
        for name in ("coords", "position", "cartesian_position"):
            if hasattr(entry, name):
                value = getattr(entry, name)
                break
    if value is None:
        raise PairMappingError(f"{description} has no Cartesian coordinates")
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise PairMappingError(f"{description} has invalid coordinates") from exc
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise PairMappingError(f"{description} must contain three finite coordinates")
    return result


def _field(entry: Any, names: tuple[str, ...]) -> Any:
    if isinstance(entry, Mapping):
        for name in names:
            if name in entry:
                return entry[name]
        return None
    for name in names:
        if hasattr(entry, name):
            return getattr(entry, name)
    return None


def _element(entry: Any, description: str) -> str:
    value = _field(entry, ("element", "species_string", "species", "label"))
    if value is None:
        raise PairMappingError(f"{description} has no species or element")
    text = str(value).strip()
    match = re.search(r"([A-Z][a-z]?)", text)
    if match is None:
        raise PairMappingError(f"{description} has an invalid species: {text!r}")
    return match.group(1)


def _integer_translation(
    displacement: np.ndarray,
    base_displacement: np.ndarray,
    lattice: np.ndarray,
    pbc: tuple[bool, bool, bool],
    tolerance: float,
) -> tuple[int, int, int] | None:
    fractional_translation = (displacement - base_displacement) @ np.linalg.inv(lattice)
    rounded = np.rint(fractional_translation).astype(int)
    if any(not periodic and value != 0 for periodic, value in zip(pbc, rounded)):
        return None
    residual = displacement - (base_displacement + rounded @ lattice)
    if float(np.linalg.norm(residual)) > tolerance:
        return None
    return tuple(int(value) for value in rounded)


def _candidate_sites(
    structure: RSStructure,
    position: np.ndarray,
    species: str,
    *,
    tolerance: float,
    site_index: int | None = None,
) -> list[tuple[int, tuple[int, int, int]]]:
    lattice = _lattice_matrix(structure)
    positions = _cartesian_positions(structure)
    candidates: list[tuple[int, tuple[int, int, int]]] = []
    indices = range(len(structure.sites)) if site_index is None else (site_index,)
    for index in indices:
        if index < 0 or index >= len(structure.sites):
            raise PairMappingError(f"reference site index {index} is out of range")
        site = structure.sites[index]
        if _element(site, f"reference site {index}") != species:
            continue
        translation = _integer_translation(
            position,
            positions[index],
            lattice,
            structure.pbc,
            tolerance,
        )
        if translation is not None:
            candidates.append((index, translation))
    return candidates


def _unique_candidate_or_fail(
    candidates: list[tuple[int, tuple[int, int, int]]],
    description: str,
) -> tuple[int, tuple[int, int, int]]:
    if not candidates:
        raise PairMappingError(f"{description} does not match a reference-cell site")
    if len(candidates) > 1:
        formatted = ", ".join(f"({index}, {translation})" for index, translation in candidates)
        raise PairMappingError(f"{description} matches ambiguously: {formatted}")
    return candidates[0]


def reference_site_for_cluster_entry(
    structure: RSStructure,
    entry: Any,
    *,
    tolerance: float = 1e-6,
    site_index: int | None = None,
) -> int:
    """Resolve a cluster center to its reference-cell site index.

    The center may be in any periodic image.  Its image translation is not
    part of the center identity and is intentionally discarded.
    """

    species = _element(entry, "cluster entry")
    candidates = _candidate_sites(
        structure,
        _position(entry, "cluster entry"),
        species,
        tolerance=tolerance,
        site_index=site_index,
    )
    return _unique_candidate_or_fail(candidates, "cluster center")[0]


def reconstruct_cartesian_displacement(
    structure: RSStructure, key: PairKey
) -> tuple[float, float, float]:
    """Reconstruct ``r_j + R·A - r_i`` from a canonical pair key."""

    if key.i >= len(structure.sites) or key.j >= len(structure.sites):
        raise PairMappingError("pair key refers to a site outside the reference structure")
    positions = _cartesian_positions(structure)
    displacement = positions[key.j] + np.asarray(key.R, dtype=float) @ _lattice_matrix(structure)
    displacement = displacement - positions[key.i]
    return tuple(float(value) for value in displacement)


def pair_key_from_cluster_neighbor(
    structure: RSStructure,
    center: Any,
    neighbor: Any,
    *,
    center_site_index: int | None = None,
    tolerance: float = 1e-6,
) -> PairKey:
    """Map one Cartesian cluster center/neighbor pair to a :class:`PairKey`."""

    center_position = _position(center, "cluster center")
    neighbor_position = _position(neighbor, "cluster neighbor")
    center_species = _element(center, "cluster center")
    neighbor_species = _element(neighbor, "cluster neighbor")
    center_candidates = _candidate_sites(
        structure,
        center_position,
        center_species,
        tolerance=tolerance,
        site_index=center_site_index,
    )
    center_index, _ = _unique_candidate_or_fail(center_candidates, "cluster center")
    positions = _cartesian_positions(structure)
    neighbor_candidates: list[tuple[int, tuple[int, int, int]]] = []
    lattice = _lattice_matrix(structure)
    displacement = neighbor_position - center_position
    for site_index in range(len(structure.sites)):
        if _element(structure.sites[site_index], f"reference site {site_index}") != neighbor_species:
            continue
        translation = _integer_translation(
            displacement,
            positions[site_index] - positions[center_index],
            lattice,
            structure.pbc,
            tolerance,
        )
        if translation is not None:
            neighbor_candidates.append((site_index, translation))
    neighbor_index, translation = _unique_candidate_or_fail(
        neighbor_candidates, "cluster neighbor"
    )
    return PairKey(center_index, neighbor_index, translation)


def map_cluster_neighbor(
    structure: RSStructure,
    cluster: Sequence[Any],
    center_index: int,
    neighbor: Any,
    *,
    center_site_index: int | None = None,
    tolerance: float = 1e-6,
    shell: int | float | str | None = None,
    shell_metadata: Mapping[str, Any] | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> ExchangePair:
    """Create an :class:`ExchangePair` for one entry from a neighbor list."""

    try:
        center = cluster[center_index]
    except (IndexError, TypeError) as exc:
        raise PairMappingError(f"cluster center index {center_index} is invalid") from exc
    key = pair_key_from_cluster_neighbor(
        structure,
        center,
        neighbor,
        center_site_index=center_site_index,
        tolerance=tolerance,
    )
    observed = _position(neighbor, "cluster neighbor") - _position(center, "cluster center")
    reconstructed = np.asarray(reconstruct_cartesian_displacement(structure, key))
    if float(np.linalg.norm(observed - reconstructed)) > tolerance:
        raise PairMappingError(
            f"pair {key.as_tuple} displacement is inconsistent with the reference lattice"
        )
    entry_index = _field(neighbor, ("index", "cluster_index"))
    pair_provenance = dict(provenance or {})
    pair_provenance.update(
        {
            "mapping": "reference_cell_fractional_translation",
            "mapping_tolerance": tolerance,
            "observed_displacement": [float(value) for value in observed],
        }
    )
    return ExchangePair(
        key=key,
        cartesian_displacement=tuple(float(value) for value in reconstructed),
        distance=float(np.linalg.norm(reconstructed)),
        species_i=structure.sites[key.i].element,
        species_j=structure.sites[key.j].element,
        length_unit=structure.length_unit,
        label_i=structure.sites[key.i].rs_label,
        label_j=structure.sites[key.j].rs_label,
        cluster_index=(int(entry_index) if entry_index is not None else None),
        shell=shell,
        shell_metadata=shell_metadata,
        provenance=pair_provenance,
    )


def map_cluster_neighbors(
    structure: RSStructure,
    cluster: Sequence[Any],
    neighbors: Mapping[int, Sequence[Any]],
    *,
    center_site_indices: Mapping[int, int] | Sequence[int] | None = None,
    tolerance: float = 1e-6,
    provenance: Mapping[str, Any] | None = None,
) -> list[ExchangePair]:
    """Map all cluster neighbors, sorted deterministically by :class:`PairKey`.

    ``center_site_indices`` is an optional explicit mapping when the caller
    already knows which reference sites are recursion centers.  If omitted,
    centers are resolved geometrically and ambiguous matches fail.
    """

    pairs: list[ExchangePair] = []
    for center_index in sorted(neighbors):
        if isinstance(center_site_indices, Mapping):
            reference_index = center_site_indices.get(center_index)
        elif center_site_indices is None:
            reference_index = None
        else:
            try:
                reference_index = center_site_indices[center_index]
            except IndexError as exc:
                raise PairMappingError(
                    f"no reference site mapping for cluster center {center_index}"
                ) from exc
        for neighbor in neighbors[center_index]:
            pairs.append(
                map_cluster_neighbor(
                    structure,
                    cluster,
                    center_index,
                    neighbor,
                    center_site_index=reference_index,
                    tolerance=tolerance,
                    provenance=provenance,
                )
            )
    pairs.sort(key=lambda pair: pair.key.sort_key)
    keys = [pair.key for pair in pairs]
    if len(set(keys)) != len(keys):
        raise PairMappingError("cluster neighbors contain duplicate canonical pair identities")
    return pairs


def legacy_displacement_mapping(pair: ExchangePair, *, one_based: bool = True) -> dict[str, Any]:
    """Serialize a canonical pair in the legacy human-readable mapping shape.

    ``R_ij`` is deliberately labeled as ``cartesian_displacement`` here.  It
    must not be consumed as the canonical lattice translation, which remains
    available as ``pair_key``.
    """

    i, j, R = pair.key.to_one_based() if one_based else pair.key.as_tuple
    return {
        "i": i,
        "j": j,
        "R_ij": list(pair.cartesian_displacement),
        "R_ij_kind": "cartesian_displacement",
        "pair_key": {"i": pair.key.i, "j": pair.key.j, "R": list(R)},
        "distance": pair.distance,
        "length_unit": pair.length_unit,
        "species_i": pair.species_i,
        "species_j": pair.species_j,
    }


def legacy_displacement_mappings(
    pairs: Sequence[ExchangePair], *, one_based: bool = True
) -> list[dict[str, Any]]:
    """Serialize canonical pairs in stable legacy displacement order."""

    return [
        legacy_displacement_mapping(pair, one_based=one_based)
        for pair in sorted(pairs, key=lambda pair: pair.key.sort_key)
    ]


def sort_pair_keys(keys: Sequence[PairKey]) -> list[PairKey]:
    """Return pair keys in the canonical deterministic order."""

    return sorted(keys, key=lambda key: key.sort_key)


# Descriptive aliases for callers that use "mapping" terminology.
map_cluster_to_pairs = map_cluster_neighbors
pair_key_from_neighbor = pair_key_from_cluster_neighbor


__all__ = [
    "ExchangePair",
    "PairKey",
    "PairMappingError",
    "legacy_displacement_mapping",
    "legacy_displacement_mappings",
    "map_cluster_neighbor",
    "map_cluster_neighbors",
    "map_cluster_to_pairs",
    "pair_key_from_cluster_neighbor",
    "pair_key_from_neighbor",
    "reconstruct_cartesian_displacement",
    "reference_site_for_cluster_entry",
    "sort_pair_keys",
]
