"""Exact periodic symmetry operations for exchange pair identities.

The crystallographic operations returned by pymatgen act on fractional
coordinates as ``x' = W x + t``.  If site ``k`` is mapped to site ``p(k)``
plus the integer translation ``n(k)``, a periodic pair transforms as::

    (i, j, R) -> (p(i), p(j), n(j) + W R - n(i))

This module keeps that identity separate from Cartesian distances and shell
labels.  Mapping a floating-point coordinate to a site and an integer shift
is the only place where a numerical tolerance is used.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from rsutil.core import PairKey


class SymmetryMappingError(ValueError):
    """A crystallographic operation cannot be mapped to exact site data."""


@dataclass(frozen=True)
class PeriodicSymmetryOperation:
    """A structure-preserving operation in exact periodic form.

    ``site_translation[k]`` is the integer ``n(k)`` in
    ``W*r[k] + t = r[p(k)] + n(k)``.  ``rotation`` is an integer matrix in
    fractional lattice coordinates.
    """

    rotation: tuple[tuple[int, int, int], ...]
    translation: tuple[float, float, float]
    site_permutation: tuple[int, ...]
    site_translation: tuple[tuple[int, int, int], ...]
    source_index: int | None = None

    def __post_init__(self) -> None:
        matrix = np.asarray(self.rotation, dtype=int)
        if matrix.shape != (3, 3):
            raise ValueError("symmetry rotation must be a 3x3 matrix")
        if len(self.translation) != 3:
            raise ValueError("symmetry translation must contain three values")
        if len(self.site_permutation) != len(self.site_translation):
            raise ValueError("site permutation and translations must have equal length")

    @property
    def rotation_matrix(self) -> np.ndarray:
        """Return the fractional-coordinate rotation as a fresh array."""

        return np.asarray(self.rotation, dtype=int)

    @property
    def translation_vector(self) -> np.ndarray:
        """Return the crystallographic fractional translation."""

        return np.asarray(self.translation, dtype=float)

    def transform(self, key: PairKey) -> PairKey:
        """Apply this operation to one canonical periodic pair key."""

        if key.i >= len(self.site_permutation) or key.j >= len(self.site_permutation):
            raise SymmetryMappingError(
                "pair key refers to a site outside the symmetry operation"
            )
        n_i = np.asarray(self.site_translation[key.i], dtype=int)
        n_j = np.asarray(self.site_translation[key.j], dtype=int)
        transformed_R = n_j + self.rotation_matrix @ np.asarray(key.R, dtype=int) - n_i
        return PairKey(
            self.site_permutation[key.i],
            self.site_permutation[key.j],
            tuple(int(value) for value in transformed_R),
        )


@dataclass(frozen=True)
class PairSymmetryReduction:
    """Representatives and available members of exact pair-key orbits."""

    representatives: tuple[PairKey, ...]
    orbits: dict[PairKey, tuple[PairKey, ...]]
    representative_for: dict[PairKey, PairKey]

    def expand(self, operations: Sequence[PeriodicSymmetryOperation]) -> frozenset[PairKey]:
        """Expand representatives through the supplied symmetry group."""

        return frozenset(
            member
            for representative in self.representatives
            for member in pair_key_orbit(representative, operations)
        )


def _as_species(site: Any) -> str:
    specie = getattr(site, "specie", None)
    if specie is not None and hasattr(specie, "symbol"):
        return str(specie.symbol)
    return str(getattr(site, "species_string", ""))


def _integer_rotation(operation: Any, tolerance: float) -> np.ndarray | None:
    matrix = np.asarray(operation.rotation_matrix, dtype=float)
    rounded = np.rint(matrix).astype(int)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        return None
    if not np.allclose(matrix, rounded, atol=tolerance, rtol=0):
        return None
    if abs(round(float(np.linalg.det(rounded)))) != 1:
        return None
    return rounded


def _site_mapping(
    structure: Structure,
    rotation: np.ndarray,
    translation: np.ndarray,
    tolerance: float,
) -> tuple[tuple[int, ...], tuple[tuple[int, int, int], ...]] | None:
    fractional = np.asarray(structure.frac_coords, dtype=float)
    species = tuple(_as_species(site) for site in structure.sites)
    permutation: list[int] = []
    shifts: list[tuple[int, int, int]] = []
    for index, position in enumerate(fractional):
        transformed = rotation @ position + translation
        candidates: list[tuple[int, np.ndarray]] = []
        for target, target_position in enumerate(fractional):
            if species[target] != species[index]:
                continue
            shift = np.rint(transformed - target_position).astype(int)
            residual = transformed - target_position - shift
            if np.linalg.norm(residual) <= tolerance:
                candidates.append((target, shift))
        if len(candidates) != 1:
            return None
        target, shift = candidates[0]
        permutation.append(target)
        shifts.append(tuple(int(value) for value in shift))
    if len(set(permutation)) != len(permutation):
        return None
    return tuple(permutation), tuple(shifts)


def _cartesian_rotation(structure: Structure, rotation: np.ndarray) -> np.ndarray | None:
    """Convert a fractional-coordinate operation to Cartesian coordinates."""

    basis = np.asarray(structure.lattice.matrix, dtype=float)
    try:
        cartesian = basis.T @ rotation @ np.linalg.inv(basis.T)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(cartesian)):
        return None
    if not np.allclose(cartesian.T @ cartesian, np.eye(3), atol=1e-7, rtol=0):
        return None
    return cartesian


def periodic_symmetry_operations(
    structure: Structure,
    analyzer: SpacegroupAnalyzer | None = None,
    *,
    soc: bool = False,
    magnetic_axis: Sequence[float] | None = None,
    tolerance: float = 1e-8,
) -> tuple[PeriodicSymmetryOperation, ...]:
    """Return exact structure-preserving operations from public pymatgen APIs.

    With ``soc=False`` all operations that map the supplied chemical
    structure to itself are retained.  With ``soc=True`` the legacy
    convention is a collinear magnetic reference with a Cartesian axis of
    ``(0, 0, 1)`` unless ``magnetic_axis`` is supplied.  The magnetization is
    an axial vector, so an operation is retained only when
    ``det(Q) * Q @ axis == axis``.  Time reversal is not assumed.  This is
    intentionally not a generic magnetic-space-group implementation.
    """

    analyzer = analyzer or SpacegroupAnalyzer(structure)
    if tolerance <= 0:
        raise ValueError("symmetry tolerance must be positive")
    axis: np.ndarray | None = None
    if soc:
        raw_axis = (0.0, 0.0, 1.0) if magnetic_axis is None else magnetic_axis
        axis = np.asarray(tuple(raw_axis), dtype=float)
        if axis.shape != (3,) or not np.all(np.isfinite(axis)):
            raise ValueError("magnetic_axis must contain three finite values")
        norm = float(np.linalg.norm(axis))
        if norm <= tolerance:
            raise ValueError("magnetic_axis must be non-zero")
        axis = axis / norm

    result: list[PeriodicSymmetryOperation] = []
    for source_index, operation in enumerate(analyzer.get_symmetry_operations()):
        rotation = _integer_rotation(operation, tolerance)
        if rotation is None:
            continue
        translation = np.asarray(operation.translation_vector, dtype=float)
        if translation.shape != (3,) or not np.all(np.isfinite(translation)):
            continue
        mapping = _site_mapping(structure, rotation, translation, tolerance)
        if mapping is None:
            # An uncertain mapping cannot be used to eliminate a calculation.
            continue
        if axis is not None:
            cartesian_rotation = _cartesian_rotation(structure, rotation)
            if cartesian_rotation is None:
                continue
            axial = float(np.linalg.det(cartesian_rotation)) * cartesian_rotation @ axis
            if not np.allclose(axial, axis, atol=1e-7, rtol=0):
                continue
        permutation, shifts = mapping
        result.append(
            PeriodicSymmetryOperation(
                rotation=tuple(tuple(int(value) for value in row) for row in rotation),
                translation=tuple(float(value) for value in translation),
                site_permutation=permutation,
                site_translation=shifts,
                source_index=source_index,
            )
        )

    # A valid analyzer always includes identity, but retaining an explicit
    # identity fallback keeps failure conservative if a third-party analyzer
    # supplies an incomplete operation list.
    if not any(
        operation.site_permutation == tuple(range(len(structure.sites)))
        and operation.rotation == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        and all(
            value == 0
            for shift in operation.site_translation
            for value in shift
        )
        for operation in result
    ):
        identity_mapping = tuple(range(len(structure.sites)))
        result.insert(
            0,
            PeriodicSymmetryOperation(
                rotation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                translation=(0.0, 0.0, 0.0),
                site_permutation=identity_mapping,
                site_translation=tuple((0, 0, 0) for _ in identity_mapping),
                source_index=None,
            ),
        )
    return tuple(result)


def pair_key_transform(
    key: PairKey, operation: PeriodicSymmetryOperation
) -> PairKey:
    """Apply one exact periodic operation to a :class:`PairKey`."""

    return operation.transform(key)


def pair_key_orbit(
    key: PairKey, operations: Sequence[PeriodicSymmetryOperation]
) -> frozenset[PairKey]:
    """Return the finite orbit generated by exact pair transformations."""

    orbit = {key}
    pending = [key]
    while pending:
        current = pending.pop()
        for operation in operations:
            transformed = operation.transform(current)
            if transformed not in orbit:
                orbit.add(transformed)
                pending.append(transformed)
    return frozenset(orbit)


def reduce_pair_keys_with_orbits(
    keys: Iterable[PairKey],
    operations: Sequence[PeriodicSymmetryOperation],
) -> PairSymmetryReduction:
    """Choose one deterministic representative for each available orbit.

    Only members present in ``keys`` are reduced.  This matters when a caller
    asks for a subset of centers: an operation that maps outside that subset
    cannot silently remove a pair that has no equivalent available job.
    """

    available = set(keys)
    remaining = set(available)
    representatives: list[PairKey] = []
    orbits: dict[PairKey, tuple[PairKey, ...]] = {}
    representative_for: dict[PairKey, PairKey] = {}
    while remaining:
        seed = min(remaining, key=lambda item: item.sort_key)
        members = tuple(sorted(pair_key_orbit(seed, operations) & available, key=lambda item: item.sort_key))
        representative = members[0]
        representatives.append(representative)
        orbits[representative] = members
        for member in members:
            representative_for[member] = representative
        remaining.difference_update(members)
    return PairSymmetryReduction(
        representatives=tuple(representatives),
        orbits=orbits,
        representative_for=representative_for,
    )


def reduce_pair_keys(
    keys: Iterable[PairKey],
    operations: Sequence[PeriodicSymmetryOperation],
) -> tuple[PairKey, ...]:
    """Return deterministic exact-orbit representatives for ``keys``."""

    return reduce_pair_keys_with_orbits(keys, operations).representatives


def allowed_public_symmetry_operations(
    analyzer: SpacegroupAnalyzer,
    structure: Structure,
    *,
    soc: bool = False,
    magnetic_axis: Sequence[float] | None = None,
    tolerance: float = 1e-8,
) -> list[Any]:
    """Compatibility view of the allowed public pymatgen operations."""

    raw = list(analyzer.get_symmetry_operations())
    allowed = {
        operation.source_index
        for operation in periodic_symmetry_operations(
            structure,
            analyzer,
            soc=soc,
            magnetic_axis=magnetic_axis,
            tolerance=tolerance,
        )
        if operation.source_index is not None
    }
    return [operation for index, operation in enumerate(raw) if index in allowed]


__all__ = [
    "PairSymmetryReduction",
    "PeriodicSymmetryOperation",
    "SymmetryMappingError",
    "allowed_public_symmetry_operations",
    "pair_key_orbit",
    "pair_key_transform",
    "periodic_symmetry_operations",
    "reduce_pair_keys",
    "reduce_pair_keys_with_orbits",
]
