"""HT-09 exact periodic and magnetic-reference symmetry coverage."""

from __future__ import annotations

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from rsutil.core import PairKey
from rsutil.exchange.symmetry import (
    pair_key_orbit,
    periodic_symmetry_operations,
    reduce_pair_keys_with_orbits,
)


def _bcc_fe() -> Structure:
    return Structure(
        Lattice.cubic(2.86),
        ["Fe", "Fe"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )


def _b2_fe_co() -> Structure:
    return Structure(
        Lattice.cubic(2.9),
        ["Fe", "Co"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )


def _fe_pt() -> Structure:
    return Structure(
        Lattice.cubic(3.0),
        ["Fe", "Pt"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )


def test_bcc_translation_bearing_pair_reduces_by_exact_periodic_orbit() -> None:
    structure = _bcc_fe()
    analyzer = SpacegroupAnalyzer(structure)
    operations = periodic_symmetry_operations(structure, analyzer)
    key = PairKey(0, 0, (1, 0, 0))
    orbit = pair_key_orbit(key, operations)

    assert len(orbit) == 12
    assert PairKey(0, 0, (-1, 0, 0)) in orbit
    assert PairKey(0, 0, (0, 1, 0)) in orbit
    reduction = reduce_pair_keys_with_orbits(orbit, operations)
    assert reduction.representatives == (min(orbit, key=lambda item: item.sort_key),)
    assert reduction.expand(operations) == orbit


def test_multisublattice_symmetry_preserves_species_and_never_merges_cross_species() -> None:
    structure = _b2_fe_co()
    analyzer = SpacegroupAnalyzer(structure)
    operations = periodic_symmetry_operations(structure, analyzer)
    keys = set(pair_key_orbit(PairKey(0, 1, (1, 0, 0)), operations))
    keys.update(pair_key_orbit(PairKey(1, 0, (0, 0, 0)), operations))
    reduction = reduce_pair_keys_with_orbits(keys, operations)

    assert len(reduction.representatives) == 2
    assert all(
        (key.i, key.j) in {(0, 1), (1, 0)} for key in reduction.representatives
    )
    assert reduction.expand(operations) == keys


def test_soc_uses_axial_z_reference_and_excludes_axis_flipping_operations() -> None:
    structure = _fe_pt()
    analyzer = SpacegroupAnalyzer(structure)
    non_soc = periodic_symmetry_operations(structure, analyzer)
    soc = periodic_symmetry_operations(structure, analyzer, soc=True)

    assert len(soc) < len(non_soc)
    assert all(
        np.allclose(
            np.linalg.det(
                np.asarray(structure.lattice.matrix).T
                @ operation.rotation_matrix
                @ np.linalg.inv(np.asarray(structure.lattice.matrix).T)
            )
            * (
                np.asarray(structure.lattice.matrix).T
                @ operation.rotation_matrix
                @ np.linalg.inv(np.asarray(structure.lattice.matrix).T)
            )
            @ np.array([0.0, 0.0, 1.0]),
            [0.0, 0.0, 1.0],
            atol=1e-7,
        )
        for operation in soc
    )
    assert any(
        not np.allclose(
            np.linalg.det(
                np.asarray(structure.lattice.matrix).T
                @ operation.rotation_matrix
                @ np.linalg.inv(np.asarray(structure.lattice.matrix).T)
            )
            * (
                np.asarray(structure.lattice.matrix).T
                @ operation.rotation_matrix
                @ np.linalg.inv(np.asarray(structure.lattice.matrix).T)
            )
            @ np.array([0.0, 0.0, 1.0]),
            [0.0, 0.0, 1.0],
            atol=1e-7,
        )
        for operation in non_soc
    )


def test_pair_reduction_does_not_change_derived_distance_or_species() -> None:
    structure = _b2_fe_co()
    analyzer = SpacegroupAnalyzer(structure)
    operations = periodic_symmetry_operations(structure, analyzer)
    keys = {
        PairKey(0, 1, (0, 0, 0)),
        PairKey(0, 1, (1, 0, 0)),
        PairKey(0, 1, (0, 1, 0)),
    }
    reduction = reduce_pair_keys_with_orbits(keys, operations)
    frac = np.asarray(structure.frac_coords)
    lattice = np.asarray(structure.lattice.matrix)

    for representative in reduction.representatives:
        displacement = (frac[representative.j] + representative.R) @ lattice - frac[
            representative.i
        ] @ lattice
        assert np.linalg.norm(displacement) > 0
        assert structure[representative.i].specie.symbol == "Fe"
        assert structure[representative.j].specie.symbol == "Co"
