"""HT-07 canonical periodic pair identity coverage."""

from __future__ import annotations

import numpy as np
import pytest

from rsutil.core import ExchangePair, PairKey, RSStructure, StructureSite
from rsutil.exchange.pairs import (
    PairMappingError,
    legacy_displacement_mapping,
    map_cluster_neighbors,
    pair_key_from_cluster_neighbor,
    reconstruct_cartesian_displacement,
)


def _bcc_fe() -> RSStructure:
    return RSStructure(
        lattice_vectors=((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5)),
        length_unit="angstrom",
        coordinate_convention="fractional",
        sites=(StructureSite(0, (0, 0, 0), "Fe", 26, rs_label="Fe1"),),
    )


def _fe_pt() -> RSStructure:
    return RSStructure(
        lattice_vectors=((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        length_unit="angstrom",
        coordinate_convention="fractional",
        sites=(
            StructureSite(0, (0, 0, 0), "Fe", 26, rs_label="Fe1"),
            StructureSite(1, (0.5, 0.5, 0.5), "Pt", 78, rs_label="Pt1"),
        ),
    )


def _fe_n_o() -> RSStructure:
    return RSStructure(
        lattice_vectors=((3, 0, 0), (0.2, 2.5, 0), (0, 0.1, 4)),
        length_unit="angstrom",
        coordinate_convention="fractional",
        sites=(
            StructureSite(0, (0, 0, 0), "Fe", 26, rs_label="Fe1"),
            StructureSite(1, (0.25, 0.25, 0), "N", 7, rs_label="N1"),
            StructureSite(2, (0.5, 0, 0.5), "O", 8, rs_label="O1"),
        ),
    )


def test_pair_key_is_hashable_reversible_and_has_explicit_index_boundaries() -> None:
    key = PairKey(2, 4, (-3, 0, 2))
    assert hash(key)
    assert key.reverse() == PairKey(4, 2, (3, 0, -2))
    assert key.to_one_based() == (3, 5, (-3, 0, 2))
    assert PairKey.from_one_based(*key.to_one_based()) == key
    assert PairKey.from_dict(key.to_dict()) == key
    with pytest.raises(ValueError):
        PairKey.from_one_based(0, 1, (0, 0, 0))


def test_bcc_fe_exact_pair_and_positive_negative_translations() -> None:
    structure = _bcc_fe()
    lattice = np.asarray(structure.lattice_vectors)
    center = {"coords": [0, 0, 0], "element": "Fe", "label": "Fe1"}
    positive = {"coords": lattice[0], "element": "Fe", "index": 10}
    negative = {"coords": -lattice[2], "element": "Fe", "index": 11}

    assert pair_key_from_cluster_neighbor(structure, center, positive) == PairKey(
        0, 0, (1, 0, 0)
    )
    pair = map_cluster_neighbors(
        structure,
        [center],
        {0: [negative, positive]},
    )
    assert [item.key for item in pair] == [PairKey(0, 0, (0, 0, -1)), PairKey(0, 0, (1, 0, 0))]
    for item in pair:
        assert np.allclose(
            item.cartesian_displacement,
            reconstruct_cartesian_displacement(structure, item.key),
        )


def test_multisite_multispecies_mapping_is_order_independent_and_crosses_boundary() -> None:
    structure = _fe_pt()
    center = {"coords": [0, 0, 0], "element": "Fe", "label": "Fe1"}
    pt = {"coords": [1, 1, 1], "element": "Pt", "index": 20}
    fe_negative = {"coords": [-2, 0, 0], "element": "Fe", "index": 21}
    first = map_cluster_neighbors(structure, [center], {0: [pt, fe_negative]})
    second = map_cluster_neighbors(structure, [center], {0: [fe_negative, pt]})
    assert [item.key for item in first] == [PairKey(0, 0, (-1, 0, 0)), PairKey(0, 1, (0, 0, 0))]
    assert [item.key for item in first] == [item.key for item in second]
    assert first[1].label_i == "Fe1"
    assert first[1].label_j == "Pt1"


def test_second_multisite_fixture_reconstructs_nonorthogonal_geometry() -> None:
    structure = _fe_n_o()
    lattice = np.asarray(structure.lattice_vectors)
    positions = np.asarray([site.position for site in structure.sites]) @ lattice
    center = {"coords": positions[1] + lattice[1], "element": "N", "label": "N1"}
    neighbor = {"coords": positions[2] + lattice[1] - lattice[2], "element": "O", "index": 4}
    key = pair_key_from_cluster_neighbor(structure, center, neighbor)
    assert key == PairKey(1, 2, (0, 0, -1))
    pair = map_cluster_neighbors(structure, [center], {0: [neighbor]})[0]
    assert np.allclose(pair.cartesian_displacement, positions[2] - positions[1] - lattice[2])


def test_ambiguous_mapping_fails_explicitly() -> None:
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="angstrom",
        coordinate_convention="fractional",
        sites=(
            StructureSite(0, (0, 0, 0), "Fe", 26),
            StructureSite(1, (0.5, 0.5, 0.5), "Fe", 26),
            StructureSite(2, (0.5, 0.5, 0.5), "Fe", 26),
        ),
    )
    with pytest.raises(PairMappingError, match="ambiguously"):
        map_cluster_neighbors(
            structure,
            [{"coords": [0, 0, 0], "element": "Fe"}],
            {0: [{"coords": [0.5, 0.5, 0.5], "element": "Fe"}]},
        )


def test_legacy_serialization_is_explicitly_cartesian() -> None:
    pair = ExchangePair(
        key=PairKey(0, 1, (-1, 0, 0)),
        cartesian_displacement=(-1.5, 0.5, 0.5),
        distance=np.sqrt(2.75),
        species_i="Fe",
        species_j="Pt",
        length_unit="angstrom",
    )
    mapping = legacy_displacement_mapping(pair)
    assert mapping["R_ij_kind"] == "cartesian_displacement"
    assert mapping["R_ij"] == [-1.5, 0.5, 0.5]
    assert mapping["length_unit"] == "angstrom"
    assert mapping["pair_key"] == {"i": 0, "j": 1, "R": [-1, 0, 0]}
