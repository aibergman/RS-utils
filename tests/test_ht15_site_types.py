"""HT-R15 source-faithful RS site/type identity controls."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rsutil.core import PairKey, RSStructure, StructureSite
from rsutil.core.hashing import hash_json
from rsutil.exchange.reference import build_rs_structure, build_structure_and_analyzer
from rsutil.results import InteractionMappingError, parse_exchange_results


def _lattice(crd, izp):
    return {
        "lattice": {
            "a": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "crd": crd,
            "izp": izp,
        }
    }


def _manifest(records: list[dict]) -> dict:
    keys = [record["key"] for record in records]
    spec = {
        "pair_keys": keys,
        "pairs": records,
        "legacy_ijpair": [[1] * len(records), [1] * len(records)],
    }
    return {
        "schema_version": "rsutil-exchange-v3",
        "spec": spec,
        "spec_hash": hash_json(spec),
        "canonical_pair_list": keys,
        "pair_keys": keys,
    }


def test_two_primitive_sites_can_share_one_rs_type() -> None:
    lattice = _lattice([[0, 0, 0], [0.5, 0.5, 0.5]], [1, 1])
    input_nml = {"atoms": {"label": ["Co"]}}

    canonical = build_rs_structure(lattice, input_nml)
    pymatgen, _analyzer = build_structure_and_analyzer(lattice, input_nml)

    assert len(canonical.sites) == 2
    assert [site.element for site in canonical.sites] == ["Co", "Co"]
    assert [site.rs_label for site in canonical.sites] == ["Co", "Co"]
    assert [site.rs_type_index for site in canonical.sites] == [1, 1]
    assert len(pymatgen) == 2
    assert pymatgen.site_properties["rs_type_index"] == [1, 1]


def test_same_chemistry_distinct_rs_types_are_not_collapsed() -> None:
    lattice = _lattice(
        [[0, 0, 0], [0.2, 0.2, 0.2], [0.4, 0.4, 0.4], [0.6, 0.6, 0.6]],
        [1, 2, 3, 4],
    )
    input_nml = {"atoms": {"label": ["Mn1", "Mn2", "Mn3", "Mn4"]}}

    canonical = build_rs_structure(lattice, input_nml)

    assert [site.element for site in canonical.sites] == ["Mn"] * 4
    assert [site.rs_label for site in canonical.sites] == ["Mn1", "Mn2", "Mn3", "Mn4"]
    assert [site.rs_type_index for site in canonical.sites] == [1, 2, 3, 4]


def test_exchange_type_columns_map_to_manifest_sites_by_type_and_geometry(tmp_path: Path) -> None:
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="angstrom",
        sites=(
            StructureSite(0, (0, 0, 0), "Co", rs_label="Co", rs_type_index=1),
            StructureSite(1, (0.5, 0.5, 0.5), "Co", rs_label="Co", rs_type_index=1),
        ),
    )
    record = {
        "key": {"i": 0, "j": 1, "R": [0, 0, 0]},
        "rs_type_i": 1,
        "rs_type_j": 1,
        "cartesian_displacement": [0.5, 0.5, 0.5],
        "distance": 3**0.5 / 2,
        "length_unit": "angstrom",
    }
    table = tmp_path / "jijso.out"
    table.write_text("1 1 0.5 0.5 0.5 2.0 0.8660254\n", encoding="utf-8")

    result = parse_exchange_results(table, manifest=_manifest([record]), structure=structure)

    assert result.pair_keys == (PairKey(0, 1, (0, 0, 0)),)

    collected = parse_exchange_results(table, manifest=_manifest([record]))
    assert collected.pair_keys == (PairKey(0, 1, (0, 0, 0)),)


def test_equal_type_equal_geometry_mapping_fails_closed(tmp_path: Path) -> None:
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="angstrom",
        sites=(
            StructureSite(0, (0, 0, 0), "Co", rs_type_index=1),
            StructureSite(1, (0, 0, 0), "Co", rs_type_index=1),
        ),
    )
    records = [
        {
            "key": {"i": 0, "j": 0, "R": [0, 0, 0]},
            "rs_type_i": 1, "rs_type_j": 1,
            "cartesian_displacement": [0, 0, 0], "distance": 0.0,
        },
        {
            "key": {"i": 1, "j": 1, "R": [0, 0, 0]},
            "rs_type_i": 1, "rs_type_j": 1,
            "cartesian_displacement": [0, 0, 0], "distance": 0.0,
        },
    ]
    table = tmp_path / "jijso.out"
    table.write_text("1 1 0 0 0 2.0 0\n", encoding="utf-8")

    with pytest.raises(InteractionMappingError, match="ambiguously"):
        parse_exchange_results(table, manifest=_manifest(records), structure=structure)


def test_structure_site_type_index_json_round_trip() -> None:
    site = StructureSite(0, (0, 0, 0), "Mn", rs_label="Mn3", rs_type_index=3)
    restored = StructureSite.from_dict(json.loads(json.dumps(site.to_dict())))
    assert restored == site
