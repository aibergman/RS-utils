"""HT-12 exchange/DMI/tensor result parsing coverage."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import f90nml

from rsutil.core import PairKey, RSStructure, StructureSite
from rsutil.core.hashing import hash_json
from rsutil.results import (
    InteractionMappingError,
    MILLI_RYDBERG,
    MalformedResultError,
    Quantity,
    RYDBERG,
    SourceUnitMismatchError,
    UNKNOWN_UNIT,
    UnsupportedResultFormatError,
    parse_exchange_results,
    write_exchange_json,
)


ROOT = Path(__file__).resolve().parents[1]


def _mnsi_fixture_structure() -> RSStructure:
    # The real export contains the four Mn RS types (1..4) as exchange
    # endpoints.  Read their actual primitive-cell positions from the
    # checked-in lattice fixture and preserve the source labels explicitly.
    lattice = f90nml.read(str(ROOT / "examples" / "export" / "MnSi" / "lattice.nml"))["lattice"]
    positions = lattice["crd"][:4]
    type_indices = tuple(int(value) for value in lattice["izp"][:4])
    labels = {index: f"Mn{index}" for index in type_indices}
    return RSStructure(
        lattice_vectors=tuple(tuple(float(value) for value in row) for row in lattice["a"]),
        length_unit="lattice",
        coordinate_convention="cartesian",
        sites=tuple(
            StructureSite(
                index,
                tuple(float(value) for value in position),
                "Mn",
                rs_label=labels[type_index],
                rs_type_index=type_index,
            )
            for index, (type_index, position) in enumerate(zip(type_indices, positions))
        ),
    )


def test_real_jij_and_dmi_fixtures_map_to_periodic_pair_keys() -> None:
    structure = _mnsi_fixture_structure()
    assert [(site.rs_label, site.rs_type_index) for site in structure.sites] == [
        ("Mn1", 1),
        ("Mn2", 2),
        ("Mn3", 3),
        ("Mn4", 4),
    ]
    result = parse_exchange_results(ROOT / "examples" / "export" / "MnSi", structure=structure)

    assert len(result.records) == 584
    periodic = result.get(PairKey(0, 0, (-2, 0, 0)))
    assert periodic is not None
    assert periodic.J is not None and periodic.J.unit is MILLI_RYDBERG
    assert periodic.J.value == pytest.approx(-0.00587)
    assert periodic.D is not None and periodic.D.unit is MILLI_RYDBERG
    assert periodic.D.value == pytest.approx((0.000054, -0.000041, 0.000138))
    assert periodic.pair is not None
    assert (periodic.pair.label_i, periodic.pair.label_j) == ("Mn1", "Mn1")
    assert (periodic.pair.rs_type_i, periodic.pair.rs_type_j) == (1, 1)
    assert periodic.channel == "so"
    assert result.validation.valid


def test_json_export_round_trips_schema_units_and_pair_identity(tmp_path: Path) -> None:
    result = parse_exchange_results(
        ROOT / "examples" / "export" / "MnSi" / "jijso.out",
        structure=_mnsi_fixture_structure(),
    )
    output = write_exchange_json(result, tmp_path / "interactions.json")
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "rsutil-exchange-result-v1"
    assert payload["records"][0]["key"] == {"i": 0, "j": 0, "R": [-2, 0, 0]}
    assert payload["records"][0]["isotropic_exchange"]["unit"] == "mRy"


def test_explicit_ry_mry_conversion_is_factor_of_1000() -> None:
    source_value = Quantity(-0.00587, MILLI_RYDBERG)

    in_ry = source_value.convert(RYDBERG)

    assert in_ry.value == pytest.approx(-0.00000587)
    assert in_ry.convert(MILLI_RYDBERG).value == pytest.approx(source_value.value)


@pytest.mark.parametrize("filename", ["jijso.out", "dijso.out", "aijso.out"])
def test_known_fable_v4_units_cannot_be_overridden_by_a_contradictory_header(
    tmp_path: Path, filename: str
) -> None:
    row = {
        "jijso.out": "1 1 0 0 0 1 0\n",
        "dijso.out": "1 1 0 0 0 1 2 3 0\n",
        "aijso.out": "1 1 0 0 0 1 2 3 4 5 6 7 8 9 0\n",
    }[filename]
    path = tmp_path / filename
    path.write_text("# unit: Ry\n" + row, encoding="utf-8")
    with pytest.raises(SourceUnitMismatchError) as error:
        parse_exchange_results(
            path,
            manifest={
                "schema_version": "rsutil-exchange-v1",
                "canonical_pair_list": [{"i": 0, "j": 0, "R": [0, 0, 0]}],
            },
        )
    assert error.value.error_code == "SOURCE_UNIT_MISMATCH"


def test_known_fable_v4_unit_header_mry_is_accepted_and_missing_header_is_mry(tmp_path: Path) -> None:
    manifest = {
        "schema_version": "rsutil-exchange-v1",
        "canonical_pair_list": [{"i": 0, "j": 0, "R": [0, 0, 0]}],
    }
    for name, header in (("jijso.out", "# unit: mRy\n"), ("jijfo.out", "")):
        path = tmp_path / name
        path.write_text(header + "1 1 0 0 0 1.25 0\n", encoding="utf-8")
        result = parse_exchange_results(path, manifest=manifest)
        assert result.interactions[0].J is not None
        assert result.interactions[0].J.unit is MILLI_RYDBERG


def test_aij_fortran_column_major_order_is_preserved_and_json_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "aijso.out"
    path.write_text(
        "1 1 0 0 0 1 2 3 4 5 6 7 8 9 0\n",
        encoding="utf-8",
    )
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="bohr",
        coordinate_convention="cartesian",
        sites=(StructureSite(0, (0, 0, 0), "Fe"),),
    )
    result = parse_exchange_results(path, structure=structure)
    record = result.interactions[0]
    assert record.anisotropic_exchange_tensor is not None
    assert record.anisotropic_exchange_tensor.value == (
        (1.0, 4.0, 7.0),
        (2.0, 5.0, 8.0),
        (3.0, 6.0, 9.0),
    )
    assert record.pair is not None and record.pair.length_unit == "bohr"
    payload = result.to_dict()["records"][0]
    assert "anisotropic_exchange_tensor" in payload
    assert "exchange_tensor" not in payload
    assert payload["pair"]["length_unit"] == "bohr"
    restored = type(result).from_json(result.to_json())
    assert restored.interactions[0].anisotropic_exchange_tensor is not None
    assert restored.interactions[0].pair is not None
    assert restored.interactions[0].pair.length_unit == "bohr"
    assert np.allclose(
        restored.interactions[0].anisotropic_exchange_tensor.value,
        record.anisotropic_exchange_tensor.value,
    )


def test_full_exchange_tensor_uses_fable_v4_j_d_a_convention(tmp_path: Path) -> None:
    files = {
        "jijso.out": "1 1 0 0 0 2 0\n",
        "dijso.out": "1 1 0 0 0 0.1 0.2 0.3 0\n",
        "aijso.out": "1 1 0 0 0 1 2 3 4 5 6 7 8 9 0\n",
    }
    paths = []
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="angstrom",
        sites=(StructureSite(0, (0, 0, 0), "Fe"),),
    )
    for name, contents in files.items():
        path = tmp_path / name
        path.write_text(contents, encoding="utf-8")
        paths.append(path)
    record = parse_exchange_results(paths, structure=structure).interactions[0]
    full = record.require_full_exchange_tensor()
    assert full.unit is MILLI_RYDBERG
    assert np.allclose(full.value, ((3.0, 4.3, 6.8), (1.7, 7.0, 8.1), (3.2, 5.9, 11.0)))


def test_full_exchange_tensor_is_unavailable_without_all_known_j_d_a_parts() -> None:
    from rsutil.results import InteractionRecord

    record = InteractionRecord(
        PairKey(0, 0, (0, 0, 0)),
        isotropic_exchange=Quantity(1.0, MILLI_RYDBERG),
    )
    assert record.full_exchange_tensor() is None
    with pytest.raises(ValueError, match="full exchange tensor is unavailable"):
        record.require_full_exchange_tensor()


def test_mismatched_manifest_fails_instead_of_guessing(tmp_path: Path) -> None:
    source = ROOT / "examples" / "export" / "MnSi" / "jijso.out"
    table = tmp_path / "jijso.out"
    table.write_text(source.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "rsutil-exchange-v1",
        "canonical_pair_list": [{"i": 0, "j": 1, "R": [1, 0, 0]}],
    }
    with pytest.raises(InteractionMappingError):
        parse_exchange_results(table, manifest=manifest, structure=_mnsi_fixture_structure())


def test_manifest_mapping_accepts_periodic_image_without_using_row_order(tmp_path: Path) -> None:
    table = tmp_path / "jijso.out"
    table.write_text("1 1 -2 0 0 -0.005870 2\n", encoding="utf-8")
    manifest = {
        "schema_version": "rsutil-exchange-v1",
        "canonical_pair_list": [{"i": 0, "j": 0, "R": [-2, 0, 0]}],
    }
    result = parse_exchange_results(table, manifest, _mnsi_fixture_structure())
    assert result.pair_keys == (PairKey(0, 0, (-2, 0, 0)),)


def test_parser_rejects_tampered_duplicate_pair_view(tmp_path: Path) -> None:
    table = tmp_path / "jijso.out"
    table.write_text("1 1 -2 0 0 -0.005870 2\n", encoding="utf-8")
    manifest = {
        "schema_version": "rsutil-exchange-v2",
        "spec": {
            "pair_keys": [{"i": 0, "j": 0, "R": [-2, 0, 0]}],
            "legacy_ijpair": [[1], [1]],
        },
        "spec_hash": hash_json({
            "pair_keys": [{"i": 0, "j": 0, "R": [-2, 0, 0]}],
            "legacy_ijpair": [[1], [1]],
        }),
        "canonical_pair_list": [{"i": 0, "j": 0, "R": [2, 0, 0]}],
    }
    with pytest.raises(MalformedResultError, match="canonical_pair_list"):
        parse_exchange_results(table, manifest=manifest, structure=_mnsi_fixture_structure())


def test_tensor_rows_require_explicit_tensor_format_and_keep_field_distinct(tmp_path: Path) -> None:
    path = tmp_path / "aijso.out"
    path.write_text(
        "# itype jtype r_ij J_ij tensor\n"
        "1 1 1 0 0 1 2 3 4 5 6 7 8 9 1\n",
        encoding="utf-8",
    )
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="lattice",
        sites=(StructureSite(0, (0, 0, 0), "Fe"),),
    )
    result = parse_exchange_results(path, structure=structure)
    record = result.get(PairKey(0, 0, (1, 0, 0)))
    assert record is not None
    assert record.J is None and record.D is None
    assert record.tensor is not None
    assert record.full_exchange_tensor() is None
    assert record.tensor.unit is MILLI_RYDBERG
    assert record.tensor.value == ((1.0, 4.0, 7.0), (2.0, 5.0, 8.0), (3.0, 6.0, 9.0))


def test_unverified_tensor_source_keeps_unit_unknown(tmp_path: Path) -> None:
    path = tmp_path / "exchange_tensor.out"
    path.write_text(
        "# itype jtype r_ij J_ij tensor\n"
        "1 1 1 0 0 1 2 3 4 5 6 7 8 9 1\n",
        encoding="utf-8",
    )
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="lattice",
        coordinate_convention="cartesian",
        sites=(StructureSite(0, (0, 0, 0), "Fe"),),
    )

    result = parse_exchange_results(path, structure=structure)
    record = result.get(PairKey(0, 0, (1, 0, 0)))

    assert record is not None and record.tensor is not None
    assert record.tensor.unit is UNKNOWN_UNIT
    assert record.provenance["energy_unit"] == "unknown"


def test_reverse_pair_validation_checks_available_directions_only(tmp_path: Path) -> None:
    path = tmp_path / "jijso.out"
    path.write_text(
        "1 1 1 0 0 0.25 1\n"
        "1 1 -1 0 0 0.25 1\n",
        encoding="utf-8",
    )
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="lattice",
        sites=(StructureSite(0, (0, 0, 0), "Fe"),),
    )
    result = parse_exchange_results(path, structure=structure)
    assert result.validation.valid
    assert not result.validation.diagnostics


def test_malformed_and_unknown_formats_fail_explicitly(tmp_path: Path) -> None:
    malformed = tmp_path / "jijso.out"
    malformed.write_text("1 1 0 0 0 not-a-number 0\n", encoding="utf-8")
    with pytest.raises(MalformedResultError):
        parse_exchange_results(malformed)

    unknown = tmp_path / "result.dat"
    unknown.write_text("1 1 0 0 0 1 1\n", encoding="utf-8")
    with pytest.raises(UnsupportedResultFormatError, match="identify supported interaction format"):
        parse_exchange_results(unknown)
