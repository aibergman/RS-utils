"""HT-12 exchange/DMI/tensor result parsing coverage."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rsutil.core import PairKey, RSStructure, StructureSite
from rsutil.results import (
    InteractionMappingError,
    MILLI_RYDBERG,
    MalformedResultError,
    Quantity,
    RYDBERG,
    UNKNOWN_UNIT,
    UnsupportedResultFormatError,
    parse_exchange_results,
    write_exchange_json,
)


ROOT = Path(__file__).resolve().parents[1]


def _mnsi_fixture_structure() -> RSStructure:
    # The positions are the four-site Cartesian reference used by the real
    # MnSi result rows; the table itself is intentionally not synthetic.
    return RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="lattice",
        coordinate_convention="cartesian",
        sites=(
            StructureSite(0, (0, 0, 0), "Mn"),
            StructureSite(1, (-0.5, 0.2238, -0.2762), "Si"),
            StructureSite(2, (0.2238, -0.2762, -0.5), "Mn"),
            StructureSite(3, (-0.2762, -0.5, 0.2238), "Si"),
        ),
    )


def test_real_jij_and_dmi_fixtures_map_to_periodic_pair_keys() -> None:
    result = parse_exchange_results(ROOT / "examples" / "export" / "MnSi", structure=_mnsi_fixture_structure())

    assert len(result.records) == 584
    periodic = result.get(PairKey(0, 0, (-2, 0, 0)))
    assert periodic is not None
    assert periodic.J is not None and periodic.J.unit is MILLI_RYDBERG
    assert periodic.J.value == pytest.approx(-0.00587)
    assert periodic.D is not None and periodic.D.unit is MILLI_RYDBERG
    assert periodic.D.value == pytest.approx((0.000054, -0.000041, 0.000138))
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
    assert record.tensor.unit is MILLI_RYDBERG
    assert record.tensor.value == ((1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0))


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
