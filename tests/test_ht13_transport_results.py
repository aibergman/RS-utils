"""HT-13 evidence-first transport/damping result contracts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from rsutil.results import (
    ConductivityResult,
    DampingResult,
    DIMENSIONLESS,
    MissingResultFileError,
    PARSER_REGISTRY,
    ParserProvenance,
    Quantity,
    RYDBERG,
    ResultAxis,
    TransportSupport,
    parse_conductivity_result,
    parse_damping_result,
)


ROOT = Path(__file__).resolve().parents[1]
AMBIGUOUS = ROOT / "tests" / "fixtures" / "transport" / "ambiguous_transport.out"


def _provenance() -> ParserProvenance:
    return ParserProvenance(source_path="fixture", content_hash="fixture-hash", format="fixture")


def test_ambiguous_fixture_is_explicitly_unsupported_and_hashed() -> None:
    result = parse_conductivity_result(AMBIGUOUS)

    assert result.status is TransportSupport.UNSUPPORTED_AMBIGUOUS
    assert not result.supported
    assert result.observable.value == "conductivity"
    assert result.provenance.content_hash == hashlib.sha256(AMBIGUOUS.read_bytes()).hexdigest()
    assert result.provenance.source_files[0].role == "transport_candidate"
    assert result.to_dict()["status"] == "UNSUPPORTED_AMBIGUOUS"


def test_damping_alias_does_not_reinterpret_an_ambiguous_file() -> None:
    result = parse_damping_result(AMBIGUOUS)

    assert result.status is TransportSupport.UNSUPPORTED_AMBIGUOUS
    assert "no values were parsed" in result.reason


@pytest.mark.parametrize("name", ["cond_total.out", "Fe_cond.out", "cond_total_orb_real.out"])
def test_fable_v4_conductivity_filenames_are_detected_but_not_guessed(tmp_path: Path, name: str) -> None:
    (tmp_path / name).write_text("0.0 1.0 2.0\n", encoding="utf-8")

    result = parse_conductivity_result(tmp_path)

    assert result.status is TransportSupport.UNSUPPORTED_AMBIGUOUS
    assert not result.supported
    assert "no values were parsed" in result.reason


@pytest.mark.parametrize("name", ["damping-energy.out", "alldampings.out"])
def test_fable_v4_damping_filenames_are_detected_but_not_guessed(tmp_path: Path, name: str) -> None:
    (tmp_path / name).write_text("# ambiguous damping output\n", encoding="utf-8")

    result = parse_damping_result(tmp_path)

    assert result.status is TransportSupport.UNSUPPORTED_AMBIGUOUS
    assert not result.supported
    assert "no values were parsed" in result.reason


def test_registry_records_source_audited_unsupported_evidence() -> None:
    evidence = {spec.observable.value: spec.evidence for spec in PARSER_REGISTRY}

    assert "conductivity.f90" in evidence["conductivity"]
    assert "exchange_dynamics.f90" in evidence["damping"]
    assert all(spec.status is TransportSupport.UNSUPPORTED_AMBIGUOUS for spec in PARSER_REGISTRY)


def test_missing_transport_file_is_not_reported_as_a_partial_result(tmp_path: Path) -> None:
    with pytest.raises(MissingResultFileError):
        parse_conductivity_result(tmp_path)


def test_conductivity_contract_checks_axes_units_and_component_order() -> None:
    result = ConductivityResult(
        values=Quantity(((1.0, 2.0, 3.0), (4.0, 5.0, 6.0)), "S/m"),
        axes=(ResultAxis("fermi_energy", (-0.5, 0.5), RYDBERG, "fermi"),),
        component_order=("xx", "xy", "xz"),
        method="future-kubo",
        channel="charge, non-SOC",
        approximation="future documented approximation",
        provenance=_provenance(),
        options={"broadening": 0.01, "soc": False},
    )

    restored = ConductivityResult.from_dict(json.loads(json.dumps(result.to_dict())))
    assert restored.values.unit.value == "S/m"
    assert restored.axes[0].role == "fermi"
    assert restored.component_order == ("xx", "xy", "xz")
    assert restored.options["soc"] is False


def test_damping_contract_checks_site_tensor_shape_and_indexing() -> None:
    result = DampingResult(
        values=Quantity((((1.0, 0.0), (0.0, 1.0)),), DIMENSIONLESS),
        axes=(ResultAxis("energy", (0.0,), RYDBERG, "energy"),),
        index_kind="site",
        index_labels=("Fe-1", "Fe-2"),
        component_order=("xx", "yy"),
        method="future-torque-correlation",
        channel="spin, SOC",
        approximation="future documented approximation",
        provenance=_provenance(),
    )

    assert DampingResult.from_dict(result.to_dict()).index_kind == "site"
    assert result.to_dict()["component_order"] == ["xx", "yy"]
