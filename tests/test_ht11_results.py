"""HT-11 coverage for versioned SCF result parsing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rsutil.exchange.fermi import parse_fermi_energy
from rsutil.results import (
    ConvergenceState,
    MalformedResultError,
    MissingResultFileError,
    RYDBERG,
    SCFResult,
    parse_run,
)


ROOT = Path(__file__).resolve().parents[1]


def test_parse_two_repository_runs_and_retain_raw_provenance() -> None:
    b2fe = parse_run(ROOT / "examples" / "B2Fe")
    fept = parse_run(ROOT / "examples" / "FePt")

    assert b2fe.convergence is ConvergenceState.NOT_CONVERGED
    assert fept.convergence is ConvergenceState.NOT_CONVERGED
    assert b2fe.fermi_energy is not None
    assert b2fe.fermi_energy.value == pytest.approx(-0.671411)
    assert b2fe.fermi_energy.unit is RYDBERG
    assert fept.fermi_energy is not None
    assert fept.fermi_energy.value == pytest.approx(-0.159456)
    assert fept.total_energy is not None
    assert fept.total_energy.unit.value == "unknown"
    assert [item.site_index for item in b2fe.site_spin_moments] == [0, 1]
    assert b2fe.site_spin_moments[0].orientation == pytest.approx((0.0, 0.0, -1.0))
    assert b2fe.charge_transfers[0].quantity.unit.value == "electron"
    assert b2fe.provenance.source_path.endswith("examples/B2Fe/report.out")
    assert b2fe.provenance.content_hash
    assert {item.role for item in b2fe.provenance.source_files} >= {
        "report",
        "run_log",
        "spin_orientation",
    }


def test_fermi_compatibility_reader_uses_structured_result() -> None:
    assert parse_fermi_energy(ROOT / "examples" / "B2Fe" / "report.out") == pytest.approx(-0.671411)
    assert parse_fermi_energy(ROOT / "examples" / "FePt" / "report.out") == pytest.approx(-0.159456)


def test_explicit_path_is_independent_of_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    result = parse_run(ROOT / "examples" / "B2Fe" / "report.out")
    assert result.fermi == pytest.approx(-0.671411)


def test_report_without_convergence_evidence_is_unknown(tmp_path: Path) -> None:
    report = (ROOT / "examples" / "B2Fe" / "report.out").read_text(encoding="utf-8")
    truncated = report[: report.index("Fermi energy:") + len("Fermi energy:  -0.671411")]
    path = tmp_path / "report.out"
    path.write_text(truncated, encoding="utf-8")

    result = parse_run(path)
    assert result.convergence is ConvergenceState.UNKNOWN
    assert result.converged is None
    assert result.validation.valid


def test_malformed_or_missing_fermi_does_not_parse_as_converged(tmp_path: Path) -> None:
    path = tmp_path / "report.out"
    path.write_text("Total energy of system: -1.0\n", encoding="utf-8")
    with pytest.raises(MissingResultFileError):
        parse_run(path)

    path.write_text("Fermi energy: not-a-number\n", encoding="utf-8")
    with pytest.raises(MalformedResultError):
        parse_run(path)


def test_scf_scalar_and_array_metadata_round_trip(tmp_path: Path) -> None:
    report = """\
Total energy of system: -1.250000
Band energy of system: -0.500000
Total spin moment: 0.0 0.0 2.0
Spin moment of atom 1: 2.0
Spin moment projections of atom 1: 0.0 0.0 2.0
Occupation at atom 1: 8.0
Up orbital occupation at atom 1: 1.0 2.0 5.0
Down orbital occupation at atom 1: 0.5 1.5 2.5
Charge transfer at atom 1: 0.0
Fermi energy: -0.125000
"""
    path = tmp_path / "report.out"
    path.write_text(report, encoding="utf-8")
    (tmp_path / "out.scf").write_text("Converged!0.0000000001\n", encoding="utf-8")

    result = parse_run(tmp_path)
    assert result.convergence is ConvergenceState.CONVERGED
    assert result.site_up_orbital_occupations[0].quantity.value == pytest.approx((1.0, 2.0, 5.0))
    assert result.site_down_orbital_occupations[0].quantity.value == pytest.approx((0.5, 1.5, 2.5))
    restored = SCFResult.from_json(result.to_json())
    assert restored.to_dict() == json.loads(result.to_json())
