"""HT-06 semantic validation, provenance, and report-schema coverage."""

from __future__ import annotations

import json
from pathlib import Path

from rsutil.core import ErrorCode, SetupReport, SetupSpec
from rsutil.setup import SetupService, setup_report, validate_setup, validate_setup_outputs
from rsutil.setup_validation import write_setup_report

from .conftest import copy_example_files, run_rsutil


def _fixture(destination: Path, example: str = "bccFe") -> None:
    names = ["lattice.nml", "input.nml"] if example != "bccFe" else ["lattice.nml", "input.nml"]
    copy_example_files(destination, example, names)


def _codes(result):
    return {diagnostic.code for diagnostic in result.diagnostics}


def test_baseline_fixtures_are_semantically_valid(tmp_path: Path) -> None:
    for example in ("bccFe", "FeN", "FePt"):
        directory = tmp_path / example
        directory.mkdir()
        _fixture(directory, example)
        validation = validate_setup(directory / "lattice.nml", directory / "input.nml")
        assert validation.valid, (example, validation.to_dict())


def test_inconsistent_ntype_is_coded_as_a_semantic_error(tmp_path: Path) -> None:
    _fixture(tmp_path, "FeN")
    path = tmp_path / "input.nml"
    path.write_text(path.read_text(encoding="utf-8").replace("ntype       = 2", "ntype       = 1"), encoding="utf-8")

    validation = validate_setup(tmp_path / "lattice.nml", path)

    assert not validation.valid
    assert ErrorCode.LABEL_COUNT_MISMATCH.value in _codes(validation)


def test_near_singular_lattice_is_rejected_before_worker_use(tmp_path: Path) -> None:
    _fixture(tmp_path)
    path = tmp_path / "lattice.nml"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        "a(:, 3) = 0.50000000, 0.50000000,-0.50000000",
        "a(:, 3) =-0.50000000, 0.50000000, 0.50000000",
    )
    path.write_text(content, encoding="utf-8")

    validation = validate_setup(path, tmp_path / "input.nml")

    assert not validation.valid
    assert ErrorCode.LATTICE_SINGULAR.value in _codes(validation)


def test_invalid_spec_and_missing_outputs_are_machine_coded(tmp_path: Path) -> None:
    _fixture(tmp_path)
    validation = validate_setup(
        tmp_path / "lattice.nml",
        tmp_path / "input.nml",
        spec={"rec_method": "not-a-method"},
    )
    assert not validation.valid
    assert ErrorCode.SETUP_SPEC_INVALID.value in _codes(validation)

    output_validation = validate_setup_outputs(tmp_path)
    assert not output_validation.valid
    assert ErrorCode.OUTPUT_MISSING.value in _codes(output_validation)


def test_bad_supercell_and_2d_options_fail_as_spec_diagnostics(tmp_path: Path) -> None:
    _fixture(tmp_path)
    for invalid_spec in (
        {"supercell": [0, 1, 1]},
        {"z_iso_factor": 2.0},
    ):
        validation = validate_setup(
            tmp_path / "lattice.nml",
            tmp_path / "input.nml",
            spec=invalid_spec,
        )
        assert not validation.valid
        assert ErrorCode.SETUP_SPEC_INVALID.value in _codes(validation)


def test_generated_recursion_and_sphere_records_are_checked(tmp_path: Path) -> None:
    _fixture(tmp_path)
    result = SetupService().run(tmp_path, tmp_path, SetupSpec(seed=42))
    assert result.success, result.to_dict()

    generated = tmp_path / "input_new.nml"
    generated.write_text(
        generated.read_text(encoding="utf-8").replace("recur = 'block'", "recur = 'invalid'"),
        encoding="utf-8",
    )
    validation = validate_setup_outputs(tmp_path, spec=SetupSpec(seed=42))
    assert not validation.valid
    assert ErrorCode.RECURSION_INVALID.value in _codes(validation)


def test_report_schema_and_scientific_hash_are_stable(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    _fixture(first)
    _fixture(second)
    spec = SetupSpec(seed=42)

    first_result = SetupService().run(first, first, spec)
    second_result = SetupService().run(second, second, spec)
    assert first_result.success and second_result.success

    first_report = setup_report(first_result)
    second_report = setup_report(second_result)
    first_payload = first_report.to_dict()
    assert first_payload["schema_version"] == "rsutil-setup-report-v1"
    assert first_payload["manifest"]["schema_version"] == "rsutil-run-v1"
    assert first_payload["result"]["semantic_validation"]["schema_version"] == "rsutil-validation-v1"
    assert first_payload["manifest"]["setup_identity_hash"] == second_report.manifest.setup_identity_hash
    assert first_payload["manifest"]["structure_hash"] == second_report.manifest.structure_hash
    assert first_payload["manifest"]["stochastic"]["reproducible"] is True
    assert first_payload["manifest"]["timestamps"]
    assert SetupReport.from_dict(first_payload).manifest.setup_identity_hash == first_report.manifest.setup_identity_hash

    report_path = tmp_path / "report.json"
    write_setup_report(first_report, report_path)
    assert json.loads(report_path.read_text(encoding="utf-8"))["schema_version"] == "rsutil-setup-report-v1"


def test_legacy_cli_can_optionally_write_json_report(tmp_path: Path) -> None:
    _fixture(tmp_path)
    result = run_rsutil(tmp_path, "-S", "--json-report", "setup-report.json")
    assert result.returncode == 0, result.stderr
    payload = json.loads((tmp_path / "setup-report.json").read_text(encoding="utf-8"))
    assert payload["schema_version"] == "rsutil-setup-report-v1"
    assert payload["result"]["success"] is True
