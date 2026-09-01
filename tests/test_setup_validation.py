"""HT-06 semantic validation, provenance, and report-schema coverage."""

from __future__ import annotations

import json
from datetime import datetime
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


def _surface_like_fixture(destination: Path) -> None:
    """Build the source-derived fable_v4 surface count fixture.

    fable_v4 ``build_clusup`` derives ntype=nbulk_bulk+nlay, nrec=nlay,
    and keeps nbas as an independent charge-distribution count.  The public
    surface example does not materialize a lattice.nml, so this is the
    smallest local fixture expressing those source semantics.
    """

    _fixture(destination)
    lattice = destination / "lattice.nml"
    lattice_text = lattice.read_text(encoding="utf-8")
    lattice_text = lattice_text.replace("nbas = 1", "nbas = 51").replace("nrec = 1", "nrec = 6")
    lattice_text = lattice_text.replace(
        "irec(1) = 1\n!\nct(1) = 3.3",
        "irec(1) = 1\nirec(2) = 2\nirec(3) = 3\nirec(4) = 4\nirec(5) = 5\nirec(6) = 6\n!\nct(1) = 3.3",
    )
    lattice.write_text(lattice_text, encoding="utf-8")

    input_path = destination / "input.nml"
    input_text = input_path.read_text(encoding="utf-8")
    input_text = input_text.replace("crystal_sym = 'file'", "crystal_sym = 'file'\nntype = 7\nnlay = 6")
    input_text = input_text.replace(
        "label(1) = 'Fe'",
        "\n".join(
            [
                "label(1) = 'Cu'",
                "label(2) = 'ES'",
                "label(3) = 'Cu-S'",
                "label(4) = 'Cu-S-1'",
                "label(5) = 'Cu-S-2'",
                "label(6) = 'Cu-S-3'",
                "label(7) = 'Cu-S-4'",
            ]
        ),
    ).replace("calctype = 'B'", "calctype = 'S'")
    input_path.write_text(input_text, encoding="utf-8")


def _impurity_like_fixture(destination: Path) -> None:
    """Build the public fable_v4 impurity label/count relationship."""

    _fixture(destination, "FeN")
    input_path = destination / "input.nml"
    input_text = input_path.read_text(encoding="utf-8")
    input_text = input_text.replace("ntype       = 2\n", "")
    input_text = input_text.replace(
        "label(1) = 'Fe1'\nlabel(2) = 'N1'",
        "label(1) = 'Co1-bulk'\nlabel(2) = 'Co2-bulk'\nlabel(3) = 'Fe-imp'",
    ).replace("calctype = 'B'", "calctype = 'I'")
    input_path.write_text(input_text, encoding="utf-8")


def test_baseline_fixtures_are_semantically_valid(tmp_path: Path) -> None:
    for example in ("bccFe", "FeN", "FePt"):
        directory = tmp_path / example
        directory.mkdir()
        _fixture(directory, example)
        validation = validate_setup(directory / "lattice.nml", directory / "input.nml")
        assert validation.valid, (example, validation.to_dict())


def test_source_derived_surface_counts_can_differ(tmp_path: Path) -> None:
    _surface_like_fixture(tmp_path)

    validation = validate_setup(tmp_path / "lattice.nml", tmp_path / "input.nml")

    assert validation.valid, validation.to_dict()
    assert not any(
        "does not match nrec" in diagnostic.message
        for diagnostic in validation.diagnostics
    )


def test_source_derived_interface_counts_are_mode_aware(tmp_path: Path) -> None:
    _surface_like_fixture(tmp_path)
    path = tmp_path / "input.nml"
    input_text = path.read_text(encoding="utf-8")
    input_text = input_text.replace("ntype = 7", "ntype = 8")
    input_text = input_text.replace("label(7) = 'Cu-S-4'", "label(7) = 'Cu-S-4'\nlabel(8) = 'Cu-S-5'")
    input_text = input_text.replace("calctype = 'S'", "calctype = 'L'")
    path.write_text(input_text, encoding="utf-8")

    validation = validate_setup(tmp_path / "lattice.nml", path)

    assert validation.valid, validation.to_dict()


def test_impurity_labels_are_not_limited_by_recursion_count(tmp_path: Path) -> None:
    _impurity_like_fixture(tmp_path)

    validation = validate_setup(tmp_path / "lattice.nml", tmp_path / "input.nml")

    assert validation.valid, validation.to_dict()


def test_recursion_array_count_is_still_rejected(tmp_path: Path) -> None:
    _surface_like_fixture(tmp_path)
    path = tmp_path / "lattice.nml"
    path.write_text(
        path.read_text(encoding="utf-8").replace("irec(6) = 6\n", ""),
        encoding="utf-8",
    )

    validation = validate_setup(path, tmp_path / "input.nml")

    assert not validation.valid
    assert ErrorCode.COUNT_MISMATCH.value in _codes(validation)


def test_surface_recursion_count_must_match_active_layers(tmp_path: Path) -> None:
    _surface_like_fixture(tmp_path)
    path = tmp_path / "lattice.nml"
    path.write_text(
        path.read_text(encoding="utf-8").replace("nrec = 6", "nrec = 5"),
        encoding="utf-8",
    )

    validation = validate_setup(path, tmp_path / "input.nml")

    assert not validation.valid
    assert ErrorCode.COUNT_MISMATCH.value in _codes(validation)


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


def test_deterministic_setup_without_seed_is_reproducible(tmp_path: Path) -> None:
    _fixture(tmp_path)
    validation = validate_setup(
        tmp_path / "lattice.nml", tmp_path / "input.nml", spec=SetupSpec()
    )
    assert validation.valid
    assert ErrorCode.UNSEEDED_SETUP.value not in _codes(validation)

    manifest = setup_report(SetupService().run(tmp_path, tmp_path, SetupSpec())).manifest
    assert manifest.stochastic["active"] is False
    assert manifest.stochastic["operations"] == []
    assert manifest.stochastic["reproducible"] is True
    assert manifest.stochastic["reason"] == "no stochastic operations requested"
    timestamps = manifest.timestamps
    assert datetime.fromisoformat(timestamps["started_at"].replace("Z", "+00:00")) <= datetime.fromisoformat(
        timestamps["finished_at"].replace("Z", "+00:00")
    )

    stochastic = SetupSpec(lattice_dist=0.1)
    stochastic_manifest = setup_report(SetupService().run(tmp_path, tmp_path, stochastic)).manifest
    assert stochastic_manifest.stochastic["active"] is True
    assert stochastic_manifest.stochastic["operations"] == ["lattice_distortion"]
    assert stochastic_manifest.stochastic["reproducible"] is False
    assert ErrorCode.UNSEEDED_SETUP.value in {
        item.code for item in validate_setup(
            tmp_path / "lattice.nml", tmp_path / "input.nml", spec=stochastic
        ).diagnostics
    }


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
