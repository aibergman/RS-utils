"""HT-04/HT-05 coverage for the path-explicit setup service."""

from __future__ import annotations

import os
from pathlib import Path

import f90nml
import pytest

from rsutil.core import (
    Diagnostic,
    DiagnosticSeverity,
    ErrorCode,
    OperationStatus,
    SetupResult,
    SetupSpec,
    ValidationResult,
    hash_file,
)
from rsutil.setup import SetupService

from .conftest import copy_example_files


def _fixture(destination: Path) -> None:
    copy_example_files(destination, "bccFe", ["Dist/lattice.nml", "Dist/input.nml"])


class _FakeRunner:
    def __init__(self, mode: str = "success") -> None:
        self.mode = mode

    def run(self, lattice_path, input_path, spec, *, workspace):
        if self.mode == "worker-failure":
            return SetupResult(
                status=OperationStatus.FAILED,
                diagnostics=(
                    Diagnostic(
                        DiagnosticSeverity.ERROR,
                        ErrorCode.WORKER_FAILED,
                        "injected worker failure",
                    ),
                ),
            )
        if self.mode == "validation-failure":
            return SetupResult(
                status=OperationStatus.FAILED,
                diagnostics=(
                    Diagnostic(
                        DiagnosticSeverity.ERROR,
                        ErrorCode.OUTPUT_VALIDATION_FAILED,
                        "injected validation failure",
                    ),
                ),
                semantic_validation=ValidationResult(valid=False),
            )
        workspace = Path(workspace)
        workspace.mkdir(parents=True, exist_ok=True)
        outputs = {
            "lattice_new.nml": "&lattice\n ntot = 1\n/\n",
            "input_new.nml": "&atoms\n label(1) = 'Fe'\n/\n",
            "radii.out": "1\n1.0\n",
        }
        for name, content in outputs.items():
            (workspace / name).write_text(content, encoding="utf-8")
        return SetupResult(
            status=OperationStatus.SUCCESS,
            candidate_files=tuple(outputs),
            semantic_validation=ValidationResult(valid=True),
        )


def test_setup_service_does_not_use_process_cwd(tmp_path: Path) -> None:
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    unrelated = tmp_path / "unrelated"
    source.mkdir()
    destination.mkdir()
    unrelated.mkdir()
    _fixture(source)
    original = {name: (source / name).read_bytes() for name in ("lattice.nml", "input.nml")}

    previous = Path.cwd()
    os.chdir(unrelated)
    try:
        result = SetupService().run(source, destination, SetupSpec())
    finally:
        os.chdir(previous)

    assert result.success, result.to_dict()
    assert result.provenance["source_dir"] == str(source.resolve())
    assert result.provenance["destination_dir"] == str(destination.resolve())
    assert (destination / "lattice_new.nml").is_file()
    assert (destination / "input_new.nml").is_file()
    assert {name: (source / name).read_bytes() for name in original} == original


def test_setup_service_supports_source_and_destination_paths_with_spaces(tmp_path: Path) -> None:
    source = tmp_path / "source with spaces"
    destination = tmp_path / "destination with spaces"
    source.mkdir()
    _fixture(source)

    result = SetupService().run(source, destination, SetupSpec(seed=42))

    assert result.success, result.to_dict()
    assert (destination / "lattice_new.nml").is_file()
    assert (destination / "input_new.nml").is_file()
    assert result.provenance["destination_dir"] == str(destination.resolve())


def test_setup_service_preserves_seeded_supercell_sqs(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    _fixture(first)
    _fixture(second)
    spec = SetupSpec(
        supercell=(2, 1, 1),
        sqs="Fe=0.5,Co=0.5",
        sqs_steps=10,
        sqs_seed=42,
    )

    first_result = SetupService().run(first, first, spec)
    second_result = SetupService().run(second, second, spec)

    assert first_result.success, first_result.to_dict()
    assert second_result.success, second_result.to_dict()
    assert f90nml.read(first / "lattice_new.nml") == f90nml.read(second / "lattice_new.nml")
    assert f90nml.read(first / "input_new.nml") == f90nml.read(second / "input_new.nml")
    assert not (first / ".rsutil_supercell_lattice.nml").exists()
    assert not (first / ".rsutil_supercell_input.nml").exists()


def test_seeded_setup_outputs_have_identical_content_hashes(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    _fixture(first)
    _fixture(second)
    spec = SetupSpec(seed=42)

    first_result = SetupService().run(first, first, spec)
    second_result = SetupService().run(second, second, spec)

    assert first_result.success, first_result.to_dict()
    assert second_result.success, second_result.to_dict()
    assert first_result.provenance["reproducible"] is True
    assert second_result.provenance["reproducible"] is True
    for name in first_result.candidate_files:
        assert hash_file(first / name) == hash_file(second / name)


def test_repeated_non_in_place_generation_replaces_only_owned_outputs(tmp_path: Path) -> None:
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.mkdir()
    _fixture(source)
    runner = _FakeRunner()

    first = SetupService(runner=runner).run(
        source, destination, SetupSpec(seed=42), collision_policy="replace"
    )
    first_bytes = {
        name: (destination / name).read_bytes() for name in first.candidate_files
    }
    second = SetupService(runner=runner).run(
        source, destination, SetupSpec(seed=42), collision_policy="replace"
    )

    assert first.success and second.success
    assert {
        name: (destination / name).read_bytes() for name in second.candidate_files
    } == first_bytes
    assert sorted(path.name for path in destination.iterdir()) == sorted(first_bytes)


@pytest.mark.parametrize("mode", ["worker-failure", "validation-failure"])
def test_failed_worker_or_validation_does_not_publish_or_modify_source(tmp_path: Path, mode: str) -> None:
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.mkdir()
    _fixture(source)
    before = {name: (source / name).read_bytes() for name in ("lattice.nml", "input.nml")}

    result = SetupService(runner=_FakeRunner(mode)).run(source, destination, SetupSpec())

    assert not result.success
    assert {name: (source / name).read_bytes() for name in before} == before
    assert not destination.exists()


def test_staging_copy_failure_does_not_publish_or_modify_source(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.mkdir()
    _fixture(source)
    before = {name: (source / name).read_bytes() for name in ("lattice.nml", "input.nml")}

    def fail_stage(*args, **kwargs):
        raise OSError("injected staging copy failure")

    monkeypatch.setattr("rsutil.setup_runner.stage_file", fail_stage)
    result = SetupService().run(source, destination, SetupSpec())

    assert not result.success
    assert result.diagnostics[0].code == ErrorCode.WORKSPACE_ERROR.value
    assert {name: (source / name).read_bytes() for name in before} == before
    assert not destination.exists()


def test_output_collision_is_explicit_and_non_destructive(tmp_path: Path) -> None:
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.mkdir()
    destination.mkdir()
    _fixture(source)
    existing = destination / "radii.out"
    existing.write_bytes(b"user-owned\n")
    before = existing.read_bytes()

    result = SetupService(runner=_FakeRunner()).run(source, destination, SetupSpec())

    assert not result.success
    assert result.diagnostics[0].code == ErrorCode.OUTPUT_COLLISION.value
    assert existing.read_bytes() == before
    assert not (destination / "lattice_new.nml").exists()


def test_in_place_commit_failure_restores_inputs_and_backups(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "source"
    source.mkdir()
    _fixture(source)
    before = {name: (source / name).read_bytes() for name in ("lattice.nml", "input.nml")}
    original_replace = os.replace
    calls = 0

    def fail_second_replace(source_path, target_path):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected commit failure")
        return original_replace(source_path, target_path)

    monkeypatch.setattr("rsutil.setup_service.os.replace", fail_second_replace)
    result = SetupService(runner=_FakeRunner()).run(
        source,
        source,
        SetupSpec(),
        in_place=True,
        backup=True,
        collision_policy="replace",
    )

    assert not result.success
    assert result.diagnostics[0].code == ErrorCode.COMMIT_FAILED.value
    assert {name: (source / name).read_bytes() for name in before} == before
    assert not (source / "lattice.nml.bak").exists()
    assert not (source / "input.nml.bak").exists()


def test_in_place_success_creates_backups_and_replaces_atomically(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    _fixture(source)
    before = {name: (source / name).read_bytes() for name in ("lattice.nml", "input.nml")}

    result = SetupService(runner=_FakeRunner()).run(
        source,
        source,
        SetupSpec(),
        in_place=True,
        backup=True,
        collision_policy="replace",
    )

    assert result.success, result.to_dict()
    assert (source / "lattice.nml.bak").read_bytes() == before["lattice.nml"]
    assert (source / "input.nml.bak").read_bytes() == before["input.nml"]
    assert (source / "lattice.nml").read_text(encoding="utf-8").startswith("&lattice")
    assert (source / "input.nml").read_text(encoding="utf-8").startswith("&atoms")
