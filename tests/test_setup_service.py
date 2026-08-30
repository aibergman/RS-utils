"""HT-04 coverage for the path-explicit setup service."""

from __future__ import annotations

import os
from pathlib import Path

import f90nml

from rsutil.core import SetupSpec
from rsutil.setup import SetupService

from .conftest import copy_example_files


def _fixture(destination: Path) -> None:
    copy_example_files(destination, "bccFe", ["Dist/lattice.nml", "Dist/input.nml"])


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
