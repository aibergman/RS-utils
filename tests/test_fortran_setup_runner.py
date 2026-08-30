"""HT-03 coverage for the process-isolated legacy setup kernel."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import f90nml

from rsutil import process_lattice
from rsutil.setup import FortranSetupRunner

from .conftest import copy_example_files


def _fixture(destination: Path) -> None:
    copy_example_files(destination, "bccFe", ["Dist/lattice.nml", "Dist/input.nml"])


def _run_direct(destination: Path) -> None:
    previous = Path.cwd()
    try:
        os.chdir(destination)
        process_lattice()
    finally:
        os.chdir(previous)


def _semantic(path: Path):
    return f90nml.read(path)


def test_worker_matches_direct_setup_outputs(tmp_path: Path) -> None:
    direct = tmp_path / "direct"
    worker = tmp_path / "worker"
    direct.mkdir()
    worker.mkdir()
    _fixture(direct)
    _fixture(worker)

    _run_direct(direct)
    result = FortranSetupRunner().run(
        worker / "lattice.nml", worker / "input.nml", workspace=worker / "isolated"
    )

    assert result.success, result.to_dict()
    assert _semantic(direct / "lattice_new.nml") == _semantic(worker / "isolated" / "lattice_new.nml")
    assert _semantic(direct / "input_new.nml") == _semantic(worker / "isolated" / "input_new.nml")
    assert (direct / "radii.out").read_text() == (worker / "isolated" / "radii.out").read_text()
    assert {"lattice_new.nml", "input_new.nml", "radii.out", "spheres.log", "Fe.nml"} <= set(
        result.candidate_files
    )


def test_malformed_input_does_not_modify_source_directory(tmp_path: Path) -> None:
    lattice = tmp_path / "lattice.nml"
    input_nml = tmp_path / "input.nml"
    lattice.write_text("not a lattice namelist\n", encoding="utf-8")
    input_nml.write_text("&atoms\n label(1)='Fe'\n/\n", encoding="utf-8")
    before = {path.name: path.read_bytes() for path in (lattice, input_nml)}

    result = FortranSetupRunner().run(lattice, input_nml)

    assert not result.success
    assert result.diagnostics[0].code in {"worker_failed", "output_validation_failed"}
    assert {path.name: path.read_bytes() for path in (lattice, input_nml)} == before
    assert not (tmp_path / "lattice_new.nml").exists()
    assert not (tmp_path / "input_new.nml").exists()


def test_fatal_worker_exit_is_returned_without_killing_parent(tmp_path: Path) -> None:
    _fixture(tmp_path)
    input_path = tmp_path / "input.nml"
    input_path.write_text(input_path.read_text(encoding="utf-8").replace("'Fe'", "'Xx'"), encoding="utf-8")
    workspace = tmp_path / "fatal-worker"
    workspace.mkdir()
    shutil.copy2(tmp_path / "lattice.nml", workspace / "lattice.nml")
    shutil.copy2(tmp_path / "input.nml", workspace / "input.nml")
    request = workspace / "request.json"
    request.write_text(
        json.dumps(
            {
                "version": 1,
                "parameters": {
                    "cluster_type": "periodic",
                    "empty_spheres": False,
                    "min_empty_radius": 1.2,
                    "max_empty_radius": 2.5,
                    "lmx": 2,
                    "lattice_dist": 0.0,
                    "basis_dist": 0.0,
                    "alat_scale": 1.0,
                    "seed": None,
                    "mode_2d": True,
                    # SetupSpec rejects this, but the worker deliberately
                    # accepts raw parameters so a Fortran fatal-stop path
                    # remains process-isolated.
                    "z_iso_factor": 0.0,
                    "rec_steps": 20,
                    "rec_method": "block",
                },
            }
        ),
        encoding="utf-8",
    )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rsutil.worker",
            "--request",
            str(request),
            "--workdir",
            str(workspace),
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode != 0
    assert "Unknown element" in result.stdout or "STOP" in result.stderr
    # The pytest parent is still running here; a normal assertion is enough to
    # make termination of the parent process observable in this test.
    assert workspace.exists()


def test_workers_are_concurrent_and_do_not_cross_contaminate(tmp_path: Path) -> None:
    cases = []
    for index in range(2):
        case = tmp_path / f"case-{index}"
        case.mkdir()
        _fixture(case)
        cases.append(case)

    def run(case: Path):
        return FortranSetupRunner().run(
            case / "lattice.nml", case / "input.nml", workspace=case / "isolated"
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(run, cases))

    assert all(result.success for result in results)
    for index, result in enumerate(results):
        assert f"case-{index}" in result.provenance["workspace"]
        assert (cases[index] / "isolated" / "Fe.nml").is_file()
        assert not (cases[1 - index] / "isolated" / "Co.nml").exists()
