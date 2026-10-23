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
from rsutil.core import SetupSpec
from rsutil.setup import FortranSetupRunner

from .conftest import copy_example_files


def _fixture(destination: Path) -> None:
    copy_example_files(destination, "bccFe", ["Dist/lattice.nml", "Dist/input.nml"])


def _fcc_nio_fixture(destination: Path) -> None:
    (destination / "lattice.nml").write_text(
        """&lattice
nbulk_bulk = 4
ntot = 4
nbas = 4
nrec = 4
a(:, 1) = 1.0, 0.5, 0.5
a(:, 2) = 0.5, 1.0, 0.5
a(:, 3) = 0.5, 0.5, 1.0
crd(:, 1) = 0.0, 0.0, 0.0
crd(:, 2) = 1.0, 1.0, 1.0
crd(:, 3) = 0.5, 0.5, 0.5
crd(:, 4) = 1.5, 1.5, 1.5
izp = 1, 2, 3, 4
no = 1, 2, 3, 4
iu = 1, 2, 3, 4
ib = 1, 2, 3, 4
irec = 1, 2, 3, 4
/
""",
        encoding="utf-8",
    )
    (destination / "input.nml").write_text(
        """&calculation
pre_processing = 'bravais'
verbose = .true.
/
&lattice
rc = 80
alat = 4.1684002458
crystal_sym = 'file'
wav = 1.2929345694
ntype = 4
ct = 3.3, 3.3, 3.3, 3.3
r2 = 10.9
pbc = .true.
b1 = .true.
b2 = .true.
b3 = .true.
n1 = 10
n2 = 10
n3 = 10
/
&atoms
database = './'
label = 'Niu', 'Nid', 'O1', 'O2'
/
&self
nstep = 35
use_kspace = .true.
cold = .true.
/
&energy
fermi = -0.187938
energy_min = -1.0
energy_max = 1.5
channels_ldos = 2000
/
&control
calctype = 'B'
nsp = 1
lld = 200
recur = 'chebyshev'
/
&mix
beta = 0.1
mixtype = 'broyden'
/
""",
        encoding="utf-8",
    )


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


def test_empty_spheres_found_on_first_pass_are_kept(tmp_path: Path) -> None:
    _fcc_nio_fixture(tmp_path)
    result = FortranSetupRunner().run(
        tmp_path / "lattice.nml",
        tmp_path / "input.nml",
        SetupSpec(empty_spheres=True),
        workspace=tmp_path / "worker",
    )

    assert result.success, result.to_dict()
    output = tmp_path / "worker"
    radii = (output / "radii.out").read_text(encoding="utf-8").splitlines()
    generated_lattice = f90nml.read(output / "lattice_new.nml")["lattice"]
    generated_input = f90nml.read(output / "input_new.nml")
    assert int(radii[0]) > 4
    assert any("/ E/" in line for line in radii[1:])
    assert generated_lattice["nrec"] == int(radii[0])
    assert generated_input["atoms"]["label"][4:] == [
        f"E{index}" for index in range(1, int(radii[0]) - 3)
    ]
    assert f90nml.read(output / "E1.nml")["element"]["atomic_number"] == 0


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


def test_worker_failure_keeps_diagnostics_before_owned_workspace_cleanup(
    tmp_path: Path, monkeypatch
) -> None:
    lattice = tmp_path / "lattice.nml"
    input_nml = tmp_path / "input.nml"
    lattice.write_text("lattice\n", encoding="utf-8")
    input_nml.write_text("input\n", encoding="utf-8")

    def fake_run(command, *, cwd, **kwargs):
        workspace = Path(cwd)
        (workspace / "spheres.log").write_text(
            "sphere line 1\nsphere line 2\n", encoding="utf-8"
        )
        (workspace / "radii.out").write_text(
            "radii line 1\nradii line 2\n", encoding="utf-8"
        )
        return subprocess.CompletedProcess(
            command,
            17,
            stdout="worker stdout\n",
            stderr="worker stderr\n",
        )

    monkeypatch.setattr("rsutil.setup_runner.subprocess.run", fake_run)
    result = FortranSetupRunner().run(lattice, input_nml)

    assert not result.success
    diagnostic = result.diagnostics[0]
    assert diagnostic.details["returncode"] == 17
    assert diagnostic.details["spheres_log_tail"] == "sphere line 1\nsphere line 2\n"
    assert diagnostic.details["radii_out_tail"] == "radii line 1\nradii line 2\n"
    assert result.stdout == "worker stdout\n"
    assert result.stderr == "worker stderr\n"
    assert not Path(result.provenance["workspace"]).exists()


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
