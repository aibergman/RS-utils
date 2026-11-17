"""Shared helpers for legacy CLI regression tests."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from functools import lru_cache
from pathlib import Path


import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
TEST_MPLCONFIGDIR = Path(tempfile.gettempdir()) / f"rsutil-pytest-mpl-{os.getpid()}"
TEST_MPLCONFIGDIR.mkdir(exist_ok=True)


def copy_example_files(destination: Path, example: str, names: list[str]) -> None:
    """Copy only the selected files from an example into an isolated directory."""

    source_root = EXAMPLES_ROOT / example
    for name in names:
        source = source_root / name
        target = destination / Path(name).name
        if not source.is_file():
            raise FileNotFoundError(f"Missing fixture source: {source}")
        shutil.copy2(source, target)


@lru_cache(maxsize=1)
def _cli_command() -> tuple[list[str], bool]:
    """Return a usable installed CLI, or the current source checkout CLI."""

    installed = shutil.which("rsutil")
    if installed:
        probe = subprocess.run(
            [installed, "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        if probe.returncode == 0:
            return [installed], True

    return [sys.executable, "-m", "rsutil.cli"], False


def run_rsutil(destination: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    """Run the installed/current rsutil CLI from an isolated working directory."""

    command, installed = _cli_command()
    environment = os.environ.copy()
    environment["MPLCONFIGDIR"] = str(TEST_MPLCONFIGDIR)
    if not installed:
        existing_pythonpath = environment.get("PYTHONPATH")
        source_path = str(REPOSITORY_ROOT)
        environment["PYTHONPATH"] = (
            source_path
            if not existing_pythonpath
            else os.pathsep.join([source_path, existing_pythonpath])
        )

    return subprocess.run(
        [*command, *arguments],
        cwd=destination,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )


def assert_cli_success(result: subprocess.CompletedProcess[str]) -> None:
    """Provide useful diagnostics when a CLI regression subprocess fails."""

    assert result.returncode == 0, (
        f"CLI failed with exit code {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "native: needs the compiled Fortran setup kernel (skipped when it is missing)"
    )


def pytest_collection_modifyitems(config, items):
    """Skip kernel-dependent tests, with a reason, when the library is unavailable."""

    import rsutil

    if rsutil.native_library_available():
        return
    skip = pytest.mark.skip(
        reason=f"Fortran kernel unavailable ({rsutil._native_error}); "
        "install a wheel or build with CMake + gfortran"
    )
    for item in items:
        if "native" in item.keywords:
            item.add_marker(skip)
