import pytest
"""Release packaging and native-artifact contract checks."""

from pathlib import Path

import rsutil


def test_package_version_is_release_candidate() -> None:
    assert rsutil.__version__ == "1.0.0"


@pytest.mark.native


def test_native_artifact_is_package_local() -> None:
    package_dir = Path(rsutil.__file__).resolve().parent

    assert rsutil._rsutil is not None
    assert rsutil.lib_path is not None
    native_path = Path(rsutil.lib_path).resolve()
    assert native_path.is_file()
    assert native_path.parent == package_dir
    assert native_path.name.startswith("rsutil")
    assert native_path.suffix in {".so", ".dylib", ".pyd", ".dll"}


def _load_setup_module():
    import importlib.util

    path = Path(__file__).resolve().parents[1] / "setup.py"
    if not path.is_file():  # tests copied away from the source tree (wheel suite)
        pytest.skip("setup.py not available outside the source tree")
    spec = importlib.util.spec_from_file_location("rsutil_build_setup", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_missing_toolchain_gives_actionable_error(monkeypatch):
    build = _load_setup_module()
    monkeypatch.delenv("FC", raising=False)
    monkeypatch.setattr(build.shutil, "which", lambda _name: None)
    with pytest.raises(RuntimeError) as excinfo:
        build._check_toolchain()
    message = str(excinfo.value)
    assert "CMake" in message and "Fortran compiler" in message
    assert "pip install rslmto-utils" in message  # points at the wheel route
    assert "apt install gfortran" in message and "brew install gcc" in message


def test_toolchain_check_passes_when_tools_exist(monkeypatch):
    build = _load_setup_module()
    monkeypatch.delenv("FC", raising=False)
    monkeypatch.setattr(build.shutil, "which", lambda name: f"/usr/bin/{name}")
    build._check_toolchain()  # must not raise
