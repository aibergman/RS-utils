"""Release packaging and native-artifact contract checks."""

from pathlib import Path

import rsutil


def test_package_version_is_release_candidate() -> None:
    assert rsutil.__version__ == "0.11.0rc1"


def test_native_artifact_is_package_local() -> None:
    package_dir = Path(rsutil.__file__).resolve().parent

    assert rsutil._rsutil is not None
    assert rsutil.lib_path is not None
    native_path = Path(rsutil.lib_path).resolve()
    assert native_path.is_file()
    assert native_path.parent == package_dir
    assert native_path.name.startswith("rsutil")
    assert native_path.suffix in {".so", ".dylib", ".pyd", ".dll"}
