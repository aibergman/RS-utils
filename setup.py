"""Build script for the native Fortran setup kernel.

Package metadata (name, version, dependencies, ...) lives in ``pyproject.toml``.
This file only builds the CMake/Fortran extension that is loaded through ctypes.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

_FORTRAN_COMPILERS = ("gfortran", "ifx", "ifort", "flang-new", "flang")


def _check_toolchain() -> None:
    """Fail early, with install hints, if CMake or a Fortran compiler is missing."""

    missing = []
    if shutil.which("cmake") is None:
        missing.append("CMake")
    compiler = os.environ.get("FC") or next(
        (name for name in _FORTRAN_COMPILERS if shutil.which(name)), None
    )
    if compiler is None or shutil.which(compiler) is None:
        missing.append("a Fortran compiler (for example gfortran)")
    if not missing:
        return
    raise RuntimeError(
        "Building the RS-utils Fortran setup kernel needs "
        + " and ".join(missing)
        + ", which could not be found.\n"
        "  - Easiest: install a pre-built wheel (Linux and macOS wheels are published "
        "on PyPI): pip install rslmto-utils\n"
        "  - Debian/Ubuntu: sudo apt install gfortran cmake\n"
        "  - macOS (Homebrew): brew install gcc cmake\n"
        "  - conda: conda install -c conda-forge gfortran cmake\n"
        "Set FC to use a different Fortran compiler."
    )


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        _check_toolchain()
        ext_fullpath = Path(self.get_ext_fullpath(ext.name)).resolve()
        extdir = ext_fullpath.parent
        extdir.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPython3_EXECUTABLE={sys.executable}",
            "-DPYTHON_BINDINGS=ON",
        ]

        debug_flags = os.environ.get("RSUTIL_DEBUG_FORTRAN_FLAGS")
        if debug_flags:
            cmake_args.append(f"-DCMAKE_Fortran_FLAGS={debug_flags}")

        if os.environ.get("RSUTIL_DEBUG_CMAKE_VERBOSE"):
            cmake_args.append("-DCMAKE_VERBOSE_MAKEFILE=ON")

        build_args = ["--config", "Release"]

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)

        # CMake emits the unqualified target name, while setuptools owns the
        # ABI-qualified destination (and its editable-install location).
        target_name = ext.name.rsplit(".", 1)[-1]
        suffixes = (".so", ".dylib", ".pyd", ".dll")
        output_dirs = (extdir, extdir / "Release", extdir / "Debug")
        candidates = [
            output_dir / f"{target_name}{suffix}"
            for output_dir in output_dirs
            for suffix in suffixes
            if (output_dir / f"{target_name}{suffix}").is_file()
        ]

        if len(candidates) != 1:
            found = ", ".join(str(path) for path in candidates) or "none"
            raise FileNotFoundError(
                f"CMake did not produce exactly one {target_name} shared library "
                f"in {extdir}; found: {found}"
            )

        shared_lib = candidates[0]
        if shared_lib != ext_fullpath:
            # A previous build may have left the ABI-qualified destination in
            # place. It is generated output, so replace it with this build's
            # artifact before allowing setuptools to expose it.
            if ext_fullpath.exists():
                ext_fullpath.unlink()
            self.move_file(shared_lib, ext_fullpath)
        if not ext_fullpath.is_file():
            raise FileNotFoundError(
                f"Expected setuptools extension path {ext_fullpath} was not created"
            )


if __name__ == "__main__":  # setuptools runs this file as __main__; tests import it
    setup(
        ext_modules=[CMakeExtension("rsutil.rsutil")],
        cmdclass={"build_ext": CMakeBuild},
    )
