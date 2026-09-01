import os
import re
import sys
import subprocess
from pathlib import Path
import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


def _package_version() -> str:
    """Read the package version without importing the package during setup."""

    version_file = Path(__file__).parent / "rsutil" / "_version.py"
    match = re.search(
        r"^__version__\s*=\s*[\"']([^\"']+)[\"']",
        version_file.read_text(),
        re.MULTILINE,
    )
    if match is None:
        raise RuntimeError(f"Could not read package version from {version_file}")
    return match.group(1)


PACKAGE_VERSION = _package_version()


RUNTIME_REQUIRES = [
    "numpy>=1.20.0",
    "pymatgen",
    "f90nml",
]

BANDS_REQUIRES = [
    "matplotlib>=3.5.0",
    "scipy>=1.7.0",
    "spglib>=2.0.0",
    "seekpath>=2.0.0",
]

TEST_REQUIRES = [
    "pytest>=7.0.0",
    "ase>=3.22",
    *BANDS_REQUIRES,
]

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        ext_fullpath = Path(self.get_ext_fullpath(ext.name)).resolve()
        extdir = ext_fullpath.parent
        extdir.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPython3_EXECUTABLE={sys.executable}",
            "-DPYTHON_BINDINGS=ON",
        ]

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


setup(
    name="rsutil",
    version=PACKAGE_VERSION,
    author="Anders Bergman",
    author_email="anders.bergman@physics.uu.se",
    description="Utilities for the RS-LMTO code including band structure and DOS analysis",
    long_description="Utilities for the RS-LMTO code including lattice processing, visualization, exchange mapping, UppASD export, band structure calculations, and density of states analysis",
    ext_modules=[CMakeExtension("rsutil.rsutil")],
    cmdclass={"build_ext": CMakeBuild},
    packages=setuptools.find_packages(exclude=["tests", "tests.*"]),
    include_package_data=True,  # ✅ Ensures package_data is included in the wheel
    zip_safe=False,
    package_data={"rsutil": ["viz/cpk_colors.json"]},
    python_requires=">=3.9",
    install_requires=RUNTIME_REQUIRES,
    extras_require={
        "ase": ["ase>=3.22"],
        "bands": BANDS_REQUIRES,
        "visualization": ["pyvista"],
        "test": TEST_REQUIRES,
        "full": [*BANDS_REQUIRES, "ase>=3.22", "pyvista"],
        "dev": [*TEST_REQUIRES, "black", "flake8"],
        "docs": ["sphinx", "myst-parser", "sphinx-rtd-theme", "sphinx-autodoc-typehints"],
    },
    entry_points={
        "console_scripts": [
            "rsutil=rsutil.cli:main"
        ],
    },
)
