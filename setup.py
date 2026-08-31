import os
import re
import sys
import subprocess
import glob
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
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPython3_EXECUTABLE={sys.executable}",
            f"-DPython3_ROOT_DIR={sys.prefix}",  # ✅ Add this to explicitly set the Python root
            "-DPYTHON_BINDINGS=ON",
        ]

        build_args = ["--config", "Release"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

        # Locate the built shared library
        shared_lib = os.path.join(extdir, "rsutil.so")
        print('Attempting to copy shared library:', shared_lib)

        if not os.path.exists(shared_lib):
            raise FileNotFoundError(f"Expected shared library {shared_lib} not found!")

        # Ensure the library is copied into the rsutil package directory
        package_dir = os.path.abspath(os.path.join(self.build_lib, "rsutil"))
        os.makedirs(package_dir, exist_ok=True)
        self.move_file(shared_lib, package_dir)


setup(
    name="rsutil",
    version=PACKAGE_VERSION,
    author="Anders Bergman",
    author_email="anders.bergman@physics.uu.se",
    description="Utilities for the RS-LMTO code including band structure and DOS analysis",
    long_description="Utilities for the RS-LMTO code including lattice processing, visualization, exchange mapping, UppASD export, band structure calculations, and density of states analysis",
    ext_modules=[CMakeExtension("rsutil")],
    cmdclass={"build_ext": CMakeBuild},
    packages=setuptools.find_packages(exclude=["tests", "tests.*"]),
    include_package_data=True,  # ✅ Ensures package_data is included in the wheel
    zip_safe=False,
    package_data={"rsutil": ["*.so", "viz/cpk_colors.json"]},
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
