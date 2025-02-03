import os
import sys
import subprocess
import glob
import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

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
            "-DPYTHON_BINDINGS=ON",
        ]
        build_args = ["--config", "Release"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

        # Locate the built shared library
        shared_lib = os.path.join(extdir, "rsutil.so")

        if not os.path.exists(shared_lib):
            raise FileNotFoundError(f"Expected shared library {shared_lib} not found!")

        # Ensure the library is copied into the rsutil package directory
        package_dir = os.path.abspath(os.path.join(self.build_lib, "rsutil"))
        os.makedirs(package_dir, exist_ok=True)
        self.copy_file(shared_lib, package_dir)

        # Explicitly add rsutil.so to package_data to ensure it's inside the wheel
        self.distribution.package_data = {
            "rsutil": ["rsutil.so"]
        }

setup(
    name="rsutil",
    version="0.9.0",
    author="Anders Bergman",
    author_email="anders.bergman@physics.uu.se",
    description="Utilities for the RS-LMTO code",
    long_description="Utilities for the RS-LMTO code",
    ext_modules=[CMakeExtension("rsutil")],
    cmdclass={"build_ext": CMakeBuild},
    packages=setuptools.find_packages(),
    include_package_data=True,  # ✅ Ensures package_data is included in the wheel
    zip_safe=False,
    package_data={"rsutil": ["rsutil.so"]},  # ✅ Explicitly include rsutil.so
    entry_points={
        "console_scripts": [
            "rsutil=rsutil.cli:main"
        ],
    },
)

