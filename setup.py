import os
import sys
import subprocess
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

setup(
    name="rsutil",
    version="0.1.0",
    author="Your Name",
    author_email="your@email.com",
    description="Fortran-based utilities for lattice processing",
    long_description="Python bindings for rsutil, a Fortran-based numerical library.",
    ext_modules=[CMakeExtension("rsutil")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)

