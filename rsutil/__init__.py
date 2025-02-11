"""
This module provides a Python wrapper for the `process_lattice` Fortran routine
from the `rsutil` shared library.

Attributes
----------
lib_filename : Optional[str]
    The filename of the `rsutil` shared library, determined dynamically.
lib_path : str
    The full path to the `rsutil` shared library.
_rsutil : ctypes.CDLL
    The loaded `rsutil` shared library.

Functions
---------
process_lattice() -> None
    Python wrapper for the `process_lattice` Fortran routine.
"""
import os
import ctypes
from typing import Optional

# Dynamically determine the correct shared library filename
lib_filename: Optional[str] = None

for filename in os.listdir(os.path.dirname(__file__)):
    if filename.startswith("rsutil") and filename.endswith(".so"):
        lib_filename = filename
        break

if lib_filename is None:
    raise ImportError("Could not find the rsutil shared library.")

# Load the shared library
lib_path: str = os.path.join(os.path.dirname(__file__), lib_filename)
_rsutil: ctypes.CDLL = ctypes.CDLL(lib_path)

# Define process_lattice() signature
_rsutil.process_lattice.argtypes = []  # Modify if needed
_rsutil.process_lattice.restype = None  # Modify if needed


def process_lattice() -> None:
    """Python wrapper for the process_lattice() Fortran routine."""
    _rsutil.process_lattice()
