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
import warnings
from typing import Optional

# Dynamically determine the correct shared library filename
lib_filename: Optional[str] = None

for filename in os.listdir(os.path.dirname(__file__)):
    if filename.startswith("rsutil") and filename.endswith(".so"):
        lib_filename = filename
        break

_rsutil: Optional[ctypes.CDLL] = None

if lib_filename is None:
    warnings.warn(
        "Could not find the rsutil shared library. "
        "The process_lattice() function will not be available. "
        "Build the package with 'python setup.py build_ext --inplace' to enable it.",
        ImportWarning
    )
else:
    # Load the shared library
    lib_path: str = os.path.join(os.path.dirname(__file__), lib_filename)
    try:
        _rsutil = ctypes.CDLL(lib_path)
    except OSError:
        warnings.warn(
            f"Could not load the rsutil shared library at {lib_path}. "
            "The process_lattice() function will not be available.",
            ImportWarning
        )

if _rsutil is not None:
    # Define process_lattice() signature with multiple arguments
    _rsutil.process_lattice.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_double, ctypes.c_double]
    _rsutil.process_lattice.restype = None


def process_lattice(cluster_type: str = "periodic", empty_spheres: bool = False,
                   min_empty_radius: float = 1.2, max_empty_radius: float = 2.5) -> None:
    """
    Python wrapper for the process_lattice() Fortran routine.

    Parameters
    ----------
    cluster_type : str, optional
        Type of cluster: 'periodic' (default) or 'free'.
        If 'periodic', sets pbc to .true.
        If 'free', sets pbc to .false.
    empty_spheres : bool, optional
        Enable empty sphere calculations (default: False).
    min_empty_radius : float, optional
        Minimum empty sphere radius (default: 1.2).
    max_empty_radius : float, optional
        Maximum empty sphere radius (default: 2.5).
    """
    if _rsutil is None:
        raise RuntimeError(
            "The rsutil shared library is not available. "
            "Build the package with 'python setup.py build_ext --inplace' to enable this function."
        )

    if cluster_type not in ["periodic", "free"]:
        raise ValueError(f"cluster_type must be 'periodic' or 'free', got '{cluster_type}'")

    # Convert empty_spheres boolean to integer (0 or 1)
    empty_spheres_int = ctypes.c_int(1 if empty_spheres else 0)
    min_radius_dbl = ctypes.c_double(min_empty_radius)
    max_radius_dbl = ctypes.c_double(max_empty_radius)

    # Call the Fortran routine with separate arguments
    _rsutil.process_lattice(cluster_type.encode('utf-8'), empty_spheres_int, min_radius_dbl, max_radius_dbl)
