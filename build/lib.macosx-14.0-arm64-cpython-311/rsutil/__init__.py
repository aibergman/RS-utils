"""
This module provides a Python wrapper for the `process_lattice` Fortran routine
from the `rsutil` shared library, along with band structure and DOS analysis tools.

Modules
-------
bands : Band structure and DOS analysis
    Provides comprehensive electronic structure analysis including band structure
    calculations, density of states, orbital projections, and spin analysis.
exchange : Exchange interaction mapping
    Tools for creating exchange interaction input files.
export : Export functionality
    Export RS-LMTO data to other formats (e.g., UppASD).
viz : Visualization tools
    3D visualization of crystal structures and data.

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

Examples
--------
Process a lattice structure:

    >>> from rsutil import process_lattice
    >>> process_lattice(cluster_type='periodic', empty_spheres=True)

Calculate band structure:

    >>> from rsutil.bands import LatticeReader, HamiltonianReader, BandStructureCalculator
    >>> lattice = LatticeReader('lattice.nml')
    >>> hamiltonian = HamiltonianReader('rs2txt.dat')
    >>> calc = BandStructureCalculator(lattice, hamiltonian)
    >>> kpoints, klabels, kdist = calc.setup_kpath(100)
    >>> bands = calc.calculate_bands(kpoints)
    >>> calc.plot_bands(kdist, bands, klabels, 'bands.png')
"""
import os
import ctypes
from typing import Optional

# Import bands module for convenient access
from . import bands

__version__ = "0.10.0"

__all__ = ['process_lattice', 'bands']

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
    if cluster_type not in ["periodic", "free"]:
        raise ValueError(f"cluster_type must be 'periodic' or 'free', got '{cluster_type}'")
    
    # Convert empty_spheres boolean to integer (0 or 1)
    empty_spheres_int = ctypes.c_int(1 if empty_spheres else 0)
    min_radius_dbl = ctypes.c_double(min_empty_radius)
    max_radius_dbl = ctypes.c_double(max_empty_radius)
    
    # Call the Fortran routine with separate arguments
    _rsutil.process_lattice(cluster_type.encode('utf-8'), empty_spheres_int, min_radius_dbl, max_radius_dbl)
