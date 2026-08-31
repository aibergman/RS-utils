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
import sys
import glob
import ctypes
import importlib
from typing import Optional

__version__ = "0.10.0"

from .setup_runner import FortranSetupRunner
from .setup_service import SetupService, run_setup
from .results import SCFResult, parse_run
from .campaign import (
    Campaign,
    CampaignEntry,
    CampaignGenerator,
    CampaignIndex,
    generate_campaign,
    load_campaign,
    validate_campaign,
)

__all__ = [
    'process_lattice',
    'bands',
    'FortranSetupRunner',
    'SetupService',
    'run_setup',
    'SCFResult',
    'parse_run',
    'Campaign',
    'CampaignEntry',
    'CampaignGenerator',
    'CampaignIndex',
    'generate_campaign',
    'load_campaign',
    'validate_campaign',
]


def __getattr__(name):
    """Load the optional, dependency-heavy bands package on first access."""

    if name == "bands":
        module = importlib.import_module(f"{__name__}.bands")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# Dynamically determine the correct shared library filename
lib_filename: Optional[str] = None

# First try: compiled shared object located inside package directory
pkg_dir = os.path.dirname(__file__)
for filename in os.listdir(pkg_dir):
    if filename.startswith("rsutil") and filename.endswith(".so"):
        lib_filename = filename
        break

# Fallback: if not found in package, search common build output directories
if lib_filename is None:
    # Look for build outputs (editable/develop builds) and wheelhouse
    search_patterns = [
        os.path.join("..", "build", "**", "rsutil*.so"),
        os.path.join("..", "build", "lib*", "**", "rsutil*.so"),
        os.path.join("..", "bdist.*", "**", "rsutil*.so"),
        os.path.join("..", "wheelhouse", "**", "rsutil*.so"),
    ]
    candidates = []
    for pat in search_patterns:
        candidates.extend(glob.glob(os.path.join(pkg_dir, pat), recursive=True))
    # If we found any candidate shared libs, pick the first (most likely)
    if candidates:
        # Use basename as lib_filename and adjust pkg_dir for loading below
        lib_full = os.path.abspath(candidates[0])
        lib_filename = os.path.basename(lib_full)
        # Prepend containing directory to sys.path so imports can find Python extensions
        lib_dir = os.path.dirname(lib_full)
        if lib_dir not in sys.path:
            sys.path.insert(0, lib_dir)
        # Also set pkg_dir so later join() will work if needed
        pkg_dir = lib_dir

# Load the shared library if found, but don't fail if not available
_rsutil = None
if lib_filename is not None:
    try:
        # lib_filename may come from package dir or discovered build dir; use pkg_dir
        lib_path: str = os.path.join(pkg_dir, lib_filename)
        _rsutil: ctypes.CDLL = ctypes.CDLL(lib_path)
        # Define process_lattice() signature with multiple arguments
        _rsutil.process_lattice.argtypes = [
            ctypes.c_char_p, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int
        ]
        _rsutil.process_lattice.restype = None
    except Exception as e:
        import warnings
        warnings.warn(f"Could not load rsutil shared library: {e}. Some Fortran functionality will be unavailable.")
        _rsutil = None
else:
    import warnings
    warnings.warn("Could not find the rsutil shared library. Building from source or ensure the package was compiled correctly.")
    _rsutil = None


def _process_lattice_unsafe(cluster_type: str = "periodic", empty_spheres: bool = False,
                            min_empty_radius: float = 1.2, max_empty_radius: float = 2.5,
                            lmx: int = 2, lattice_dist: float = 0.0, basis_dist: float = 0.0,
                            alat_scale: float = 1.0, seed: Optional[int] = None,
                            mode_2d: bool = False, z_iso_factor: Optional[float] = None,
                            rec_steps: int = 20, rec_method: str = "block") -> None:
    """Invoke the legacy ctypes entry point in the current process.

    This function is deliberately explicit about being unsafe: the Fortran
    call graph contains ``STOP``/``ERROR STOP`` paths that can terminate the
    Python interpreter.  The HT setup worker is the only new caller of this
    function.  Keep the ctypes argument conversion here identical to the
    historical public wrapper.
    """
    if _rsutil is None:
        raise RuntimeError("The rsutil shared library is not available. Please rebuild the package with: pip install .")

    # Convert empty_spheres boolean to integer (0 or 1)
    empty_spheres_int = ctypes.c_int(1 if empty_spheres else 0)
    min_radius_dbl = ctypes.c_double(min_empty_radius)
    max_radius_dbl = ctypes.c_double(max_empty_radius)

    lattice_dist_dbl = ctypes.c_double(lattice_dist)
    basis_dist_dbl = ctypes.c_double(basis_dist)
    alat_scale_dbl = ctypes.c_double(alat_scale)
    seed_flag = ctypes.c_int(1 if seed is not None else 0)
    seed_value = ctypes.c_int(0 if seed is None else int(seed))
    mode_2d_flag = ctypes.c_int(1 if mode_2d else 0)
    z_iso_factor_dbl = ctypes.c_double(-1.0 if z_iso_factor is None else float(z_iso_factor))
    rec_steps_int = ctypes.c_int(int(rec_steps))
    rec_method_map = {"block": 0, "chebyshev": 1, "lanczos": 2}
    rec_method_flag = ctypes.c_int(rec_method_map.get(str(rec_method).lower(), 0))

    # Keep the low-level call signature and argument order unchanged.
    _rsutil.process_lattice(
        cluster_type.encode('utf-8'),
        empty_spheres_int,
        min_radius_dbl,
        max_radius_dbl,
        ctypes.c_int(lmx),
        lattice_dist_dbl,
        basis_dist_dbl,
        alat_scale_dbl,
        seed_flag,
        seed_value,
        mode_2d_flag,
        z_iso_factor_dbl,
        rec_steps_int,
        rec_method_flag,
    )


def process_lattice(cluster_type: str = "periodic", empty_spheres: bool = False,
                   min_empty_radius: float = 1.2, max_empty_radius: float = 2.5,
                   lmx: int = 2, lattice_dist: float = 0.0, basis_dist: float = 0.0,
                   alat_scale: float = 1.0, seed: Optional[int] = None,
                   mode_2d: bool = False, z_iso_factor: Optional[float] = None,
                   rec_steps: int = 20, rec_method: str = "block") -> None:
    """Compatibility wrapper for the historical in-process setup call."""
    if cluster_type not in ["periodic", "free"]:
        raise ValueError(f"cluster_type must be 'periodic' or 'free', got '{cluster_type}'")
    _process_lattice_unsafe(
        cluster_type=cluster_type,
        empty_spheres=empty_spheres,
        min_empty_radius=min_empty_radius,
        max_empty_radius=max_empty_radius,
        lmx=lmx,
        lattice_dist=lattice_dist,
        basis_dist=basis_dist,
        alat_scale=alat_scale,
        seed=seed,
        mode_2d=mode_2d,
        z_iso_factor=z_iso_factor,
        rec_steps=rec_steps,
        rec_method=rec_method,
    )
