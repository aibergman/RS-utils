#!/usr/bin/env python3
"""
DEPRECATED: Use rsutil.bands.readers instead.

This module is kept for backward compatibility.
All functionality has been moved to rsutil/bands/readers/ subpackage.
"""

import warnings

warnings.warn(
    "Importing from rsutil.bands.readers module is deprecated. "
    "Please import from rsutil.bands.readers instead:\n"
    "  from rsutil.bands.readers import LatticeReader, HamiltonianReader\n"
    "  from rsutil.bands.readers import load_lattice, load_hamiltonian",
    DeprecationWarning,
    stacklevel=2
)

# Re-export for backward compatibility
from .readers import (
    LatticeReader,
    HamiltonianReader,
    BaseLatticeReader,
    BaseHamiltonianReader,
    load_lattice,
    load_hamiltonian,
    detect_lattice_format,
    detect_hamiltonian_format,
)

__all__ = [
    'LatticeReader',
    'HamiltonianReader',
    'BaseLatticeReader',
    'BaseHamiltonianReader',
    'load_lattice',
    'load_hamiltonian',
    'detect_lattice_format',
    'detect_hamiltonian_format',
]
