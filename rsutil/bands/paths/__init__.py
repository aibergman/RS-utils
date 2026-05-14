"""
K-path generation and management module.

This module provides comprehensive k-path generation tools for band structure
calculations, including automatic high-symmetry path detection via seekpath,
custom k-path creation, and coordinate system conversions.

Developed by Ramon Cardias with integration into rsutil package.

Functions
---------
seekpath_kpath
    Generate high-symmetry k-paths using seekpath library.

build_custom_kpath
    Create custom k-paths from user-specified high-symmetry points.

get_reciprocal_lattice
    Compute reciprocal lattice vectors.

cartesian_to_fractional
    Convert k-points from Cartesian to fractional coordinates.

fractional_to_cartesian
    Convert k-points from fractional to Cartesian coordinates.

apply_prime_convention
    Apply RS-LMTO prime convention for symmetry-related k-points.

normalize_label
    Normalize and format k-point labels (e.g., Γ → G).

tex_format_label
    Format k-point labels for LaTeX output.

allocate_points_per_segment
    Distribute k-points evenly across path segments.

Typical Usage
-------------
Automatic high-symmetry path:

    >>> from rsutil.bands.paths import seekpath_kpath
    >>> from rsutil.bands import load_lattice
    >>> lattice = load_lattice('lattice.nml')
    >>> kpoints, klabels, kdist = seekpath_kpath(lattice, num_points=200)

Custom k-path:

    >>> from rsutil.bands.paths import build_custom_kpath
    >>> custom_path = build_custom_kpath(lattice, ['G', 'X', 'M', 'G'])
    >>> kpts, labels, distances = custom_path.generate(npoints=150)
"""

from .coordinates import (
    get_reciprocal_lattice,
    cartesian_to_fractional,
    fractional_to_cartesian,
)
from .symmetry import apply_prime_convention
from .labels import normalize_label, tex_format_label
from .point_allocation import allocate_points_per_segment
from .seekpath_interface import seekpath_kpath
from .custom_paths import build_custom_kpath

__all__ = [
    'seekpath_kpath',
    'build_custom_kpath',
    'get_reciprocal_lattice',
    'cartesian_to_fractional',
    'fractional_to_cartesian',
    'apply_prime_convention',
    'normalize_label',
    'tex_format_label',
    'allocate_points_per_segment',
]
