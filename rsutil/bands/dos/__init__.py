"""
Density of States (DOS) calculation module.

This module provides DOS calculation and analysis for TB electronic structures
using the tetrahedron integration method with optional Gaussian smearing.

Classes
-------
DOSCalculator
    Main class for computing density of states from regular k-point mesh.
    Supports tetrahedron method and Gaussian broadening.

TetrahedronDOS
    Tetrahedron integration kernel for accurate DOS near band edges.

FullBZCalculator
    Full Brillouin zone k-point mesh generation with symmetry reduction.

DOSMomentsCalculator
    Calculate moments of DOS and determine Fermi level properties.

Functions
---------
analyze_dos_moments
    Quick analysis of DOS moments (<E>, width, etc.).

Typical Usage
-------------
Basic DOS calculation:

    >>> from rsutil.bands.dos import DOSCalculator, FullBZCalculator
    >>> from rsutil.bands import load_lattice, load_hamiltonian
    >>> lattice = load_lattice('lattice.nml')
    >>> hamiltonian = load_hamiltonian('rs2txt.dat')
    >>> bz_calc = FullBZCalculator(lattice, hamiltonian)
    >>> dos_calc = DOSCalculator(bz_calc)
    >>> dos_data = dos_calc.calculate_dos(mesh_size=(8, 8, 8))
    >>> dos_calc.plot_dos(dos_data, fermi_energy=0.0, output='dos.png')

With DOS moments analysis:

    >>> from rsutil.bands.dos import DOSMomentsCalculator, analyze_dos_moments
    >>> moments_calc = DOSMomentsCalculator(dos_calc, dos_data)
    >>> moments = moments_calc.calculate_dos_moments()
    >>> analysis = analyze_dos_moments(dos_data, lattice)
"""

from .dos import DOSCalculator, FullBZCalculator, TetrahedronDOS
from .dos_moments import DOSMomentsCalculator, analyze_dos_moments

__all__ = [
    'DOSCalculator',
    'FullBZCalculator',
    'TetrahedronDOS',
    'DOSMomentsCalculator',
    'analyze_dos_moments',
]
