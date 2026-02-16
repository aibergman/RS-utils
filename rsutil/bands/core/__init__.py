"""
Core band structure calculation module.

This module provides the fundamental band structure calculation engine
for TB-LMTO electronic structure data.

Classes
-------
BandStructureCalculator
    Core class for calculating band energies along arbitrary k-paths.
    Supports eigenvector storage for downstream analysis (spin texture,
    effective mass, berry curvature, etc.).

Typical Usage
-------------
Basic band structure calculation:

    >>> from rsutil.bands.core import BandStructureCalculator
    >>> from rsutil.bands import load_lattice, load_hamiltonian
    >>> lattice = load_lattice('lattice.nml')
    >>> hamiltonian = load_hamiltonian('rs2txt.dat')
    >>> calc = BandStructureCalculator(lattice, hamiltonian)
    >>> kpoints, klabels, kdist = calc.setup_kpath(npoints=100)
    >>> bands = calc.calculate_bands(kpoints)
    >>> calc.plot_bands(kdist, bands, klabels, 'bands.png')

With eigenvector storage for downstream analysis:

    >>> calc.calculate_bands(kpoints, store_eigenvectors=True)
    >>> # Now eigenvectors are available for spin/curvature/mass analysis
    >>> from rsutil.bands.analysis.spin import SpinTextureCalculator
    >>> spin_calc = SpinTextureCalculator(calc)
    >>> spin_texture = spin_calc.calculate_spin_texture()
"""

from .band_structure import BandStructureCalculator
from .hamiltonian import fourier_transform_hamiltonian, HamiltonianDumper
from .plotting import plot_bands, plot_bands_collinear

__all__ = [
    'BandStructureCalculator',
    'fourier_transform_hamiltonian',
    'HamiltonianDumper',
    'plot_bands',
    'plot_bands_collinear',
]
