"""
Band structure and DOS analysis for RS-LMTO outputs.

This module provides comprehensive electronic structure analysis including:
- Band structure calculations along high-symmetry paths
- Density of states (DOS) with tetrahedron method
- Orbital and site projections
- DOS moments and Fermi level determination
- Non-collinear spin texture analysis

Examples
--------
Basic band structure calculation:

    >>> from rsutil.bands import LatticeReader, HamiltonianReader, BandStructureCalculator
    >>> lattice = LatticeReader('lattice.nml')
    >>> hamiltonian = HamiltonianReader('rs2txt.dat')
    >>> calc = BandStructureCalculator(lattice, hamiltonian)
    >>> kpoints, klabels, kdist = calc.setup_kpath(npoints=100)
    >>> bands = calc.calculate_bands(kpoints)
    >>> calc.plot_bands(kdist, bands, klabels, 'bands.png')

DOS calculation:

    >>> from rsutil.bands import DOSCalculator, FullBZCalculator
    >>> bz_calc = FullBZCalculator(lattice, hamiltonian)
    >>> dos_calc = DOSCalculator(bz_calc)
    >>> dos_data = dos_calc.calculate_dos(mesh_size=(8, 8, 8))
    >>> dos_calc.plot_dos(dos_data, fermi_energy=0.0, output_file='dos.png')
"""

from .readers import LatticeReader, HamiltonianReader
from .band_structure import BandStructureCalculator
from .dos import DOSCalculator, FullBZCalculator, TetrahedronDOS
from .dos_moments import DOSMomentsCalculator
from .projections import ProjectedBandCalculator
from .spin_analysis import SpinAnalyzer
from .band_analysis import BandAnalyzer

__version__ = "1.0.0"

__all__ = [
    "LatticeReader",
    "HamiltonianReader",
    "BandStructureCalculator",
    "DOSCalculator",
    "FullBZCalculator",
    "TetrahedronDOS",
    "DOSMomentsCalculator",
    "ProjectedBandCalculator",
    "SpinAnalyzer",
    "BandAnalyzer",
]
