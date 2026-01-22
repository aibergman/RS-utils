"""
Band structure and DOS analysis for RS-LMTO outputs.

This module provides comprehensive electronic structure analysis including:
- Band structure calculations along high-symmetry paths
- Density of states (DOS) with tetrahedron method
- Orbital and site projections
- DOS moments and Fermi level determination
- Non-collinear spin texture analysis
- Support for multiple Hamiltonian formats (RS-LMTO, PAOFLOW, KITE)

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

Using format detection:

    >>> from rsutil.bands import load_lattice, load_hamiltonian, BandStructureCalculator
    >>> lattice = load_lattice('lattice.nml')  # Auto-detects RS-LMTO
    >>> hamiltonian = load_hamiltonian('rs2txt.dat')
    >>> calc = BandStructureCalculator(lattice, hamiltonian)

DOS calculation:

    >>> from rsutil.bands import DOSCalculator, FullBZCalculator
    >>> bz_calc = FullBZCalculator(lattice, hamiltonian)
    >>> dos_calc = DOSCalculator(bz_calc)
    >>> dos_data = dos_calc.calculate_dos(mesh_size=(8, 8, 8))
    >>> dos_calc.plot_dos(dos_data, fermi_energy=0.0, output_file='dos.png')
"""

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
from .band_structure import BandStructureCalculator
from .dos import DOSCalculator, FullBZCalculator, TetrahedronDOS
from .dos_moments import DOSMomentsCalculator
from .projections import ProjectedBandCalculator
from .spin_analysis import SpinAnalyzer
from .band_analysis import BandAnalyzer
from .projection_utils import (
    guess_basis_order,
    degeneracy_groups,
    gauge_fix_evecs,
    atom_orbital_indices,
    parse_atom_selection,
    get_projection_vrange,
    compute_quantity_from_sexp,
)

__version__ = "1.0.0"

__all__ = [
    # Readers and factories
    "LatticeReader",
    "HamiltonianReader",
    "BaseLatticeReader",
    "BaseHamiltonianReader",
    "load_lattice",
    "load_hamiltonian",
    "detect_lattice_format",
    "detect_hamiltonian_format",
    # Calculators
    "BandStructureCalculator",
    "DOSCalculator",
    "FullBZCalculator",
    "TetrahedronDOS",
    "DOSMomentsCalculator",
    "ProjectedBandCalculator",
    "SpinAnalyzer",
    "BandAnalyzer",
    # Projection utilities
    "guess_basis_order",
    "degeneracy_groups",
    "gauge_fix_evecs",
    "atom_orbital_indices",
    "parse_atom_selection",
    "get_projection_vrange",
    "compute_quantity_from_sexp",
]

