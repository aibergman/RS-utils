"""
Band structure and DOS analysis for RS-LMTO outputs.

This module provides comprehensive electronic structure analysis including:
- Band structure calculations along high-symmetry paths
- Density of states (DOS) with tetrahedron method
- Orbital and site projections
- DOS moments and Fermi level determination
- Non-collinear spin texture analysis
- Support for multiple Hamiltonian formats (RS-LMTO, PAOFLOW, KITE)
- Advanced k-path generation with seekpath integration
- Interactive k-path selection and custom path creation
- Prime convention for symmetry-related k-points
- Degenerate band stabilization strategies
- 3D Brillouin zone visualization

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

Using seekpath k-path with interactive selection:

    >>> from rsutil.bands import BandStructureCalculator
    >>> calc = BandStructureCalculator(lattice, hamiltonian)
    >>> kpoints, klabels, kdist = calc.setup_kpath_seekpath(num_points=200, interactive=True)

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

Authors
-------
Original implementation: rslmto_official
Advanced k-path features (seekpath, interactive selection, prime convention): Ramon Cardias
Package integration: Anders Bergman
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
from .dos_moments import DOSMomentsCalculator, analyze_dos_moments
from .projections import ProjectedBandCalculator
from .spin_analysis import SpinAnalyzer
from .band_analysis import BandAnalyzer
from .magnetic_moments import MagneticMomentCalculator
from .effective_mass import EffectiveMassCalculator
from .altermagnetic import AltermagneticAnalyzer
from .berry_curvature import BerryCurvatureCalculator
from .structure_analysis import StructureAnalyzer, analyze_structure
from .complete_analysis import plot_combined_analysis
from .quick_analysis import check_files, run_analysis
from .spin_texture import SpinTextureCalculator
from .fat_bands import FatBandsPlotter
from .projection_utils import (
    guess_basis_order,
    degeneracy_groups,
    gauge_fix_evecs,
    atom_orbital_indices,
    parse_atom_selection,
    get_projection_vrange,
    compute_quantity_from_sexp,
    stabilize_degenerate_projections,
)
from .kpath_utils import (
    get_reciprocal_lattice,
    cartesian_to_fractional,
    fractional_to_cartesian,
    apply_prime_convention,
    normalize_label,
    tex_format_label,
    allocate_points_per_segment,
    seekpath_kpath,
    build_custom_kpath,
)
from .viz_utils import (
    plot_brillouin_zone,
    plot_bands_with_projections,
    export_kpath,
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
    "MagneticMomentCalculator",
    "EffectiveMassCalculator",
    "AltermagneticAnalyzer",
    "BerryCurvatureCalculator",
    "StructureAnalyzer",
    "SpinTextureCalculator",
    "FatBandsPlotter",
    # Analysis utilities
    "analyze_dos_moments",
    "analyze_structure",
    "plot_combined_analysis",
    "check_files",
    "run_analysis",
    # Projection utilities
    "guess_basis_order",
    "degeneracy_groups",
    "gauge_fix_evecs",
    "atom_orbital_indices",
    "parse_atom_selection",
    "get_projection_vrange",
    "compute_quantity_from_sexp",
    "stabilize_degenerate_projections",
    # K-path utilities (Ramon Cardias)
    "get_reciprocal_lattice",
    "cartesian_to_fractional",
    "fractional_to_cartesian",
    "apply_prime_convention",
    "normalize_label",
    "tex_format_label",
    "allocate_points_per_segment",
    "seekpath_kpath",
    "build_custom_kpath",
    # Visualization utilities (Ramon Cardias)
    "plot_brillouin_zone",
    "plot_bands_with_projections",
    "export_kpath",
]

