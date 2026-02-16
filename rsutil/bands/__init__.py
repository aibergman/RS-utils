"""
Band structure and DOS analysis for RS-LMTO outputs.

This module provides comprehensive electronic structure analysis including:

**Core Functionality**
- Band structure calculations along high-symmetry paths
- Density of states (DOS) with tetrahedron method
- Orbital and site projections with gauge fixing
- DOS moments and Fermi level determination

**Advanced Analysis**
- Spin texture and magnetic moment analysis (two approaches)
- Berry phase calculations and topological properties
- Effective mass tensors
- Altermagnetic material detection

**Visualization & Utilities**
- Fat bands with spin/orbital resolution
- 3D Brillouin zone visualization
- Advanced k-path generation with seekpath integration
- Interactive k-path selection and custom path creation
- Prime convention for symmetry-related k-points

**Data Input**
- Support for multiple Hamiltonian formats (RS-LMTO, PAOFLOW, KITE)
- Automatic format detection
- Flexible lattice input

Module Organization
-------------------
The package is organized into focused submodules:

- **core**: Basic band structure calculation
- **dos**: Density of states computation
- **paths**: K-path generation (seekpath, custom)
- **projections**: Orbital and spin projections
- **analysis**: High-level analysis tools
  - spin/: Spin analysis (two complementary approaches)
  - topological/: Berry phase, Chern numbers, altermagnetic
  - magnetic/: Magnetic moments, effective mass
- **visualization**: Plotting and export utilities

All submodules are re-exported at the package level for convenience,
maintaining backward compatibility with existing code.

Examples
--------
Basic band structure calculation:

    >>> from rsutil.bands import load_lattice, load_hamiltonian, BandStructureCalculator
    >>> lattice = load_lattice('lattice.nml')
    >>> hamiltonian = load_hamiltonian('rs2txt.dat')
    >>> calc = BandStructureCalculator(lattice, hamiltonian)
    >>> kpoints, klabels, kdist = calc.setup_kpath(npoints=100)
    >>> energies = calc.calculate_bands(kpoints)
    >>> calc.plot_bands(kdist, energies, klabels, 'bands.png')

DOS calculation:

    >>> from rsutil.bands import DOSCalculator, FullBZCalculator
    >>> bz_calc = FullBZCalculator(lattice, hamiltonian)
    >>> dos_calc = DOSCalculator(bz_calc)
    >>> dos_data = dos_calc.calculate_dos(mesh_size=(8, 8, 8))
    >>> dos_calc.plot_dos(dos_data, fermi_energy=0.0)

Spin texture analysis (eigenvector-based):

    >>> from rsutil.bands.analysis.spin import SpinTextureCalculator
    >>> calc.calculate_bands(kpoints, store_eigenvectors=True)
    >>> spin = SpinTextureCalculator(calc)
    >>> rashba = spin.extract_rashba_parameters(band_idx=3)

Topological analysis:

    >>> from rsutil.bands.analysis.topological import BerryCurvatureCalculator
    >>> berry = BerryCurvatureCalculator(calc)
    >>> chern = berry.calculate_chern_number(band_range=(0, 10))

Authors & Credits
-----------------
Original implementation: rslmto_official
Advanced k-path features (seekpath, interactive selection, prime convention): Ramon Cardias
Recent enhancements & package integration: Anders Bergman
"""

# Import from readers (unchanged location)
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

# Import from core module
from .core import BandStructureCalculator

# Import from dos module
from .dos import (
    DOSCalculator,
    FullBZCalculator,
    TetrahedronDOS,
    DOSMomentsCalculator,
    analyze_dos_moments,
)

# Import from projections module
from .projections import (
    ProjectedBandCalculator,
    guess_basis_order,
    degeneracy_groups,
    gauge_fix_evecs,
    atom_orbital_indices,
    parse_atom_selection,
    get_projection_vrange,
    compute_quantity_from_sexp,
    stabilize_degenerate_projections,
)

# Import from paths module (k-path utilities)
from .paths import (
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

# Import from analysis module  
from .analysis import (
    BandAnalyzer,
    StructureAnalyzer,
    analyze_structure,
    plot_combined_analysis,
    check_files,
    run_analysis,
)

# Import from analysis.spin submodule
from .analysis.spin import (
    SpinAnalyzer,
    SpinTextureCalculator,
)

# Import from analysis.topological submodule
from .analysis.topological import (
    BerryCurvatureCalculator,
    AltermagneticAnalyzer,
)

# Import from analysis.magnetic submodule
from .analysis.magnetic import (
    MagneticMomentCalculator,
    EffectiveMassCalculator,
)

# Import from visualization module
from .visualization import (
    FatBandsPlotter,
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

