"""
Spin analysis module.

Two complementary approaches for spin analysis:

1. **SpinMomentAnalyzer** (Hamiltonian-based):
   - Input: ProjectedBandCalculator with H(k) decomposition
   - Output: Spin moments on atoms/orbitals
   - Use case: Magnetic structure analysis, skyrmions

    (SpinAnalyzer is kept as a backward-compatible alias.)

2. **SpinTextureCalculator** (Eigenvector-based):
   - Input: BandStructureCalculator with eigenvectors
   - Output: Spin texture in k-space
   - Use case: Rashba splitting, topological properties

Classes
-------
SpinMomentAnalyzer
    Hamiltonian block-based spin analysis.
    Computes magnetic moments and spin anisotropy energy.

SpinAnalyzer
    Backward-compatible alias for SpinMomentAnalyzer.

SpinTextureCalculator
    Eigenvector-based spin texture calculation.
    Detects Rashba splitting and spin-momentum locking.

Typical Usage
-------------
Magnetic structure (using SpinMomentAnalyzer):

    >>> from rsutil.bands.analysis.spin import SpinMomentAnalyzer
    >>> from rsutil.bands.projections import ProjectedBandCalculator
    >>> proj_calc = ProjectedBandCalculator(lattice, hamiltonian)
    >>> proj_calc.calculate_projected_bands(kpoints)
    >>> spin = SpinMomentAnalyzer(proj_calc)
    >>> moments = spin.calculate_spin_moments()
    >>> spin.plot_spin_texture()

K-space spin structure (using SpinTextureCalculator):

    >>> from rsutil.bands.analysis.spin import SpinTextureCalculator
    >>> from rsutil.bands.core import BandStructureCalculator
    >>> band_calc = BandStructureCalculator(lattice, hamiltonian)
    >>> band_calc.calculate_bands(kpoints, store_eigenvectors=True)
    >>> spin_tex = SpinTextureCalculator(band_calc)
    >>> rashba = spin_tex.extract_rashba_parameters(band_idx=3)
    >>> spin_tex.plot_fermi_surface_spin_texture()
"""

from .spin_analysis import SpinAnalyzer, SpinMomentAnalyzer
from .spin_texture import SpinTextureCalculator
from .plotting import plot_spin_texture, plot_orbital_magnetism
from .texture_plotting import (
    plot_spin_texture_2d,
    plot_fermi_surface_spin_texture,
)

__all__ = [
    'SpinMomentAnalyzer',
    'SpinAnalyzer',
    'SpinTextureCalculator',
    'plot_spin_texture',
    'plot_orbital_magnetism',
    'plot_spin_texture_2d',
    'plot_fermi_surface_spin_texture',
]
