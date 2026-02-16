"""
Electronic structure analysis tools.

This subpackage provides high-level analysis of band structures and
electronic properties, organized by specialty:

- spin/: Spin texture and magnetism analysis
- topological/: Berry phase and topological properties
- magnetic/: Magnetic moments and effective mass

Additionally provides generic analysis tools and quick-access functions.

Classes
-------
BandAnalyzer
    Generic band structure analysis and property extraction.

StructureAnalyzer
    Crystal structure and symmetry-based band analysis.

Functions
---------
analyze_structure
    Quick structural analysis of band data.

plot_combined_analysis
    Create multi-panel combined analysis plots.

check_files
    Verify input file availability and format.

run_analysis
    Quick-start analysis workflow.

Subpackages
-----------
spin
    Spin texture analysis and magnetic moment calculations.
    - SpinMomentAnalyzer: Hamiltonian-based spin analysis
    - SpinAnalyzer: Backward-compatible alias
    - SpinTextureCalculator: Eigenvector-based spin analysis

topological
    Topological property calculations.
    - BerryCurvatureCalculator: Berry curvature and Chern numbers
    - AltermagneticAnalyzer: Altermagnetic material detection

magnetic
    Magnetic property calculations.
    - MagneticMomentCalculator: Magnetic moment analysis
    - EffectiveMassCalculator: Effective mass calculations

Typical Usage
-------------
Generic band analysis:

    >>> from rsutil.bands.analysis import BandAnalyzer
    >>> analyzer = BandAnalyzer(band_calc)
    >>> properties = analyzer.analyze()

Spin analysis (two complementary approaches):

    >>> from rsutil.bands.analysis.spin import SpinMomentAnalyzer, SpinTextureCalculator
    >>> # Hamiltonian-based (site-resolved moments)
    >>> spin_mom = SpinMomentAnalyzer(proj_calc)
    >>> moments = spin_mom.calculate_spin_moments()
    >>> # Eigenvector-based (k-space spin distribution)
    >>> spin_tex = SpinTextureCalculator(band_calc)
    >>> texture = spin_tex.calculate_spin_texture()

Topological properties:

    >>> from rsutil.bands.analysis.topological import BerryCurvatureCalculator
    >>> berry = BerryCurvatureCalculator(band_calc)
    >>> chern = berry.calculate_chern_number()
"""

from .band_analysis import BandAnalyzer
from .structure_analysis import StructureAnalyzer, analyze_structure
from .complete_analysis import plot_combined_analysis
from .quick_analysis import check_files, run_analysis

__all__ = [
    'BandAnalyzer',
    'StructureAnalyzer',
    'analyze_structure',
    'plot_combined_analysis',
    'check_files',
    'run_analysis',
]
