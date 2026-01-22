"""
rs_band - Real-space to Band structure analysis for TB-LMTO

A comprehensive Python package for analyzing TB-LMTO electronic structure data,
including band structure calculations, density of states, orbital projections,
and non-collinear spin analysis.
"""

from .band_structure_calculator import BandStructureCalculator, LatticeReader, HamiltonianReader
from .projected_bands import ProjectedBandCalculator
from .spin_analysis import SpinAnalyzer
from .full_bz_dos import FullBZCalculator, DOSCalculator, TetrahedronDOS
from .band_analysis import BandAnalyzer

__version__ = "1.0.0"
__author__ = "TB-LMTO Analysis Suite"
__email__ = "contact@example.com"

__all__ = [
    "BandStructureCalculator",
    "LatticeReader", 
    "HamiltonianReader",
    "ProjectedBandCalculator",
    "SpinAnalyzer",
    "FullBZCalculator",
    "DOSCalculator",
    "TetrahedronDOS",
    "BandAnalyzer"
]
