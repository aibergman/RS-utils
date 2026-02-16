"""
Magnetic property analysis module.

Calculate magnetic properties of electronic structures including
magnetic moments and effective carrier mass.

Classes
-------
MagneticMomentCalculator
    Calculate various magnetic moments and susceptibilities.

EffectiveMassCalculator
    Calculate effective masses from band curvature.

Typical Usage
-------------
Magnetic moment analysis:

    >>> from rsutil.bands.analysis.magnetic import MagneticMomentCalculator
    >>> mag = MagneticMomentCalculator(band_calc)
    >>> moments = mag.calculate_magnetic_moments()
    >>> susceptibility = mag.calculate_susceptibility()

Effective mass calculations:

    >>> from rsutil.bands.analysis.magnetic import EffectiveMassCalculator
    >>> eff_mass = EffectiveMassCalculator(band_calc)
    >>> masses = eff_mass.calculate_effective_masses(band_ids=[0, 1, 2])
    >>> eff_mass.plot_effective_mass_tensors()
"""

from .magnetic_moments import MagneticMomentCalculator
from .effective_mass import EffectiveMassCalculator

__all__ = [
    'MagneticMomentCalculator',
    'EffectiveMassCalculator',
]
