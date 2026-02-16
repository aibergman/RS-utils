"""
Topological property analysis module.

Calculate and analyze topological properties of band structures including
Berry phase effects, Chern numbers, and topological material classification.

Classes
-------
BerryCurvatureCalculator
    Calculate Berry curvature and related topological properties.
    - Chern numbers
    - Anomalous Hall conductivity
    - Orbital magnetization
    - Topological band identification

AltermagneticAnalyzer
    Detect and analyze altermagnetic materials and d-wave magnetism.

Typical Usage
-------------
Berry curvature analysis:

    >>> from rsutil.bands.analysis.topological import BerryCurvatureCalculator
    >>> berry = BerryCurvatureCalculator(band_calc)
    >>> chern = berry.calculate_chern_number(band_range=(0, 10))
    >>> berry.plot_berry_curvature(band_idx=5)

Altermagnetic analysis:

    >>> from rsutil.bands.analysis.topological import AltermagneticAnalyzer
    >>> alt = AltermagneticAnalyzer(band_calc)
    >>> is_altermagnetic = alt.detect_altermagnetic()
    >>> alt.plot_altermagnetic_signature()
"""

from .berry_curvature import BerryCurvatureCalculator
from .altermagnetic import AltermagneticAnalyzer

__all__ = [
    'BerryCurvatureCalculator',
    'AltermagneticAnalyzer',
]
