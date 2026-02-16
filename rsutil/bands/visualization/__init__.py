"""
Visualization and plotting module.

Provides comprehensive visualization tools for band structures, DOS,
projections, and other electronic structure properties.

Functions
---------
plot_brillouin_zone
    Visualize the first Brillouin zone with high-symmetry points.

plot_bands_with_projections
    Plot band structure with color-coded projections.

export_kpath
    Export k-path data in various formats.

Classes
-------
FatBandsPlotter
    Create fat bands plots with orbital/spin-resolved line widths.
    - Spin-resolved bands (spin up/down coloring)
    - Orbital character (line width or color intensity)
    - Fermi surface spin texture

Typical Usage
-------------
Basic Brillouin zone visualization:

    >>> from rsutil.bands.visualization import plot_brillouin_zone
    >>> plot_brillouin_zone(lattice, output='bz.png')

Fat bands (orbital character):

    >>> from rsutil.bands.visualization import FatBandsPlotter
    >>> fat = FatBandsPlotter(band_calc, proj_calc)
    >>> fat.plot_fat_bands_color(kpt_labels, output='fat_bands.png')

Spin-resolved fat bands:

    >>> fat.plot_spin_resolved_fat_bands(kpt_labels, energies, projections)
"""

from .fat_bands import FatBandsPlotter
from .viz_utils import (
    plot_brillouin_zone,
    plot_bands_with_projections,
    export_kpath,
)

__all__ = [
    'FatBandsPlotter',
    'plot_brillouin_zone',
    'plot_bands_with_projections',
    'export_kpath',
]
