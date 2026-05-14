"""
Band structure visualization utilities.

Provides functions for plotting band structures with publication-quality styling,
supporting both non-collinear and collinear spin-polarized calculations.

Author: Ramon Cardias
"""

import numpy as np
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator


def plot_bands(
    kdistances: np.ndarray,
    band_energies: np.ndarray,
    klabels: List[str],
    output_file: Optional[str] = None,
    energy_range: Optional[Tuple[float, float]] = None,
) -> None:
    """Plot the band structure with publication-quality styling.
    
    Parameters
    ----------
    kdistances : np.ndarray
        K-point distances along path
    band_energies : np.ndarray
        Band energies (nk × nbands)
    klabels : List[str]
        High-symmetry point labels
    output_file : str, optional
        Output file path. If None, displays plot interactively.
    energy_range : tuple, optional
        (E_min, E_max) for y-axis limits
        
    Author
    ------
    Ramon Cardias
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot all bands with clean styling
    for band in range(band_energies.shape[1]):
        ax.plot(
            kdistances,
            band_energies[:, band],
            color="#2E86AB",
            linewidth=1.5,
            alpha=0.8,
        )

    # Set up k-point labels with explicit tick control
    tick_positions = []
    tick_labels = []
    for i, label in enumerate(klabels):
        if label:
            tick_positions.append(kdistances[i])
            tick_labels.append(label)

    if tick_positions:
        ax.xaxis.set_major_locator(FixedLocator(tick_positions))
        ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))
        ax.xaxis.set_minor_locator(NullLocator())

        # Add thin vertical lines at symmetry points
        for pos in tick_positions[1:-1]:
            ax.axvline(
                pos,
                color="black",
                linestyle="-",
                linewidth=0.5,
                alpha=0.5,
                zorder=0,
            )

    # Set axis limits and styling
    ax.set_xlim(kdistances[0], kdistances[-1])
    ax.margins(x=0)

    # Compute tight y-limits
    if energy_range:
        target_ylim = tuple(energy_range)
    else:
        e_min = float(np.min(band_energies))
        e_max = float(np.max(band_energies))
        e_range = e_max - e_min
        e_pad = 0.02 * e_range
        target_ylim = (e_min - e_pad, e_max + e_pad)

    ax.set_ylim(*target_ylim)
    ax.autoscale(False)

    ax.set_xlabel("k-path", fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_title("TB-LMTO Band Structure", fontsize=14)

    # Clean appearance
    ax.grid(False)
    ax.minorticks_off()
    ax.set_axisbelow(True)

    if output_file:
        fig.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Band structure saved to {output_file}")
        plt.close()
    else:
        plt.show()


def plot_bands_collinear(
    kdistances: np.ndarray,
    band_energies_up: np.ndarray,
    band_energies_down: np.ndarray,
    klabels: List[str],
    output_file: Optional[str] = None,
    energy_range: Optional[Tuple[float, float]] = None,
) -> None:
    """Plot spin-polarized band structure with publication-quality styling.
    
    Spin-up bands shown in red, spin-down in blue. Exchange splitting
    is visible as energy difference between channels.
    
    Parameters
    ----------
    kdistances : np.ndarray
        Distances along k-path for x-axis
    band_energies_up : np.ndarray
        Spin-up band energies (nk × nbands)
    band_energies_down : np.ndarray
        Spin-down band energies (nk × nbands)
    klabels : List[str]
        Labels for high-symmetry k-points
    output_file : str, optional
        Output filename for saving the plot. If None, displays interactively.
    energy_range : tuple, optional
        (E_min, E_max) for y-axis limits
        
    Author
    ------
    Ramon Cardias
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot spin-up bands in red
    for band in range(band_energies_up.shape[1]):
        ax.plot(
            kdistances,
            band_energies_up[:, band],
            color="#E63946",
            linewidth=1.5,
            alpha=0.75,
            label="Spin ↑" if band == 0 else "",
        )

    # Plot spin-down bands in blue
    for band in range(band_energies_down.shape[1]):
        ax.plot(
            kdistances,
            band_energies_down[:, band],
            color="#457B9D",
            linewidth=1.5,
            alpha=0.75,
            label="Spin ↓" if band == 0 else "",
        )

    # Set up k-point labels with explicit tick control
    tick_positions = []
    tick_labels = []
    for i, label in enumerate(klabels):
        if label:
            tick_positions.append(kdistances[i])
            tick_labels.append(label)

    if tick_positions:
        ax.xaxis.set_major_locator(FixedLocator(tick_positions))
        ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))
        ax.xaxis.set_minor_locator(NullLocator())

        # Add thin vertical lines at symmetry points
        for pos in tick_positions[1:-1]:
            ax.axvline(
                pos,
                color="black",
                linestyle="-",
                linewidth=0.5,
                alpha=0.5,
                zorder=0,
            )

    # Set axis limits and styling
    ax.set_xlim(kdistances[0], kdistances[-1])
    ax.margins(x=0)

    # Compute tight y-limits from both spin channels
    if energy_range:
        target_ylim = tuple(energy_range)
    else:
        e_min = float(min(np.min(band_energies_up), np.min(band_energies_down)))
        e_max = float(max(np.max(band_energies_up), np.max(band_energies_down)))
        e_range = e_max - e_min
        e_pad = 0.02 * e_range
        target_ylim = (e_min - e_pad, e_max + e_pad)

    ax.set_ylim(*target_ylim)
    ax.autoscale(False)

    ax.set_xlabel("k-path", fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_title("TB-LMTO Spin-Polarized Band Structure", fontsize=14)
    ax.legend(loc="best", frameon=False)

    # Clean appearance
    ax.grid(False)
    ax.minorticks_off()
    ax.set_axisbelow(True)

    if output_file:
        fig.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Spin-polarized band structure saved to {output_file}")
        plt.close()
    else:
        plt.show()
