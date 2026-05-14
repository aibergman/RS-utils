"""
Spin analysis plotting utilities.

Provides plotting helpers for spin textures and orbital-resolved magnetism.

Author: Ramon Cardias
"""

from typing import Dict, List, Optional

import numpy as np
import matplotlib.pyplot as plt


def plot_spin_texture(
    kdistances: np.ndarray,
    band_energies: np.ndarray,
    spin_moments: Dict,
    klabels: List[str],
    band_indices: Optional[List[int]] = None,
    output_file: Optional[str] = None,
) -> None:
    """Plot spin texture for selected bands."""
    if band_indices is None:
        # Find bands near Fermi level (assuming E_F = 0)
        fermi_bands = []
        for band in range(band_energies.shape[1]):
            if np.any(np.abs(band_energies[:, band]) < 0.5):
                fermi_bands.append(band)
        band_indices = fermi_bands[:5]

    fig, axes = plt.subplots(len(band_indices), 1, figsize=(12, 3 * len(band_indices)))
    if len(band_indices) == 1:
        axes = [axes]

    for idx, band in enumerate(band_indices):
        ax = axes[idx]

        # Plot band structure
        ax.plot(kdistances, band_energies[:, band], "k-", linewidth=2, label="Band")

        # Color by spin polarization
        spin_z = spin_moments["z"][:, band]
        scatter = ax.scatter(
            kdistances, band_energies[:, band], c=spin_z, cmap="RdBu", s=50, alpha=0.8
        )

        # Set up k-point labels
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:
                tick_positions.append(kdistances[i])
                tick_labels.append(label)

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)

        for pos in tick_positions:
            ax.axvline(x=pos, color="gray", linestyle="--", alpha=0.3)

        ax.set_ylabel("Energy (eV)")
        ax.set_title(f"Band {band + 1} - Spin Texture (z-component)")
        ax.grid(True, alpha=0.3)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Spin z-component")

    axes[-1].set_xlabel("k-path")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Spin texture plot saved to {output_file}")

    plt.show()


def plot_orbital_magnetism(
    kdistances: np.ndarray,
    band_energies: np.ndarray,
    projections: Dict,
    klabels: List[str],
    output_file: Optional[str] = None,
) -> None:
    """Plot orbital-resolved magnetic moments."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    orbital_groups = ["s", "p", "d"]
    colors = ["red", "blue", "green"]

    for orbital_idx, (orbital, color) in enumerate(zip(orbital_groups, colors)):
        ax = axes[orbital_idx]

        # Calculate orbital-resolved spin polarization
        spin_up_proj = projections["orbital"][orbital]
        # Simplified: use spin-up projections as orbital magnetism proxy
        spin_polarization = spin_up_proj

        # Plot bands colored by orbital magnetism
        for band in range(band_energies.shape[1]):
            weights = spin_polarization[:, band]
            if np.max(weights) > 0.05:
                sizes = weights * 50
                ax.scatter(
                    kdistances,
                    band_energies[:, band],
                    s=sizes,
                    c=color,
                    alpha=0.6,
                    label=f"{orbital}" if band == 0 else "",
                )

        # Set up k-point labels
        tick_positions = []
        tick_labels = []
        for i, label in enumerate(klabels):
            if label:
                tick_positions.append(kdistances[i])
                tick_labels.append(label)

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)

        for pos in tick_positions:
            ax.axvline(x=pos, color="gray", linestyle="--", alpha=0.3)

        ax.set_ylabel("Energy (eV)")
        ax.set_title(f"{orbital.upper()}-orbital Magnetic Contributions")
        ax.grid(True, alpha=0.3)
        ax.legend()

    axes[-1].set_xlabel("k-path")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Orbital magnetism plot saved to {output_file}")

    plt.show()
