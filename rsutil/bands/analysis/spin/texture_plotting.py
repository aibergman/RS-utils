"""
Spin texture plotting utilities.

Provides 2D and 3D visualization helpers for spin texture analysis.

Author: Ramon Cardias
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt


def plot_spin_texture_2d(
    kpoints: np.ndarray,
    spin_texture: np.ndarray,
    spin_polarization: np.ndarray,
    band_idx: int,
    component: str = "z",
    plane: str = "xy",
    plane_position: float = 0.0,
    output_file: Optional[str] = None,
) -> None:
    """Plot 2D spin texture with arrows showing spin direction."""
    # Select k-points near the specified plane
    plane_map = {"xy": 2, "xz": 1, "yz": 0}
    perp_idx = plane_map[plane]

    plane_mask = np.abs(kpoints[:, perp_idx] - plane_position) < 0.05

    if not np.any(plane_mask):
        raise ValueError(f"No k-points found near {plane} plane at {plane_position}")

    # Get coordinates for plotting
    axes_map = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}
    ax1, ax2 = axes_map[plane]

    kx = kpoints[plane_mask, ax1]
    ky = kpoints[plane_mask, ax2]

    spin_vecs = spin_texture[plane_mask, band_idx, :]

    # Determine color
    comp_map = {"x": 0, "y": 1, "z": 2}
    if component == "magnitude":
        colors = spin_polarization[plane_mask, band_idx]
    else:
        colors = spin_vecs[:, comp_map[component]]

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot arrows
    scale = 0.02
    quiver = ax.quiver(
        kx,
        ky,
        spin_vecs[:, ax1],
        spin_vecs[:, ax2],
        colors,
        cmap="RdBu_r",
        scale=scale,
        scale_units="xy",
        width=0.005,
        headwidth=3,
        headlength=4,
    )

    plt.colorbar(quiver, ax=ax, label=f"<sigma_{component}>")

    ax.set_xlabel(f"k{['x', 'y', 'z'][ax1]} (crystal coords)")
    ax.set_ylabel(f"k{['x', 'y', 'z'][ax2]} (crystal coords)")
    ax.set_title(f"Spin Texture - Band {band_idx} - {plane} plane")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved spin texture plot to {output_file}")
    else:
        plt.show()

    plt.close()


def plot_fermi_surface_spin_texture(
    kpoints: np.ndarray,
    energies: np.ndarray,
    spin_texture: np.ndarray,
    fermi_energy: float,
    band_idx: int,
    energy_window: float = 0.05,
    output_file: Optional[str] = None,
) -> None:
    """Plot 3D Fermi surface colored by spin polarization."""
    # Find k-points on Fermi surface
    band_energies = energies[:, band_idx]
    on_fermi = np.abs(band_energies - fermi_energy) < energy_window

    if not np.any(on_fermi):
        print(f"Warning: No states found within {energy_window} eV of Fermi level")
        return

    # Get Fermi surface k-points and spins
    k_fermi = kpoints[on_fermi]
    spin_fermi = spin_texture[on_fermi, band_idx, :]

    # Calculate spin polarization for coloring
    spin_mag = np.linalg.norm(spin_fermi, axis=1)

    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Plot Fermi surface points
    scatter = ax.scatter(
        k_fermi[:, 0],
        k_fermi[:, 1],
        k_fermi[:, 2],
        c=spin_mag,
        cmap="viridis",
        s=50,
        alpha=0.7,
    )

    # Add spin direction arrows (subsample for clarity)
    step = max(1, len(k_fermi) // 50)
    ax.quiver(
        k_fermi[::step, 0],
        k_fermi[::step, 1],
        k_fermi[::step, 2],
        spin_fermi[::step, 0],
        spin_fermi[::step, 1],
        spin_fermi[::step, 2],
        length=0.1,
        normalize=True,
        color="red",
        alpha=0.6,
        arrow_length_ratio=0.3,
    )

    ax.set_xlabel("kx (crystal coords)")
    ax.set_ylabel("ky (crystal coords)")
    ax.set_zlabel("kz (crystal coords)")
    ax.set_title(f"Fermi Surface Spin Texture - Band {band_idx}")

    plt.colorbar(scatter, ax=ax, label="|<sigma>|", shrink=0.6)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved 3D Fermi surface plot to {output_file}")
    else:
        plt.show()

    plt.close()
