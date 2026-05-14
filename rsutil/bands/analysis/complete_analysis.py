#!/usr/bin/env python3
"""
Complete TB-LMTO Analysis Suite

This script provides a unified interface for:
1. Line band structure calculations
2. Full BZ DOS calculations
3. Projected band analysis
4. Spin structure analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

from ..readers import LatticeReader, HamiltonianReader
from ..core import BandStructureCalculator
from ..projections import ProjectedBandCalculator
from ..dos import FullBZCalculator, DOSCalculator
from .structure_analysis import StructureAnalyzer


def plot_combined_analysis(line_data: dict, dos_data: dict, output_file: str = None):
    """Create a comprehensive 4-panel plot with all analyses"""

    fig = plt.figure(figsize=(16, 12))

    # Create custom grid layout
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Panel 1: Band structure (large, left)
    ax1 = fig.add_subplot(gs[:2, :2])

    # Plot bands
    kdistances = line_data["kdistances"]
    band_energies = line_data["band_energies"]
    klabels = line_data["klabels"]

    for band in range(band_energies.shape[1]):
        ax1.plot(kdistances, band_energies[:, band], "b-", linewidth=1)

    # Set up k-point labels
    tick_positions = []
    tick_labels = []
    for i, label in enumerate(klabels):
        if label:
            tick_positions.append(kdistances[i])
            tick_labels.append(label)

    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels)

    for pos in tick_positions:
        ax1.axvline(x=pos, color="k", linestyle="--", alpha=0.3)

    ax1.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax1.set_xlabel("k-path")
    ax1.set_ylabel("Energy (eV)")
    ax1.set_title("Band Structure")
    ax1.grid(True, alpha=0.3)

    # Panel 2: Total DOS (top right)
    ax2 = fig.add_subplot(gs[0, 2])
    energies = dos_data["energies"]
    ax2.plot(dos_data["total_dos"], energies, "k-", linewidth=2)
    ax2.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax2.set_xlabel("DOS (states/eV)")
    ax2.set_ylabel("Energy (eV)")
    ax2.set_title("Total DOS")
    ax2.grid(True, alpha=0.3)

    # Match y-axis with band structure
    y_min = min(np.min(band_energies), np.min(energies))
    y_max = max(np.max(band_energies), np.max(energies))
    ax1.set_ylim(y_min, y_max)
    ax2.set_ylim(y_min, y_max)

    # Panel 3: Orbital projected DOS (middle right)
    ax3 = fig.add_subplot(gs[1, 2])
    colors = ["red", "blue", "green"]
    for i, (orbital, color) in enumerate(zip(["s", "p", "d"], colors)):
        if orbital in dos_data["orbital_dos"]:
            ax3.plot(
                dos_data["orbital_dos"][orbital],
                energies,
                color=color,
                linewidth=2,
                label=f"{orbital}",
            )
    ax3.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax3.set_xlabel("DOS (states/eV)")
    ax3.set_ylabel("Energy (eV)")
    ax3.set_title("Orbital DOS")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(y_min, y_max)

    # Panel 4: Spin DOS (bottom right)
    ax4 = fig.add_subplot(gs[2, 2])
    if "up" in dos_data["spin_dos"] and "down" in dos_data["spin_dos"]:
        ax4.plot(dos_data["spin_dos"]["up"], energies, "r-", linewidth=2, label="↑")
        ax4.plot(-dos_data["spin_dos"]["down"], energies, "b-", linewidth=2, label="↓")
        ax4.axvline(x=0, color="k", linestyle="-", alpha=0.3)
    ax4.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax4.set_xlabel("DOS (states/eV)")
    ax4.set_ylabel("Energy (eV)")
    ax4.set_title("Spin DOS")
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(y_min, y_max)

    # Panel 5: Projected bands (bottom left)
    ax5 = fig.add_subplot(gs[2, :2])

    # Plot background bands
    for band in range(band_energies.shape[1]):
        ax5.plot(kdistances, band_energies[:, band], "gray", alpha=0.3, linewidth=0.5)

    # Overlay orbital projections if available
    if "projections" in line_data and line_data["projections"] is not None:
        projections = line_data["projections"]
        orbital_colors = {"s": "red", "p": "blue", "d": "green"}

        for band in range(band_energies.shape[1]):
            for orbital_group, color in orbital_colors.items():
                if orbital_group in projections["orbital"]:
                    weights = projections["orbital"][orbital_group][:, band]
                    mask = weights > 0.1
                    if np.any(mask):
                        sizes = weights[mask] * 30
                        ax5.scatter(
                            kdistances[mask],
                            band_energies[mask, band],
                            s=sizes,
                            c=color,
                            alpha=0.7,
                        )

    ax5.set_xticks(tick_positions)
    ax5.set_xticklabels(tick_labels)
    for pos in tick_positions:
        ax5.axvline(x=pos, color="k", linestyle="--", alpha=0.3)
    ax5.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax5.set_xlabel("k-path")
    ax5.set_ylabel("Energy (eV)")
    ax5.set_title("Orbital Projected Bands")
    ax5.grid(True, alpha=0.3)
    ax5.set_ylim(y_min, y_max)

    plt.suptitle("TB-LMTO Electronic Structure Analysis", fontsize=16)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Combined analysis saved to {output_file}")

    plt.show()


def main():
    """Main analysis function"""

    parser = argparse.ArgumentParser(
        description="Complete TB-LMTO electronic structure analysis"
    )
    parser.add_argument(
        "--lattice", default="lattice.nml", help="Lattice file (default: lattice.nml)"
    )
    parser.add_argument(
        "--hamiltonian",
        default="rs2txt.dat",
        help="Hamiltonian file (default: rs2txt.dat)",
    )
    parser.add_argument(
        "--output-prefix",
        default="complete_analysis",
        help="Output file prefix (default: complete_analysis)",
    )

    # Line band structure options
    parser.add_argument(
        "--nk-line",
        type=int,
        default=100,
        help="Number of k-points for line (default: 100)",
    )

    # DOS options
    parser.add_argument(
        "--dos-mesh",
        nargs=3,
        type=int,
        default=[6, 6, 6],
        help="k-point mesh for DOS (default: 6 6 6)",
    )
    parser.add_argument(
        "--energy-range", nargs=2, type=float, help="Energy range for DOS (emin emax)"
    )
    parser.add_argument(
        "--num-energies",
        type=int,
        default=500,
        help="Number of energy points for DOS (default: 500)",
    )
    parser.add_argument(
        "--fermi-energy", type=float, default=0.0, help="Fermi energy (default: 0.0)"
    )

    # Analysis options
    parser.add_argument(
        "--skip-dos", action="store_true", help="Skip DOS calculation (faster)"
    )
    parser.add_argument(
        "--skip-projections", action="store_true", help="Skip projection calculations"
    )
    parser.add_argument(
        "--no-symmetry", action="store_true", help="Disable symmetry reduction for DOS"
    )

    args = parser.parse_args()

    # Check if files exist
    if not Path(args.lattice).exists():
        print(f"Error: Lattice file {args.lattice} not found")
        return

    if not Path(args.hamiltonian).exists():
        print(f"Error: Hamiltonian file {args.hamiltonian} not found")
        return

    print("=" * 60)
    print("TB-LMTO Complete Electronic Structure Analysis")
    print("=" * 60)

    print("\n1. Reading input files...")
    lattice_reader = LatticeReader(args.lattice)
    hamiltonian_reader = HamiltonianReader(args.hamiltonian)

    # Structure analysis will be performed later when generating the k-path
    # (StructureAnalyzer is created by the band-structure code), so we avoid
    # instantiating it here to prevent duplicate output.

    print(f"   Crystal structure: {lattice_reader.nbas} atoms")
    print(f"   Lattice vectors:\n{lattice_reader.lattice_vectors}")

    print("\n2. Calculating line band structure...")
    if args.skip_projections:
        calculator = BandStructureCalculator(lattice_reader, hamiltonian_reader)
        kpoints, klabels, kdistances = calculator.setup_kpath(args.nk_line)
        band_energies = calculator.calculate_bands(kpoints)
        projections = None
    else:
        calculator = ProjectedBandCalculator(lattice_reader, hamiltonian_reader)
        kpoints, klabels, kdistances = calculator.setup_kpath(args.nk_line)
        band_energies, projections = calculator.calculate_projected_bands(kpoints)

    line_data = {
        "kdistances": kdistances,
        "band_energies": band_energies,
        "klabels": klabels,
        "kpoints": kpoints,
        "projections": projections,
    }

    # Find band gap
    gap_info = ""
    if projections:
        # Simple gap estimate
        fermi_bands = []
        for band in range(band_energies.shape[1]):
            if np.any(np.abs(band_energies[:, band]) < 1.0):
                fermi_bands.append(band)

        if len(fermi_bands) > 0:
            vb_max = np.max(band_energies[:, : fermi_bands[-1] // 2])
            cb_min = np.min(band_energies[:, fermi_bands[-1] // 2 :])
            gap = cb_min - vb_max
            gap_info = f"   Estimated band gap: {gap:.3f} eV"

    print(f"   Calculated {band_energies.shape[1]} bands over {len(kpoints)} k-points")
    if gap_info:
        print(gap_info)

    dos_data = None
    if not args.skip_dos:
        print(f"\n3. Calculating full BZ DOS (mesh: {args.dos_mesh})...")

        full_bz_calc = FullBZCalculator(lattice_reader, hamiltonian_reader)
        dos_calc = DOSCalculator(full_bz_calc)

        energy_range = tuple(args.energy_range) if args.energy_range else None

        dos_data = dos_calc.calculate_dos(
            mesh_size=tuple(args.dos_mesh),
            energy_range=energy_range,
            num_energies=args.num_energies,
            use_symmetry=not args.no_symmetry,
            use_tetrahedron=True,
        )

        print(f"   DOS calculated over {len(dos_data['kpoints'])} k-points")
        print(
            f"   Energy range: {dos_data['energies'][0]:.2f} to {dos_data['energies'][-1]:.2f} eV"
        )

    print("\n4. Creating plots...")

    # Individual plots
    calculator_for_plot = (
        calculator
        if hasattr(calculator, "plot_bands")
        else BandStructureCalculator(lattice_reader, hamiltonian_reader)
    )
    calculator_for_plot.plot_bands(
        kdistances, band_energies, klabels, f"{args.output_prefix}_bands.png"
    )

    if projections and not args.skip_projections:
        calculator.plot_projection_analysis(
            kdistances,
            band_energies,
            projections,
            klabels,
            output_file=f"{args.output_prefix}_projections.png",
        )

    if dos_data:
        dos_calc.plot_dos(
            dos_data,
            fermi_energy=args.fermi_energy,
            output_file=f"{args.output_prefix}_dos.png",
        )

        # Combined plot
        plot_combined_analysis(
            line_data, dos_data, f"{args.output_prefix}_complete.png"
        )

        # Save DOS data
        output_data = {
            "energies": dos_data["energies"],
            "total_dos": dos_data["total_dos"],
            "orbital_dos": dos_data["orbital_dos"],
            "site_dos": dos_data["site_dos"],
            "spin_dos": dos_data["spin_dos"],
        }
        np.savez(f"{args.output_prefix}_dos_data.npz", **output_data)
        print(f"   DOS data saved to {args.output_prefix}_dos_data.npz")

    print("\n5. Analysis summary:")
    print(f"   - Band structure: {args.output_prefix}_bands.png")
    if projections:
        print(f"   - Projected bands: {args.output_prefix}_projections.png")
    if dos_data:
        print(f"   - DOS analysis: {args.output_prefix}_dos.png")
        print(f"   - Complete analysis: {args.output_prefix}_complete.png")

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
