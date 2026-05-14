"""
Command-line interface for spin analysis.

Runs spin structure analysis and plotting for non-collinear TB-LMTO inputs.

Author: Ramon Cardias
"""

import argparse
import numpy as np

from ...projections import ProjectedBandCalculator
from ...readers import LatticeReader, HamiltonianReader
from .spin_analysis import SpinMomentAnalyzer


def main() -> None:
    """Main function for spin analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze non-collinear TB-LMTO spin structure"
    )
    parser.add_argument(
        "--lattice",
        default="lattice.nml",
        help="Lattice file (default: lattice.nml)",
    )
    parser.add_argument(
        "--hamiltonian",
        default="rs2txt.dat",
        help="Hamiltonian file (default: rs2txt.dat)",
    )
    parser.add_argument(
        "--output-prefix",
        default="spin_analysis",
        help="Output file prefix (default: spin_analysis)",
    )
    parser.add_argument(
        "--nk",
        type=int,
        default=50,
        help="Number of k-points (default: 50)",
    )

    args = parser.parse_args()

    print("Setting up projected band calculator...")
    lattice_reader = LatticeReader(args.lattice)
    hamiltonian_reader = HamiltonianReader(args.hamiltonian)
    calculator = ProjectedBandCalculator(lattice_reader, hamiltonian_reader)

    print("Generating k-path...")
    kpoints, klabels, kdistances = calculator.setup_kpath(args.nk)

    print("Calculating projected bands...")
    band_energies, projections = calculator.calculate_projected_bands(kpoints)

    print("Analyzing spin structure...")
    spin_analyzer = SpinMomentAnalyzer(calculator)

    # Calculate spin moments
    spin_moments = spin_analyzer.calculate_spin_moments(
        band_energies, projections, kpoints
    )

    # Plot spin texture
    spin_analyzer.plot_spin_texture(
        kdistances,
        band_energies,
        spin_moments,
        klabels,
        output_file=f"{args.output_prefix}_texture.png",
    )

    # Plot orbital magnetism
    spin_analyzer.plot_orbital_magnetism(
        kdistances,
        band_energies,
        projections,
        klabels,
        output_file=f"{args.output_prefix}_orbital_mag.png",
    )

    # Analyze magnetic anisotropy
    print("Analyzing magnetic anisotropy...")
    anisotropy = spin_analyzer.analyze_magnetic_anisotropy(kpoints)

    print("Magnetic anisotropy analysis:")
    print(f"Average H0 trace: {np.mean(anisotropy['H0_trace']):.6f}")
    print(f"Average Hz trace: {np.mean(anisotropy['Hz_trace']):.6f}")
    print(f"Average Hx trace: {np.mean(anisotropy['Hx_trace']):.6f}")
    print(f"Average Hy trace: {np.mean(anisotropy['Hy_trace']):.6f}")

    print("Done!")


if __name__ == "__main__":
    main()
