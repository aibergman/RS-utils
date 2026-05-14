#!/usr/bin/env python3
"""
TB-LMTO Band Structure Calculator

This program reads TB-LMTO Hamiltonian data and lattice information to:
1. Parse the real-space Hamiltonian from rs2txt.dat
2. Read crystal structure from lattice.nml
3. Use spglib to determine crystal symmetry and generate k-path
4. Fourier transform the Hamiltonian to k-space
5. Diagonalize and calculate the band structure
6. Plot the results

Authors: Original code by rslmto_official
         KPath enhancements by Ramon Cardias
"""

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..readers import HamiltonianReader, LatticeReader
from .band_calculations import (
    extract_spin_blocks as extract_spin_blocks_fn,
    calculate_bands as calculate_bands_fn,
    calculate_bands_collinear as calculate_bands_collinear_fn,
)
from .hamiltonian import fourier_transform_hamiltonian, HamiltonianDumper
from .plotting import plot_bands, plot_bands_collinear
from .kpath_generator import (
    setup_kpath,
    setup_kpath_seekpath,
    _setup_kpath_spglib,
    _setup_kpath_pymatgen,
    _setup_kpath_manual,
    _generate_kpath,
)


class BandStructureCalculator:
    """Calculate band structure using TB-LMTO Hamiltonian

    Supports both non-collinear (default) and collinear spin-polarized calculations.

    Parameters
    ----------
    lattice_reader : LatticeReader
        Reader for lattice structure
    hamiltonian_reader : HamiltonianReader
        Reader for Hamiltonian data
    collinear : bool, optional
        If True, separate spin-up and spin-down bands for collinear magnetism.
        If False (default), use full non-collinear Hamiltonian.

    Attributes
    ----------
    collinear : bool
        Whether to use collinear spin treatment
    """

    def __init__(
        self,
        lattice_reader: LatticeReader,
        hamiltonian_reader: HamiltonianReader,
        collinear: bool = False,
    ):
        self.lattice = lattice_reader
        self.hamiltonian = hamiltonian_reader
        self.collinear = collinear
        self.kpath = None
        self.band_energies = None
        self.band_energies_up = None
        self.band_energies_down = None
        # Eigenvector storage (Phase 3)
        self.eigenvectors = None
        self.eigenvectors_up = None
        self.eigenvectors_down = None

    def setup_kpath(
        self, num_points: int = 100
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-point path using structure analysis.
        
        Delegates to kpath_generator.setup_kpath() for implementation.
        """
        return setup_kpath(self.lattice, num_points)

    def setup_kpath_seekpath(
        self, num_points: int = 100, interactive: bool = False
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-point path using seekpath with optional interactive selection.
        
        Delegates to kpath_generator.setup_kpath_seekpath() for implementation.
        """
        return setup_kpath_seekpath(self.lattice, num_points, interactive)

    def _setup_kpath_spglib(
        self, num_points: int
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-path using spglib.
        
        Delegates to kpath_generator._setup_kpath_spglib() for implementation.
        """
        return _setup_kpath_spglib(self.lattice, num_points)

    def _setup_kpath_pymatgen(
        self, num_points: int
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Setup k-path using pymatgen.
        
        Delegates to kpath_generator._setup_kpath_pymatgen() for implementation.
        """
        return _setup_kpath_pymatgen(self.lattice, num_points)

    def _setup_kpath_manual(
        self, num_points: int
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Fallback manual k-path setup.
        
        Delegates to kpath_generator._setup_kpath_manual() for implementation.
        """
        return _setup_kpath_manual(self.lattice, num_points)

    def _generate_kpath(
        self, special_points: Dict, path: List[str], num_points: int
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Generate k-point path from special points.
        
        Delegates to kpath_generator._generate_kpath() for implementation.
        """
        return _generate_kpath(special_points, path, num_points)

    def extract_spin_blocks(self, H_full: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Extract spin-up and spin-down blocks from full Hamiltonian.

        Delegates to band_calculations.extract_spin_blocks() for implementation.
        """
        return extract_spin_blocks_fn(H_full, self.lattice.nbas)

    def fourier_transform_hamiltonian(self, k_point: np.ndarray) -> np.ndarray:
        """Fourier transform the real-space Hamiltonian to k-space.
        
        Delegates to module function fourier_transform_hamiltonian() for implementation.
        """
        return fourier_transform_hamiltonian(k_point, self.lattice, self.hamiltonian)

    def dump_hamiltonian_at_points(self, output_prefix: str = "H_k"):
        """Dump Hamiltonian matrices at Gamma and H points for BCC structure.
        
        Delegates to HamiltonianDumper for implementation.
        """
        dumper = HamiltonianDumper(self.lattice, self.hamiltonian)
        return dumper.dump_at_points(output_prefix)

    def dump_hamiltonian_at_custom_points(
        self, k_points: List[Tuple[str, np.ndarray]], output_prefix: str = "H_k"
    ):
        """Dump Hamiltonian matrices at custom k-points.
        
        Delegates to HamiltonianDumper for implementation.
        """
        dumper = HamiltonianDumper(self.lattice, self.hamiltonian)
        return dumper.dump_at_custom_points(k_points, output_prefix)



    def calculate_bands(
        self, kpoints: np.ndarray, store_eigenvectors: bool = False
    ) -> np.ndarray:
        """Calculate band structure for given k-points.

        Delegates to band_calculations.calculate_bands() for implementation.
        """
        if self.collinear:
            raise ValueError(
                "For collinear calculations, use calculate_bands_collinear() "
                "which returns separate up and down bands"
            )

        band_energies, eigenvectors = calculate_bands_fn(
            kpoints,
            self.lattice,
            self.hamiltonian,
            self.fourier_transform_hamiltonian,
            store_eigenvectors=store_eigenvectors,
        )

        if store_eigenvectors:
            self.eigenvectors = eigenvectors

        return band_energies

    def calculate_bands_collinear(
        self, kpoints: np.ndarray, store_eigenvectors: bool = False
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate spin-polarized band structure for collinear magnetism.

        Delegates to band_calculations.calculate_bands_collinear() for implementation.
        """
        if not self.collinear:
            raise ValueError(
                "calculate_bands_collinear() requires collinear=True. "
                "Use calculate_bands() for non-collinear calculations."
            )

        band_energies_up, band_energies_down, eigenvectors_up, eigenvectors_down = (
            calculate_bands_collinear_fn(
                kpoints,
                self.lattice,
                self.fourier_transform_hamiltonian,
                store_eigenvectors=store_eigenvectors,
            )
        )

        self.band_energies_up = band_energies_up
        self.band_energies_down = band_energies_down

        if store_eigenvectors:
            self.eigenvectors_up = eigenvectors_up
            self.eigenvectors_down = eigenvectors_down

        return band_energies_up, band_energies_down

    def plot_bands(
        self,
        kdistances: np.ndarray,
        band_energies: np.ndarray,
        klabels: List[str],
        output_file: Optional[str] = None,
        energy_range: Optional[Tuple[float, float]] = None,
    ):
        """Plot the band structure with publication-quality styling.
        
        Delegates to plotting.plot_bands() for implementation.
        """
        return plot_bands(kdistances, band_energies, klabels, output_file, energy_range)

    def plot_bands_collinear(
        self,
        kdistances: np.ndarray,
        band_energies_up: np.ndarray,
        band_energies_down: np.ndarray,
        klabels: List[str],
        output_file: Optional[str] = None,
        energy_range: Optional[Tuple[float, float]] = None,
    ):
        """Plot spin-polarized band structure with publication-quality styling.
        
        Delegates to plotting.plot_bands_collinear() for implementation.
        """
        return plot_bands_collinear(kdistances, band_energies_up, band_energies_down, klabels, output_file, energy_range)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Calculate TB-LMTO band structure")
    parser.add_argument(
        "--lattice", default="lattice.nml", help="Lattice file (default: lattice.nml)"
    )
    parser.add_argument(
        "--hamiltonian",
        default="rs2txt.dat",
        help="Hamiltonian file (default: rs2txt.dat)",
    )
    parser.add_argument(
        "--output",
        default="band_structure.png",
        help="Output plot file (default: band_structure.png)",
    )
    parser.add_argument(
        "--nk", type=int, default=100, help="Number of k-points (default: 100)"
    )
    parser.add_argument(
        "--dump-hamiltonian",
        action="store_true",
        help="Dump H(k) matrices at Gamma and H points",
    )

    args = parser.parse_args()

    # Check if files exist
    if not Path(args.lattice).exists():
        print(f"Error: Lattice file {args.lattice} not found")
        return

    if not Path(args.hamiltonian).exists():
        print(f"Error: Hamiltonian file {args.hamiltonian} not found")
        return

    print("Reading lattice information...")
    lattice_reader = LatticeReader(args.lattice)

    print("Reading Hamiltonian data...")
    hamiltonian_reader = HamiltonianReader(args.hamiltonian)

    print("Setting up band structure calculator...")
    calculator = BandStructureCalculator(lattice_reader, hamiltonian_reader)

    print(f"Generating k-path with {args.nk} points...")
    kpoints, klabels, kdistances = calculator.setup_kpath(args.nk)

    print("Calculating band structure...")
    band_energies = calculator.calculate_bands(kpoints)

    if args.dump_hamiltonian:
        print("Dumping Hamiltonian matrices...")
        calculator.dump_hamiltonian_at_points()

    print("Plotting results...")
    calculator.plot_bands(kdistances, band_energies, klabels, args.output)

    print("Done!")


if __name__ == "__main__":
    main()
