#!/usr/bin/env python3
"""
Example: Basic band structure calculation with rs_band

This example demonstrates how to use rs_band to calculate 
band structure from TB-LMTO data.
"""

import rs_band
import numpy as np
import matplotlib.pyplot as plt

def main():
    """Run basic band structure example"""
    
    # Check if data files exist
    try:
        # Load lattice and Hamiltonian data
        print("Loading TB-LMTO data...")
        lattice = rs_band.LatticeReader('lattice.nml')
        hamiltonian = rs_band.HamiltonianReader('rs2txt.dat')
        
        print(f"Crystal structure: {lattice.nbas} atoms")
        print(f"Lattice vectors:\n{lattice.lattice_vectors}")
        
        # Calculate band structure
        print("Calculating band structure...")
        calculator = rs_band.BandStructureCalculator(lattice, hamiltonian)
        
        # Generate k-path
        kpoints, klabels, kdistances = calculator.setup_kpath(100)
        print(f"Generated k-path with {len(kpoints)} points")
        
        # Calculate bands
        band_energies = calculator.calculate_bands(kpoints)
        print(f"Calculated {band_energies.shape[1]} bands")
        
        # Plot results
        calculator.plot_bands(kdistances, band_energies, klabels, 'example_bands.png')
        print("Band structure saved to example_bands.png")
        
        # Optional: Calculate projections
        print("Calculating orbital projections...")
        proj_calculator = rs_band.ProjectedBandCalculator(lattice, hamiltonian)
        band_energies, projections = proj_calculator.calculate_projected_bands(kpoints)
        
        # Plot projected bands
        proj_calculator.plot_projected_bands(
            kdistances, band_energies, projections, klabels,
            projection_type='orbital', output_file='example_projections.png'
        )
        print("Projected bands saved to example_projections.png")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure lattice.nml and rs2txt.dat are in the current directory")
        print("You can find example files in the rs_band examples directory")

if __name__ == "__main__":
    main()
