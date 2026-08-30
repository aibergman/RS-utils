"""
This module provides a command-line interface (CLI) for the rsutil package.
It allows users to process lattice structures, visualize data, create
exchange input files, and perform band structure/DOS analysis through
various command-line arguments and subcommands.

Functions
---------
main()
    Sets up the command-line argument parser and handles the execution of
    different functionalities based on the provided arguments.

Usage
-----
Existing functionality (flags):
    rsutil --setup
    rsutil --lattice --data <data_directory>
    rsutil --clust --data <data_directory>
    rsutil --exchange
    rsutil --export

New subcommand functionality:
    rsutil bands [options]          # Band structure calculation
    rsutil dos [options]            # DOS calculation
    rsutil moments [options]        # DOS moments analysis
    rsutil projections [options]    # Projected bands
    rsutil complete [options]       # Complete analysis

Flags (existing functionality)
-------------------------------
--setup, -S
    Process lattice.
--lattice, -L
    Run visualizer for lattice input.
--clust, -C
    Run visualizer for clust file.
--exchange, -J
    Create exchange input files.
--export, -E
    Export to UppASD format.
--data
    Directory for data files to be visualized (default: current directory).
"""

import argparse
import os
import numpy as np
import sys
from pathlib import Path
from .viz.rs_visualizer import rs_visualizer
from .exchange.exchange_map import exchange_map
from .export.uppasd import export_from_files
from .export.lmto47 import export_ctrl
from .core.contracts import Diagnostic, DiagnosticSeverity, ErrorCode, OperationStatus, SetupResult, SetupSpec
from .setup_service import SetupService
from .setup_service import validate_setup_inputs as service_validate_setup_inputs
from .setup_service import validate_setup_parameters as service_validate_setup_parameters
from .setup_validation import write_setup_report


def _make_common_parser():
    """Create common parser for band/DOS subcommands"""
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument(
        "--lattice", default="lattice.nml", help="Lattice file (default: lattice.nml)"
    )
    p.add_argument(
        "--hamiltonian",
        default="rs2txt.dat",
        help="Hamiltonian file (default: rs2txt.dat)",
    )
    p.add_argument(
        "--lattice-format",
        choices=["rs-lmto", "paoflow", "kite"],
        default=None,
        help="Lattice file format (auto-detect if not specified)",
    )
    p.add_argument(
        "--hamiltonian-format",
        choices=["rs-lmto", "paoflow", "kite"],
        default=None,
        help="Hamiltonian file format (auto-detect if not specified)",
    )
    return p


def run_bands_subcommand(args):
    """Run band structure calculation, projections, altermagnetic, or Berry curvature analysis"""
    from .bands import load_lattice, load_hamiltonian, BandStructureCalculator

    print(f"Reading lattice from {args.lattice}...")
    lattice = load_lattice(args.lattice, format=args.lattice_format)

    print(f"Reading Hamiltonian from {args.hamiltonian}...")
    hamiltonian = load_hamiltonian(args.hamiltonian, format=args.hamiltonian_format)

    # Helper function to setup k-path based on mode
    def setup_kpath_from_args(calc, args):
        """Setup k-path based on --kpath-mode and related arguments"""
        kpath_mode = getattr(args, "kpath_mode", "auto")
        
        if kpath_mode == "seekpath":
            print(f"Setting up k-path with seekpath (symprec={args.seekpath_symprec}, distance={args.seekpath_distance})...")
            try:
                kpoints, klabels, kdistances = calc.setup_kpath_seekpath(
                    num_points=args.nk,
                    interactive=False
                )
                return kpoints, klabels, kdistances
            except ImportError:
                print("Warning: seekpath not available, falling back to auto mode")
                return calc.setup_kpath(args.nk)
        
        elif kpath_mode == "custom":
            if not args.custom_path:
                print("Error: --custom-path required with --kpath-mode custom")
                raise ValueError("--custom-path not specified")
            
            from .bands import seekpath_kpath, build_custom_kpath
            
            print(f"Setting up seekpath data for custom path...")
            kpath_data = seekpath_kpath(
                lattice_vectors=lattice.lattice_vectors,
                atomic_positions=lattice.atomic_positions,
                atomic_numbers=lattice.atomic_numbers,
                num_points=args.nk,
                symprec=args.seekpath_symprec,
                reference_distance=args.seekpath_distance
            )
            
            # Show available symmetry points
            print("\nAvailable symmetry points:")
            point_coords = kpath_data['point_coords']
            available_points = sorted(point_coords.keys())
            for point in available_points:
                coords = point_coords[point]
                print(f"  {point:6s}: ({coords[0]:7.4f}, {coords[1]:7.4f}, {coords[2]:7.4f})")
            
            labels = args.custom_path.split('-')
            print(f"\nBuilding custom k-path: {' -> '.join(labels)}")
            kpoints_frac, xvals, xticks, xticklabels, metadata = build_custom_kpath(
                point_coords=kpath_data['point_coords'],
                reciprocal_lattice=kpath_data['reciprocal_lattice'],
                labels_list=labels,
                total_kpoints=args.nk
            )
            
            # keep k-points in fractional coordinates for calculate_bands()
            kpoints = kpoints_frac
            kdistances = xvals
            klabels_list = xticklabels
            
            # Store for later use (export, plotting)
            calc._kpath_data = {
                'metadata': metadata,
                'reciprocal_lattice': kpath_data['reciprocal_lattice'],
                'point_coords': kpath_data['point_coords']
            }
            
            return kpoints, klabels_list, kdistances
        
        elif kpath_mode == "interactive":
            print("Interactive k-path selection mode...")
            try:
                kpoints, klabels, kdistances = calc.setup_kpath_seekpath(
                    num_points=args.nk,
                    interactive=True
                )
                return kpoints, klabels, kdistances
            except ImportError:
                print("Warning: seekpath not available for interactive mode")
                return calc.setup_kpath(args.nk)
        
        else:  # auto mode (default)
            print(f"Generating k-path with {args.nk} points (auto mode)...")
            return calc.setup_kpath(args.nk)

    # Helper function to handle BZ plotting and k-path export
    def handle_kpath_output(args, lattice, kpoints, klabels, kdistances, kpath_data=None):
        """Plot Brillouin zone and export k-path if requested"""
        if args.plot_bz:
            try:
                from .bands import seekpath_kpath, plot_brillouin_zone
                
                # Get BZ data if not provided
                if kpath_data is None:
                    kpath_data = seekpath_kpath(
                        lattice_vectors=lattice.lattice_vectors,
                        atomic_positions=lattice.atomic_positions,
                        atomic_numbers=lattice.atomic_numbers,
                        symprec=args.seekpath_symprec,
                        reference_distance=args.seekpath_distance
                    )
                
                print(f"Plotting Brillouin zone to {args.output_prefix}_bz.png...")
                fig, ax = plot_brillouin_zone(
                    reciprocal_lattice=kpath_data['reciprocal_lattice'],
                    point_coords=kpath_data['point_coords'],
                    label_hs_points=True,
                    label_primes=args.prime_points,
                    figsize=(10, 10)
                )
                fig.savefig(f"{args.output_prefix}_bz.png", dpi=300, bbox_inches='tight')
                print("Brillouin zone plot saved.")
            except ImportError as e:
                print(f"Warning: Could not plot Brillouin zone: {e}")
        
        if args.export_kpath:
            from .bands import export_kpath
            
            print(f"Exporting k-path to {args.export_kpath}...")
            export_kpath(
                filename=args.export_kpath,
                kpath_x=kdistances,
                kpoints_abs=kpoints,
                kpoint_labels=klabels,
                header='x_linearcoord  kx  ky  kz'
            )
            print(f"K-path exported to {args.export_kpath}")

    # Check for altermagnetic analysis
    if args.altermagnetic:
        from .bands.altermagnetic import AltermagneticAnalyzer

        # Initialize band structure calculator with collinear support
        print("Setting up band structure calculator (collinear)...")
        calc = BandStructureCalculator(lattice, hamiltonian, collinear=True)

        kpoints, klabels, kdistances = setup_kpath_from_args(calc, args)

        print("Calculating spin-polarized bands...")
        bands_up, bands_down = calc.calculate_bands_collinear(kpoints)

        # Create altermagnetic analyzer
        print("Running altermagnetic analysis...")
        analyzer = AltermagneticAnalyzer(
            bands_up, bands_down, kpoints, kdistances, fermi_energy=args.fermi_energy
        )

        # Get identification result
        result = analyzer.identify_altermagnetic_bands(
            magnetization_threshold=args.mag_threshold
        )

        # Print summary
        print()
        analyzer.print_altermagnetic_summary()

        # Handle BZ plotting and k-path export
        handle_kpath_output(args, lattice, kpoints, klabels, kdistances)

        # Generate plots if requested
        if args.plot_splitting:
            output_file = f"{args.output_prefix}_splitting.png"
            print(f"\nPlotting spin splitting map to {output_file}...")
            analyzer.plot_spin_splitting_map(output_file)

        if args.plot_character:
            output_file = f"{args.output_prefix}_character.png"
            print(f"Plotting altermagnetic character to {output_file}...")
            analyzer.plot_altermagnetic_character(output_file)

        # Save band structure plot
        if args.plot_bands:
            output_file = f"{args.output_prefix}_bands.png"
            print(f"Plotting collinear bands to {output_file}...")
            calc.plot_bands_collinear(
                kdistances, bands_up, bands_down, klabels, output_file
            )

        print("\nDone!")

    # Check for Berry curvature analysis
    elif args.berry:
        from .bands.berry_curvature import BerryCurvatureCalculator

        print("Setting up band structure calculator...")
        calc = BandStructureCalculator(lattice, hamiltonian)

        kpoints, klabels, kdistances = setup_kpath_from_args(calc, args)

        print("Calculating band structure with eigenvector storage...")
        energies = calc.calculate_bands(kpoints, store_eigenvectors=True)

        # Use real eigenvectors from Hamiltonian diagonalization (Phase 3)
        eigenvectors = calc.eigenvectors
        print(f"Eigenvector shape: {eigenvectors.shape}")

        # Initialize Berry curvature calculator
        print("Computing Berry curvature...")
        a = lattice.lattice_vectors
        reciprocal_lattice = (
            2 * np.pi * np.linalg.inv(a.T)
        )  # Reciprocal lattice vectors as columns

        berry_calc = BerryCurvatureCalculator(
            energies,
            eigenvectors,
            kpoints,
            reciprocal_lattice,
            fermi_energy=args.fermi_energy,
        )

        # Identify topological bands
        print("\nAnalyzing topological properties...")
        topology = berry_calc.identify_topological_bands(chern_tolerance=args.chern_tol)

        print("\nTopological Band Classification:")
        print(f"  Trivial bands: {topology['trivial_bands']}")
        print(f"  Non-trivial bands: {topology['nontrivial_bands']}")

        # Print band-resolved properties
        if args.num_bands:
            print(
                f"\nBand-Resolved Properties (first {min(args.num_bands, nbands)} bands):"
            )
            berry_calc.print_berry_summary(num_bands=min(args.num_bands, nbands))

        # Calculate Chern numbers for selected bands
        if args.chern_bands:
            print("\nChern Numbers:")
            for band in args.chern_bands:
                if band < nbands:
                    chern = berry_calc.calculate_chern_number(band)
                    print(f"  Band {band}: C = {chern:.4f}")

        # Calculate anomalous Hall conductivity
        if args.ahc:
            print("\nAnomalous Hall Conductivity:")
            for band in range(min(5, nbands)):
                sigma_xy = berry_calc.calculate_anomalous_hall_conductivity(
                    band, temperature=args.temperature
                )
                print(f"  Band {band}: σ_xy = {sigma_xy:.4e} (e²/h)")

        # Calculate orbital magnetization
        if args.orbital_mag:
            print("\nOrbital Magnetization:")
            for band in range(min(5, nbands)):
                m_orb = berry_calc.calculate_orbital_magnetization(
                    band, temperature=args.temperature
                )
                print(f"  Band {band}: M_orb = {m_orb:.4e} μ_B")

        # Handle BZ plotting and k-path export
        handle_kpath_output(args, lattice, kpoints, klabels, kdistances)

        # Generate plots if requested
        if args.plot_berry:
            band = min(0, args.plot_berry - 1)  # Default to band 0
            output_file = f"{args.output_prefix}_berry_band{band}.png"
            print(f"\nPlotting Berry curvature for band {band} to {output_file}...")
            berry_calc.plot_berry_curvature_map(band, output_file)

        if args.plot_chern:
            output_file = f"{args.output_prefix}_chern.png"
            print(f"Plotting Chern numbers to {output_file}...")
            berry_calc.plot_chern_number_histogram(output_file)

        # Save band structure plot
        if args.plot_bands:
            output_file = f"{args.output_prefix}_bands.png"
            print(f"Plotting band structure to {output_file}...")
            calc.plot_bands(kdistances, energies, klabels, output_file)

        print("\nDone!")

    # Check for spin texture analysis (Phase 4)
    elif args.spin_texture:
        from .bands.spin_texture import SpinTextureCalculator

        print("Setting up band structure calculator with eigenvector storage...")
        calc = BandStructureCalculator(lattice, hamiltonian)

        kpoints, klabels, kdistances = setup_kpath_from_args(calc, args)

        print("Calculating band structure with eigenvector storage...")
        energies = calc.calculate_bands(kpoints, store_eigenvectors=True)

        # Use real eigenvectors from Hamiltonian diagonalization (Phase 4)
        eigenvectors = calc.eigenvectors
        print(f"Eigenvector shape: {eigenvectors.shape}")

        # Initialize spin texture calculator
        print("Computing spin texture...")
        a = lattice.lattice_vectors
        reciprocal_lattice = (
            2 * np.pi * np.linalg.inv(a.T)
        )  # Reciprocal lattice vectors as columns

        spin_calc = SpinTextureCalculator(
            energies,
            eigenvectors,
            kpoints,
            reciprocal_lattice,
            fermi_energy=args.fermi_energy,
        )

        spin_texture = spin_calc.calculate_spin_texture()
        print(f"Spin texture shape: {spin_texture.shape}")
        print(f"Average spin polarization: {np.mean(spin_calc.spin_polarization):.4f}")

        # Analyze specific band if requested
        if args.spin_band is not None:
            band_idx = args.spin_band
            print(f"\nAnalyzing spin texture for band {band_idx}...")

            # Spin-momentum locking
            locking = spin_calc.identify_spin_momentum_locking(
                band_idx, energy_window=0.1
            )
            print(f"  Locking type: {locking['locking_type']}")
            print(f"  Locking strength: {locking['locking_strength']:.3f}")
            print(f"  Dominant spin orientation: {locking['spin_orientation']}")

            # Plot 2D spin texture
            output_file = f"{args.output_prefix}_spin_texture_band{band_idx}.png"
            print(f"Plotting 2D spin texture to {output_file}...")
            spin_calc.plot_spin_texture_2d(
                band_idx,
                component=args.spin_component,
                plane=args.spin_plane,
                output_file=output_file,
            )

            # Plot 3D Fermi surface if requested
            if args.plot_spin_fermi:
                output_file = f"{args.output_prefix}_fermi_spin_band{band_idx}.png"
                print(f"Plotting 3D Fermi surface to {output_file}...")
                spin_calc.plot_fermi_surface_spin_texture(
                    band_idx, output_file=output_file
                )

            # Rashba analysis if requested
            if args.rashba_analysis:
                print(f"\nExtracting Rashba parameters for band {band_idx}...")
                try:
                    rashba_params = spin_calc.extract_rashba_parameters(band_idx)
                    print(f"  αR = {rashba_params['alpha_R']:.4e} eV·Å")
                    print(f"  E0 = {rashba_params['E_0']:.4f} eV")
                    print(f"  Fit quality (R²) = {rashba_params['fit_quality']:.4f}")
                except Exception as e:
                    print(f"  Rashba analysis failed: {e}")

        # Handle BZ plotting and k-path export
        handle_kpath_output(args, lattice, kpoints, klabels, kdistances)

        print("\nDone!")

    # Check for fat bands visualization (Phase 4)
    elif args.fat_bands:
        from .bands.fat_bands import FatBandsPlotter

        print("Setting up band structure calculator with eigenvector storage...")
        calc = BandStructureCalculator(lattice, hamiltonian)

        kpoints, klabels, kdistances = setup_kpath_from_args(calc, args)

        print("Calculating band structure with eigenvector storage...")
        energies = calc.calculate_bands(kpoints, store_eigenvectors=True)
        eigenvectors = calc.eigenvectors

        # Initialize fat bands plotter
        print("Preparing fat bands visualization...")
        fat_plotter = FatBandsPlotter(
            energies, eigenvectors, kdistances, klabels, fermi_energy=args.fermi_energy
        )

        # Use energy range only if explicitly provided
        energy_range = tuple(args.energy_range) if args.energy_range else None

        # Plot based on type
        if args.fat_type == "spin":
            output_file = f"{args.output_prefix}_fat_bands_spin.png"
            print(f"Plotting spin-resolved fat bands to {output_file}...")
            fat_plotter.plot_spin_resolved_fat_bands(
                energy_range=energy_range, output_file=output_file
            )

        elif args.fat_type == "linewidth" and args.fat_orbitals:
            projections = fat_plotter.calculate_orbital_projections(args.fat_orbitals)
            output_file = f"{args.output_prefix}_fat_bands_orbitals.png"
            print(f"Plotting orbital fat bands (linewidth) to {output_file}...")
            fat_plotter.plot_fat_bands_linewidth(
                projections,
                scale_factor=args.fat_scale,
                energy_range=energy_range,
                output_file=output_file,
                title=f"Fat Bands - Orbitals {args.fat_orbitals}",
            )

        elif args.fat_type == "color" and args.fat_orbitals:
            projections = fat_plotter.calculate_orbital_projections(args.fat_orbitals)
            output_file = f"{args.output_prefix}_fat_bands_orbitals_color.png"
            print(f"Plotting orbital fat bands (color) to {output_file}...")
            fat_plotter.plot_fat_bands_color(
                projections,
                energy_range=energy_range,
                output_file=output_file,
                title=f"Fat Bands - Orbitals {args.fat_orbitals}",
            )

        elif args.fat_orbital_groups:
            import json

            orbital_groups = json.loads(args.fat_orbital_groups)
            output_file = f"{args.output_prefix}_fat_bands_groups.png"
            print(f"Plotting orbital character fat bands to {output_file}...")
            fat_plotter.plot_orbital_character(
                orbital_groups, energy_range=energy_range, output_file=output_file
            )

        else:
            print("Default: Plotting spin-resolved fat bands...")
            output_file = f"{args.output_prefix}_fat_bands_spin.png"
            fat_plotter.plot_spin_resolved_fat_bands(
                energy_range=energy_range, output_file=output_file
            )

        # Handle BZ plotting and k-path export
        handle_kpath_output(args, lattice, kpoints, klabels, kdistances)

        print("\nDone!")

    # Regular projections or band structure
    elif args.projections:
        # Run projected bands
        from .bands.projections import ProjectedBandCalculator

        print("Setting up projected band structure calculator...")
        # Create calculator with specified spin projection mode
        proj_calc = ProjectedBandCalculator(
            lattice,
            hamiltonian,
            spin_projection_mode=getattr(args, "projections_mode", "expectation"),
        )

        kpoints, klabels, kdistances = setup_kpath_from_args(proj_calc, args)

        print("Calculating projected band structure...")
        band_energies, projections = proj_calc.calculate_projected_bands(kpoints)

        projection_type = args.projection
        output_file = f"{args.output_prefix}_{projection_type}.png"
        print(
            f"Plotting {projection_type} projections (mode: {proj_calc.spin_projection_mode}) to {output_file}..."
        )
        proj_calc.plot_projected_bands(
            kdistances,
            band_energies,
            projections,
            klabels,
            projection_type,
            output_file,
        )

        # Handle BZ plotting and k-path export
        handle_kpath_output(args, lattice, kpoints, klabels, kdistances)

        print("Done!")
    else:
        # Run regular band structure
        print("Setting up band structure calculator...")
        calc = BandStructureCalculator(lattice, hamiltonian)

        kpoints, klabels, kdistances = setup_kpath_from_args(calc, args)

        print("Calculating band structure...")
        band_energies = calc.calculate_bands(kpoints)

        # Handle BZ plotting and k-path export
        handle_kpath_output(args, lattice, kpoints, klabels, kdistances)

        print(f"Plotting results to {args.output}...")
        calc.plot_bands(kdistances, band_energies, klabels, args.output)

        print("Done!")


def run_dos_subcommand(args):
    """Run DOS calculation or DOS moments analysis"""
    from .bands import load_lattice, load_hamiltonian, FullBZCalculator, DOSCalculator

    if args.moments:
        # Run DOS moments analysis
        from .bands.dos_moments import DOSMomentsCalculator

        lattice = load_lattice(args.lattice, format=args.lattice_format)
        hamiltonian = load_hamiltonian(args.hamiltonian, format=args.hamiltonian_format)
        bz_calc = FullBZCalculator(lattice, hamiltonian)
        dos_calc = DOSCalculator(bz_calc)

        print("Calculating DOS for moments analysis...")
        use_tetrahedron = not getattr(args, "no_tetrahedron", False)

        dos_data = dos_calc.calculate_dos(
            mesh_size=tuple(args.mesh),
            energy_range=tuple(args.energy_range) if args.energy_range else None,
            num_energies=args.num_energies,
            use_symmetry=not args.no_symmetry,
            use_tetrahedron=use_tetrahedron,
            gaussian_smearing=0.1,
        )

        # Perform moments analysis
        moments_calc = DOSMomentsCalculator(
            dos_data["energies"],
            dos_data["total_dos"],
            dos_data.get("orbital_dos"),
            temperature=args.temperature,
        )

        if args.electrons:
            print(f"Finding Fermi level for {args.electrons} electrons...")
            energy_range = tuple(args.energy_range) if args.energy_range else None
            fermi_level = moments_calc.find_fermi_level(args.electrons, energy_range)
            print(f"Fermi level: {fermi_level:.6f} eV")
        else:
            print("Auto-detecting electron count...")
            # Count electrons at zero energy
            total_electrons = moments_calc.total_dos.sum() * (
                moments_calc.energies[1] - moments_calc.energies[0]
            )
            fermi_level = moments_calc.find_fermi_level(total_electrons)
            print(f"Fermi level: {fermi_level:.6f} eV")

        # Calculate moments
        print("\nCalculating DOS moments...")
        moments = moments_calc.calculate_moments(fermi_level)

        # Print results
        print("\nDOS Moments:")
        if "total" in moments:
            print(f"\nTotal DOS:")
            total_moments = moments["total"]
            print(f"  Occupation (electrons):      {total_moments['occupation']:.4f}")
            print(
                f"  Avg Energy (rel to EF):      {total_moments['avg_energy_rel_ef']:.4f} eV"
            )
            print(
                f"  Energy Spread (width):       {total_moments['band_width']:.4f} eV"
            )
            print(
                f"  Band Center:                 {total_moments['band_center']:.4f} eV"
            )

        # Save results
        moments_calc.save_moments(f"{args.output_prefix}_moments.npz")
        print(f"\nResults saved to {args.output_prefix}_moments.npz")
        print("Done!")
    else:
        # Run regular DOS calculation
        lattice = load_lattice(args.lattice, format=args.lattice_format)
        hamiltonian = load_hamiltonian(args.hamiltonian, format=args.hamiltonian_format)
        bz_calc = FullBZCalculator(lattice, hamiltonian)
        dos_calc = DOSCalculator(bz_calc)

        # If electrons are specified, use Fermi level determination
        if args.electrons is not None:
            print(
                f"Calculating DOS with Fermi level determination for {args.electrons} electrons..."
            )
            energy_range_val = tuple(args.energy_range) if args.energy_range else None
            dos_data, moments_calc = dos_calc.calculate_dos_with_moments(
                mesh_size=tuple(args.mesh),
                target_electrons=args.electrons,
                temperature=getattr(args, "temperature", 300.0),
                energy_range=energy_range_val,
                num_energies=args.num_energies,
                use_symmetry=not args.no_symmetry,
                use_tetrahedron=not args.no_tetrahedron,
                gaussian_smearing=getattr(args, "gaussian_smearing", 0.1),
                output_prefix=args.output_prefix,
            )

            fermi_level = moments_calc.fermi_level
            if fermi_level is not None:
                print(f"Fermi level determined: {fermi_level:.6f} eV")
                dos_calc.plot_dos(
                    dos_data,
                    fermi_energy=fermi_level,
                    output_file=f"{args.output_prefix}.png",
                )
            else:
                print("Warning: Could not determine Fermi level")
                dos_calc.plot_dos(
                    dos_data, fermi_energy=0.0, output_file=f"{args.output_prefix}.png"
                )
        else:
            # Simple DOS calculation
            print("Calculating DOS...")
            dos_data = dos_calc.calculate_dos(
                mesh_size=tuple(args.mesh),
                energy_range=tuple(args.energy_range) if args.energy_range else None,
                num_energies=args.num_energies,
                use_symmetry=not args.no_symmetry,
                use_tetrahedron=not args.no_tetrahedron,
                gaussian_smearing=getattr(args, "gaussian_smearing", 0.1),
            )

            fermi_energy = getattr(args, "fermi_energy", 0.0)
            dos_calc.plot_dos(
                dos_data,
                fermi_energy=fermi_energy,
                output_file=f"{args.output_prefix}.png",
            )

        print("Done!")


def create_bands_parser(subparsers):
    """Create bands subcommand parser"""
    parser = subparsers.add_parser(
        "bands",
        parents=[_make_common_parser()],
        help="Calculate band structure and projections",
    )
    parser.add_argument(
        "--nk",
        type=int,
        default=100,
        help="Number of k-points along path (default: 100)",
    )
    parser.add_argument(
        "--output",
        default="band_structure.png",
        help="Output plot file (default: band_structure.png)",
    )
    parser.add_argument(
        "--projections",
        action="store_true",
        help="Calculate projected bands instead of total bands",
    )
    parser.add_argument(
        "--projection",
        choices=["orbital", "site", "spin", "all"],
        default="all",
        help="Type of projection (default: all)",
    )
    parser.add_argument(
        "--projections-mode",
        choices=["expectation", "weight_up", "weight_down", "polarization"],
        default="expectation",
        help=(
            "Spin projection visualization mode (for --projections --projection spin): "
            "expectation: <S_z> in [-0.5,+0.5]; "
            "weight_up: P_up in [0,1]; "
            "weight_down: P_down in [0,1]; "
            "polarization: signed in [-1,+1] (default: expectation)"
        ),
    )
    parser.add_argument(
        "--output-prefix", default="bands", help="Output file prefix (default: bands)"
    )
    parser.add_argument(
        "--energy-range",
        nargs=2,
        type=float,
        help="Y-axis energy range (emin emax) for band plots",
    )
    # Altermagnetic analysis options
    parser.add_argument(
        "--altermagnetic",
        action="store_true",
        help="Analyze altermagnetic materials (spin splitting without net magnetization)",
    )
    parser.add_argument(
        "--mag-threshold",
        type=float,
        default=0.1,
        help="Magnetization threshold for altermagnetic detection (default: 0.1)",
    )
    parser.add_argument(
        "--plot-splitting",
        action="store_true",
        help="Plot spin splitting map in k-space (with --altermagnetic)",
    )
    parser.add_argument(
        "--plot-character",
        action="store_true",
        help="Plot altermagnetic character, magnetization vs splitting (with --altermagnetic)",
    )

    # Berry curvature options
    parser.add_argument(
        "--berry",
        action="store_true",
        help="Calculate Berry curvature and topological properties",
    )
    parser.add_argument(
        "--fermi-energy",
        type=float,
        default=0.0,
        help="Fermi energy in eV (default: 0.0)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.0,
        help="Temperature in K for thermal effects (default: 0.0 K)",
    )
    parser.add_argument(
        "--chern-bands",
        nargs="+",
        type=int,
        help="Band indices for Chern number calculation (space-separated)",
    )
    parser.add_argument(
        "--chern-tol",
        type=float,
        default=0.1,
        help="Tolerance for Chern number quantization (default: 0.1)",
    )
    parser.add_argument(
        "--ahc",
        action="store_true",
        help="Calculate anomalous Hall conductivity for first 5 bands (with --berry)",
    )
    parser.add_argument(
        "--orbital-mag",
        action="store_true",
        help="Calculate orbital magnetization for first 5 bands (with --berry)",
    )
    parser.add_argument(
        "--num-bands", type=int, help="Number of bands to analyze (default: all)"
    )
    parser.add_argument(
        "--plot-berry",
        type=int,
        help="Plot Berry curvature for specified band (with --berry)",
    )
    parser.add_argument(
        "--plot-chern",
        action="store_true",
        help="Plot Chern numbers histogram (with --berry)",
    )
    parser.add_argument(
        "--plot-bands",
        action="store_true",
        help="Plot band structure (with --altermagnetic or --berry)",
    )

    # Advanced K-path options (seekpath, custom paths, prime convention)
    parser.add_argument(
        "--kpath-mode",
        choices=["auto", "seekpath", "custom", "interactive"],
        default="auto",
        help=(
            "K-path generation mode: "
            "auto: Simple path from pymatgen (fast); "
            "seekpath: High-symmetry path with explicit points; "
            "custom: User-specified path (requires --custom-path); "
            "interactive: Prompt user to choose (default: auto)"
        ),
    )
    parser.add_argument(
        "--seekpath-symprec",
        type=float,
        default=1e-5,
        help="Symmetry tolerance for seekpath (default: 1e-5)",
    )
    parser.add_argument(
        "--seekpath-distance",
        type=float,
        default=0.025,
        help="Reference distance for seekpath BZ sampling (default: 0.025 Å⁻¹)",
    )
    parser.add_argument(
        "--custom-path",
        type=str,
        help="Custom k-path labels, e.g. 'K-Gamma-K' or 'X-Gamma-X-M' (requires --kpath-mode custom)",
    )
    parser.add_argument(
        "--prime-points",
        action="store_true",
        help="Apply prime convention to k-points (K'=(-K_x,-K_y,K_z); useful for custom paths)",
    )
    parser.add_argument(
        "--stabilize-degenerate",
        choices=["none", "block_average", "pair_cancel"],
        default="block_average",
        help="Stabilization mode for degenerate bands in projections (default: block_average)",
    )
    parser.add_argument(
        "--plot-bz",
        action="store_true",
        help="Plot Brillouin zone polyhedron with high-symmetry points (requires pymatgen)",
    )
    parser.add_argument(
        "--export-kpath",
        type=str,
        help="Export k-path to file for reproducibility (e.g., 'kpath.dat')",
    )

    # Spin texture options (Phase 4)
    parser.add_argument(
        "--spin-texture",
        action="store_true",
        help="Calculate and visualize spin texture",
    )
    parser.add_argument(
        "--spin-band",
        type=int,
        help="Band index for spin texture analysis (with --spin-texture)",
    )
    parser.add_argument(
        "--spin-plane",
        choices=["xy", "xz", "yz"],
        default="xy",
        help="k-space plane for 2D spin texture plot (default: xy)",
    )
    parser.add_argument(
        "--spin-component",
        choices=["x", "y", "z", "magnitude"],
        default="z",
        help="Spin component for coloring (default: z)",
    )
    parser.add_argument(
        "--rashba-analysis",
        action="store_true",
        help="Extract Rashba/Dresselhaus parameters (with --spin-texture)",
    )
    parser.add_argument(
        "--plot-spin-fermi",
        action="store_true",
        help="Plot 3D Fermi surface with spin texture (with --spin-texture)",
    )

    # Fat bands options (Phase 4)
    parser.add_argument(
        "--fat-bands",
        action="store_true",
        help="Plot fat bands showing orbital/spin projections",
    )
    parser.add_argument(
        "--fat-type",
        choices=["linewidth", "color", "spin"],
        default="spin",
        help="Fat bands visualization type (default: spin)",
    )
    parser.add_argument(
        "--fat-orbitals",
        nargs="+",
        type=int,
        help="Orbital indices for fat bands projection",
    )
    parser.add_argument(
        "--fat-orbital-groups",
        type=str,
        help='Orbital groups as JSON string, e.g., \'{"s":[0,9],"p":[1,2,3,10,11,12]}\'',
    )
    parser.add_argument(
        "--fat-scale",
        type=float,
        default=5.0,
        help="Scale factor for fat bands linewidth (default: 5.0)",
    )

    parser.set_defaults(func=run_bands_subcommand)
    return parser


def create_dos_parser(subparsers):
    """Create DOS subcommand parser"""
    parser = subparsers.add_parser(
        "dos",
        parents=[_make_common_parser()],
        help="Calculate density of states or DOS moments",
    )
    parser.add_argument(
        "--mesh",
        nargs=3,
        type=int,
        default=[6, 6, 6],
        help="k-point mesh size (default: 6 6 6)",
    )
    parser.add_argument(
        "--energy-range", nargs=2, type=float, help="Energy range (emin emax)"
    )
    parser.add_argument(
        "--num-energies",
        type=int,
        default=1000,
        help="Number of energy points (default: 1000)",
    )
    parser.add_argument(
        "--fermi-energy",
        type=float,
        default=0.0,
        help="Fermi energy for plots (default: 0.0)",
    )
    parser.add_argument(
        "--electrons",
        type=float,
        help="Number of electrons (auto-determine Fermi level)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperature in K (default: 300.0)",
    )
    parser.add_argument(
        "--no-symmetry", action="store_true", help="Disable symmetry reduction"
    )
    parser.add_argument(
        "--no-tetrahedron",
        action="store_true",
        help="Use Gaussian broadening instead of tetrahedron",
    )
    parser.add_argument(
        "--gaussian-smearing",
        type=float,
        default=0.1,
        help="Gaussian smearing width in eV (default: 0.1)",
    )
    parser.add_argument(
        "--output-prefix", default="dos", help="Output file prefix (default: dos)"
    )
    parser.add_argument(
        "--moments",
        action="store_true",
        help="Calculate DOS moments analysis instead of DOS plot",
    )
    parser.set_defaults(func=run_dos_subcommand)
    return parser


def validate_setup_inputs():
    """
    Validate that required input files exist for --setup workflow.
    
    Returns
    -------
    bool
        True if all required files exist, False otherwise
    """
    required_files = ['lattice.nml', 'input.nml']
    missing_files = []
    
    for filename in required_files:
        if not Path(filename).exists():
            missing_files.append(filename)
    
    if missing_files:
        print(f"ERROR: Required input files not found: {', '.join(missing_files)}", file=sys.stderr)
        print("\nThe --setup workflow requires the following files in the current directory:", file=sys.stderr)
        print("  - lattice.nml   : Lattice structure data (lattice_t namelist)", file=sys.stderr)
        print("  - input.nml     : Atomic positions and parameters (&lattice, &atoms namelists)", file=sys.stderr)
        print("\nPlease ensure these files exist before running --setup.", file=sys.stderr)
        return False
    
    return True


def validate_setup_parameters(min_radius, max_radius, lmx, lattice_dist, basis_dist, alat_scale,
                              mode_2d, z_iso_factor, rec_steps):
    """
    Validate setup parameter ranges.
    
    Parameters
    ----------
    min_radius : float
        Minimum empty sphere radius
    max_radius : float
        Maximum empty sphere radius
    
    Returns
    -------
    bool
        True if parameters are valid, False otherwise
    """
    if min_radius <= 0:
        print(f"ERROR: --min_empty_radius must be positive, got {min_radius}", file=sys.stderr)
        return False
    
    if max_radius <= 0:
        print(f"ERROR: --max_empty_radius must be positive, got {max_radius}", file=sys.stderr)
        return False
    
    if min_radius >= max_radius:
        print(f"ERROR: --min_empty_radius ({min_radius}) must be less than --max_empty_radius ({max_radius})", file=sys.stderr)
        return False
    
    if max_radius > 10.0:
        print(f"WARNING: --max_empty_radius ({max_radius}) is very large. Typical values are 1.5-3.0.", file=sys.stderr)
    
    if min_radius < 0.5:
        print(f"WARNING: --min_empty_radius ({min_radius}) is very small. Typical values are 1.0-1.5.", file=sys.stderr)

    if lmx < 0 or lmx > 3:
        print(f"ERROR: --lmx must be in [0,3], got {lmx}", file=sys.stderr)
        return False

    if lattice_dist < 0:
        print(f"ERROR: --lattice-dist must be >= 0, got {lattice_dist}", file=sys.stderr)
        return False

    if basis_dist < 0:
        print(f"ERROR: --basis-dist must be >= 0, got {basis_dist}", file=sys.stderr)
        return False

    if alat_scale <= 0:
        print(f"ERROR: --alat-scale must be > 0, got {alat_scale}", file=sys.stderr)
        return False

    if z_iso_factor is not None and z_iso_factor <= 0:
        print(f"ERROR: --z-iso-factor must be > 0, got {z_iso_factor}", file=sys.stderr)
        return False

    if z_iso_factor is not None and not mode_2d:
        print("ERROR: --z-iso-factor requires --2d", file=sys.stderr)
        return False

    if rec_steps <= 0:
        print(f"ERROR: --rec-steps must be > 0, got {rec_steps}", file=sys.stderr)
        return False
    
    return True


def validate_sqs_parameters(sqs, sqs_shells, sqs_steps):
    if sqs is None:
        return True
    try:
        parse_sqs_spec(sqs)
    except Exception as exc:
        print(f"ERROR: invalid --sqs specification: {exc}", file=sys.stderr)
        return False
    if sqs_shells < 1:
        print(f"ERROR: --sqs-shells must be >= 1, got {sqs_shells}", file=sys.stderr)
        return False
    if sqs_steps < 0:
        print(f"ERROR: --sqs-steps must be >= 0, got {sqs_steps}", file=sys.stderr)
        return False
    return True


def validate_output_files():
    """
    Validate that output files were created and are readable.
    
    Returns
    -------
    dict
        Dictionary with validation results for each file
    """
    validation_results = {}
    expected_files = ['lattice_new.nml', 'input_new.nml']
    
    for filename in expected_files:
        filepath = Path(filename)
        result = {
            'exists': filepath.exists(),
            'readable': False,
            'size': 0,
            'valid_format': False
        }
        
        if result['exists']:
            result['size'] = filepath.stat().st_size
            try:
                with open(filepath, 'r') as f:
                    content = f.read()
                    result['readable'] = True
                    # Basic format check - should contain namelist markers
                    result['valid_format'] = '&' in content and '/' in content
            except Exception:
                pass
        
        validation_results[filename] = result
    
    return validation_results


def print_validation_report(validation_results):
    """Print validation report for output files."""
    print("\nValidation Report:")
    print("-" * 70)
    
    all_valid = True
    for filename, result in validation_results.items():
        status = "✓ PASS" if result['exists'] and result['readable'] and result['valid_format'] else "✗ FAIL"
        print(f"{status} {filename}")
        
        if not result['exists']:
            print(f"     ERROR: File not created")
            all_valid = False
        elif not result['readable']:
            print(f"     ERROR: File not readable")
            all_valid = False
        elif not result['valid_format']:
            print(f"     WARNING: File may not be valid namelist format")
        else:
            print(f"     Size: {result['size']} bytes")
    
    print("-" * 70)
    if all_valid:
        print("All validation checks passed ✓")
    else:
        print("Some validation checks failed ✗")
    
    return all_valid


def create_backup_files():
    """Create backup copies of input files."""
    input_files = ['lattice.nml', 'input.nml']
    backed_up = []
    
    for filename in input_files:
        if Path(filename).exists():
            backup_path = Path(f"{filename}.bak")
            try:
                import shutil
                shutil.copy2(filename, backup_path)
                backed_up.append(str(backup_path))
            except Exception as e:
                print(f"WARNING: Could not create backup of {filename}: {e}", file=sys.stderr)
    
    if backed_up:
        print(f"Created backups: {', '.join(backed_up)}")
    
    return backed_up


def apply_in_place_updates():
    """Apply in-place updates by renaming _new files to originals."""
    updates = [
        ('lattice_new.nml', 'lattice.nml'),
        ('input_new.nml', 'input.nml')
    ]
    
    updated = []
    for new_file, orig_file in updates:
        if Path(new_file).exists():
            # Create backup first
            if Path(orig_file).exists():
                backup_path = Path(f"{orig_file}.bak")
                try:
                    import shutil
                    shutil.copy2(orig_file, backup_path)
                    print(f"Backup created: {backup_path}")
                except Exception as e:
                    print(f"WARNING: Could not create backup of {orig_file}: {e}", file=sys.stderr)
            
            # Replace original with new
            try:
                import shutil
                shutil.move(new_file, orig_file)
                updated.append(orig_file)
                print(f"Updated: {orig_file}")
            except Exception as e:
                print(f"ERROR: Could not update {orig_file}: {e}", file=sys.stderr)
    
    return updated


def parse_sqs_spec(sqs_spec: str):
    """
    Parse SQS composition spec like "Fe=0.5,Ni=0.3,Co=0.2".

    Returns
    -------
    tuple[list[str], np.ndarray]
        Species labels and normalized concentrations.
    """
    parts = [p.strip() for p in sqs_spec.split(",") if p.strip()]
    if not parts:
        raise ValueError("Empty --sqs specification")

    species = []
    conc = []
    for part in parts:
        if "=" not in part:
            raise ValueError(f"Invalid --sqs token '{part}', expected Label=fraction")
        label, value = [x.strip() for x in part.split("=", 1)]
        if not label:
            raise ValueError(f"Invalid --sqs token '{part}', missing label")
        frac = float(value)
        if frac < 0:
            raise ValueError(f"Invalid --sqs fraction for {label}: {frac} (must be >= 0)")
        species.append(label)
        conc.append(frac)

    conc = np.array(conc, dtype=float)
    s = conc.sum()
    if s <= 0:
        raise ValueError("Sum of --sqs fractions must be > 0")
    conc /= s
    return species, conc


def _allocate_species_counts(concentrations: np.ndarray, n_sites: int):
    """Allocate integer species counts matching target concentrations."""
    raw = concentrations * n_sites
    counts = np.floor(raw).astype(int)
    remainder = n_sites - counts.sum()
    if remainder > 0:
        order = np.argsort(raw - counts)[::-1]
        for i in range(remainder):
            counts[order[i % len(order)]] += 1
    elif remainder < 0:
        order = np.argsort(raw - counts)
        for i in range(-remainder):
            idx = order[i % len(order)]
            if counts[idx] > 0:
                counts[idx] -= 1
    return counts


def _build_shell_pairs(positions: np.ndarray, n_shells: int, distance_tol: float = 1e-5):
    """Build pair lists grouped by increasing distance shells (i<j)."""
    n = positions.shape[0]
    pairs = []
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(positions[i] - positions[j]))
            pairs.append((i, j))
            dists.append(d)

    if not pairs:
        return []

    order = np.argsort(dists)
    shells = []
    current = [pairs[order[0]]]
    dref = dists[order[0]]
    for idx in order[1:]:
        d = dists[idx]
        if abs(d - dref) <= distance_tol:
            current.append(pairs[idx])
        else:
            shells.append(current)
            if len(shells) >= n_shells:
                break
            current = [pairs[idx]]
            dref = d
    if len(shells) < n_shells and current:
        shells.append(current)
    return shells[:n_shells]


def _sqs_objective(assignment: np.ndarray, n_species: int, shells, concentrations: np.ndarray):
    """Simple pair-correlation mismatch objective for unordered pair probabilities."""
    target = np.zeros((n_species, n_species), dtype=float)
    for a in range(n_species):
        for b in range(a, n_species):
            if a == b:
                target[a, b] = concentrations[a] * concentrations[a]
            else:
                target[a, b] = 2.0 * concentrations[a] * concentrations[b]

    total = 0.0
    for shell_pairs in shells:
        obs = np.zeros((n_species, n_species), dtype=float)
        for i, j in shell_pairs:
            a = int(assignment[i])
            b = int(assignment[j])
            if a <= b:
                obs[a, b] += 1.0
            else:
                obs[b, a] += 1.0
        obs /= max(len(shell_pairs), 1)
        diff = obs - target
        total += float(np.sum(diff * diff))
    return total


def build_sqs_assignment(positions: np.ndarray, species: list[str], concentrations: np.ndarray,
                         n_shells: int = 3, n_steps: int = 5000, seed: int | None = None):
    """
    Build a simple SQS-like assignment by MC swaps minimizing pair-shell mismatch.
    """
    n_sites = positions.shape[0]
    n_species = len(species)
    counts = _allocate_species_counts(concentrations, n_sites)
    if counts.sum() != n_sites:
        raise RuntimeError("Internal error: SQS counts do not sum to n_sites")

    rng = np.random.default_rng(seed)
    assignment = np.empty(n_sites, dtype=int)
    start = 0
    for ispec, c in enumerate(counts):
        assignment[start:start + c] = ispec
        start += c
    rng.shuffle(assignment)

    shells = _build_shell_pairs(positions, n_shells=n_shells)
    if not shells:
        raise ValueError("Could not construct SQS pair shells from structure")

    best = assignment.copy()
    best_obj = _sqs_objective(best, n_species, shells, concentrations)
    curr = best.copy()
    curr_obj = best_obj

    for _ in range(max(n_steps, 0)):
        i, j = rng.integers(0, n_sites, size=2)
        if i == j or curr[i] == curr[j]:
            continue
        curr[i], curr[j] = curr[j], curr[i]
        new_obj = _sqs_objective(curr, n_species, shells, concentrations)
        if new_obj <= curr_obj:
            curr_obj = new_obj
            if new_obj < best_obj:
                best_obj = new_obj
                best = curr.copy()
        else:
            curr[i], curr[j] = curr[j], curr[i]

    return best, counts, best_obj


def create_supercell(
    nx, ny, nz,
    lattice_path: str = "lattice.nml",
    input_path: str = "input.nml",
    output_lattice_path: str = "lattice.nml",
    output_input_path: str = "input.nml",
    in_place: bool = True,
    sqs_spec: str | None = None,
    sqs_shells: int = 3,
    sqs_steps: int = 5000,
    sqs_seed: int | None = None,
):
    """
    Create a supercell by repeating the unit cell nx×ny×nz times.
    
    Reads lattice.nml, creates supercell, and writes updated file.
    Coordinates in lattice.nml are Cartesian (normalized, lattice parameter = 1).
    The supercell uses broken symmetry indexing (1..N) for index arrays, expands
    input.nml labels into unique per-atom species, updates ntype/ct accordingly,
    and creates per-species .nml files by copying from the original species files.
    
    Parameters
    ----------
    nx, ny, nz : int
        Repetitions along each lattice vector
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    import re
    import shutil
    import numpy as np
    
    print(f"\nCreating {nx}×{ny}×{nz} supercell...")
    
    lattice_file = Path(lattice_path)
    if not lattice_file.exists():
        print("ERROR: lattice.nml not found", file=sys.stderr)
        return False
    
    def set_or_insert_lattice_scalar(text, key, value):
        """Set key in &lattice block, or insert it if missing."""
        block = re.search(r'&lattice\b(.*?)/', text, flags=re.IGNORECASE | re.DOTALL)
        if not block:
            return text
        btxt = block.group(0)
        # Replace existing assignment in lattice block only
        new_btxt, nrep = re.subn(
            rf'(^\s*{key}\s*=\s*).*$',
            rf'\g<1>{value}',
            btxt,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        if nrep == 0:
            # Insert before terminating "/"
            new_btxt = re.sub(r'/\s*$', f"{key} = {value}\n/", btxt, flags=re.MULTILINE)
        return text[:block.start()] + new_btxt + text[block.end():]

    try:
        # === Process lattice.nml ===
        with open(lattice_file, 'r') as f:
            content = f.read()
        
        # Extract lattice vectors (a(:,1), a(:,2), a(:,3))
        lattice_vectors = []
        for i in [1, 2, 3]:
            pattern = rf'a\s*\(\s*:\s*,\s*{i}\s*\)\s*=\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)'
            match = re.search(pattern, content)
            if match:
                vec = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                lattice_vectors.append(vec)
        
        if len(lattice_vectors) != 3:
            print(f"ERROR: Could not parse all 3 lattice vectors (found {len(lattice_vectors)})", file=sys.stderr)
            return False
        
        a_vec, b_vec, c_vec = np.array(lattice_vectors)
        
        # Extract atomic positions (Cartesian coordinates, normalized)
        positions = []
        crd_pattern = r'crd\s*\(\s*:\s*,\s*(\d+)\s*\)\s*=\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)'
        for match in re.finditer(crd_pattern, content):
            atom_idx = int(match.group(1))
            x, y, z = float(match.group(2)), float(match.group(3)), float(match.group(4))
            positions.append(np.array([x, y, z]))
        
        if not positions:
            print("ERROR: Could not parse atomic positions from lattice.nml", file=sys.stderr)
            return False
        
        n_atoms_orig = len(positions)
        print(f"  Original unit cell: {n_atoms_orig} atoms")
        
        # Create supercell positions (Cartesian)
        # For each repetition (ix, iy, iz), shift by ix*a + iy*b + iz*c
        supercell_positions = []
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for pos in positions:
                        new_pos = pos + ix * a_vec + iy * b_vec + iz * c_vec
                        supercell_positions.append(new_pos)
        
        n_atoms_super = len(supercell_positions)
        print(f"  Supercell: {n_atoms_super} atoms ({nx}×{ny}×{nz} = {nx*ny*nz} unit cells)")
        
        # Update nbas and related counts
        content = re.sub(r'(nbulk_bulk\s*=\s*)\d+', f'nbulk_bulk = {n_atoms_super}', content, flags=re.IGNORECASE)
        content = re.sub(r'(ntot\s*=\s*)\d+', f'ntot = {n_atoms_super}', content, flags=re.IGNORECASE)
        content = re.sub(r'(nbas\s*=\s*)\d+', f'nbas = {n_atoms_super}', content, flags=re.IGNORECASE)
        content = re.sub(r'(nrec\s*=\s*)\d+', f'nrec = {n_atoms_super}', content, flags=re.IGNORECASE)
        
        # Scale lattice vectors
        content = re.sub(
            r'(a\s*\(\s*:\s*,\s*1\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
            lambda m: f"{m.group(1)}{float(m.group(2))*nx:.10f}, {float(m.group(3))*nx:.10f}, {float(m.group(4))*nx:.10f}",
            content,
            flags=re.IGNORECASE
        )
        
        content = re.sub(
            r'(a\s*\(\s*:\s*,\s*2\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
            lambda m: f"{m.group(1)}{float(m.group(2))*ny:.10f}, {float(m.group(3))*ny:.10f}, {float(m.group(4))*ny:.10f}",
            content,
            flags=re.IGNORECASE
        )
        
        content = re.sub(
            r'(a\s*\(\s*:\s*,\s*3\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
            lambda m: f"{m.group(1)}{float(m.group(2))*nz:.10f}, {float(m.group(3))*nz:.10f}, {float(m.group(4))*nz:.10f}",
            content,
            flags=re.IGNORECASE
        )
        
        # Replace atomic positions
        crd_start = re.search(r'crd\s*\(\s*:\s*,\s*1\s*\)', content)
        if not crd_start:
            print("ERROR: Could not find coordinate section", file=sys.stderr)
            return False
        
        # Find end of crd section
        crd_section_start = crd_start.start()
        lines = content[crd_section_start:].split('\n')
        crd_lines_count = sum(1 for line in lines if 'crd(' in line)
        
        temp_pos = crd_section_start
        for i in range(crd_lines_count):
            temp_pos = content.find('\n', temp_pos) + 1
        crd_section_end = temp_pos
        
        # Build new coordinate section
        new_crd_section = ""
        for i, pos in enumerate(supercell_positions, start=1):
            new_crd_section += f"crd(:, {i}) = {pos[0]:.10f}, {pos[1]:.10f}, {pos[2]:.10f}\n"
        
        new_content = content[:crd_section_start] + new_crd_section + content[crd_section_end:]
        
        # Update arrays that depend on atom count (broken symmetry: 1..N)
        for array_name in ['izp', 'no', 'iu', 'ib', 'irec']:
            array_start = re.search(rf'{array_name}\s*\(\s*1\s*\)', new_content)
            if array_start:
                array_section_start = array_start.start()
                lines = new_content[array_section_start:].split('\n')
                array_lines = sum(1 for line in lines if f'{array_name}(' in line)
                
                temp_pos = array_section_start
                for i in range(array_lines):
                    temp_pos = new_content.find('\n', temp_pos) + 1
                array_section_end = temp_pos
                
                new_array_section = ""
                for i in range(1, n_atoms_super + 1):
                    new_array_section += f"{array_name}({i}) = {i}\n"
                
                new_content = new_content[:array_section_start] + new_array_section + new_content[array_section_end:]
        
        output_lattice_file = Path(output_lattice_path)
        if in_place:
            backup_file = Path('lattice_original.nml.bak')
            shutil.copy2(lattice_file, backup_file)
            print(f"  ✓ Backup lattice file: {backup_file}")

        with open(output_lattice_file, 'w') as f:
            f.write(new_content)

        if in_place:
            print(f"  ✓ Updated {output_lattice_file} with supercell")
        else:
            print(f"  ✓ Wrote staged supercell lattice to {output_lattice_file}")
        
        # === Update input.nml to replicate labels ===
        input_file = Path(input_path)
        if input_file.exists():
            with open(input_file, 'r') as f:
                input_content = f.read()
            
            # Extract labels
            label_pattern = r"label\s*\(\s*(\d+)\s*\)\s*=\s*['\"]([^'\"]+)['\"]"
            orig_labels = {}
            for match in re.finditer(label_pattern, input_content):
                idx = int(match.group(1))
                label = match.group(2)
                if idx <= n_atoms_orig:
                    orig_labels[idx] = label
            
            # Extract ct values (match original label ordering)
            ct_pattern = r'ct\s*\(\s*(\d+)\s*\)\s*=\s*([\-\d.eE+]+)'
            orig_ct = {}
            for match in re.finditer(ct_pattern, input_content):
                idx = int(match.group(1))
                val = match.group(2)
                orig_ct[idx] = val
            
            if orig_labels:
                # Replicate labels for supercell (unique labels, broken symmetry)
                orig_label_keys = sorted(orig_labels.keys())
                label_meta = []
                for key in orig_label_keys:
                    label = orig_labels[key]
                    match = re.match(r'^(.*?)(\d+)$', label)
                    if match:
                        base = match.group(1)
                        start_num = int(match.group(2))
                    else:
                        base = label
                        start_num = 0
                    label_meta.append((base, start_num, label))
                
                orig_ct_values = []
                if orig_ct:
                    for key in orig_label_keys:
                        if key in orig_ct:
                            orig_ct_values.append(orig_ct[key])
                
                label_counts = {}
                super_labels = []
                super_ct = []
                source_labels = []

                if sqs_spec:
                    sqs_species, sqs_conc = parse_sqs_spec(sqs_spec)
                    assign, counts, sqs_obj = build_sqs_assignment(
                        positions=np.array(supercell_positions, dtype=float),
                        species=sqs_species,
                        concentrations=sqs_conc,
                        n_shells=sqs_shells,
                        n_steps=sqs_steps,
                        seed=sqs_seed,
                    )
                    print(f"  SQS enabled: {sqs_spec}")
                    print(f"  SQS shells={sqs_shells}, steps={sqs_steps}, seed={sqs_seed}")
                    print(f"  SQS objective={sqs_obj:.6e}")

                    for isite, ispec in enumerate(assign):
                        base = sqs_species[int(ispec)]
                        if base not in label_counts:
                            label_counts[base] = 0
                        label_counts[base] += 1
                        super_labels.append(f"{base}{label_counts[base]}")
                        source_labels.append(base)
                        if orig_ct_values:
                            super_ct.append(orig_ct_values[0])
                else:
                    for _ in range(nx * ny * nz):
                        for idx, key in enumerate(orig_label_keys):
                            base, start_num, src_label = label_meta[idx]
                            if base not in label_counts:
                                label_counts[base] = start_num - 1 if start_num > 0 else 0
                            label_counts[base] += 1
                            super_labels.append(f"{base}{label_counts[base]}")
                            source_labels.append(src_label)
                            if orig_ct_values:
                                super_ct.append(orig_ct_values[idx])
                new_ntype = len(super_labels)
                
                # Find label section and replace
                label_start = re.search(r"label\s*\(\s*1\s*\)", input_content)
                if label_start:
                    # Find end of label section
                    label_section_start = label_start.start()
                    lines = input_content[label_section_start:].split('\n')
                    label_lines = sum(1 for line in lines if 'label(' in line)
                    
                    temp_pos = label_section_start
                    for i in range(label_lines):
                        temp_pos = input_content.find('\n', temp_pos) + 1
                    label_section_end = temp_pos
                    
                    # Build new label section
                    new_label_section = ""
                    for i, label in enumerate(super_labels, start=1):
                        new_label_section += f"label({i}) = '{label}'\n"
                    
                    input_content = input_content[:label_section_start] + new_label_section + input_content[label_section_end:]
                    
                    # Ensure ntype in &lattice matches expanded labels
                    input_content = set_or_insert_lattice_scalar(input_content, "ntype", str(new_ntype))
                    
                    # Expand ct definitions to match ntype.
                    # Supports both ct(:)=... and ct(i)=... input styles.
                    if super_ct:
                        ct_start = re.search(r'ct\s*\(\s*(:|\d+)\s*\)', input_content, flags=re.IGNORECASE)
                        if ct_start:
                            ct_section_start = ct_start.start()
                            lines = input_content[ct_section_start:].split('\n')
                            ct_lines = sum(1 for line in lines if re.search(r'\bct\s*\(', line, flags=re.IGNORECASE))
                            
                            temp_pos = ct_section_start
                            for i in range(ct_lines):
                                temp_pos = input_content.find('\n', temp_pos) + 1
                            ct_section_end = temp_pos
                            
                            new_ct_section = ""
                            for i, val in enumerate(super_ct, start=1):
                                new_ct_section += f"ct({i}) = {val}\n"
                            
                            input_content = input_content[:ct_section_start] + new_ct_section + input_content[ct_section_end:]
                        else:
                            # No ct definition existed; insert per-type ct lines in &lattice
                            insert_ct = "".join([f"ct({i}) = {val}\n" for i, val in enumerate(super_ct, start=1)])
                            block = re.search(r'&lattice\b(.*?)/', input_content, flags=re.IGNORECASE | re.DOTALL)
                            if block:
                                btxt = block.group(0)
                                new_btxt = re.sub(r'/\s*$', insert_ct + "/", btxt, flags=re.MULTILINE)
                                input_content = input_content[:block.start()] + new_btxt + input_content[block.end():]
                    
                    output_input_file = Path(output_input_path)
                    if in_place:
                        input_backup = Path('input_original.nml.bak')
                        shutil.copy2(input_file, input_backup)
                        print(f"  ✓ Backup input file: {input_backup}")

                    with open(output_input_file, 'w') as f:
                        f.write(input_content)

                    if in_place:
                        print(f"  ✓ Updated {output_input_file} with replicated labels")
                    else:
                        print(f"  ✓ Wrote staged supercell input to {output_input_file}")
                    
                    # Create per-species .nml files for unique labels only for in-place mode
                    if in_place:
                        for new_label, src_label in zip(super_labels, source_labels):
                            src_file = Path(f"{src_label}.nml")
                            dest_file = Path(f"{new_label}.nml")
                            if src_file.exists() and not dest_file.exists():
                                shutil.copy2(src_file, dest_file)
        
        print(f"  ✓ Supercell creation complete")
        
        return True
        
    except Exception as e:
        print(f"ERROR creating supercell: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False


def print_summary_diff():
    """Print a summary of changes between old and new files."""
    print("\nSummary of Changes:")
    print("=" * 70)

    old_lattice = Path('lattice.nml')
    new_lattice = Path('lattice_new.nml')
    if old_lattice.exists() and new_lattice.exists():
        print(f"\nlattice.nml:")
        print(f"  Old size: {old_lattice.stat().st_size} bytes")
        print(f"  New size: {new_lattice.stat().st_size} bytes")
        print(f"  Change:   {new_lattice.stat().st_size - old_lattice.stat().st_size:+d} bytes")
    old_input = Path('input.nml')
    new_input = Path('input_new.nml')
    if old_input.exists() and new_input.exists():
        print(f"\ninput.nml:")
        print(f"  Old size: {old_input.stat().st_size} bytes")
        print(f"  New size: {new_input.stat().st_size} bytes")
        print(f"  Change:   {new_input.stat().st_size - old_input.stat().st_size:+d} bytes")
    print("=" * 70)


def _run_setup_cli(args):
    """Translate legacy ``-S`` arguments, present results, and call the service."""
    print("=" * 70)
    print("RS-utils Lattice Setup")
    print("=" * 70)
    rec_steps = args.rec_steps
    if rec_steps is None:
        rec_steps = 200 if args.rec_method == "chebyshev" else 20
    try:
        spec = SetupSpec(
            cluster_type=args.cluster,
            empty_spheres=args.empty_spheres,
            min_empty_radius=args.min_empty_radius,
            max_empty_radius=args.max_empty_radius,
            lmx=args.lmx,
            rec_method=args.rec_method,
            rec_steps=rec_steps,
            lattice_dist=args.lattice_dist,
            basis_dist=args.basis_dist,
            alat_scale=args.alat_scale,
            seed=args.seed,
            mode_2d=args.mode_2d,
            z_iso_factor=args.z_iso_factor,
            supercell=tuple(args.supercell or (1, 1, 1)),
            sqs=args.sqs,
            sqs_shells=args.sqs_shells,
            sqs_steps=args.sqs_steps,
            sqs_seed=args.sqs_seed,
        )
        service_validate_setup_parameters(spec)
        service_validate_setup_inputs(Path.cwd())
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if args.json_report:
            try:
                failed_result = SetupResult(
                    status=OperationStatus.FAILED,
                    diagnostics=(
                        Diagnostic(
                            DiagnosticSeverity.ERROR,
                            ErrorCode.MISSING_INPUT if isinstance(exc, FileNotFoundError) else ErrorCode.SETUP_SPEC_INVALID,
                            str(exc),
                        ),
                    ),
                )
                write_setup_report(failed_result, args.json_report)
            except OSError as report_exc:
                print(f"ERROR: could not write JSON setup report: {report_exc}", file=sys.stderr)
        return 1

    print("\nConfiguration:")
    print(f"  Cluster type:      {spec.cluster_type}")
    print(f"  Empty spheres:     {'enabled' if spec.empty_spheres else 'disabled'}")
    if spec.empty_spheres:
        print(f"  Min radius:        {spec.min_empty_radius}")
        print(f"  Max radius:        {spec.max_empty_radius}")
    print(f"  LMAX override:     {spec.lmx}")
    print(f"  Rec method:        {spec.rec_method}")
    print(f"  Rec steps:         {spec.rec_steps}")
    print(f"  Lattice dist:      {spec.lattice_dist} alat")
    print(f"  Basis dist:        {spec.basis_dist} alat")
    print(f"  Alat scale:        {spec.alat_scale}")
    print(f"  2D mode:           {'enabled' if spec.mode_2d else 'disabled'}")
    if spec.mode_2d:
        rule = "(zmax-zmin) + 1.0*alat" if spec.z_iso_factor is None else f"isotropic (factor={spec.z_iso_factor})"
        print(f"  2D z rule:         {rule}")
    if spec.seed is not None:
        print(f"  Seed:              {spec.seed}")
    if spec.sqs is not None:
        print(f"  SQS spec:          {spec.sqs}")
        print(f"  SQS shells:        {spec.sqs_shells}")
        print(f"  SQS steps:         {spec.sqs_steps}")
        if spec.sqs_seed is not None:
            print(f"  SQS seed:          {spec.sqs_seed}")
    if args.supercell:
        print(f"  Supercell:         {'×'.join(str(value) for value in spec.supercell)}")
    for label, enabled in (("Verbose mode", args.verbose), ("Validation", args.validate), ("Summary", args.summary), ("Backup", args.backup), ("In-place update", args.in_place)):
        if enabled:
            print(f"  {label}:           enabled")
    print()
    if args.backup or args.in_place:
        print("Creating backup files...\n")
    result = SetupService().run(
        Path.cwd(),
        Path.cwd(),
        spec,
        in_place=args.in_place,
        backup=args.backup or args.in_place,
        # Preserve the established interactive behavior for owned cwd output
        # names while the library API defaults to collision-safe ``error``.
        collision_policy="replace",
    )
    if args.json_report:
        try:
            write_setup_report(result, args.json_report)
        except OSError as exc:
            print(f"ERROR: could not write JSON setup report: {exc}", file=sys.stderr)
            return 1
    if not result.success:
        for diagnostic in result.diagnostics:
            print(f"{diagnostic.severity.value.upper()}: {diagnostic.message}", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="" if result.stderr.endswith("\n") else "\n")
        return 1
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if args.verbose:
        print("Starting lattice processing...")
        print("Lattice processing completed.")
    print("\n" + "=" * 70)
    print("SUCCESS: Lattice processing completed")
    print("=" * 70)
    print("\nGenerated files:")
    for filename in result.candidate_files:
        path = Path(filename)
        if path.is_file():
            print(f"  ✓ {filename} ({path.stat().st_size} bytes)")
    if args.validate:
        print_validation_report(validate_output_files())
    if args.summary:
        print_summary_diff()
    if args.in_place:
        print("\nApplying in-place updates...")
        print("\nUpdated 2 file(s) in place.")
        print("Original files backed up with .bak extension.")
    else:
        print("\nNext steps:")
        print("  1. Review the generated files")
        if not args.validate:
            print("     (run with --validate to check file integrity)")
        if not args.summary:
            print("     (run with --summary to see detailed changes)")
        print("  2. Rename them to replace original files:")
        print("     mv lattice_new.nml lattice.nml")
        print("     mv input_new.nml input.nml")
        print("     (or use --in-place flag to do this automatically)")
        print("  3. Continue with exchange calculations:")
        print("     rsutil --exchange --cutoff 2.5")
    return 0
    
    # Compare lattice files
    old_lattice = Path('lattice.nml')
    new_lattice = Path('lattice_new.nml')
    
    if old_lattice.exists() and new_lattice.exists():
        old_size = old_lattice.stat().st_size
        new_size = new_lattice.stat().st_size
        print(f"\nlattice.nml:")
        print(f"  Old size: {old_size} bytes")
        print(f"  New size: {new_size} bytes")
        print(f"  Change:   {new_size - old_size:+d} bytes")
    
    # Compare input files  
    old_input = Path('input.nml')
    new_input = Path('input_new.nml')
    
    if old_input.exists() and new_input.exists():
        old_size = old_input.stat().st_size
        new_size = new_input.stat().st_size
        print(f"\ninput.nml:")
        print(f"  Old size: {old_size} bytes")
        print(f"  New size: {new_size} bytes")
        print(f"  Change:   {new_size - old_size:+d} bytes")
        
        # Try to extract some key info
        try:
            with open(new_input, 'r') as f:
                content = f.read()
                if 'nbas' in content:
                    # Try to find nbas value
                    for line in content.split('\n'):
                        if 'nbas' in line.lower() and '=' in line:
                            print(f"  {line.strip()}")
                            break
        except Exception:
            pass
    
    print("=" * 70)


def main():
    """
    Command-line interface for rsutil.

    This function sets up the command-line argument parser and handles the
    execution of different functionalities based on the provided arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(
        description="RS Tools Command Line Interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Existing functionality (flags)
  rsutil --setup                            # Process lattice
  rsutil --lattice --data examples/bccFe    # Visualize lattice
  rsutil --exchange --cutoff 1.5            # Create exchange maps
  rsutil --export --jij jijso.out           # Export to UppASD
  
  # Band structure & DOS (subcommands)
  rsutil bands --lattice lattice.nml --nk 200            # Band structure
  rsutil bands --projections --projection orbital        # Orbital projections
  rsutil dos --mesh 8 8 8 --electrons 8.0                # DOS with Fermi level
  rsutil dos --mesh 10 10 10 --electrons 8.0 --moments   # DOS moments analysis
  
  # Band-Related Properties (Phase 2 - flags on bands subcommand)
  rsutil bands --altermagnetic --nk 100 --plot-splitting --plot-bands
                                      # Detect altermagnetic materials
  rsutil bands --berry --nk 150 --ahc --orbital-mag --plot-chern
                                      # Topological properties analysis
  rsutil bands --berry --nk 150 --chern-bands 0 1 2 --plot-berry 0
                                      # Calculate specific Chern numbers
        """,
    )

    # Create subparsers for new band/DOS functionality
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # Add band/DOS subcommands
    create_bands_parser(subparsers)
    create_dos_parser(subparsers)

    # Existing flags for backward compatibility
    parser.add_argument(
        "--setup", "-S", action="store_true", default=False, help="Process lattice"
    )
    parser.add_argument(
        "--lattice",
        "-L",
        action="store_true",
        default=False,
        help="Run visualizer for lattice input",
    )
    parser.add_argument(
        "--clust",
        "-C",
        action="store_true",
        default=False,
        help="Run visualizer for clust file",
    )
    parser.add_argument(
        "--exchange",
        "-J",
        action="store_true",
        default=False,
        help="Create exchange input files",
    )
    parser.add_argument(
        "--data",
        type=str,
        default=".",
        help="Directory for data files to be visualized",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default="1.0",
        help="Cutoff for exchange interactions, in lattice units",
    )
    parser.add_argument(
        "--unique",
        action="store_true",
        default=False,
        help="Find only unique neighbors per shell (one representative per unique species-vector combination)",
    )
    parser.add_argument(
        "--shell-width",
        type=float,
        default=0.1,
        help="Width of shells for grouping neighbors when using --unique (default: 0.1)",
    )
    parser.add_argument(
        "--soc",
        action="store_true",
        default=False,
        help="Account for spin-orbit coupling in symmetry considerations",
    )
    parser.add_argument(
        "--rij-out",
        type=str,
        default=None,
        help="Write the legacy R_ij neighbor mapping to this file",
    )
    parser.add_argument(
        "--rij-format",
        choices=["json", "txt"],
        default="json",
        help="Format for --rij-out (default: json)",
    )
    parser.add_argument(
        "--export",
        "-E",
        action="store_true",
        default=False,
        help="Export RS-LMTO inputs to UppASD format",
    )
    parser.add_argument(
        "--jij",
        type=str,
        default=None,
        help="Path to jijso.out file to include in export",
    )
    parser.add_argument(
        "--dmi",
        type=str,
        default=None,
        help="Path to dijso.out / dmi file to include in export (DM vectors)",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="uppasd_out",
        help="Output directory for exported UppASD files (relative to --data)",
    )
    parser.add_argument(
        "--lmto47",
        action="store_true",
        default=False,
        help="Export RS-LMTO lattice/input to LMTO47 CTRL format",
    )
    parser.add_argument(
        "--nspin",
        type=int,
        default=2,
        help="NSPIN value for LMTO47 export (default: 2)",
    )
    parser.add_argument(
        "--lmx",
        type=int,
        default=2,
        help="Maximum L channel override (default: 2; spd)",
    )
    parser.add_argument(
        "--rec-method",
        type=str,
        choices=["block", "chebyshev", "lanczos"],
        default="block",
        help="Recursion method for &control recur in setup output (default: block)",
    )
    parser.add_argument(
        "--rec-steps",
        type=int,
        default=None,
        help="Recursion steps for &control lld in setup output (default: 20, or 200 for chebyshev)",
    )
    parser.add_argument(
        "--lattice-dist",
        type=float,
        default=0.0,
        help="Std.dev of Gaussian lattice-vector distortion in units of alat (default: 0.0)",
    )
    parser.add_argument(
        "--basis-dist",
        type=float,
        default=0.0,
        help="Std.dev of Gaussian basis-vector distortion in units of alat (default: 0.0)",
    )
    parser.add_argument(
        "--alat-scale",
        type=float,
        default=1.0,
        help="Scale factor applied to alat before setup/export transforms (default: 1.0)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible lattice/basis distortions (default: random)",
    )
    parser.add_argument(
        "--2d",
        dest="mode_2d",
        action="store_true",
        default=False,
        help="Enable 2D setup mode and derive z-axis length before WS/ASA processing",
    )
    parser.add_argument(
        "--z-iso-factor",
        type=float,
        default=None,
        help="Isotropic-density surface margin factor for 2D z-axis derivation (requires --2d)",
    )
    parser.add_argument(
        "--sqs",
        type=str,
        default=None,
        help='Enable SQS assignment on all sites, e.g. --sqs "Fe=0.5,Ni=0.5"',
    )
    parser.add_argument(
        "--sqs-shells",
        type=int,
        default=3,
        help="Number of pair distance shells in SQS objective (default: 3)",
    )
    parser.add_argument(
        "--sqs-steps",
        type=int,
        default=5000,
        help="Number of MC swap steps for SQS optimization (default: 5000)",
    )
    parser.add_argument(
        "--sqs-seed",
        type=int,
        default=None,
        help="Random seed for SQS optimizer (default: random)",
    )
    parser.add_argument(
        "--cluster",
        type=str,
        default="periodic",
        choices=["periodic", "free"],
        help="Cluster type: periodic (default) or free",
    )
    parser.add_argument(
        "--empty_spheres",
        action="store_true",
        default=False,
        help="Enable empty sphere calculations",
    )
    parser.add_argument(
        "--min_empty_radius",
        type=float,
        default=1.2,
        help="Minimum empty sphere radius (default: 1.2)",
    )
    parser.add_argument(
        "--max_empty_radius",
        type=float,
        default=2.5,
        help="Maximum empty sphere radius (default: 2.5)",
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true', default=False,
        help='Show detailed progress during processing'
    )
    parser.add_argument(
        '--summary', action='store_true', default=False,
        help='Print summary of changes after processing'
    )
    parser.add_argument(
        '--validate', action='store_true', default=False,
        help='Validate output files after generation'
    )
    parser.add_argument(
        '--json-report', type=str, default=None, metavar='PATH',
        help='Write a versioned machine-readable setup report to PATH'
    )
    parser.add_argument(
        '--backup', action='store_true', default=False,
        help='Create backup files before any operations'
    )
    parser.add_argument(
        '--in-place', action='store_true', default=False,
        help='Update input files in place (creates .bak backups)'
    )
    parser.add_argument(
        '--supercell', nargs=3, type=int, metavar=('NX', 'NY', 'NZ'),
        help='Create supercell by repeating unit cell NX×NY×NZ times'
    )

    args = parser.parse_args()

    # Handle subcommands first (new band/DOS functionality)
    if args.subcommand:
        if hasattr(args, "func"):
            return args.func(args)
        else:
            parser.print_help()
            return

    # Handle existing flag-based functionality (backward compatibility)
    # Process lattice if --setup flag is provided
    if args.setup:
        return _run_setup_cli(args)
        # Legacy implementation retained below only as historical context;
        # all setup execution now goes through SetupService.
        print("=" * 70)
        print("RS-utils Lattice Setup")
        print("=" * 70)
        
        rec_steps = args.rec_steps
        if rec_steps is None:
            rec_steps = 200 if args.rec_method == "chebyshev" else 20

        # Validate input files exist
        if not validate_setup_inputs():
            sys.exit(1)
        
        # Validate parameter ranges
        if not validate_setup_parameters(
            args.min_empty_radius,
            args.max_empty_radius,
            args.lmx,
            args.lattice_dist,
            args.basis_dist,
            args.alat_scale,
            args.mode_2d,
            args.z_iso_factor,
            rec_steps,
        ):
            sys.exit(1)
        if not validate_sqs_parameters(args.sqs, args.sqs_shells, args.sqs_steps):
            sys.exit(1)
        if args.sqs is not None and args.supercell is None:
            print("ERROR: --sqs currently requires --supercell", file=sys.stderr)
            sys.exit(1)
        
        # Display configuration
        print(f"\nConfiguration:")
        print(f"  Cluster type:      {args.cluster}")
        print(f"  Empty spheres:     {'enabled' if args.empty_spheres else 'disabled'}")
        if args.empty_spheres:
            print(f"  Min radius:        {args.min_empty_radius}")
            print(f"  Max radius:        {args.max_empty_radius}")
        print(f"  LMAX override:     {args.lmx}")
        print(f"  Rec method:        {args.rec_method}")
        print(f"  Rec steps:         {rec_steps}")
        print(f"  Lattice dist:      {args.lattice_dist} alat")
        print(f"  Basis dist:        {args.basis_dist} alat")
        print(f"  Alat scale:        {args.alat_scale}")
        print(f"  2D mode:           {'enabled' if args.mode_2d else 'disabled'}")
        if args.mode_2d:
            if args.z_iso_factor is None:
                print("  2D z rule:         (zmax-zmin) + 1.0*alat")
            else:
                print(f"  2D z rule:         isotropic (factor={args.z_iso_factor})")
        if args.seed is not None:
            print(f"  Seed:              {args.seed}")
        if args.sqs is not None:
            print(f"  SQS spec:          {args.sqs}")
            print(f"  SQS shells:        {args.sqs_shells}")
            print(f"  SQS steps:         {args.sqs_steps}")
            if args.sqs_seed is not None:
                print(f"  SQS seed:          {args.sqs_seed}")
        if args.supercell:
            nx, ny, nz = args.supercell
            print(f"  Supercell:         {nx}×{ny}×{nz}")
        if args.verbose:
            print(f"  Verbose mode:      enabled")
        if args.validate:
            print(f"  Validation:        enabled")
        if args.summary:
            print(f"  Summary:           enabled")
        if args.backup:
            print(f"  Backup:            enabled")
        if args.in_place:
            print(f"  In-place update:   enabled")
        print()
        
        # Create backups if requested
        if args.backup or args.in_place:
            print("Creating backup files...")
            create_backup_files()
            print()
        
        # Create supercell if requested (do this BEFORE processing)
        staged_supercell = False
        staged_lattice = ".rsutil_supercell_lattice.nml"
        staged_input = ".rsutil_supercell_input.nml"
        if args.supercell:
            nx, ny, nz = args.supercell
            if nx < 1 or ny < 1 or nz < 1:
                print("ERROR: Supercell dimensions must be positive integers", file=sys.stderr)
                sys.exit(1)

            if args.in_place:
                ok_sc = create_supercell(
                    nx, ny, nz, in_place=True,
                    sqs_spec=args.sqs, sqs_shells=args.sqs_shells,
                    sqs_steps=args.sqs_steps, sqs_seed=args.sqs_seed,
                )
            else:
                ok_sc = create_supercell(
                    nx, ny, nz,
                    lattice_path="lattice.nml",
                    input_path="input.nml",
                    output_lattice_path=staged_lattice,
                    output_input_path=staged_input,
                    in_place=False,
                    sqs_spec=args.sqs, sqs_shells=args.sqs_shells,
                    sqs_steps=args.sqs_steps, sqs_seed=args.sqs_seed,
                )
                staged_supercell = ok_sc

            if not ok_sc:
                print("ERROR: Failed to create supercell", file=sys.stderr)
                sys.exit(1)
            print()
        
        # Process lattice with error handling
        try:
            if args.verbose:
                print("Starting lattice processing...")
            
            if staged_supercell:
                import shutil
                orig_lat_tmp = ".rsutil_original_lattice.nml"
                orig_inp_tmp = ".rsutil_original_input.nml"
                shutil.move("lattice.nml", orig_lat_tmp)
                shutil.move("input.nml", orig_inp_tmp)
                shutil.copy2(staged_lattice, "lattice.nml")
                shutil.copy2(staged_input, "input.nml")
                try:
                    process_lattice(
                        cluster_type=args.cluster,
                        empty_spheres=args.empty_spheres,
                        min_empty_radius=args.min_empty_radius,
                        max_empty_radius=args.max_empty_radius,
                        lmx=args.lmx,
                        lattice_dist=args.lattice_dist,
                        basis_dist=args.basis_dist,
                        alat_scale=args.alat_scale,
                        seed=args.seed,
                        mode_2d=args.mode_2d,
                        z_iso_factor=args.z_iso_factor,
                        rec_steps=rec_steps,
                        rec_method=args.rec_method,
                    )
                finally:
                    shutil.move(orig_lat_tmp, "lattice.nml")
                    shutil.move(orig_inp_tmp, "input.nml")
                    if Path(staged_lattice).exists():
                        Path(staged_lattice).unlink()
                    if Path(staged_input).exists():
                        Path(staged_input).unlink()
            else:
                process_lattice(
                    cluster_type=args.cluster,
                    empty_spheres=args.empty_spheres,
                    min_empty_radius=args.min_empty_radius,
                    max_empty_radius=args.max_empty_radius,
                    lmx=args.lmx,
                    lattice_dist=args.lattice_dist,
                    basis_dist=args.basis_dist,
                    alat_scale=args.alat_scale,
                    seed=args.seed,
                    mode_2d=args.mode_2d,
                    z_iso_factor=args.z_iso_factor,
                    rec_steps=rec_steps,
                    rec_method=args.rec_method,
                )
            
            if args.verbose:
                print("Lattice processing completed.")
            
            # Check if output files were created
            output_files = ['lattice_new.nml', 'input_new.nml']
            created_files = [f for f in output_files if Path(f).exists()]
            
            if created_files:
                print("\n" + "=" * 70)
                print("SUCCESS: Lattice processing completed")
                print("=" * 70)
                print("\nGenerated files:")
                for filename in created_files:
                    file_size = Path(filename).stat().st_size
                    print(f"  ✓ {filename} ({file_size} bytes)")
                
                # Run validation if requested
                if args.validate:
                    validation_results = validate_output_files()
                    print_validation_report(validation_results)
                
                # Print summary if requested
                if args.summary:
                    print_summary_diff()
                
                # Apply in-place updates if requested
                if args.in_place:
                    print("\nApplying in-place updates...")
                    updated = apply_in_place_updates()
                    if updated:
                        print(f"\nUpdated {len(updated)} file(s) in place.")
                        print("Original files backed up with .bak extension.")
                        print("\nNext steps:")
                        print("  Continue with exchange calculations:")
                        print("     rsutil --exchange --cutoff 2.5")
                    else:
                        print("No files were updated.", file=sys.stderr)
                else:
                    # Standard workflow message
                    print("\nNext steps:")
                    print("  1. Review the generated files")
                    if not args.validate:
                        print("     (run with --validate to check file integrity)")
                    if not args.summary:
                        print("     (run with --summary to see detailed changes)")
                    print("  2. Rename them to replace original files:")
                    print("     mv lattice_new.nml lattice.nml")
                    print("     mv input_new.nml input.nml")
                    print("     (or use --in-place flag to do this automatically)")
                    print("  3. Continue with exchange calculations:")
                    print("     rsutil --exchange --cutoff 2.5")
            else:
                print("\nWARNING: Processing completed but expected output files were not found.", file=sys.stderr)
                print("Check for error messages above.", file=sys.stderr)
                
        except FileNotFoundError as e:
            print(f"\nERROR: File not found during processing: {e}", file=sys.stderr)
            print("This may indicate missing data files or incorrect file paths.", file=sys.stderr)
            sys.exit(1)
        except ValueError as e:
            print(f"\nERROR: Invalid parameter or data: {e}", file=sys.stderr)
            sys.exit(1)
        except RuntimeError as e:
            print(f"\nERROR: Processing failed: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"\nERROR: Unexpected error during lattice processing:", file=sys.stderr)
            print(f"  {type(e).__name__}: {e}", file=sys.stderr)
            print("\nPlease check your input files and try again.", file=sys.stderr)
            print("If the problem persists, please report this error.", file=sys.stderr)
            sys.exit(1)
    # Run visualizer for lattice or clust file if --lattice or --clust flag is provided
    elif args.lattice or args.clust:
        rs_visualizer(args.data, clust=args.clust)
    # Create exchange input files if --exchange flag is provided
    elif args.exchange:
        exchange_map(
            cutoff=args.cutoff,
            unique=args.unique,
            shell_width=args.shell_width,
            soc=args.soc,
            rij_out=args.rij_out,
            rij_fmt=args.rij_format,
        )
    # Export to LMTO47 CTRL if requested
    elif args.lmto47:
        if args.lmx < 0 or args.lmx > 3:
            print(f"ERROR: --lmx must be in [0,3], got {args.lmx}", file=sys.stderr)
            sys.exit(1)
        if args.lattice_dist < 0 or args.basis_dist < 0:
            print("ERROR: --lattice-dist and --basis-dist must be >= 0", file=sys.stderr)
            sys.exit(1)
        if args.alat_scale <= 0:
            print(f"ERROR: --alat-scale must be > 0, got {args.alat_scale}", file=sys.stderr)
            sys.exit(1)
        if args.z_iso_factor is not None and args.z_iso_factor <= 0:
            print(f"ERROR: --z-iso-factor must be > 0, got {args.z_iso_factor}", file=sys.stderr)
            sys.exit(1)
        if args.z_iso_factor is not None and not args.mode_2d:
            print("ERROR: --z-iso-factor requires --2d", file=sys.stderr)
        if not validate_sqs_parameters(args.sqs, args.sqs_shells, args.sqs_steps):
            sys.exit(1)
        if args.sqs is not None and args.supercell is None:
            print("ERROR: --sqs currently requires --supercell", file=sys.stderr)
            sys.exit(1)
        lattice_path = os.path.join(args.data, "lattice.nml")
        input_path = os.path.join(args.data, "input.nml")
        out_dir = os.path.join(args.data, args.out) if args.data else args.out

        staged_supercell = False
        staged_lattice = os.path.join(args.data, ".rsutil_supercell_lattice.nml")
        staged_input = os.path.join(args.data, ".rsutil_supercell_input.nml")

        if args.supercell:
            nx, ny, nz = args.supercell
            if nx < 1 or ny < 1 or nz < 1:
                print("ERROR: Supercell dimensions must be positive integers", file=sys.stderr)
                sys.exit(1)

            if args.in_place:
                ok_sc = create_supercell(
                    nx, ny, nz,
                    lattice_path=lattice_path,
                    input_path=input_path,
                    output_lattice_path=lattice_path,
                    output_input_path=input_path,
                    in_place=True,
                    sqs_spec=args.sqs, sqs_shells=args.sqs_shells,
                    sqs_steps=args.sqs_steps, sqs_seed=args.sqs_seed,
                )
            else:
                ok_sc = create_supercell(
                    nx, ny, nz,
                    lattice_path=lattice_path,
                    input_path=input_path,
                    output_lattice_path=staged_lattice,
                    output_input_path=staged_input,
                    in_place=False,
                    sqs_spec=args.sqs, sqs_shells=args.sqs_shells,
                    sqs_steps=args.sqs_steps, sqs_seed=args.sqs_seed,
                )
                staged_supercell = ok_sc

            if not ok_sc:
                print("ERROR: Failed to create supercell for LMTO47 export", file=sys.stderr)
                sys.exit(1)

        export_lattice = staged_lattice if staged_supercell else lattice_path
        export_input = staged_input if staged_supercell else input_path
        try:
            ctrl_file = export_ctrl(
                source=("rs_lmto", export_lattice, export_input),
                out_dir=out_dir,
                header="RS-tools export",
                nspin=args.nspin,
                lmx=args.lmx,
                lattice_dist=args.lattice_dist,
                basis_dist=args.basis_dist,
                alat_scale=args.alat_scale,
                seed=args.seed,
                mode_2d=args.mode_2d,
                z_iso_factor=args.z_iso_factor,
            )
            print(f"LMTO47 CTRL written to: {ctrl_file}")
        finally:
            if staged_supercell:
                if Path(staged_lattice).exists():
                    Path(staged_lattice).unlink()
                if Path(staged_input).exists():
                    Path(staged_input).unlink()
    # Export to UppASD if requested
    elif args.export:
        # default lattice path in data dir
        lattice_path = os.path.join(args.data, "lattice.nml")
        out_dir = os.path.join(args.data, args.out) if args.data else args.out
        export_from_files(
            lattice_nml=lattice_path,
            input_nml=None,
            clust=None,
            out_dir=out_dir,
            jij_file=args.jij,
            dij_file=args.dmi,
        )
    # Print help message if no valid flag is provided
    else:
        parser.print_help()


if __name__ == "__main__":
    raise SystemExit(main())
