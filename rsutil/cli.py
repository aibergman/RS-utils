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
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer
from .exchange.exchange_map import exchange_map
from .export.uppasd import export_from_files


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
  rsutil bands --altermagnetic --nk 150 --  rsutil bands --altermagnetic --nk 150 --  rsutil bands    rsutil bands --altermae altermagneti  rsutil bands --til ban  rsutil bands --alt--ahc --orbital-mag --plot-chern
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
        process_lattice(
            cluster_type=args.cluster,
            empty_spheres=args.empty_spheres,
            min_empty_radius=args.min_empty_radius,
            max_empty_radius=args.max_empty_radius,
        )
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
        )
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
    main()
