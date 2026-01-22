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
import sys
import os
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer
from .exchange.exchange_map import exchange_map
from .export.uppasd import export_from_files


def _make_common_parser():
    """Create common parser for band/DOS subcommands"""
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument('--lattice', default='lattice.nml', help='Lattice file (default: lattice.nml)')
    p.add_argument('--hamiltonian', default='rs2txt.dat', help='Hamiltonian file (default: rs2txt.dat)')
    p.add_argument('--lattice-format', choices=['rs-lmto', 'paoflow', 'kite'], 
                   default=None, help='Lattice file format (auto-detect if not specified)')
    p.add_argument('--hamiltonian-format', choices=['rs-lmto', 'paoflow', 'kite'],
                   default=None, help='Hamiltonian file format (auto-detect if not specified)')
    return p


def run_bands_subcommand(args):
    """Run band structure calculation or projections"""
    from .bands import load_lattice, load_hamiltonian, BandStructureCalculator
    
    print(f"Reading lattice from {args.lattice}...")
    lattice = load_lattice(args.lattice, format=args.lattice_format)
    
    print(f"Reading Hamiltonian from {args.hamiltonian}...")
    hamiltonian = load_hamiltonian(args.hamiltonian, format=args.hamiltonian_format)
    
    if args.projections:
        # Run projected bands
        from .bands.projections import ProjectedBandCalculator
        
        print("Setting up projected band structure calculator...")
        # Create calculator with specified spin projection mode
        proj_calc = ProjectedBandCalculator(
            lattice, hamiltonian,
            spin_projection_mode=getattr(args, 'projections_mode', 'expectation')
        )
        
        print(f"Generating k-path with {args.nk} points...")
        kpoints, klabels, kdistances = proj_calc.setup_kpath(args.nk)
        
        print("Calculating projected band structure...")
        band_energies, projections = proj_calc.calculate_projected_bands(kpoints)
        
        projection_type = args.projection
        output_file = f"{args.output_prefix}_{projection_type}.png"
        print(f"Plotting {projection_type} projections (mode: {proj_calc.spin_projection_mode}) to {output_file}...")
        proj_calc.plot_projected_bands(kdistances, band_energies, projections, 
                                       klabels, projection_type, output_file)
        
        print("Done!")
    else:
        # Run regular band structure
        print("Setting up band structure calculator...")
        calc = BandStructureCalculator(lattice, hamiltonian)
        
        print(f"Generating k-path with {args.nk} points...")
        kpoints, klabels, kdistances = calc.setup_kpath(args.nk)
        
        print("Calculating band structure...")
        band_energies = calc.calculate_bands(kpoints)
        
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
        use_tetrahedron = not getattr(args, 'no_tetrahedron', False)
        
        dos_data = dos_calc.calculate_dos(
            mesh_size=tuple(args.mesh),
            energy_range=tuple(args.energy_range) if args.energy_range else None,
            num_energies=args.num_energies,
            use_symmetry=not args.no_symmetry,
            use_tetrahedron=use_tetrahedron,
            gaussian_smearing=0.1
        )
        
        # Perform moments analysis
        moments_calc = DOSMomentsCalculator(
            dos_data['energies'], 
            dos_data['total_dos'],
            dos_data.get('orbital_dos'),
            temperature=args.temperature
        )
        
        if args.electrons:
            print(f"Finding Fermi level for {args.electrons} electrons...")
            energy_range = tuple(args.energy_range) if args.energy_range else None
            fermi_level = moments_calc.find_fermi_level(args.electrons, energy_range)
            print(f"Fermi level: {fermi_level:.6f} eV")
        else:
            print("Auto-detecting electron count...")
            # Count electrons at zero energy
            total_electrons = moments_calc.total_dos.sum() * (moments_calc.energies[1] - moments_calc.energies[0]) 
            fermi_level = moments_calc.find_fermi_level(total_electrons)
            print(f"Fermi level: {fermi_level:.6f} eV")
        
        # Calculate moments
        print("\nCalculating DOS moments...")
        moments = moments_calc.calculate_moments(fermi_level)
        
        # Print results
        print("\nDOS Moments:")
        if 'total' in moments:
            print(f"\nTotal DOS:")
            total_moments = moments['total']
            print(f"  Occupation (electrons):      {total_moments['occupation']:.4f}")
            print(f"  Avg Energy (rel to EF):      {total_moments['avg_energy_rel_ef']:.4f} eV")
            print(f"  Energy Spread (width):       {total_moments['band_width']:.4f} eV")
            print(f"  Band Center:                 {total_moments['band_center']:.4f} eV")
        
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
            print(f"Calculating DOS with Fermi level determination for {args.electrons} electrons...")
            energy_range_val = tuple(args.energy_range) if args.energy_range else None
            dos_data, moments_calc = dos_calc.calculate_dos_with_moments(
                mesh_size=tuple(args.mesh),
                target_electrons=args.electrons,
                temperature=getattr(args, 'temperature', 300.0),
                energy_range=energy_range_val,
                num_energies=args.num_energies,
                use_symmetry=not args.no_symmetry,
                use_tetrahedron=not args.no_tetrahedron,
                gaussian_smearing=getattr(args, 'gaussian_smearing', 0.1),
                output_prefix=args.output_prefix
            )
            
            fermi_level = moments_calc.fermi_level
            if fermi_level is not None:
                print(f"Fermi level determined: {fermi_level:.6f} eV")
                dos_calc.plot_dos(dos_data, fermi_energy=fermi_level, 
                                output_file=f"{args.output_prefix}.png")
            else:
                print("Warning: Could not determine Fermi level")
                dos_calc.plot_dos(dos_data, fermi_energy=0.0,
                                output_file=f"{args.output_prefix}.png")
        else:
            # Simple DOS calculation
            print("Calculating DOS...")
            dos_data = dos_calc.calculate_dos(
                mesh_size=tuple(args.mesh),
                energy_range=tuple(args.energy_range) if args.energy_range else None,
                num_energies=args.num_energies,
                use_symmetry=not args.no_symmetry,
                use_tetrahedron=not args.no_tetrahedron,
                gaussian_smearing=getattr(args, 'gaussian_smearing', 0.1)
            )
            
            fermi_energy = getattr(args, 'fermi_energy', 0.0)
            dos_calc.plot_dos(dos_data, fermi_energy=fermi_energy,
                            output_file=f"{args.output_prefix}.png")
        
        print("Done!")


def create_bands_parser(subparsers):
    """Create bands subcommand parser"""
    parser = subparsers.add_parser('bands', parents=[_make_common_parser()],
                                   help='Calculate band structure and projections')
    parser.add_argument('--nk', type=int, default=100,
                       help='Number of k-points along path (default: 100)')
    parser.add_argument('--output', default='band_structure.png',
                       help='Output plot file (default: band_structure.png)')
    parser.add_argument('--projections', action='store_true',
                       help='Calculate projected bands instead of total bands')
    parser.add_argument('--projection', choices=['orbital', 'site', 'spin', 'all'],
                       default='all', help='Type of projection (default: all)')
    parser.add_argument('--projections-mode', 
                       choices=['expectation', 'weight_up', 'weight_down', 'polarization'],
                       default='expectation',
                       help=('Spin projection visualization mode (for --projections --projection spin): '
                            'expectation: <S_z> in [-0.5,+0.5]; '
                            'weight_up: P_up in [0,1]; '
                            'weight_down: P_down in [0,1]; '
                            'polarization: signed in [-1,+1] (default: expectation)'))
    parser.add_argument('--output-prefix', default='projected',
                       help='Output file prefix for projections (default: projected)')
    parser.set_defaults(func=run_bands_subcommand)
    return parser


def create_dos_parser(subparsers):
    """Create DOS subcommand parser"""
    parser = subparsers.add_parser('dos', parents=[_make_common_parser()],
                                   help='Calculate density of states or DOS moments')
    parser.add_argument('--mesh', nargs=3, type=int, default=[6, 6, 6],
                       help='k-point mesh size (default: 6 6 6)')
    parser.add_argument('--energy-range', nargs=2, type=float,
                       help='Energy range (emin emax)')
    parser.add_argument('--num-energies', type=int, default=1000,
                       help='Number of energy points (default: 1000)')
    parser.add_argument('--fermi-energy', type=float, default=0.0,
                       help='Fermi energy for plots (default: 0.0)')
    parser.add_argument('--electrons', type=float,
                       help='Number of electrons (auto-determine Fermi level)')
    parser.add_argument('--temperature', type=float, default=300.0,
                       help='Temperature in K (default: 300.0)')
    parser.add_argument('--no-symmetry', action='store_true',
                       help='Disable symmetry reduction')
    parser.add_argument('--no-tetrahedron', action='store_true',
                       help='Use Gaussian broadening instead of tetrahedron')
    parser.add_argument('--gaussian-smearing', type=float, default=0.1,
                       help='Gaussian smearing width in eV (default: 0.1)')
    parser.add_argument('--output-prefix', default='dos',
                       help='Output file prefix (default: dos)')
    parser.add_argument('--moments', action='store_true',
                       help='Calculate DOS moments analysis instead of DOS plot')
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
  
  # New band/DOS functionality (subcommands)
  rsutil bands --lattice lattice.nml --nk 200            # Band structure
  rsutil bands --projections --projection orbital        # Orbital projections
  rsutil dos --mesh 8 8 8 --electrons 8.0                # DOS with Fermi level
  rsutil dos --mesh 10 10 10 --electrons 8.0 --moments   # DOS moments analysis
        """
    )
    
    # Create subparsers for new band/DOS functionality
    subparsers = parser.add_subparsers(dest='subcommand', help='Available subcommands')
    
    # Add band/DOS subcommands
    create_bands_parser(subparsers)
    create_dos_parser(subparsers)
    
    # Existing flags for backward compatibility
    parser.add_argument(
        '--setup', '-S', action='store_true', default=False,
        help='Process lattice'
    )
    parser.add_argument(
        '--lattice', '-L', action='store_true', default=False,
        help='Run visualizer for lattice input'
    )
    parser.add_argument(
        '--clust', '-C', action='store_true', default=False,
        help='Run visualizer for clust file'
    )
    parser.add_argument(
        '--exchange', '-J', action='store_true', default=False,
        help='Create exchange input files'
    )
    parser.add_argument(
        '--data', type=str, default='.',
        help='Directory for data files to be visualized'
    )
    parser.add_argument(
        '--cutoff', type=float, default='1.0',
        help='Cutoff for exchange interactions, in lattice units'
    )
    parser.add_argument(
        '--unique', action='store_true', default=False,
        help='Find only unique neighbors per shell (one representative per unique species-vector combination)'
    )
    parser.add_argument(
        '--shell-width', type=float, default=0.1,
        help='Width of shells for grouping neighbors when using --unique (default: 0.1)'
    )
    parser.add_argument(
        '--soc', action='store_true', default=False,
        help='Account for spin-orbit coupling in symmetry considerations'
    )
    parser.add_argument(
        '--export', '-E', action='store_true', default=False,
        help='Export RS-LMTO inputs to UppASD format'
    )
    parser.add_argument(
        '--jij', type=str, default=None,
        help='Path to jijso.out file to include in export'
    )
    parser.add_argument(
        '--dmi', type=str, default=None,
        help='Path to dijso.out / dmi file to include in export (DM vectors)'
    )
    parser.add_argument(
        '--out', type=str, default='uppasd_out',
        help='Output directory for exported UppASD files (relative to --data)'
    )
    parser.add_argument(
        '--cluster', type=str, default='periodic', choices=['periodic', 'free'],
        help='Cluster type: periodic (default) or free'
    )
    parser.add_argument(
        '--empty_spheres', action='store_true', default=False,
        help='Enable empty sphere calculations'
    )
    parser.add_argument(
        '--min_empty_radius', type=float, default=1.2,
        help='Minimum empty sphere radius (default: 1.2)'
    )
    parser.add_argument(
        '--max_empty_radius', type=float, default=2.5,
        help='Maximum empty sphere radius (default: 2.5)'
    )

    args = parser.parse_args()

    # Handle subcommands first (new band/DOS functionality)
    if args.subcommand:
        if hasattr(args, 'func'):
            return args.func(args)
        else:
            parser.print_help()
            return

    # Handle existing flag-based functionality (backward compatibility)
    # Process lattice if --setup flag is provided
    if args.setup:
        process_lattice(cluster_type=args.cluster, empty_spheres=args.empty_spheres, 
                       min_empty_radius=args.min_empty_radius, max_empty_radius=args.max_empty_radius)
    # Run visualizer for lattice or clust file if --lattice or --clust flag is provided
    elif args.lattice or args.clust:
        rs_visualizer(args.data, clust=args.clust)
    # Create exchange input files if --exchange flag is provided
    elif args.exchange:
        exchange_map(cutoff=args.cutoff, unique=args.unique, shell_width=args.shell_width, soc=args.soc)
    # Export to UppASD if requested
    elif args.export:
        # default lattice path in data dir
        lattice_path = os.path.join(args.data, 'lattice.nml')
        out_dir = os.path.join(args.data, args.out) if args.data else args.out
        export_from_files(lattice_nml=lattice_path, input_nml=None, clust=None,
                          out_dir=out_dir, jij_file=args.jij, dij_file=args.dmi)
    # Print help message if no valid flag is provided
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
