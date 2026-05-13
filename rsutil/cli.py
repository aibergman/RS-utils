"""
This module provides a command-line interface (CLI) for the rsutil package.
It allows users to process lattice structures, visualize data, and create
exchange input files through various command-line arguments.

Functions
---------
    Sets up the command-line argument parser and handles the execution of
    different functionalities based on the provided arguments.

Usage
-----
To use this CLI, run the script with the desired flags:
    python cli.py --setup
    python cli.py --lattice --data <data_directory>
    python cli.py --clust --data <data_directory>
    python cli.py --exchange

Flags
-----
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
--lmto47, -T
    Export to LMTO47 CTRL file.
--data
    Directory for data files to be visualized (default: current directory).
"""

import argparse
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer
from .exchange.exchange_map import exchange_map
from .export.uppasd import export_from_files
from .export.lmto47 import export_ctrl
import os


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
        description="RS Tools Command Line Interface"
    )
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
        '--rij-out', type=str, default=None,
        help='Write R_ij mappings to this file (JSON)'
    )
    parser.add_argument(
        '--rij-format', type=str, default='json', choices=['json', 'txt'],
        help='Format for R_ij output (json or txt; default: json)'
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
        '--lmto47', '-T', action='store_true', default=False,
        help='Export to LMTO47 CTRL file'
    )
    parser.add_argument(
        '--nspin', type=int, default=2,
        help='Number of spin channels for LMTO47 (default: 2 for spin-polarized)'
    )
    parser.add_argument(
        '--lmx', type=int, default=None,
        help='Maximum L quantum number for all atoms in LMTO47 (default: None = use LMTO47 default)'
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

    # Process lattice if --setup flag is provided
    if args.setup:
        process_lattice(cluster_type=args.cluster, empty_spheres=args.empty_spheres, 
                       min_empty_radius=args.min_empty_radius, max_empty_radius=args.max_empty_radius)
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
    # Export to UppASD if requested
    elif args.export:
        # default lattice path in data dir
        lattice_path = os.path.join(args.data, 'lattice.nml')
        out_dir = os.path.join(args.data, args.out) if args.data else args.out
        export_from_files(lattice_nml=lattice_path, input_nml=None, clust=None,
                          out_dir=out_dir, jij_file=args.jij, dij_file=args.dmi)
    # Export to LMTO47 CTRL if requested
    elif args.lmto47:
        lattice_path = os.path.join(args.data, 'lattice.nml')
        input_path = os.path.join(args.data, 'input.nml')
        out_dir = os.path.join(args.data, args.out) if args.data else args.out
        ctrl_path = export_ctrl(("rs_lmto", lattice_path, input_path), out_dir=out_dir,
                               nspin=args.nspin, lmx=args.lmx)
        print(f"CTRL file written to: {ctrl_path}")
    # Print help message if no valid flag is provided
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
