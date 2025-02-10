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
--data
    Directory for data files to be visualized (default: current directory).
"""

import argparse
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer
from .exchange.exchange_map import exchange_map


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

    args = parser.parse_args()

    # Process lattice if --setup flag is provided
    if args.setup:
        process_lattice()
    # Run visualizer for lattice or clust file if --lattice or --clust flag is provided
    elif args.lattice or args.clust:
        rs_visualizer(args.data, clust=args.clust)
    # Create exchange input files if --exchange flag is provided
    elif args.exchange:
        exchange_map(cutoff=args.cutoff)
    # Print help message if no valid flag is provided
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
