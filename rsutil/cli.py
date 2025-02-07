import argparse
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer


def main():
    """Command-line interface for rsutil."""
    parser = argparse.ArgumentParser(description="RS Tools Command Line Interface")
    parser.add_argument('--lattice', action='store_true', default=False, help='Process lattice')
    parser.add_argument('--visualize', action='store_true', default=False, help='Run visualizer')

    args = parser.parse_args()

    if args.lattice:
        process_lattice()
    elif args.visualize:
        rs_visualizer()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
