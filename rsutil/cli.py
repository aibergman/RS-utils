import argparse
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer


def main():
    """Command-line interface for rsutil."""
    parser = argparse.ArgumentParser(description="RS Tools Command Line Interface")
    parser.add_argument('--setup', '-S', action='store_true', default=False, help='Process lattice')
    parser.add_argument('--lattice', '-L', action='store_true', default=False, help='Run visualizer for lattice input')
    parser.add_argument('--clust', '-C', action='store_true', default=False, help='Run visualizer for clust file')
    parser.add_argument('--data', type=str, default='.', help='Directory for data files to be visualized')

    args = parser.parse_args()

    if args.setup:
        process_lattice()
    elif args.lattice or args.clust:
        rs_visualizer(args.data, clust=args.clust)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
