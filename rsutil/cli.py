import click
from . import process_lattice

@click.command()
def cli():
    """Runs the Fortran process_lattice function"""
    process_lattice()

if __name__ == "__main__":
    cli()

