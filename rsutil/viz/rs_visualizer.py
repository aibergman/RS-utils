"""
This module provides classes and functions for visualizing atomic structures
using pyvista and pymatgen. It includes classes for representing atoms, bonds,
and spins, as well as methods for reading structure files and rendering the
structures interactively.

Classes:
    Atom: Represents an atom for visualization.
    Bond: Represents a bond (rendered as a cylinder).
    Spin: Represents a magnetic moment as an arrow.
    Structure: Contains all the structural information.
    InputReader: Provides methods to read structure files.
    Visualizer: Visualizes a Structure using pyvista.

Functions:
    get_element_color(element): Returns a standard color for a given element
    symbol.
"""
from inputreader import InputReader
from visualizer import Visualizer


# ----------------------------
# 5. Example usage
# ----------------------------
if __name__ == "__main__":
    # Read a POSCAR file using pymatgen
    # structure = InputReader.read_poscar(
    #     "/Users/andersb/Jobb/RS_tools_python/POSCAR"
    # )
    structure = InputReader.read_rs_lmto(
        "/Users/andersb/Jobb/RS_tools_python/examples/Mn3Sn/lattice.nml",
        "/Users/andersb/Jobb/RS_tools_python/examples/Mn3Sn/input.nml",
    )  # <-- update with your file path

    # Optionally, compute the conventional cell using pymatgen's symmetry routines.
    try:
        structure.compute_conventional_cell()
    except (ValueError, AttributeError) as e:
        print("Conventional cell computation failed:", e)

    if structure.unit_cell is not None:
        print("Unit Cell:")
        print(structure.unit_cell)
        print(len(structure.atoms_uc))

    if structure.conventional_cell is not None:
        print("Conventional Cell:")
        print(len(structure.atoms_cc))

    # Create a Visualizer instance with the structure.
    viz = Visualizer(structure=structure)
    # Uncomment toggles as desired:

    viz.render()
