from itertools import product
import numpy as np
import pyvista as pv
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class InputReader:
    """
    A class to read and parse POSCAR files.

    Attributes:
        filename (str): The path to the POSCAR file.

    Methods:
        read_poscar():
            Reads the POSCAR file and extracts lattice vectors, atom positions,
            and atomic symbols.
    """
    def __init__(self, filename):
        self.filename = filename

    def read_poscar(self):
        with open(self.filename, "r", encoding="utf-8") as f:
            lines = f.readlines()

        scale = float(lines[1].strip())
        lattice_vectors = np.array(
            [list(map(float, line.split())) for line in lines[2:5]]
        ) * scale
        atom_types = lines[5].split()
        atom_counts = list(map(int, lines[6].split()))
        total_atoms = sum(atom_counts)
        atom_positions = np.array(
            [list(map(float, line.split())) for line in lines[8: 8 + total_atoms]]
        )
        atomic_symbols = []
        for atom_type, count in zip(atom_types, atom_counts):
            atomic_symbols.extend([atom_type] * count)
        return (lattice_vectors, atom_positions, atomic_symbols)


class Atom:
    def __init__(self, position, atom_type, color, radius=1.0):
        self.position = position
        self.atom_type = atom_type
        self.color = color
        self.radius = radius


class AtomicStructureVisualizer:
    def __init__(self, poscar_file):
        reader = InputReader(poscar_file)
        self.plotter = pv.Plotter()
        self.unit_cell = reader.read_poscar()
        self.lattice, self.positions, self.atomic_symbols = self.unit_cell
        self.atom_colors = self._generate_atom_colors()
        self.unit_cell = self._create_atoms(self.positions, self.atomic_symbols)
        self.find_conventional_cell()

    def _generate_atom_colors(self):
        color_palette = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        ]
        atom_colors = {}
        for i, atom in enumerate(set(self.atomic_symbols)):
            color = color_palette[i % len(color_palette)]
            atom_colors[atom] = color
        return atom_colors

    def _create_atoms(self, positions, atomic_symbols):
        atoms = []
        for pos, atom_symbol in zip(positions, atomic_symbols):
            color = self.atom_colors[atom_symbol]
            atoms.append(Atom(pos, atom_symbol, color))
        return atoms

    def find_conventional_cell(self):
        structure = Structure(self.lattice, self.atomic_symbols, self.positions)
        sga = SpacegroupAnalyzer(structure)
        conventional_structure = sga.get_conventional_standard_structure()

        self.lattice = conventional_structure.lattice.matrix
        self.positions = conventional_structure.cart_coords
        self.atomic_symbols = [site.species_string for site in conventional_structure]

        # Expand to capture all atoms at the boundaries
        fractional_positions = conventional_structure.frac_coords
        expanded_positions = []
        expanded_symbols = []

        shifts = np.array(list(product([-1, 0, 1], repeat=3)))
        for shift in shifts:
            for pos, symbol in zip(fractional_positions, self.atomic_symbols):
                new_pos = pos + shift
                if np.all((new_pos >= -0.01) & (new_pos < 1.01)):
                    expanded_positions.append(new_pos)
                    expanded_symbols.append(symbol)

        expanded_positions = np.dot(expanded_positions, self.lattice)  # Convert to Cartesian
        self.conventional_cell = self._create_atoms(expanded_positions, expanded_symbols)

    def show_atoms(self, atoms):
        for atom in atoms:
            sphere = pv.Sphere(radius=atom.radius, center=atom.position)
            print(atom.color)
            self.plotter.add_mesh(sphere, color=atom.color, smooth_shading=True)

    def visualize(self, show_bonds=False, show_cell=False, use_conventional_cell=False, bond_threshold=2.0):
        atoms_to_show = self.conventional_cell if use_conventional_cell else self.unit_cell
        self.show_atoms(atoms_to_show)
        self.plotter.show()


# Example usage
if __name__ == "__main__":
    visualizer = AtomicStructureVisualizer("POSCAR")
    visualizer.visualize(show_bonds=False, show_cell=True, use_conventional_cell=True)
