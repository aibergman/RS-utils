import numpy as np
import pyvista as pv
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class AtomicStructureVisualizer:
    def __init__(self, poscar_file):
        self.lattice, self.atom_types, self.atom_counts, self.positions = (
            self.read_poscar(poscar_file)
        )
        self.plotter = pv.Plotter()
        self.plotter.enable_ssao()
        self.atom_colors = self._generate_atom_colors()
        self.sphere_radius = 0.5

    def read_poscar(self, filename):
        with open(filename, "r") as f:
            lines = f.readlines()

        scale = float(lines[1].strip())
        lattice_vectors = np.array(
            [list(map(float, line.split())) for line in lines[2:5]]
        )  # * scale
        atom_types = lines[5].split()
        print(atom_types)
        atom_counts = list(map(int, lines[6].split()))
        total_atoms = sum(atom_counts)
        atom_positions = np.array(
            [list(map(float, line.split())) for line in lines[8 : 8 + total_atoms]]
        )
        return lattice_vectors, atom_types, atom_counts, atom_positions

    def _generate_atom_colors(self):
        color_palette = list(pv.plotting.colors.hexcolors.values())  # Use hex colors
        return {
            atom: color_palette[i % len(color_palette)]
            for i, atom in enumerate(self.atom_types)
        }

    def show_atoms(self):
        for i, pos in enumerate(self.positions):
            cumulative_counts = np.cumsum(self.atom_counts)
            atom_index = next(
                (j for j, v in enumerate(cumulative_counts) if v > i),
                len(cumulative_counts) - 1,
            )
            atom_type = self.atom_types[atom_index]
            color = self.atom_colors[atom_type]
            sphere = pv.Sphere(radius=self.sphere_radius, center=pos)
            self.plotter.add_mesh(sphere, color=color, smooth_shading=True)

    def show_bonds(self, bond_threshold=2.0):
        for i, pos_i in enumerate(self.positions):
            for j, pos_j in enumerate(self.positions[i + 1 :], start=i + 1):
                dist = np.linalg.norm(pos_i - pos_j)
                if dist < bond_threshold:
                    line = pv.Line(pos_i, pos_j)
                    self.plotter.add_mesh(line, color="gray", line_width=3)

    def show_cell(self):
        origin = np.zeros(3)
        points = [
            origin,
            self.lattice[0],
            self.lattice[0] + self.lattice[1],
            self.lattice[1],
            self.lattice[2],
            self.lattice[0] + self.lattice[2],
            self.lattice[0] + self.lattice[1] + self.lattice[2],
            self.lattice[1] + self.lattice[2],
        ]

        edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
        ]

        lines = pv.PolyData()
        lines.points = np.array(points)
        lines.lines = np.hstack([[2, edge[0], edge[1]] for edge in edges])

        self.plotter.add_mesh(lines, color="black", line_width=1)

    def find_conventional_cell(self, repeat_atoms=False):
        # Create a pymatgen Structure object
        structure = Structure(self.lattice, self.atom_types, self.positions)

        # Use SpacegroupAnalyzer to find the conventional cell
        sga = SpacegroupAnalyzer(structure)
        conventional_structure = sga.get_conventional_standard_structure()

        # Update the lattice and positions with the conventional cell
        self.lattice = conventional_structure.lattice.matrix
        self.positions = conventional_structure.frac_coords

        if repeat_atoms:
            # Repeat atoms to include those at the border
            repeated_structure = conventional_structure * (4, 4, 4)
            repeated_structure.translate_sites(
                range(len(repeated_structure)), [-2, -2, -2]
            )
            self.positions = repeated_structure.cart_coords
            self.atom_types = [site.species_string for site in repeated_structure]
            self.atom_counts = [
                self.atom_types.count(atom) for atom in set(self.atom_types)
            ]

            # Remove atoms outside the conventional cell with a tolerance of +- 5%
            tolerance = 0.05
            within_cell_indices = [
                i
                for i, pos in enumerate(self.positions)
                if np.all((pos >= -tolerance) & (pos <= 1 + tolerance))
            ]
            self.positions = self.positions[within_cell_indices]
            self.atom_types = [self.atom_types[i] for i in within_cell_indices]
            self.atom_counts = [
                self.atom_types.count(atom) for atom in set(self.atom_types)
            ]

    def visualize(
        self,
        show_bonds=False,
        show_cell=False,
        use_conventional_cell=False,
        bond_threshold=2.0,
    ):
        if use_conventional_cell:
            self.find_conventional_cell(repeat_atoms=True)

        self.show_atoms()
        if show_bonds:
            self.show_bonds(bond_threshold)
        if show_cell:
            self.show_cell()
        self.plotter.show()


# Example usage
if __name__ == "__main__":
    visualizer = AtomicStructureVisualizer("POSCAR")
    visualizer.visualize(show_bonds=False, show_cell=True, use_conventional_cell=True)
