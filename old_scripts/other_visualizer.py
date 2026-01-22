import numpy as np
import pyvista as pv
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

class Atom:
    def __init__(self, position, species, atom_type, radius=0.2):
        self.position = np.array(position)
        self.species = species
        self.atom_type = atom_type
        self.radius = radius

class Bond:
    def __init__(self, start, end, radius=0.1):
        self.start = np.array(start)
        self.end = np.array(end)
        self.radius = radius

class Spin:
    def __init__(self, position, direction, magnitude=1.0, base_radius=0.05, top_radius=0.1):
        self.position = np.array(position)
        self.direction = np.array(direction)
        self.magnitude = magnitude
        self.base_radius = base_radius
        self.top_radius = top_radius

class Structure:
    def __init__(self, lattice, atom_data, bond_cutoff=2.5):
        self.lattice = lattice
        self.atom_data = atom_data  # List of tuples: (position, species, type)
        self.bond_cutoff = bond_cutoff
        self.atoms = [Atom(pos, species, atom_type) for pos, species, atom_type in atom_data]
        self.spins = [Spin(pos, [0, 0, 0]) for pos, _, _ in atom_data]  # Default zero spins
        self.bonds = self._calculate_bonds()

    def _calculate_bonds(self):
        bonds = []
        for i, atom1 in enumerate(self.atoms):
            for j, atom2 in enumerate(self.atoms):
                if i < j:
                    dist = np.linalg.norm(atom1.position - atom2.position)
                    if dist < self.bond_cutoff:
                        bonds.append(Bond(atom1.position, atom2.position))
        return bonds

    def find_conventional_cell(self):
        structure = Structure(self.lattice, [a.species for a in self.atoms], [a.position for a in self.atoms])
        sga = SpacegroupAnalyzer(structure)
        conventional_structure = sga.get_conventional_standard_structure()
        self.lattice = conventional_structure.lattice.matrix
        self.atoms = [Atom(pos, species, "") for pos, species in zip(conventional_structure.frac_coords, conventional_structure.species)]

class InputReader:
    @staticmethod
    def read_poscar(filename):
        with open(filename, "r") as f:
            lines = f.readlines()
        scale = float(lines[1].strip())
        lattice = np.array([list(map(float, line.split())) for line in lines[2:5]]) * scale
        atom_types = lines[5].split()
        atom_counts = list(map(int, lines[6].split()))
        positions = np.array([list(map(float, line.split())) for line in lines[8:8+sum(atom_counts)]])
        atom_data = [(pos, atom_types[i], atom_types[i]) for i, pos in enumerate(positions)]
        return Structure(lattice, atom_data)

class Visualizer:
    def __init__(self, structure):
        self.structure = structure
        self.plotter = pv.Plotter()
        self.show_bonds = False
        self.show_spins = False
        self.show_extents = False
        self.use_conventional_cell = False

    def render_atoms(self):
        for atom in self.structure.atoms:
            sphere = pv.Sphere(radius=atom.radius, center=atom.position)
            self.plotter.add_mesh(sphere, color="blue", smooth_shading=True)

    def render_bonds(self):
        if not self.show_bonds:
            return
        for bond in self.structure.bonds:
            cylinder = pv.Cylinder(center=(bond.start + bond.end) / 2, direction=bond.end - bond.start, radius=bond.radius)
            self.plotter.add_mesh(cylinder, color="gray")

    def render_spins(self):
        if not self.show_spins:
            return
        for spin in self.structure.spins:
            arrow = pv.Arrow(start=spin.position, direction=spin.direction, scale=spin.magnitude)
            self.plotter.add_mesh(arrow, color="red")

    def visualize(self):
        self.plotter.clear()
        self.render_atoms()
        self.render_bonds()
        self.render_spins()
        self.plotter.show()

# Example usage
if __name__ == "__main__":
    structure = InputReader.read_poscar("/Users/andersb/Jobb/RS_tools_python/POSCAR")
    visualizer = Visualizer(structure)
    visualizer.visualize()
