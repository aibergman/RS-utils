import numpy as np
import pyvista as pv
import re
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element

class InputReader:
    def __init__(self, filename):
        self.filename = filename

    def read_poscar(self):
        with open(self.filename, "r", encoding="utf-8") as f:
            lines = f.readlines()

        scale = float(lines[1].strip())
        lattice_vectors = np.array(
            [list(map(float, line.split())) for line in lines[2:5]]
        )
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

    def read_rs_lmto(self, lattice_file, input_file):
        """
        Reads RS-LMTO input files (lattice.nml and input.nml) and returns a Structure object.
        """
        # Read lattice.nml
        with open(lattice_file, "r", encoding="utf-8") as f:
            lines = f.readlines()

        # Extract lattice vectors
        lattice_vectors = []
        basis_vectors = []
        extracting_lattice = False
        extracting_basis = False
        for line in lines:
            if "a(:," in line:
                extracting_lattice = True
                parts = line.split("=")[-1].strip().split(",")
                lattice_vectors.append([float(p) for p in parts])
            elif extracting_lattice and "!" in line:
                extracting_lattice = False
            elif "crd(:," in line:
                extracting_basis = True
                parts = line.split("=")[-1].strip().split(",")
                basis_vectors.append([float(p) for p in parts])
            elif extracting_basis and "!" in line:
                extracting_basis = False
        
        lattice_vectors = np.array(lattice_vectors).T  # Ensure correct shape
        basis_vectors = np.array(basis_vectors)
        
        # Read input.nml to extract species labels
        species_labels = []
        with open(input_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        for line in lines:
            match = re.search(r"label\(\d+\)\s*=\s*'([^']+)'", line)
            if match:
                species = re.sub(r'[-0-9]', '', match.group(1))  # Clean label
                species_labels.append(species)
        
        # Assign species to basis vectors
        species_per_site = species_labels[: len(basis_vectors)]
        
        return Structure(lattice_vectors, basis_vectors, species_per_site)

class Atom:
    def __init__(self, position, atom_type, color, radius=1.0):  # Increased radius
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
        self._setup_lighting()  # Ensure lighting is properly configured

    def _setup_lighting(self):
        self.plotter.add_light(pv.Light(position=(10, 10, 10), focal_point=(0, 0, 0), intensity=1.5))
        self.plotter.enable_anti_aliasing()
        self.plotter.enable_shadows()

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

    def show_atoms(self, atoms):
        for atom in atoms:
            sphere = pv.Sphere(radius=atom.radius, center=atom.position)
            sphere.compute_normals(cell_normals=True, point_normals=True, auto_orient_normals=True)
            self.plotter.add_mesh(
                sphere,
                color=atom.color,
                smooth_shading=True,
                specular=0.5,
                specular_power=10,
                lighting=True
            )

    def visualize(self, show_bonds=False, show_cell=False, use_conventional_cell=False, bond_threshold=2.0):
        atoms_to_show = self.conventional_cell if use_conventional_cell else self.unit_cell
        self.show_atoms(atoms_to_show)
        self.plotter.show()

# Example usage
if __name__ == "__main__":
    visualizer = AtomicStructureVisualizer("POSCAR")
    visualizer.visualize(show_bonds=False, show_cell=True, use_conventional_cell=True)
