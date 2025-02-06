from glyphs import Atom
from helpers import get_element_color
from structure import Structure
from pymatgen.io.vasp import Poscar


# ----------------------------
# 3. Input Reader Class
# ----------------------------
class InputReader:
    """
    A utility class for reading different file formats and converting them
    into Structure objects.

    Methods:
        read_xyz(file_path):
            Reads a simple .xyz file and returns a Structure object.

        read_custom(file_path):
            Stub for reading a custom file format. Implement your own file
            parser here.

        read_poscar(file_path):
            Reads a POSCAR file using pymatgen and returns a Structure object.
    """
    @staticmethod
    def read_xyz(file_path):
        """
        Reads a simple .xyz file and returns a Structure.
        Assumes the first line is the number of atoms and the second line is a comment.
        Each subsequent line: species x y z
        """
        atoms = []
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        num_atoms = int(lines[0].strip())
        # Comment line (unused here)
        _ = lines[1].strip()
        for line in lines[2: 2 + num_atoms]:
            parts = line.split()
            species = parts[0]
            position = [float(x) for x in parts[1:4]]
            atoms.append(
                Atom(
                    position=position,
                    species=species,
                    atom_type=species,
                    color=get_element_color(species),
                )
            )
        return Structure(atoms=atoms)

    @staticmethod
    def read_custom(file_path):
        """
        Stub for reading a custom file format.
        Implement your own file parser here.
        """
        pass

    @staticmethod
    def read_poscar(file_path):
        """
        Reads a POSCAR file using pymatgen and returns a Structure.
        """
        poscar = Poscar.from_file(file_path)
        pg_struct = poscar.structure
        atoms = []
        for site in pg_struct:
            atoms.append(
                Atom(
                    position=site.coords,
                    species=site.species_string,
                    atom_type=site.species_string,
                    color=get_element_color(site.species_string),
                )
            )
        return Structure(atoms=atoms, unit_cell=pg_struct.lattice.matrix)
