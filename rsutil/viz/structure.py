from itertools import product
import numpy as np
from pymatgen.core import Structure as PymatgenStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from glyphs import Atom, Bond, Spin
from helpers import get_element_color


class Structure:
    """
    Represents a molecular structure with atoms, bonds, and spins.

    Attributes:
        atomsc (list): List of Atom objects.
        bonds (list): List of Bond objects.
        spins (list): List of Spin objects.
        atoms_cc (list): List of Atom objects in the conventional cell.
        bonds_cc (list): List of Bond objects in the conventional cell.
        spins_cc (list): List of Spin objects in the conventional cell.
        unit_cell (np.ndarray): Lattice vectors for the unit cell.
        conventional_cell (np.ndarray): Lattice vectors for the
            conventional cell.

    Methods:
        get_atom_spin_objects():
            Returns a list of tuples (Atom, Spin). If a spin is not
            provided for a given atom, a default (zero magnitude) spin
            is created.
        get_bond_objects():
            Returns a list of Bond objects. If the bonds list is empty,
            automatically generates bonds based on the distance cutoff.
        compute_conventional_cell(tol=0.01):
            Compute the conventional cell of the structure using
            pymatgen's SpacegroupAnalyzer.
    """

    def __init__(
        self,
        atoms,
        bonds=None,
        spins=None,
        unit_cell=None,
        conventional_cell=None,
        bond_cutoff=3.0,
    ):
        """
        Contains all the structural information.

        Args:
            atoms (list): List of Atom objects.
            bonds (list): List of Bond objects (can be empty; bonds may be generated).
            spins (list): List of Spin objects (optional; default spin will be provided if missing).
            unit_cell (3x3 array-like): Lattice vectors for the unit cell.
            conventional_cell (3x3 array-like): Lattice vectors for the conventional cell.
            bond_cutoff (float): Distance below which atoms are bonded.
        """
        self.atoms_uc = atoms
        self.bonds_uc = bonds if bonds is not None else []
        self.spins_uc = spins if spins is not None else []
        self.atoms_cc = []
        self.bonds_cc = []
        self.spins_cc = []
        self.unit_cell = np.array(unit_cell) if unit_cell is not None else None
        self.conventional_cell = (
            np.array(conventional_cell) if conventional_cell is not None else None
        )
        self.bond_cutoff = bond_cutoff

    def _create_pymatgen_structure(self):
        """
        Create a pymatgen Structure object from the current unit cell and atom list.

        Returns:
            PymatgenStructure: A pymatgen Structure object.
        """
        lattice = self.unit_cell
        species = [atom.species for atom in self.atoms_uc]
        positions = [np.dot(np.linalg.inv(lattice).T, atom.position) for atom in self.atoms_uc]
        return PymatgenStructure(lattice, species, positions)

    def get_atom_spin_objects(self, conventional=False):
        """
        Returns a list of tuples (Atom, Spin). If a spin is not provided
        for a given atom, a default (zero magnitude) spin is created.

        Args:
            conventional (bool): If True, use atoms_cc and spins_cc; otherwise, use atoms_uc and spins_uc.
        """
        atom_spin_list = []
        atoms = self.atoms_cc if conventional else self.atoms_uc
        spins = self.spins_cc if conventional else self.spins_uc

        for atom in atoms:
            # Look for a spin whose position matches the atom (within tolerance)
            spin = next(
                (s for s in spins if np.allclose(s.position, atom.position)), None
            )
            if spin is None:
                spin = Spin(position=atom.position, direction=(0, 0, 0), magnitude=0.0)
            atom_spin_list.append((atom, spin))
        return atom_spin_list

    def get_bond_objects(self, conventional=False):
        """
        Returns a list of Bond objects. If the bonds list is empty, automatically
        generates bonds based on the distance cutoff.

        Args:
            conventional (bool): If True, use bonds_cc; otherwise, use bonds_uc.
        """
        if conventional:
            if not self.bonds_cc:
                n = len(self.atoms_cc)
                for i in range(n):
                    for j in range(i + 1, n):
                        pos_i = self.atoms_cc[i].position
                        pos_j = self.atoms_cc[j].position
                        if np.linalg.norm(pos_i - pos_j) <= self.bond_cutoff:
                            self.bonds_cc.append(Bond(start=pos_i, end=pos_j))
            return self.bonds_cc
        else:
            if not self.bonds_uc:
                n = len(self.atoms_uc)
                for i in range(n):
                    for j in range(i + 1, n):
                        pos_i = self.atoms_uc[i].position
                        pos_j = self.atoms_uc[j].position
                        if np.linalg.norm(pos_i - pos_j) <= self.bond_cutoff:
                            self.bonds_uc.append(Bond(start=pos_i, end=pos_j))
            return self.bonds_uc

    def compute_conventional_cell(self, tol=0.01):
        """
        Compute the conventional cell of the structure using pymatgen's SpacegroupAnalyzer.
        This method creates the following new attributes:
          - self.conventional_cell : a 3x3 array representing the conventional cell lattice.
          - self.atoms_cc     : a list of Atom objects for the conventional cell,
                                including those at the boundaries.
          - self.bonds_cc     : a list of Bond objects computed from atoms_cc based on self.bond_cutoff.
          - self.spins_cc     : a list of Spin objects (default zero magnitude) corresponding to atoms_cc.

        The procedure is as follows:
          1. Build a pymatgen Structure from the current unit cell and atom list.
          2. Use SpacegroupAnalyzer to get the conventional standard structure.
          3. Retrieve the fractional coordinates and species from the conventional structure.
          4. Expand the fractional positions using shifts in [-1, 0, 1]³ to capture boundary atoms.
          5. Convert the expanded fractional coordinates into Cartesian coordinates using the new lattice.
          6. Create the new Atom, Bond, and Spin lists accordingly.

        Args:
            tol (float): Tolerance for including atoms that lie on the boundaries (default 0.01).
        """

        # Build a pymatgen structure from current unit cell and atoms.
        pg_struct = self._create_pymatgen_structure()
        print("pg_struct"+'-'*50)
        print(pg_struct)

        # Use SpacegroupAnalyzer to get the conventional standard structure.
        sga = SpacegroupAnalyzer(pg_struct)
        print("spacegroup_struct"+'-'*50)
        print(sga.get_conventional_standard_structure())
        print(sga.get_space_group_number())
        print(sga.get_space_group_symbol())
        conv_struct = sga.get_conventional_standard_structure()

        # Save the conventional cell lattice.
        self.conventional_cell = conv_struct.lattice.matrix

        # Retrieve fractional coordinates and corresponding species.
        frac_coords = conv_struct.frac_coords
        symbols = [site.species_string for site in conv_struct]

        # Expand fractional coordinates to capture boundary atoms.
        expanded_frac = []
        expanded_symbols = []
        # Generate shifts in each direction: [-1, 0, 1] for x, y, z.
        shifts = np.array(list(product([-1, 0, 1], repeat=3)))
        for shift in shifts:
            for pos, symbol in zip(frac_coords, symbols):
                new_pos = pos + shift
                # Only keep positions that (within a tolerance) lie in [0, 1) in all directions.
                if np.all(new_pos >= -tol) and np.all(new_pos < 1 + tol):
                    expanded_frac.append(new_pos)
                    expanded_symbols.append(symbol)

        # Convert the expanded fractional positions to Cartesian coordinates.
        expanded_cart = np.dot(expanded_frac, self.conventional_cell)

        # Create a new list of Atom objects for the conventional cell.
        self.atoms_cc = []
        for pos, symbol in zip(expanded_cart, expanded_symbols):
            self.atoms_cc.append(
                Atom(
                    position=pos,
                    species=symbol,
                    atom_type=symbol,
                )
            )

        # Compute bonds for the conventional cell.
        self.bonds_cc = []
        n = len(self.atoms_cc)
        for i in range(n):
            for j in range(i + 1, n):
                pos_i = self.atoms_cc[i].position
                pos_j = self.atoms_cc[j].position
                if np.linalg.norm(pos_i - pos_j) <= self.bond_cutoff:
                    self.bonds_cc.append(Bond(start=pos_i, end=pos_j))

        # Create a default spin list (all spins with zero magnitude).
        self.spins_cc = []
        for atom in self.atoms_cc:
            self.spins_cc.append(
                Spin(position=atom.position, direction=(0, 0, 0), magnitude=0.0)
            )
