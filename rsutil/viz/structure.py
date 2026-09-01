"""
Module for representing and manipulating molecular structures with atoms,
bonds, and spins. This module provides a `Structure` class that allows for
the creation, manipulation, and analysis of molecular structures, including
the computation of conventional cells using pymatgen's SpacegroupAnalyzer.

Classes:
    Structure: Represents a molecular structure with atoms, bonds, and spins.

Attributes:
    atoms_uc (list): List of Atom objects in the unit cell.
    bonds_uc (list): List of Bond objects in the unit cell.
    spins_uc (list): List of Spin objects in the unit cell.
    conventional_cell (np.ndarray): Lattice vectors for the conventional cell.

Methods:
    __init__(self, atoms, bonds=None, spins=None, unit_cell=None,
             conventional_cell=None, bond_cutoff=3.0):
        Initializes the Structure object with atoms, bonds, spins, and lattice
        vectors.

    _create_pymatgen_structure(self) -> PymatgenStructure:
        Creates a pymatgen Structure object from the current unit cell and atom
        list.

    get_atom_spin_objects(self, conventional=False) -> list[tuple[Atom, Spin]]:
        Returns a list of tuples (Atom, Spin). If a spin is not provided for a
        given atom, a default (zero magnitude) spin is created.

    get_bond_objects(self, conventional=False) -> list[Bond]:
        Returns a list of Bond objects. If the bonds list is empty, generates
        bonds based on the distance cutoff.

    compute_conventional_cell(self, tol=0.01) -> None:
        Computes the conventional cell of the structure using pymatgen's
        SpacegroupAnalyzer. This method creates new attributes for the
        conventional cell lattice, atoms, bonds, and spins.
"""
from itertools import product
import numpy as np

from pymatgen.core import Structure as PymatgenStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .glyphs import Atom, Bond, Spin


class Structure:
    """
    Represents a molecular structure with atoms, bonds, and spins.

    Attributes
    ----------
    atoms_uc : list
        List of Atom objects in the unit cell.
    bonds_uc : list
        List of Bond objects in the unit cell.
    spins_uc : list
        List of Spin objects in the unit cell.
    atoms_cc : list
        List of Atom objects in the conventional cell.
    bonds_cc : list
        List of Bond objects in the conventional cell.
    spins_cc : list
        List of Spin objects in the conventional cell.
    unit_cell : np.ndarray
        Lattice vectors for the unit cell.
    conventional_cell : np.ndarray
        Lattice vectors for the conventional cell.

    """

    def __init__(
        self,
        atoms: list,
        bonds: list = None,
        spins: list = None,
        unit_cell: np.ndarray = None,
        conventional_cell: np.ndarray = None,
        bond_cutoff: float = 3.0,
    ):
        """
        Initialize the Structure object.

        Args:
            atoms (list): List of Atom objects.
            bonds (list): List of Bond objects.
            spins (list): List of Spin objects.
            unit_cell (array-like): Lattice vectors for the unit cell.
            conventional_cell (array-like): Lattice vectors for the conventional cell.
            bond_cutoff (float): Distance below which atoms are bonded.
        """
        self.atoms_uc: list = atoms
        self.bonds_uc: list = bonds if bonds is not None else []
        self.spins_uc: list = spins if spins is not None else []
        self.atoms_cc: list = []
        self.bonds_cc: list = []
        self.spins_cc: list = []
        self.unit_cell: np.ndarray = np.array(unit_cell) if unit_cell is not None else None
        self.conventional_cell: np.ndarray = (
            np.array(conventional_cell) if conventional_cell is not None else None
        )
        self.bond_cutoff: float = bond_cutoff

    def _create_pymatgen_structure(self) -> PymatgenStructure:
        """
        Create a pymatgen Structure from the current unit cell and atoms.

        Returns:
            PymatgenStructure: A pymatgen Structure object.
        """
        lattice: np.ndarray = self.unit_cell.T
        species: list = [atom.species for atom in self.atoms_uc]
        positions: list = [
            np.dot(np.linalg.inv(lattice).T, atom.position) for atom in self.atoms_uc
        ]

        return PymatgenStructure(lattice, species, positions)

    def get_atom_spin_objects(self, conventional: bool = False) -> list[tuple[Atom, Spin]]:
        """
        Get list of (Atom, Spin) tuples. Create default spin if not provided.

        Args:
            conventional (bool): Use conventional cell if True.
        """
        atom_spin_list: list[tuple[Atom, Spin]] = []
        atoms: list[Atom] = self.atoms_cc if conventional else self.atoms_uc
        spins: list[Spin] = self.spins_cc if conventional else self.spins_uc

        for atom in atoms:
            # Look for a spin whose position matches the atom (within tolerance)
            spin: Spin = next(
                (s for s in spins if np.allclose(s.position, atom.position)), None
            )
            if spin is None:
                spin = Spin(position=atom.position, direction=np.array([0, 0, 0]), magnitude=0.0)
            atom_spin_list.append((atom, spin))
        return atom_spin_list

    def get_bond_objects(self, conventional: bool = False) -> list[Bond]:
        """
        Get list of Bond objects. Generate bonds if list is empty.

        Args:
            conventional (bool): Use conventional cell if True.
        """
        if conventional:
            if not self.bonds_cc:
                n: int = len(self.atoms_cc)
                for i in range(n):
                    for j in range(i + 1, n):
                        pos_i: np.ndarray = self.atoms_cc[i].position
                        pos_j: np.ndarray = self.atoms_cc[j].position
                        if np.linalg.norm(pos_i - pos_j) <= self.bond_cutoff:
                            self.bonds_cc.append(Bond(start=pos_i, end=pos_j))
            return self.bonds_cc
        else:
            if not self.bonds_uc:
                n: int = len(self.atoms_uc)
                for i in range(n):
                    for j in range(i + 1, n):
                        pos_i: np.ndarray = self.atoms_uc[i].position
                        pos_j: np.ndarray = self.atoms_uc[j].position
                        if np.linalg.norm(pos_i - pos_j) <= self.bond_cutoff:
                            self.bonds_uc.append(Bond(start=pos_i, end=pos_j))
            return self.bonds_uc

    def compute_conventional_cell(self, tol: float = 0.01) -> None:
        """
        Compute the conventional cell of the structure using pymatgen's 
        SpacegroupAnalyzer. This method creates the following new attributes:
        
        - self.conventional_cell : a 3x3 array representing the conventional 
          cell lattice.
        - self.atoms_cc : a list of Atom objects for the conventional cell, 
          including those at the boundaries.
        - self.bonds_cc : a list of Bond objects computed from atoms_cc based 
          on self.bond_cutoff.
        - self.spins_cc : a list of Spin objects (default zero magnitude) 
          corresponding to atoms_cc.

        The procedure is as follows:
        
        1. Build a pymatgen Structure from the current unit cell and atom list.
        2. Use SpacegroupAnalyzer to get the conventional standard structure.
        3. Retrieve the fractional coordinates and species from the 
           conventional structure.
        4. Expand the fractional positions using shifts in [-1, 0, 1]³ to 
           capture boundary atoms.
        5. Convert the expanded fractional coordinates into Cartesian 
           coordinates using the new lattice.
        6. Create the new Atom, Bond, and Spin lists accordingly.

        Args:
            tol (float): Tolerance for including atoms that lie on the 
                 boundaries (default 0.01).
        """

        # Build a pymatgen structure from current unit cell and atoms.
        pg_struct: PymatgenStructure = self._create_pymatgen_structure()
        print("--Unit cell" + "-" * 39)
        print(pg_struct.lattice)
        print(pg_struct)
        # Use SpacegroupAnalyzer to get the conventional standard structure.
        sga: SpacegroupAnalyzer = SpacegroupAnalyzer(pg_struct)
        conv_struct: PymatgenStructure = sga.get_conventional_standard_structure()
        print("--Conventional cell" + "-" * 31)
        print(conv_struct.lattice)
        print(conv_struct.lattice.matrix)
        print(conv_struct.frac_coords)
        print(conv_struct)
        print("-" * 50)
        print("Symmetry:", sga.get_space_group_number(), sga.get_space_group_symbol())
        print("-" * 50)

        # Save the conventional cell lattice.
        self.conventional_cell: np.ndarray = conv_struct.lattice.matrix.T

        # Retrieve fractional coordinates and corresponding species.
        frac_coords: np.ndarray = conv_struct.frac_coords
        symbols: list[str] = [site.species_string for site in conv_struct]

        # Expand fractional coordinates to capture boundary atoms.
        expanded_frac: list[np.ndarray] = []
        expanded_symbols: list[str] = []
        # Generate shifts in each direction: [-1, 0, 1] for x, y, z.
        shifts: np.ndarray = np.array(list(product([-1, 0, 1], repeat=3)))
        for shift in shifts:
            for pos, symbol in zip(frac_coords, symbols):
                new_pos = pos + shift
                # Only keep positions that (within a tolerance) lie in [0, 1) in all directions.
                if np.all(new_pos >= -tol) and np.all(new_pos < 1 + tol):
                    expanded_frac.append(new_pos)
                    expanded_symbols.append(symbol)

        # Convert the expanded fractional positions to Cartesian coordinates.
        expanded_cart: np.ndarray = np.dot(expanded_frac, self.conventional_cell.T)

        # Create a new list of Atom objects for the conventional cell.
        self.atoms_cc: list[Atom] = []
        for pos, symbol in zip(expanded_cart, expanded_symbols):
            self.atoms_cc.append(
                Atom(
                    position=np.array(pos),
                    species=symbol,
                    atom_type=symbol,
                )
            )

        # # Compute bonds for the conventional cell.
        # self.bonds_cc = []
        # n = len(self.atoms_cc)
        # for i in range(n):
        #     for j in range(i + 1, n):
        #         pos_i = self.atoms_cc[i].position
        #         pos_j = self.atoms_cc[j].position
        #         if np.linalg.norm(pos_i - pos_j) <= self.bond_cutoff:
        #             self.bonds_cc.append(Bond(start=pos_i, end=pos_j))

        # Create a default spin list (all spins with zero magnitude).
        self.spins_cc: list[Spin] = []
        for atom in self.atoms_cc:
            self.spins_cc.append(
                Spin(position=np.array(atom.position), direction=np.array([0, 0, 0]), magnitude=0.0)
            )
