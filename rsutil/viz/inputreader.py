import re
import f90nml
import numpy as np
from pymatgen.io.vasp import Poscar

from .glyphs import Atom
from .structure import Structure

# ----------------------------
# 3. Input Reader Class
# ----------------------------
class InputReader:
    """
    read_rs_lmto(lattice_file, input_file):
        Parses lattice vectors and basis coordinates from lattice.nml.
        Parses species labels from input.nml.
        Returns a Structure object with the list of Atom objects and the unit cell.
    read_rs_lmto_nml(lattice_file, input_file):
        Parses lattice vectors and basis coordinates from lattice.nml.
        Parses species labels and scaling factor from input.nml.
        Returns a Structure object with the list of Atom objects and the scaled unit cell.
    read_rs_clust(clust_file, input_file):
        Parses atom positions and species from the cluster file.
        Parses species labels from input.nml.
        Returns a Structure object built from the parsed cluster file and input file.
    A utility class for reading different file formats and converting them
    into Structure objects.
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
            cart_coords = pg_struct.lattice.get_cartesian_coords(site.frac_coords)
            atoms.append(
                Atom(
                    position=cart_coords,
                    species=site.species_string,
                    atom_type=site.species_string,
                )
            )
        return Structure(atoms=atoms, unit_cell=pg_struct.lattice.matrix)
    
    @staticmethod
    def read_rs_lmto(lattice_file, input_file):
        """
        Reads RS-LMTO input from a lattice.nml file and an input.nml file.
        
        The lattice.nml file should contain a &lattice block with lines like:
          a(:, 1) = 1.00000000, 0.00000000, 0.00000000
          a(:, 2) = 0.00000000, 1.00000000, 0.00000000
          a(:, 3) = 0.00000000, 0.00000000, 1.00000000
          ...
          crd(:, 1) = 0.0, 0.0, 0.0
          crd(:, 2) = 0.5, 0.5, 0.5
          ...
        
        The input.nml file should contain an &atoms block with lines like:
          label(1) = 'Fe'
          label(2) = 'Co'
        
        Species labels are scrubbed (only alphabetic characters are kept).
        
        Returns:
            A Structure object with the list of Atom objects and the unit cell.
        """
        alat = 1.0 # default lattice parameter
        # --- Parse lattice.nml ---
        with open(lattice_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        # Remove comments (lines starting with "!") and blank lines.
        lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('!')]
        
        in_lattice_block = False
        lattice_vectors = []  # to store vectors from lines with "a(:, i)"
        basis_coords = []     # to store basis atom coordinates from "crd(:, i)"
        for line in lines:
            # Start of lattice block:
            if line.lower().startswith("&lattice"):
                in_lattice_block = True
                continue
            if line.startswith("/") or line.startswith("&end"):
                in_lattice_block = False
            if in_lattice_block:
                # Look for lattice vector lines
                if re.search(r"^a\(:,\s*\d+\)\s*=", line, re.IGNORECASE):
                    # e.g. "a(:, 1) = 1.00000000, 0.00000000, 0.00000000"
                    parts = line.split("=")
                    if len(parts) >= 2:
                        vec_str = parts[1].strip()
                        # Split by comma and convert to float
                        try:
                            vec = [float(val) for val in vec_str.split(',')]
                        except Exception as e:
                            raise ValueError(f"Error parsing lattice vector in line: {line}\n{e}")
                        lattice_vectors.append(vec)
                # Look for basis coordinate lines
                elif re.search(r"^crd\(:,\s*\d+\)\s*=", line, re.IGNORECASE):
                    # e.g. "crd(:, 1) = 0.0, 0.0, 0.0"
                    parts = line.split("=")
                    if len(parts) >= 2:
                        coord_str = parts[1].strip()
                        try:
                            coord = [float(val) for val in coord_str.split(',')]
                        except Exception as e:
                            raise ValueError(f"Error parsing basis coordinate in line: {line}\n{e}")
                        basis_coords.append(coord)
        # In our RS-LMTO file, the lines for "a" are given as columns.
        # For a 3x3 cell, we expect 3 lines. Convert to a NumPy array and transpose.
        if len(lattice_vectors) != 3:
            raise ValueError("Expected 3 lattice vector lines in lattice.nml")
        lattice_array = np.array(lattice_vectors)  # shape (3, 3) where each row is a vector from a(:, i)
        # Transpose so that each row of unit_cell is one lattice vector.
        unit_cell = lattice_array.T

        # --- Parse input.nml for species labels ---
        with open(input_file, 'r') as f:
            input_lines = f.readlines()
        input_lines = [line.strip() for line in input_lines if line.strip() and not line.strip().startswith('!')]
        in_atoms_block = False
        in_lattice_block = False
        species_labels = []
        for line in input_lines:
            if line.lower().startswith("&lattice"):
                in_lattice_block = True
                continue
            if line.startswith("/") or line.lower().startswith("&end"):
                in_lattice_block = False
            if in_lattice_block:
                # Look for the line starting with "alat"
                if line.lower().startswith("alat"):
                    parts = line.split("=")
                    if len(parts) >= 2:
                        try:
                            alat = float(parts[1].strip())
                        except ValueError as e:
                            raise ValueError(f"Error parsing 'alat' value in line: {line}\n{e}")
                    
            if line.lower().startswith("&atoms"):
                in_atoms_block = True
                continue
            if line.startswith("/") or line.lower().startswith("&end"):
                in_atoms_block = False
            if in_atoms_block:
                # Look for lines starting with "label("
                if re.search(r"^label\(\s*\d+\s*\)\s*=", line, re.IGNORECASE):
                    parts = line.split("=")
                    if len(parts) >= 2:
                        label_str = parts[1].strip()
                        # Remove any quotes
                        label_clean = label_str.strip("'\"")
                        # Remove hyphens and numbers (keep only letters)
                        label_clean = ''.join(ch for ch in label_clean if ch.isalpha())
                        species_labels.append(label_clean)
        if not species_labels:
            raise ValueError("No species labels found in input.nml &atoms block.")

        # Scale with alat
        unit_cell *= alat
        # Scale basis coordinates with alat
        basis_coords = [np.array(coord) * alat for coord in basis_coords]
        
        # --- Create Atom objects ---
        if len(basis_coords) != len(species_labels):
            # If there is only one species label, apply it to all atoms.
            if len(species_labels) == 1:
                species_labels = species_labels * len(basis_coords)
            else:
                raise ValueError("Number of basis coordinates does not match number of species labels.")
        atoms = []
        for coord, sp in zip(basis_coords, species_labels):
            atoms.append(Atom(position=coord, species=sp, atom_type=sp))
        
        # --- Create and return the Structure ---
        return Structure(atoms=atoms, unit_cell=unit_cell)
    
    @staticmethod
    def read_rs_lmto_nml(lattice_file, input_file):
        """
        Reads RS-LMTO input using f90nml.

        Parameters:
            lattice_file (str): Path to the lattice.nml file.
            input_file (str): Path to the input.nml file containing the &atoms and &lattice blocks.

        Returns:
            Structure: A Structure object built from the parsed namelists.
        """
        # Read the namelist files.
        lat_nml = f90nml.read(lattice_file)
        inp_nml = f90nml.read(input_file)

        # --- Parse lattice.nml ---
        # Assume the &lattice block in lattice.nml contains the following keys:
        #   "a(:, 1)", "a(:, 2)", "a(:, 3)" for the lattice vectors.
        #   "crd(:, 1)", "crd(:, 2)", "crd(:, 3)" for the basis (atom) coordinates.
        lattice_data = lat_nml.get('lattice', {})

        try:
            a1 = lattice_data["a(:, 1)"]
            a2 = lattice_data["a(:, 2)"]
            a3 = lattice_data["a(:, 3)"]
        except KeyError:
            raise ValueError("Could not find lattice vectors with keys 'a(:, 1)', 'a(:, 2)', 'a(:, 3)' in lattice.nml.")

        # Build a 3x3 lattice matrix.
        # The file appears to list the vectors as columns, so we form an array and transpose.
        lattice_array = np.array([a1, a2, a3]) # shape (3, 3)
        unit_cell = lattice_array  # Now each row is a lattice vector. Note lack of transpose!

        # Get basis (atom) coordinates.
        try:
            crd1 = lattice_data["crd(:, 1)"]
            crd2 = lattice_data["crd(:, 2)"]
            crd3 = lattice_data["crd(:, 3)"]
        except KeyError:
            raise ValueError("Could not find basis coordinates with keys 'crd(:, 1)', etc. in lattice.nml.")

        # Assume all three coordinate lists are of equal length.
        n_atoms = len(crd1)
        basis_coords = []
        for i in range(n_atoms):
            coord = [crd1[i], crd2[i], crd3[i]]
            basis_coords.append(coord)

        # --- Parse input.nml for species labels and scaling factor ---
        # The &atoms block should provide the species labels.
        atoms_data = inp_nml.get('atoms', {})
        species_labels = []
        i = 1
        while f"label({i})" in atoms_data:
            raw_label = atoms_data[f"label({i})"]
            # Remove any quotes and keep only alphabetic characters.
            clean_label = ''.join(ch for ch in raw_label.strip("'\"") if ch.isalpha())
            species_labels.append(clean_label)
            i += 1
        if not species_labels:
            raise ValueError("No species labels found in the &atoms block of input.nml.")

        # Read the scaling factor 'alat' from the &lattice block in input.nml.
        inp_lattice = inp_nml.get('lattice', {})
        alat = inp_lattice.get('alat', 1.0)  # default scaling factor is 1.0 if not provided

        # --- Scale the lattice and basis coordinates ---
        unit_cell *= alat
        basis_coords = [np.array(coord) * alat for coord in basis_coords]

        # Check the number of species. If there is only one label and more than one atom,
        # apply the single label to all atoms.
        if len(species_labels) == 1 and n_atoms > 1:
            species_labels = species_labels * n_atoms
        elif len(species_labels) != n_atoms:
            raise ValueError(f"Mismatch: {n_atoms} basis coordinates but {len(species_labels)} species labels.")

        # --- Create Atom objects ---
        atoms = []
        for coord, sp in zip(basis_coords, species_labels):
            atoms.append(Atom(position=coord, species=sp, atom_type=sp))

        # Return a Structure with the list of Atom objects and the scaled unit cell.
        return Structure(atoms=atoms, unit_cell=unit_cell)
    
    @staticmethod
    def read_rs_clust(clust_file, input_file):
        """
        Reads a cluster file and an input file to create a Structure object.

        Parameters:
            clust_file (str): Path to the cluster file.
            input_file (str): Path to the input.nml file containing the &atoms block.

        Returns:
            Structure: A Structure object built from the parsed cluster file and input file.
        """
        # Example of clust_file: first line gives number of atoms
        # Followin lines contain atom positions and species for two atoms per row.
        #   II =   8940
        #  0.00000000    0.00000000    0.00000000   1   1    0.50000000    0.50000000    0.50000000   2   2
        # -10.00000000    0.00000000    0.00000000   1   1   -9.00000000   -4.00000000   -1.00000000   1   1
        #  -9.00000000   -4.00000000    0.00000000   1   1   -9.00000000   -4.00000000    1.00000000   1   1
        # Parse input.nml for species labels
        inp_nml = f90nml.read(input_file)
        atoms_data = inp_nml.get('atoms', {})
        # Remove hyphens and digits from species labels
        species_lookup = [''.join(filter(str.isalpha, label)) for label in atoms_data['label']]
        # Create Atom objects with species inferred from species_type
        atoms = []
        with open(clust_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        _ = int(lines[0].strip().split('=')[1].strip())  # num_atoms is not used
        for line in lines[1:]:
            parts = line.split()
            for i in range(0, len(parts), 10):
                position1 = [float(parts[i]), float(parts[i+1]), float(parts[i+2])]
                species_type1 = int(parts[i+3])
                species1 = species_lookup[species_type1 - 1]
                position2 = [float(parts[i+5]), float(parts[i+6]), float(parts[i+7])]
                species_type2 = int(parts[i+8])
                species2 = species_lookup[species_type2 - 1]
                atoms.append(Atom(position=position1, species=species1, atom_type=species1))
                atoms.append(Atom(position=position2, species=species2, atom_type=species2))
        
        # If the number of atoms is odd, remove the last element from atoms
        if len(atoms) % 2 != 0:
            atoms.pop()
        # Determine the largest absolute distance in positions to create a cubic unit cell
        max_distance = max(np.linalg.norm(atom.position) for atom in atoms)
        unit_cell = np.eye(3) * (2 * max_distance)
        # Create and return the Structure
        return Structure(atoms=atoms, unit_cell=None)