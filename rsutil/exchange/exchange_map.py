"""
This module provides functions to read and process cluster and namelist files,
create structures and analyzers, and determine neighbors and matching indices
for atoms in a cluster. It also includes functionality to update and store
input namelists with post-processing information.

Functions
---------
read_rs_clust(clust_file, input_file)

read_namelists(input_file, lattice_file, verbose=False)
    Reads input and lattice namelist files.

create_structure_and_analyzer(lattice_nml, input_nml)
    Creates a Structure object and SpacegroupAnalyzer from lattice and input
    NML files.

get_neighbors_list(cluster, sites, cutoff=2.0, verbose=False)
    Gets the neighbors list for each atom in the cluster.

find_matching_indices(cluster, neighbors)
    Finds the indices of the atoms in the cluster that match the neighbor list.

store_jij_namelist(input_nml, matching_indices, atom_index)
    Updates the input namelist with post-processing information and matching
    indices.

exchange_map(verbose=False)
    Main function to read namelists, create structure and analyzer, get
    neighbors, and symmetry reduced neighbors.
"""
import os
import shutil

import f90nml
import numpy as np
from pymatgen.core import Lattice, Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def read_rs_clust(clust_file, input_file):
    """
    Reads a cluster file and an input file to create a Structure object.

    Parameters:
        clust_file (str): Path to the cluster file.
        input_file (str): Path to the input.nml file containing the &atoms block.

    Returns:
        Structure: A Structure object built from the parsed cluster file and input file.
    """
    # Parse input.nml for species labels
    inp_nml = f90nml.read(input_file)
    atoms_data = inp_nml.get("atoms", {})
    # Remove hyphens and digits from species labels
    species_lookup = [
        "".join(filter(str.isalpha, label)) for label in atoms_data["label"]
    ]
    # Create Atom objects with species inferred from species_type
    atoms = []
    with open(clust_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
    _ = int(lines[0].strip().split("=")[1].strip())  # num_atoms is not used
    positions = []
    atoms = []
    for line in lines[1:]:
        parts = line.split()
        for i in range(0, len(parts), 10):
            position1 = [float(parts[i]), float(parts[i + 1]), float(parts[i + 2])]
            species_type1 = int(parts[i + 3])
            species1 = species_lookup[species_type1 - 1]
            position2 = [float(parts[i + 5]), float(parts[i + 6]), float(parts[i + 7])]
            species_type2 = int(parts[i + 8])
            species2 = species_lookup[species_type2 - 1]
            positions.append(position1)
            positions.append(position2)
            atoms.append(species1)
            atoms.append(species2)

    # If the number of atoms is odd, remove the last element from atoms
    if len(atoms) % 2 != 0:
        atoms.pop()
        # Create a pymatgen Molecule object

    cluster = Molecule(atoms, positions)
    return cluster


def read_namelists(input_file: str, lattice_file: str, verbose=False) -> tuple:
    """
    Read input and lattice namelist files.

    Parameters
    ----------
    input_file : str
        Path to the input namelist file.
    lattice_file : str
        Path to the lattice namelist file.
    verbose : bool, optional
        Flag to control printouts, by default False.

    Returns
    -------
    tuple
        A tuple containing the parsed input and lattice namelists.
    """
    input_nml = f90nml.read(input_file)
    lattice_nml = f90nml.read(lattice_file)

    if verbose:
        print("Input Namelist:", input_nml)
        print("Lattice Namelist:", lattice_nml)

    return input_nml, lattice_nml


def create_structure_and_analyzer(lattice_nml, input_nml):
    """
    Create a Structure object and SpacegroupAnalyzer from lattice and input NML files.

    Parameters
    ----------
    lattice_nml : dict
        Parsed lattice NML file.
    input_nml : dict
        Parsed input NML file.

    Returns
    -------
    tuple
        A tuple containing the Structure object and SpacegroupAnalyzer.
    """
    lattice = np.array(lattice_nml["lattice"]["a"])
    basis = np.array(lattice_nml["lattice"]["crd"])
    atoms = input_nml["atoms"]["label"]

    lattice_matrix = Lattice(lattice)
    structure = Structure(lattice_matrix, atoms, basis)
    spacegroup_analyzer = SpacegroupAnalyzer(structure)

    rec_sites = [int(site) for site in lattice_nml["lattice"]["irec"]]
    return structure, spacegroup_analyzer, rec_sites


def get_neighbors_list(cluster, sites, cutoff=2.0, verbose=False):
    """
    Get the neighbors list for each atom in the cluster.

    Parameters
    ----------
    cluster : Molecule
        The cluster object containing atomic sites.
    sites : list
        List of atom labels.
    cutoff : float, optional
        Cutoff distance to consider neighbors, by default 2.0.
    verbose : bool, optional
        Flag to control printouts, by default False.

    Returns
    -------
    dict
        Dictionary of neighbors for each atom, sorted by distance.
    """
    neighbors = {}
    for iatom, atom in enumerate(sites):
        neighbors[iatom] = cluster.get_neighbors(cluster[iatom], cutoff)
        neighbors[iatom] = sorted(
            neighbors[iatom], key=lambda neighbor: neighbor.nn_distance
        )

    if verbose:
        print("Neighbors:")
        for iatom, neighbor_list in neighbors.items():
            print(f"Neighbors for {iatom}:")
            for neighbor in neighbor_list:
                print(
                    f"  {neighbor.species_string} at {neighbor.coords} "
                    f"with distance {neighbor.nn_distance}"
                )

    return neighbors


# Not working yet
# def get_symmetry_reduced_neighbors(neighbors, spacegroup_analyzer, verbose=False):
#     """
#     Get the symmetry reduced neighbors list for each atom in the cluster.
#
#     Parameters
#     ----------
#     neighbors : dict
#         Dictionary of neighbors for each atom.
#     spacegroup_analyzer : SpacegroupAnalyzer
#         The spacegroup analyzer object for the structure.
#     verbose : bool, optional
#         Flag to control printouts, by default False.
#
#     Returns
#     -------
#     dict
#         Dictionary of symmetry reduced neighbors for each atom.
#     """
#     symmetry_reduced_neighbors = {}
#     symm_ops = spacegroup_analyzer.get_symmetry_operations()
#
#     for atom, neighbor_list in neighbors.items():
#         reduced_neighbors = {}
#         for neighbor in neighbor_list:
#             for symm_op in symm_ops:
#                 symm_neighbor = symm_op.operate(neighbor.coords)
#                 distance = np.linalg.norm(symm_neighbor - cluster[atom].coords)
#                 distance = round(distance, 5)  # Round to avoid floating point issues
#                 symm_key = (distance, tuple(np.round(symm_neighbor, 5)))
#                 if symm_key not in reduced_neighbors:
#                     reduced_neighbors[symm_key] = symm_neighbor
#
#         symmetry_reduced_neighbors[atom] = [coords for _, coords in reduced_neighbors.values()]
#
#     if verbose:
#         print("Symmetry Reduced Neighbors:")
#         for atom, neighbor_list in symmetry_reduced_neighbors.items():
#             print(f"Symmetry reduced neighbors for {atom}:")
#             for neighbor in neighbor_list:
#                 print(f"  {neighbor}")
#
#     return symmetry_reduced_neighbors


def find_matching_indices(cluster, neighbors):
    """
    Find the indices of the atoms in the cluster that match the neighbor list.

    Parameters
    ----------
    cluster : Molecule
        The cluster object containing atomic sites.
    neighbors : dict
        Dictionary of neighbors for each atom.

    Returns
    -------
    dict
        Dictionary of atom indices that match the neighbor list.
    """
    matching_indices = {}
    for atom, neighbor_list in neighbors.items():
        atom_indices = []
        for neighbor in neighbor_list:
            for i, site in enumerate(cluster):
                if site.species_string == neighbor.species_string and np.allclose(
                    site.coords, neighbor.coords, atol=1e-5
                ):
                    atom_indices.append(i + 1)
                    break
        matching_indices[atom] = atom_indices

    return matching_indices


def store_jij_namelist(input_nml, matching_indices, atom_index):
    """
    Update the input namelist with post-processing information and matching indices.

    Parameters
    ----------
    input_nml : dict
        Parsed input namelist.
    matching_indices : dict
        Dictionary of atom indices that match the neighbor list.
    atom_index : int
        Index of the atom in rec_sites.

    Returns
    -------
    dict
        Updated input namelist.
    """
    new_input_nml = input_nml.copy()

    # Edit and write input.nml
    new_input_nml["calculation"]["post_processing"] = "exchange"
    try:
        del new_input_nml["calculation"]["pre_processing"]
    except KeyError:
        pass

    # Update njij and ijpair
    njij = len(matching_indices[atom_index - 1])
    new_input_nml["lattice"]["njij"] = njij
    new_input_nml["lattice"]["ijpair"] = matching_indices[atom_index - 1]

    iij = [atom_index] * njij
    jij = [i for i in matching_indices[atom_index - 1]]
    ijpair = [iij, jij]

    new_input_nml["lattice"]["ijpair"] = ijpair

    # Create directory Atom_i
    dir_name = f"Atom_{atom_index}"
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    # Copy lattice.nml and converged atom files to the new directory
    shutil.copy("lattice.nml", os.path.join(dir_name, "lattice.nml"))
    for file in os.listdir():
        if file.endswith("_out.nml"):
            shutil.copy(file, os.path.join(dir_name, file))

    # Save new_input_nml to input.nml in the new directory
    f90nml.write(new_input_nml, os.path.join(dir_name, "input.nml"), force=True)

    print(f"Updated input_nml saved to {os.path.join(dir_name, 'input.nml')}")
    print(f"lattice.nml copied to {os.path.join(dir_name, 'lattice.nml')}")

    return new_input_nml


# Create lists for input and lattice NML
# input_nml_list = [np.array(input_nml[key]) for key in input_nml]
# lattice_nml_list = [np.array(lattice_nml[key]) for key in lattice_nml]
# atoms_list = [atom for atom in atoms]
#
def exchange_map(cutoff=1.0, verbose=False):
    """
    Main function to read namelists, create structure and analyzer,
    get neighbors, and symmetry reduced neighbors.

    Parameters
    ----------
    verbose : bool, optional
        Flag to control printouts, by default False.
    """
    input_nml, lattice_nml = read_namelists("input.nml", "lattice.nml", verbose=verbose)
    _, _, rec_sites = create_structure_and_analyzer(lattice_nml, input_nml)
    # structure, spacegroup_analyzer, rec_sites = (
    # create_structure_and_analyzer(lattice_nml, input_nml)
    # )

    if verbose:
        print("Recursion sites:", rec_sites)

    cluster = read_rs_clust("clust", "input.nml")

    neighbors = get_neighbors_list(cluster, rec_sites, cutoff=cutoff, verbose=verbose)

    if verbose:
        print(neighbors)

    matching_indices = find_matching_indices(cluster, neighbors)

    if verbose:
        print("Matching Indices:", matching_indices)

    for atom in rec_sites:
        store_jij_namelist(input_nml, matching_indices, atom)


if __name__ == "__main__":
    exchange_map()
