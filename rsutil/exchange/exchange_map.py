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


def read_rs_clust(clust_file: str, input_file: str) -> Molecule:
    """
    Reads a cluster file and an input file to create a Structure object.

    Parameters:
        clust_file (str): Path to the cluster file.
        input_file (str): Path to the input.nml file containing the &atoms block.

    Returns:
        Molecule: A Molecule object built from the parsed cluster file and input file.
    """
    # Parse input.nml for species labels
    inp_nml: dict = f90nml.read(input_file)
    atoms_data: dict = inp_nml.get("atoms", {})
    # Remove hyphens and digits from species labels
    species_lookup: list[str] = [
        "".join(filter(str.isalpha, label)) for label in atoms_data["label"]
    ]
    # Create Atom objects with species inferred from species_type
    atoms: list[str] = []
    with open(clust_file, "r", encoding="utf-8") as f:
        lines: list[str] = f.readlines()
    _ = int(lines[0].strip().split("=")[1].strip())  # num_atoms is not used
    positions: list[list[float]] = []
    for line in lines[1:]:
        parts: list[str] = line.split()
        for i in range(0, len(parts), 10):
            position1: list[float] = [float(parts[i]), float(parts[i + 1]), float(parts[i + 2])]
            species_type1: int = int(parts[i + 3])
            species1: str = species_lookup[species_type1 - 1]
            position2: list[float] = [float(parts[i + 5]), float(parts[i + 6]), float(parts[i + 7])]
            species_type2: int = int(parts[i + 8])
            species2: str = species_lookup[species_type2 - 1]
            positions.append(position1)
            positions.append(position2)
            atoms.append(species1)
            atoms.append(species2)

    # If the number of atoms is odd, remove the last element from atoms
    if len(atoms) % 2 != 0:
        atoms.pop()

    # Create a pymatgen Molecule object
    cluster: Molecule = Molecule(atoms, positions)
    return cluster


def read_namelists(input_file: str, lattice_file: str, verbose: bool = False) -> tuple[dict, dict]:
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
    input_nml: dict = f90nml.read(input_file)
    lattice_nml: dict = f90nml.read(lattice_file)

    if verbose:
        print("Input Namelist:", input_nml)
        print("Lattice Namelist:", lattice_nml)

    return input_nml, lattice_nml


def create_structure_and_analyzer(
    lattice_nml: dict, input_nml: dict, cartesian: bool = True
) -> tuple[Structure, SpacegroupAnalyzer, list[int]]:
    """
    Create a Structure object and SpacegroupAnalyzer from lattice and input NML files.

    Parameters
    ----------
    lattice_nml : dict
        Parsed lattice NML file.
    input_nml : dict
        Parsed input NML file.
    cartesian : bool, optional
        If True (default), treat basis vectors as Cartesian coordinates and convert
        to fractional. If False, treat as fractional coordinates.

    Returns
    -------
    tuple
        A tuple containing the Structure object and SpacegroupAnalyzer.
    """
    lattice = np.array(lattice_nml["lattice"]["a"])
    basis = np.array(lattice_nml["lattice"]["crd"])
    atoms = input_nml["atoms"]["label"]

    lattice_matrix = Lattice(lattice)
    
    if cartesian:
        # Convert Cartesian coordinates to fractional
        basis_fractional = lattice_matrix.get_fractional_coords(basis)
        structure = Structure(lattice_matrix, atoms, basis_fractional)
    else:
        # Use fractional coordinates as-is
        structure = Structure(lattice_matrix, atoms, basis)
    
    spacegroup_analyzer = SpacegroupAnalyzer(structure)

    rec_sites = [int(site) for site in lattice_nml["lattice"]["irec"]]
    return structure, spacegroup_analyzer, rec_sites


def get_neighbors_list(
    cluster: Molecule, sites: list[int], cutoff: float = 2.0,
        verbose: bool = False) -> dict[int, list]:
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
    neighbors: dict[int, list] = {}
    for iatom, _ in enumerate(sites):
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


def get_unique_neighbors_per_shell(
    cluster: Molecule, sites: list[int], cutoff: float = 2.0,
        shell_width: float = 0.1, verbose: bool = False, spacegroup_analyzer: SpacegroupAnalyzer = None) -> dict[int, list]:
    """
    Get unique neighbors per shell for each atom in the cluster.
    
    For each neighbor shell (defined by distance), this function returns only one 
    representative neighbor per unique (species, symmetry-equivalent vector) combination.
    Symmetry-equivalent vectors (e.g., [1,2,2] and [-1,-2,-2]) are grouped together.
    This is useful for identifying structurally distinct neighbor environments.

    Parameters
    ----------
    cluster : Molecule
        The cluster object containing atomic sites.
    sites : list
        List of atom labels.
    cutoff : float, optional
        Cutoff distance to consider neighbors, by default 2.0.
    shell_width : float, optional
        Width of each shell to define neighbor groups, by default 0.1.
    verbose : bool, optional
        Flag to control printouts, by default False.
    spacegroup_analyzer : SpacegroupAnalyzer, optional
        The spacegroup analyzer for applying symmetry operations, by default None.

    Returns
    -------
    dict
        Dictionary of unique neighbors per shell for each atom.
    """
    unique_neighbors: dict[int, list] = {}
    
    def are_vectors_symmetry_equivalent(v1: np.ndarray, v2: np.ndarray, 
                                       symm_ops, tolerance: float = 1e-4) -> bool:
        """Check if two vectors are related by a symmetry operation."""
        for symm_op in symm_ops:
            rotated = symm_op.rotation_matrix @ v1
            if np.allclose(rotated, v2, atol=tolerance):
                return True
        return False
    
    for iatom, _ in enumerate(sites):
        center = cluster[iatom]
        all_neighbors = cluster.get_neighbors(center, cutoff)
        
        # Group neighbors by distance shells
        shell_neighbors = {}
        for neighbor in all_neighbors:
            distance = neighbor.nn_distance
            # Find which shell this neighbor belongs to
            shell_idx = round(distance / shell_width)
            shell_distance = shell_idx * shell_width
            
            if shell_distance not in shell_neighbors:
                shell_neighbors[shell_distance] = []
            shell_neighbors[shell_distance].append(neighbor)
        
        # For each shell, keep only one representative per (species, symmetry-equivalent vector) combination
        unique_list = []
        seen_combinations = []  # List of (species, representative_vector) tuples
        symm_ops = spacegroup_analyzer.get_symmetry_operations() if spacegroup_analyzer else []
        
        for shell_distance in sorted(shell_neighbors.keys()):
            for neighbor in shell_neighbors[shell_distance]:
                vector = neighbor.coords - center.coords
                species = neighbor.species_string
                
                # Check if this (species, vector) is symmetry-equivalent to any seen combination
                is_new = True
                if spacegroup_analyzer is not None:
                    for seen_species, seen_vector in seen_combinations:
                        if seen_species == species:
                            if are_vectors_symmetry_equivalent(vector, seen_vector, symm_ops):
                                is_new = False
                                break
                else:
                    # Fallback without symmetry: exact vector match
                    for seen_species, seen_vector in seen_combinations:
                        if seen_species == species and np.allclose(vector, seen_vector, atol=1e-5):
                            is_new = False
                            break
                
                if is_new:
                    seen_combinations.append((species, vector))
                    unique_list.append(neighbor)
        
        # Sort by distance
        unique_neighbors[iatom] = sorted(
            unique_list, key=lambda neighbor: neighbor.nn_distance
        )
        
        if verbose:
            print(f"Unique neighbors for atom {iatom}:")
            for neighbor in unique_neighbors[iatom]:
                vector = neighbor.coords - center.coords
                print(
                    f"  {neighbor.species_string} at {neighbor.coords} "
                    f"(vector: {vector}), distance: {neighbor.nn_distance:.4f}"
                )
    
    return unique_neighbors


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


def get_fermi_energy(report_file: str = "report.out") -> float:
    """
    Extract the Fermi energy from the report file.

    Parameters
    ----------
    report_file : str, optional
        Path to the report file, by default "report.out".

    Returns
    -------
    float
        The Fermi energy value.
    """
    fermi_energy: float = None
    with open(report_file, "r", encoding="utf-8") as file:
        for line in file:
            if "Fermi energy:" in line:
                fermi_energy = float(line.split()[2])
                break

    if fermi_energy is None:
        raise ValueError("Fermi energy not found in the report file.")

    return fermi_energy


def find_matching_indices(cluster: Molecule, neighbors: dict) -> dict:
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
    matching_indices: dict[int, list[int]] = {}
    for atom, neighbor_list in neighbors.items():
        atom_indices: list[int] = []
        for neighbor in neighbor_list:
            atom_indices.append(neighbor.index + 1)
            # for i, site in enumerate(cluster):
            #     if site.species_string == neighbor.species_string and np.allclose(
            #         site.coords, neighbor.coords, atol=1e-5
            #     ):
            #         atom_indices.append(i + 1)
            #         break
        matching_indices[atom] = atom_indices

    return matching_indices


def store_jij_namelist(input_nml: dict, matching_indices: dict, atom_index: int) -> dict:
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
    new_input_nml: dict = input_nml.copy()

    # Update the fermi energy in the energy namelist
    fermi_energy: float = get_fermi_energy()
    new_input_nml["energy"]["fermi"] = fermi_energy

    # Edit and write input.nml
    new_input_nml["calculation"]["post_processing"] = "exchange"
    try:
        del new_input_nml["calculation"]["pre_processing"]
    except KeyError:
        pass

    # Store all lattice entries and rebuild to control order
    lattice_section = new_input_nml["lattice"]
    ct_value = lattice_section.get("ct", None)
    r2_value = lattice_section.get("r2", None)

    # Build a new lattice section with controlled ordering
    new_lattice = f90nml.Namelist()
    
    # Copy all existing entries except njij, ijpair, ct, and r2
    for key, value in lattice_section.items():
        if key not in ["njij", "ijpair", "ct", "r2"]:
            new_lattice[key] = value
    
    # Add njij and ijpair in the desired position
    njij: int = len(matching_indices[atom_index - 1])
    new_lattice["njij"] = njij
    
    iij: list[int] = [atom_index] * njij
    jij: list[int] = [i for i in matching_indices[atom_index - 1]]
    ijpair: list[list[int]] = [iij, jij]
    new_lattice["ijpair"] = ijpair
    
    # Add ct and r2 after njij and ijpair
    if ct_value is not None:
        new_lattice["ct"] = ct_value
    if r2_value is not None:
        new_lattice["r2"] = r2_value
    
    # Replace the lattice section
    new_input_nml["lattice"] = new_lattice

    # Create directory Atom_i
    dir_name: str = f"Atom_{atom_index}"
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
def exchange_map(cutoff: float = 1.0, verbose: bool = False, unique: bool = False, 
                  shell_width: float = 0.1) -> None:
    """
    Main function to read namelists, create structure and analyzer,
    get neighbors, and symmetry reduced neighbors.

    Parameters
    ----------
    cutoff : float, optional
        Cutoff distance to consider neighbors, by default 1.0.
    verbose : bool, optional
        Flag to control printouts, by default False.
    unique : bool, optional
        If True, find only unique neighbors per shell (one representative per 
        unique species-vector combination). Symmetry-equivalent vectors are grouped.
        by default False.
    shell_width : float, optional
        Width of shells for grouping neighbors (used with unique=True), 
        by default 0.1.
    """
    input_nml, lattice_nml = read_namelists("input.nml", "lattice.nml", verbose=verbose)
    structure, spacegroup_analyzer, rec_sites = create_structure_and_analyzer(lattice_nml, input_nml)

    if verbose:
        print("Recursion sites:", rec_sites)
        print(f"Detected spacegroup: {spacegroup_analyzer.get_space_group_symbol()} "
              f"(#{spacegroup_analyzer.get_space_group_number()})")
        print(f"Number of symmetry operations: {len(spacegroup_analyzer.get_symmetry_operations())}")

    cluster = read_rs_clust("clust", "input.nml")

    if unique:
        neighbors = get_unique_neighbors_per_shell(
            cluster, rec_sites, cutoff=cutoff, shell_width=shell_width, 
            verbose=verbose, spacegroup_analyzer=spacegroup_analyzer
        )
    else:
        neighbors = get_neighbors_list(cluster, rec_sites, cutoff=cutoff, verbose=verbose)

    if verbose:
        print(neighbors)

    matching_indices = find_matching_indices(cluster, neighbors)

    if verbose:
        print("Matching Indices:", matching_indices)

    for atom in rec_sites:
        store_jij_namelist(input_nml, matching_indices, atom)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate exchange interaction mapping from cluster files."
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.0,
        help="Cutoff distance to consider neighbors (default: 1.0)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "--unique",
        action="store_true",
        help="Find only unique neighbors per shell (one representative per unique species-vector combination)"
    )
    parser.add_argument(
        "--shell-width",
        type=float,
        default=0.1,
        help="Width of shells for grouping neighbors when using --unique (default: 0.1)"
    )
    
    args = parser.parse_args()
    
    exchange_map(
        cutoff=args.cutoff,
        verbose=args.verbose,
        unique=args.unique,
        shell_width=args.shell_width
    )
