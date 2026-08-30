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
import json

from .pairs import (
    ExchangePair,
    PairKey,
    PairMappingError,
    legacy_displacement_mapping,
    legacy_displacement_mappings,
    map_cluster_neighbor,
    map_cluster_neighbors,
    pair_key_from_cluster_neighbor,
    reconstruct_cartesian_displacement,
)


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
            position1: list[float] = [
                float(parts[i]),
                float(parts[i + 1]),
                float(parts[i + 2]),
            ]
            species_type1: int = int(parts[i + 3])
            species1: str = species_lookup[species_type1 - 1]
            position2: list[float] = [
                float(parts[i + 5]),
                float(parts[i + 6]),
                float(parts[i + 7]),
            ]
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


def read_namelists(
    input_file: str, lattice_file: str, verbose: bool = False
) -> tuple[dict, dict]:
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
    lattice_nml: dict, input_nml: dict, cartesian: bool = True, soc: bool = False
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
    soc : bool, optional
        If True, account for spin-orbit coupling in symmetry considerations.
        SOC breaks time-reversal symmetry, reducing the spacegroup. Default False.

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

    # With SOC, we need to remove the horizontal mirror plane (z-inversion)
    # If a symmetry operation exists with [xrow; yrow; 0 0 z], remove its counterpart
    # with [xrow; yrow; 0 0 -z] (z-mirror operation)
    if soc:
        symm_ops = spacegroup_analyzer.get_symmetry_operations()
        filtered_ops = []
        seen_xy_planes = set()

        # First pass: identify which xy-plane operations we have
        for op in symm_ops:
            rot = op.rotation_matrix
            # Check if this is a z-preserving operation (z-component unchanged)
            # i.e., [a b 0; c d 0; 0 0 ±1]
            if np.allclose(rot[0:2, 2], 0) and np.allclose(rot[2, 0:2], 0):
                # Get the xy-plane part and z-component
                xy_part = tuple(np.round(rot[0:2, 0:2], 5).flatten())
                z_component = np.sign(rot[2, 2])
                key = (xy_part, z_component)
                seen_xy_planes.add(key)

        # Second pass: keep operations, but remove z-mirrors where z-preserving exists
        for op in symm_ops:
            rot = op.rotation_matrix
            # Check if this is a z-mirroring operation
            if np.allclose(rot[0:2, 2], 0) and np.allclose(rot[2, 0:2], 0):
                xy_part = tuple(np.round(rot[0:2, 0:2], 5).flatten())
                z_component = np.sign(rot[2, 2])

                # Keep this operation if:
                # - z_component is +1 (z-preserving), OR
                # - z_component is -1 (z-mirror) AND there's no +1 counterpart
                if z_component > 0 or (xy_part, 1) not in seen_xy_planes:
                    filtered_ops.append(op)
            else:
                # Non-z-operations: keep all
                filtered_ops.append(op)

        spacegroup_analyzer._symmetry_operations = filtered_ops

    rec_sites = [int(site) for site in lattice_nml["lattice"]["irec"]]
    return structure, spacegroup_analyzer, rec_sites


def get_neighbors_list(
    cluster: Molecule, sites: list[int], cutoff: float = 2.0, verbose: bool = False
) -> dict[int, list]:
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
        shell_width: float = 0.1, verbose: bool = False, spacegroup_analyzer: SpacegroupAnalyzer = None,
        structure: Structure = None) -> dict[int, list]:
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
    
    def are_vectors_symmetry_equivalent(v1_cart: np.ndarray, v2_cart: np.ndarray,
                                       symm_ops, tol: float = 1e-4) -> bool:
        """Check if two Cartesian vectors are related by a symmetry operation.

        Symmetry operations from pymatgen are given in fractional coordinates,
        so convert Cartesian vectors to fractional using the provided
        `structure` lattice before applying rotations.
        """
        if structure is None:
            # Fallback: compare Cartesian vectors directly
            for symm_op in symm_ops:
                rotated = symm_op.rotation_matrix @ v1_cart
                if np.allclose(rotated, v2_cart, atol=tol):
                    return True
            return False

        # Convert absolute Cartesian coordinates to fractional relative vectors
        # We convert via absolute coords because pymatgen lattice methods expect
        # absolute positions; use an arbitrary origin shift that cancels when
        # subtracting.
        # Use center at origin: get fractional coords of the endpoints then subtract.
        # Here v1_cart and v2_cart are difference vectors; create dummy points
        # and convert by adding to origin.
        frac_v1 = structure.lattice.get_fractional_coords([v1_cart])[0]
        frac_v2 = structure.lattice.get_fractional_coords([v2_cart])[0]

        for symm_op in symm_ops:
            rotated = symm_op.rotation_matrix @ frac_v1
            # account for periodic images by checking equivalence modulo 1
            diff = np.mod(rotated - frac_v2 + 0.5, 1.0) - 0.5
            if np.allclose(diff, 0.0, atol=tol):
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
        symm_ops = (
            spacegroup_analyzer.get_symmetry_operations() if spacegroup_analyzer else []
        )

        for shell_distance in sorted(shell_neighbors.keys()):
            # Track seen combinations per shell to allow distinct shells to retain
            # their own representatives while deduplicating within the shell.
            shell_seen_combinations = []
            for neighbor in shell_neighbors[shell_distance]:
                vector = neighbor.coords - center.coords
                species = neighbor.species_string

                # Check if this (species, vector) is symmetry-equivalent to any seen combination
                is_new = True
                if spacegroup_analyzer is not None:
                    for seen_species, seen_vector in shell_seen_combinations:
                        if seen_species == species:
                            if are_vectors_symmetry_equivalent(
                                vector, seen_vector, symm_ops
                            ):
                                is_new = False
                                break
                else:
                    # Fallback without symmetry: exact vector match
                    for seen_species, seen_vector in shell_seen_combinations:
                        if seen_species == species and np.allclose(vector, seen_vector, atol=1e-5):
                            is_new = False
                            break

                if is_new:
                    shell_seen_combinations.append((species, vector))
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


def get_Rij_mappings(cluster: Molecule, neighbors: dict, one_based: bool = True) -> list:
    """
    Construct a flat list of R_ij records for each center atom from neighbor lists.

    Each record has fields: i, j, R_ij (list), distance, species_i, species_j.
    ``R_ij`` is a Cartesian displacement for legacy compatibility, and is
    explicitly labeled with ``R_ij_kind``. It is not a periodic lattice
    translation; canonical periodic identities are provided by
    :mod:`rsutil.exchange.pairs`.

    Parameters
    ----------
    cluster : Molecule
        Cluster/molecule containing sites.
    neighbors : dict
        Output from `get_neighbors_list` or `get_unique_neighbors_per_shell`.
    one_based : bool, optional
        If True, return atom indices in 1-based Fortran-style indexing (default True).

    Returns
    -------
    list
        List of records, each a dict: {"i": int, "j": int, "R_ij": [x,y,z], "distance": float, "species_i": str, "species_j": str}
    """
    records: list = []
    for i_atom, neighbor_list in neighbors.items():
        center = cluster[i_atom]
        species_i = getattr(center, "species_string", str(center.species))
        i_idx = (i_atom + 1) if one_based else i_atom
        for neighbor in neighbor_list:
            j_idx = (neighbor.index + 1) if one_based else neighbor.index
            R = (neighbor.coords - center.coords).tolist()
            rec = {
                "i": int(i_idx),
                "j": int(j_idx),
                "R_ij": [float(x) for x in R],
                "R_ij_kind": "cartesian_displacement",
                "distance": float(neighbor.nn_distance),
                "species_i": str(species_i),
                "species_j": str(neighbor.species_string),
            }
            records.append(rec)
    return records


def write_Rij_mappings(mappings: dict, filename: str, fmt: str = "json") -> None:
    """
    Write R_ij mappings to disk in supported formats. Currently supports JSON.

    Parameters
    ----------
    mappings : dict
        Output from `get_Rij_mappings`.
    filename : str
        Path to write the mappings to.
    fmt : str, optional
        Output format ('json' supported), by default 'json'.
    """
    fmt_l = fmt.lower()
    # Normalize mappings input: accept dict mapping->list or flat list of records
    if isinstance(mappings, dict):
        # convert dict-of-lists to flat list
        records = []
        for key in sorted(mappings.keys()):
            for entry in mappings[key]:
                rec = {
                    "i": int(key),
                    "j": int(entry.get("j")),
                    "R_ij": [float(x) for x in entry.get("R_ij", [])],
                    "distance": float(entry.get("distance", 0.0)),
                    "species_i": str(entry.get("species_i", "")),
                    "species_j": str(entry.get("species", "")),
                }
                records.append(rec)
    else:
        records = list(mappings)

    if fmt_l == "json":
        with open(filename, "w", encoding="utf-8") as fh:
            json.dump(records, fh, indent=2)
    elif fmt_l == "txt":
        # Human readable plain-text format
        with open(filename, "w", encoding="utf-8") as fh:
            # fixed-width columns: i(4) j(4) Rx(12.8f) Ry(12.8f) Rz(12.8f) dist(12.8f) species_i(10) species_j(10)
            fh.write("# R_ij is Cartesian displacement; canonical identity is a PairKey\n")
            fh.write("# i j R_x R_y R_z distance species_i species_j\n")
            fmt_line = "{:4d} {:4d} {:12.8f} {:12.8f} {:12.8f} {:12.8f} {:10s} {:10s}\n"
            for rec in records:
                i = int(rec.get("i"))
                j = int(rec.get("j"))
                rx, ry, rz = rec.get("R_ij", [0.0, 0.0, 0.0])
                dist = float(rec.get("distance", 0.0))
                si = str(rec.get("species_i", ""))[:10]
                sj = str(rec.get("species_j", ""))[:10]
                fh.write(fmt_line.format(i, j, rx, ry, rz, dist, si, sj))
    else:
        raise ValueError(f"Unsupported format: {fmt}")


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


def store_jij_namelist(
    input_nml: dict, matching_indices: dict, atom_index: int
) -> dict:
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
                  shell_width: float = 0.1, soc: bool = False,
                  rij_out: str | None = None, rij_fmt: str = "json") -> None:
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
    soc : bool, optional
        If True, account for spin-orbit coupling in symmetry considerations.
        SOC breaks time-reversal symmetry, reducing effective symmetry. Default False.
    """
    input_nml, lattice_nml = read_namelists("input.nml", "lattice.nml", verbose=verbose)
    structure, spacegroup_analyzer, rec_sites = create_structure_and_analyzer(
        lattice_nml, input_nml, soc=soc
    )

    if verbose:
        print("Recursion sites:", rec_sites)
        print(
            f"Detected spacegroup: {spacegroup_analyzer.get_space_group_symbol()} "
            f"(#{spacegroup_analyzer.get_space_group_number()})"
        )
        symm_ops = spacegroup_analyzer.get_symmetry_operations()
        n_ops = len(symm_ops) if hasattr(symm_ops, "__len__") else len(list(symm_ops))
        print(f"Number of symmetry operations: {n_ops}")
        if soc:
            print("Spin-orbit coupling enabled: only proper rotations included")

    cluster = read_rs_clust("clust", "input.nml")

    if unique:
        neighbors = get_unique_neighbors_per_shell(
            cluster, rec_sites, cutoff=cutoff, shell_width=shell_width,
            verbose=verbose, spacegroup_analyzer=spacegroup_analyzer, structure=structure
        )
    else:
        neighbors = get_neighbors_list(
            cluster, rec_sites, cutoff=cutoff, verbose=verbose
        )

    if verbose:
        print(neighbors)

    matching_indices = find_matching_indices(cluster, neighbors)

    if verbose:
        print("Matching Indices:", matching_indices)

    for atom in rec_sites:
        store_jij_namelist(input_nml, matching_indices, atom)

    # Optionally write R_ij mappings for the neighbor network
    if rij_out:
        try:
            mappings = get_Rij_mappings(cluster, neighbors, one_based=True)
            write_Rij_mappings(mappings, rij_out, fmt=rij_fmt)
            if verbose:
                print(f"Wrote R_ij mappings to {rij_out} (format={rij_fmt})")
        except Exception as exc:
            raise RuntimeError(
                f"Failed to write R_ij mappings to {rij_out}: {exc}"
            ) from exc


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate exchange interaction mapping from cluster files."
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.0,
        help="Cutoff distance to consider neighbors (default: 1.0)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument(
        "--unique",
        action="store_true",
        help="Find only unique neighbors per shell (one representative per unique species-vector combination)",
    )
    parser.add_argument(
        "--shell-width",
        type=float,
        default=0.1,
        help="Width of shells for grouping neighbors when using --unique (default: 0.1)",
    )
    parser.add_argument(
        "--soc",
        action="store_true",
        default=False,
        help="Account for spin-orbit coupling in symmetry considerations (breaks time-reversal symmetry)",
    )

    args = parser.parse_args()

    exchange_map(
        cutoff=args.cutoff,
        verbose=args.verbose,
        unique=args.unique,
        shell_width=args.shell_width,
        soc=args.soc,
    )
