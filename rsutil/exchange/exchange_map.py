"""Legacy exchange helpers backed by the path-explicit exchange service."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import f90nml
import numpy as np
from pymatgen.core import Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .fermi import parse_fermi_energy
from .jobs import build_exchange_input
from .neighbors import discover_default_neighbors, discover_legacy_unique_neighbors
from .reference import build_structure_and_analyzer, read_namelists as _read_namelists, read_rs_clust as _read_rs_clust
from .service import ExchangeService
from .writer import write_Rij_mappings


def _as_list(value: Any) -> list[Any]:
    return list(value) if isinstance(value, (list, tuple, np.ndarray)) else [value]


def read_rs_clust(clust_file: str = "clust", input_file: str = "input.nml") -> Molecule:
    """Compatibility wrapper for explicit cluster loading."""

    return _read_rs_clust(clust_file, input_file)


def read_namelists(
    input_file: str = "input.nml", lattice_file: str = "lattice.nml", verbose: bool = False
) -> tuple[dict, dict]:
    """Compatibility wrapper for explicit namelist loading."""

    return _read_namelists(input_file, lattice_file, verbose=verbose)


def create_structure_and_analyzer(
    lattice_nml: dict,
    input_nml: dict,
    cartesian: bool = True,
    soc: bool = False,
) -> tuple[Structure, SpacegroupAnalyzer, list[int]]:
    """Build the legacy structure/analyzer tuple.

    ``soc`` remains accepted for compatibility; SOC reduction is now applied
    by the isolated neighbor-reduction hook rather than private analyzer state.
    """

    structure, analyzer = build_structure_and_analyzer(
        lattice_nml, input_nml, cartesian=cartesian
    )
    rec_sites = [int(value) for value in _as_list(lattice_nml["lattice"]["irec"])]
    return structure, analyzer, rec_sites


def get_neighbors_list(
    cluster: Molecule, sites: list[int], cutoff: float = 2.0, verbose: bool = False
) -> dict[int, list]:
    """Compatibility wrapper for default non-unique neighbor discovery."""

    return discover_default_neighbors(cluster, sites, cutoff=cutoff, verbose=verbose)


def get_unique_neighbors_per_shell(
    cluster: Molecule,
    sites: list[int],
    cutoff: float = 2.0,
    shell_width: float = 0.1,
    verbose: bool = False,
    spacegroup_analyzer: SpacegroupAnalyzer | None = None,
    structure: Structure | None = None,
    magnetic_reference: str | None = None,
) -> dict[int, list]:
    """Compatibility wrapper for the pre-HT-09 unique-neighbor reducer."""

    return discover_legacy_unique_neighbors(
        cluster,
        sites,
        cutoff=cutoff,
        shell_width=shell_width,
        analyzer=spacegroup_analyzer,
        structure=structure,
        magnetic_reference=magnetic_reference,
        verbose=verbose,
    )


def get_Rij_mappings(cluster: Molecule, neighbors: dict, one_based: bool = True) -> list:
    """Return the historical Cartesian-displacement mapping shape."""

    records: list[dict[str, Any]] = []
    for center_index, neighbor_list in neighbors.items():
        center = cluster[center_index]
        i_index = center_index + 1 if one_based else center_index
        for neighbor in neighbor_list:
            j_index = neighbor.index + 1 if one_based else neighbor.index
            displacement = (neighbor.coords - center.coords).tolist()
            records.append(
                {
                    "i": int(i_index),
                    "j": int(j_index),
                    "R_ij": [float(value) for value in displacement],
                    "R_ij_kind": "cartesian_displacement",
                    "distance": float(neighbor.nn_distance),
                    "length_unit": "angstrom",
                    "species_i": str(center.species_string),
                    "species_j": str(neighbor.species_string),
                }
            )
    return records


def get_fermi_energy(report_file: str = "report.out") -> float:
    """Compatibility wrapper for structured Fermi-energy parsing."""

    return parse_fermi_energy(report_file)


def find_matching_indices(cluster: Molecule, neighbors: dict) -> dict:
    """Return legacy one-based cluster-neighbor indices."""

    return {
        int(center): [int(neighbor.index) + 1 for neighbor in neighbor_list]
        for center, neighbor_list in neighbors.items()
    }


def store_jij_namelist(
    input_nml: dict,
    matching_indices: dict,
    atom_index: int,
    *,
    report_file: str = "report.out",
    lattice_file: str = "lattice.nml",
    output_dir: str | Path | None = None,
    source_dir: str | Path = ".",
) -> dict:
    """Compatibility wrapper for writing one legacy exchange input."""

    root = Path(source_dir).expanduser().resolve()
    report_path = Path(report_file).expanduser()
    if not report_path.is_absolute():
        report_path = root / report_path
    lattice_path = Path(lattice_file).expanduser()
    if not lattice_path.is_absolute():
        lattice_path = root / lattice_path
    destination = Path(output_dir).expanduser() if output_dir is not None else root / f"Atom_{atom_index}"
    if not destination.is_absolute():
        destination = root / destination
    generated = build_exchange_input(
        input_nml,
        atom_index=atom_index,
        matching_indices=matching_indices[atom_index - 1],
        fermi_energy=parse_fermi_energy(report_path),
    )
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(lattice_path, destination / "lattice.nml")
    for output_file in sorted(root.glob("*_out.nml")):
        shutil.copy2(output_file, destination / output_file.name)
    f90nml.write(generated, str(destination / "input.nml"), force=True)
    print(f"Updated input_nml saved to {destination / 'input.nml'}")
    print(f"lattice.nml copied to {destination / 'lattice.nml'}")
    return generated


def exchange_map(
    cutoff: float = 1.0,
    verbose: bool = False,
    unique: bool = False,
    shell_width: float = 0.1,
    soc: bool | None = None,
    full_pairs: bool = False,
    magnetic_axis: tuple[float, float, float] | None = None,
    magnetic_reference: str | None = None,
    max_pairs_per_job: int = 0,
    rij_out: str | None = None,
    rij_fmt: str = "json",
    source_dir: str | Path = ".",
    collision_policy: str = "error",
    force_legacy_collision: bool = False,
) -> None:
    """Legacy entry point delegating to :class:`ExchangeService`."""

    ExchangeService(source_dir).generate(
        cutoff=cutoff,
        verbose=verbose,
        unique=unique,
        shell_width=shell_width,
        soc=soc,
        full_pairs=full_pairs,
        magnetic_axis=magnetic_axis,
        magnetic_reference=magnetic_reference,
        max_pairs_per_job=max_pairs_per_job,
        rij_out=rij_out,
        rij_fmt=rij_fmt,
        collision_policy=collision_policy,
        force_legacy_collision=force_legacy_collision,
    )


__all__ = [
    "create_structure_and_analyzer",
    "exchange_map",
    "find_matching_indices",
    "get_Rij_mappings",
    "get_fermi_energy",
    "get_neighbors_list",
    "get_unique_neighbors_per_shell",
    "read_namelists",
    "read_rs_clust",
    "store_jij_namelist",
    "write_Rij_mappings",
]
