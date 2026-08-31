"""Path-explicit loading of a converged RS-LMTO exchange calculation."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import f90nml
import numpy as np
from pymatgen.core import Lattice, Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from rsutil.core import RSStructure, StructureSite

from .fermi import parse_fermi_energy
from .symmetry import NspSemantics, nsp_semantics


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, (list, tuple, np.ndarray)):
        return list(value)
    return [value]


def _element(label: str) -> str:
    match = re.search(r"([A-Z][a-z]?)", str(label))
    if match is None:
        raise ValueError(f"could not determine chemical element from label {label!r}")
    return match.group(1)


def read_rs_clust(clust_path: str | Path, input_path: str | Path) -> Molecule:
    """Read an RS ``clust`` file without relying on process CWD."""

    input_nml = f90nml.read(str(input_path))
    labels = [_element(label) for label in _as_list(input_nml["atoms"]["label"])]
    atoms: list[str] = []
    positions: list[list[float]] = []
    lines = Path(clust_path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"cluster file is empty: {clust_path}")
    for line_number, line in enumerate(lines[1:], start=2):
        parts = line.split()
        if not parts:
            continue
        if len(parts) % 10:
            raise ValueError(
                f"cluster line {line_number} has {len(parts)} fields; expected a multiple of 10"
            )
        for offset in range(0, len(parts), 10):
            first = [float(parts[offset]), float(parts[offset + 1]), float(parts[offset + 2])]
            second = [
                float(parts[offset + 5]),
                float(parts[offset + 6]),
                float(parts[offset + 7]),
            ]
            first_type = int(parts[offset + 3])
            second_type = int(parts[offset + 8])
            try:
                atoms.extend((labels[first_type - 1], labels[second_type - 1]))
            except IndexError as exc:
                raise ValueError(f"invalid species type in cluster line {line_number}") from exc
            positions.extend((first, second))
    if len(atoms) % 2:
        atoms.pop()
        positions.pop()
    return Molecule(atoms, positions)


def read_namelists(
    input_path: str | Path, lattice_path: str | Path, *, verbose: bool = False
) -> tuple[dict, dict]:
    """Read the reference input and lattice namelists from explicit paths."""

    input_nml = f90nml.read(str(input_path))
    lattice_nml = f90nml.read(str(lattice_path))
    if verbose:
        print("Input Namelist:", input_nml)
        print("Lattice Namelist:", lattice_nml)
    return input_nml, lattice_nml


def build_structure_and_analyzer(
    lattice_nml: dict, input_nml: dict, *, cartesian: bool = True
) -> tuple[Structure, SpacegroupAnalyzer]:
    """Build the pymatgen structure from RS primitive sites and ``izp``.

    ``crd`` is indexed by primitive-cell site while ``atoms.label`` is
    indexed by RS electronic type.  ``izp`` is the source-faithful bridge
    between those two tables.
    """

    lattice = np.asarray(lattice_nml["lattice"]["a"], dtype=float)
    basis = np.asarray(lattice_nml["lattice"]["crd"], dtype=float)
    labels, type_indices = _rs_site_types(lattice_nml, input_nml)
    lattice_matrix = Lattice(lattice)
    if cartesian:
        basis = lattice_matrix.get_fractional_coords(basis)
    structure = Structure(
        lattice_matrix,
        [_element(labels[type_index - 1]) for type_index in type_indices],
        basis,
        site_properties={
            "rs_label": [labels[type_index - 1] for type_index in type_indices],
            "rs_type_index": type_indices,
        },
    )
    return structure, SpacegroupAnalyzer(structure)


def build_rs_structure(lattice_nml: dict, input_nml: dict) -> RSStructure:
    """Build the dependency-light canonical structure for pair mapping."""

    lattice = np.asarray(lattice_nml["lattice"]["a"], dtype=float)
    basis = np.asarray(lattice_nml["lattice"]["crd"], dtype=float)
    labels, type_indices = _rs_site_types(lattice_nml, input_nml)
    sites = tuple(
        StructureSite(
            index=index,
            position=tuple(float(value) for value in position),
            element=_element(labels[type_index - 1]),
            rs_label=labels[type_index - 1],
            rs_type_index=type_index,
        )
        for index, (type_index, position) in enumerate(zip(type_indices, basis, strict=True))
    )
    return RSStructure(
        lattice_vectors=tuple(tuple(float(value) for value in row) for row in lattice),
        length_unit="lattice",
        coordinate_convention="cartesian",
        sites=sites,
    )


def _rs_site_types(lattice_nml: dict, input_nml: dict) -> tuple[list[str], tuple[int, ...]]:
    """Return source RS labels and the one-based type for each basis site."""

    try:
        labels = [str(value) for value in _as_list(input_nml["atoms"]["label"])]
        type_values = tuple(int(value) for value in _as_list(lattice_nml["lattice"]["izp"]))
        basis = _as_list(lattice_nml["lattice"]["crd"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("RS reference requires atoms.label, lattice.crd, and lattice.izp") from exc
    if not labels:
        raise ValueError("RS reference must define at least one atoms.label RS type")
    if len(type_values) != len(basis):
        raise ValueError(
            "lattice.izp must contain exactly one RS type index for each lattice.crd site"
        )
    if any(type_index < 1 or type_index > len(labels) for type_index in type_values):
        raise ValueError("lattice.izp contains an RS type index outside atoms.label")
    return labels, type_values


def _recursion_sites(lattice_nml: dict) -> tuple[int, ...]:
    values = _as_list(lattice_nml["lattice"]["irec"])
    return tuple(int(value) for value in values)


def _nsp(input_nml: dict) -> int | None:
    """Read ``control.nsp`` without inventing a default for old inputs."""

    try:
        value = input_nml["control"].get("nsp")
    except (AttributeError, KeyError, TypeError):
        return None
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"reference control.nsp must be an integer, got {value!r}") from exc


@dataclass(frozen=True)
class ReferenceCalculation:
    """All source data required to plan exchange jobs."""

    source_dir: Path
    input_path: Path
    lattice_path: Path
    cluster_path: Path
    report_path: Path
    input_nml: dict
    lattice_nml: dict
    structure: Structure
    analyzer: SpacegroupAnalyzer
    canonical_structure: RSStructure
    cluster: Molecule
    recursion_sites: tuple[int, ...]
    fermi_energy: float
    nsp: int | None = None

    @property
    def nsp_semantics(self) -> NspSemantics:
        """Return the centralized source-derived ``nsp`` interpretation."""

        return nsp_semantics(self.nsp)

    @property
    def is_collinear(self) -> bool | None:
        """Whether the source calculation uses a collinear representation."""

        return self.nsp_semantics.is_collinear

    @property
    def is_noncollinear(self) -> bool | None:
        """Whether the source calculation uses a noncollinear representation."""

        return self.nsp_semantics.is_noncollinear

    @property
    def has_soc(self) -> bool | None:
        """Whether the source calculation includes spin-orbit coupling."""

        return self.nsp_semantics.has_soc


def read_reference_calculation(
    source_dir: str | Path = ".",
    *,
    input_path: str | Path | None = None,
    lattice_path: str | Path | None = None,
    cluster_path: str | Path | None = None,
    report_path: str | Path | None = None,
    verbose: bool = False,
) -> ReferenceCalculation:
    """Load a reference calculation using only explicit filesystem paths."""

    root = Path(source_dir).expanduser().resolve()

    def resolve(path: str | Path | None, default: str) -> Path:
        candidate = Path(default if path is None else path).expanduser()
        return candidate if candidate.is_absolute() else root / candidate

    input_file = resolve(input_path, "input.nml")
    lattice_file = resolve(lattice_path, "lattice.nml")
    cluster_file = resolve(cluster_path, "clust")
    report_file = resolve(report_path, "report.out")
    input_nml, lattice_nml = read_namelists(input_file, lattice_file, verbose=verbose)
    structure, analyzer = build_structure_and_analyzer(lattice_nml, input_nml)
    return ReferenceCalculation(
        source_dir=root,
        input_path=input_file,
        lattice_path=lattice_file,
        cluster_path=cluster_file,
        report_path=report_file,
        input_nml=input_nml,
        lattice_nml=lattice_nml,
        structure=structure,
        analyzer=analyzer,
        canonical_structure=build_rs_structure(lattice_nml, input_nml),
        cluster=read_rs_clust(cluster_file, input_file),
        recursion_sites=_recursion_sites(lattice_nml),
        fermi_energy=parse_fermi_energy(report_file),
        nsp=_nsp(input_nml),
    )


__all__ = [
    "ReferenceCalculation",
    "build_rs_structure",
    "build_structure_and_analyzer",
    "read_namelists",
    "read_reference_calculation",
    "read_rs_clust",
]
