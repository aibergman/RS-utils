"""Prepare RS-LMTO-ASA surface calculation input from a converged bulk run.

This module only prepares the surface input bundle.  The RS-LMTO executable's
``buildsurf`` preprocessor remains responsible for generating its cluster and
surface maps when the calculation starts.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import tempfile
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import f90nml
import numpy as np

from .core.contracts import SetupSpec
from .setup_runner import FortranSetupRunner


# P_l = n_l + 1/2 - atan(l)/pi for a zero-potential free-electron
# partial wave on the lowest radial branch (n_l = l + 1).
_EMPTY_SPHERE_FREE_ELECTRON_PL = tuple(
    l + 1.5 - math.atan(l) / math.pi for l in range(3)
)


@dataclass(frozen=True)
class SurfaceSpec:
    """User-controlled surface setup parameters."""

    miller: tuple[int, ...] = (0, 0, 1)
    nlay: int = 6
    empty_spheres: int = 1
    empty_sphere_charge: float = 0.2
    charge_sigma: float = 0.5
    crystal_sym: str | None = None
    orientation: str = "preserve"
    n1: int | None = 32
    n2: int | None = 32
    n3: int | None = 32
    rc: int | None = 100
    self_steps: int = 15
    mixing_beta: float = 0.01

    def __post_init__(self) -> None:
        if len(self.miller) not in {3, 4} or any(isinstance(x, bool) or int(x) != x for x in self.miller):
            raise ValueError("miller must contain three integers (or four for hcp Miller-Bravais indices)")
        if not any(self.miller):
            raise ValueError("miller indices cannot all be zero")
        if self.nlay < 1:
            raise ValueError("nlay must be at least 1")
        if self.empty_spheres < 0 or self.empty_spheres > self.nlay:
            raise ValueError("empty_spheres must be between zero and nlay")
        if not math.isfinite(self.empty_sphere_charge) or self.empty_sphere_charge < 0:
            raise ValueError("empty_sphere_charge must be finite and non-negative")
        if not math.isfinite(self.charge_sigma) or self.charge_sigma <= 0:
            raise ValueError("charge_sigma must be finite and positive")
        if self.crystal_sym is not None and self.crystal_sym.lower() not in {"fcc", "bcc", "hcp", "file"}:
            raise ValueError("crystal_sym must be fcc, bcc, hcp, or file")
        if self.orientation not in {"preserve", "global-z"}:
            raise ValueError("orientation must be 'preserve' or 'global-z'")
        for name in ("n1", "n2", "n3"):
            value = getattr(self, name)
            if value is not None and value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.rc is not None and self.rc < 1:
            raise ValueError("rc must be a positive integer")
        if self.self_steps < 1:
            raise ValueError("self_steps must be a positive integer")
        if not math.isfinite(self.mixing_beta) or not 0 < self.mixing_beta <= 1:
            raise ValueError("mixing_beta must be in (0, 1]")


def _find_key(group: Mapping[str, Any], name: str) -> str | None:
    return next((str(key) for key in group if str(key).lower() == name.lower()), None)


def _set(group: f90nml.Namelist, name: str, value: Any) -> None:
    key = _find_key(group, name) or name
    group[key] = value


def _write_lattice_namelist(path: Path, document: f90nml.Namelist, basis: np.ndarray) -> None:
    """Write a lattice namelist with each ``crd`` entry as a 3-vector.

    f90nml squeezes a 3x1 array to 1D and writes ``crd(:,1)``, ``crd(:,2)``,
    and ``crd(:,3)``. The Fortran lattice reader expects one three-vector per
    basis site, so restore the intended column assignments after serialization.
    """
    document.write(str(path), force=True)
    lines = path.read_text(encoding="utf-8").splitlines()
    crd_line = re.compile(r"^\s*crd\s*\(\s*:\s*,\s*\d+\s*\)\s*=", re.IGNORECASE)
    first = next((index for index, line in enumerate(lines) if crd_line.match(line)), None)
    if first is None:
        raise ValueError(f"lattice namelist {path} has no serialized crd array")
    last = first
    while last + 1 < len(lines) and crd_line.match(lines[last + 1]):
        last += 1
    replacement = [
        f"    crd(:,{site}) = " + ", ".join(f"{float(value):.16g}" for value in coordinate)
        for site, coordinate in enumerate(np.asarray(basis, dtype=float), start=1)
    ]
    lines = lines[:first] + replacement + lines[last + 1:]

    # The Fortran reader first reads this namelist before allocating its
    # allocatable arrays. It therefore needs nrec/nbulk_bulk ahead of a/crd;
    # alphabetical f90nml output places those dimensions too late.
    group_start = next(
        (index for index, line in enumerate(lines) if re.match(r"^\s*&lattice\b", line, re.IGNORECASE)),
        None,
    )
    if group_start is None:
        raise ValueError(f"lattice namelist {path} has no &lattice group")
    dimensions = re.compile(r"^\s*(?:nbulk_bulk|ntot|nbas|ntype|nrec)\s*=", re.IGNORECASE)
    dimension_lines: list[str] = []
    remaining: list[str] = []
    in_lattice = False
    for index, line in enumerate(lines):
        if index == group_start:
            in_lattice = True
            remaining.append(line)
            continue
        if in_lattice and line.strip() == "/":
            in_lattice = False
            remaining.append(line)
            continue
        if in_lattice and dimensions.match(line):
            dimension_lines.append(line)
        else:
            remaining.append(line)
    remaining[group_start + 1:group_start + 1] = dimension_lines
    path.write_text("\n".join(remaining) + "\n", encoding="utf-8")


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, (list, tuple)):
        return list(value)
    return [] if value is None else [value]


def _reduce_miller(miller: tuple[int, ...]) -> tuple[int, ...]:
    divisor = math.gcd(*(abs(value) for value in miller))
    return tuple(int(value // divisor) for value in miller)


def _rotation_to_global_z(normal: np.ndarray) -> np.ndarray:
    """Return a proper rotation taking a Cartesian normal onto positive Z."""
    vector = np.asarray(normal, dtype=float)
    vector /= np.linalg.norm(vector)
    target = np.array([0.0, 0.0, 1.0])
    cross = np.cross(vector, target)
    sine = np.linalg.norm(cross)
    cosine = float(np.dot(vector, target))
    if sine < 1.0e-14:
        if cosine > 0:
            return np.eye(3)
        return np.diag([1.0, -1.0, -1.0])
    skew = np.array(
        [[0.0, -cross[2], cross[1]], [cross[2], 0.0, -cross[0]], [-cross[1], cross[0], 0.0]]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / (sine * sine))


def _builtin_cell(symmetry: str, c_over_a: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return RS-LMTO primitive cell, basis sites, and site type indices."""
    if symmetry == "fcc":
        return np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0.]]), np.zeros((1, 3)), np.array([1])
    if symmetry == "bcc":
        return np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]), np.zeros((1, 3)), np.array([1])
    if symmetry == "hcp":
        cell = np.array([[1, 0, 0], [-0.5, math.sqrt(3) / 2, 0], [0, 0, c_over_a]], dtype=float)
        basis = np.array([[0, 0, 0], [0, 1 / math.sqrt(3), c_over_a / 2]], dtype=float)
        return cell, basis, np.array([1, 2])
    raise ValueError(f"cannot construct a lattice for crystal_sym={symmetry!r}")


def _rotate_file_lattice(
    source_file: Path | None,
    symmetry: str,
    c_over_a: float,
    hkl: tuple[int, ...],
    bulk_ntype: int,
) -> tuple[f90nml.Namelist, np.ndarray]:
    """Read or build a bulk lattice and rewrite it in a surface-aligned basis."""
    if source_file is not None and source_file.is_file():
        lattice_doc = f90nml.read(source_file)
        if not lattice_doc:
            raise ValueError(f"empty lattice namelist: {source_file}")
        lattice_key = _find_key(lattice_doc, "lattice")
        if lattice_key is None:
            raise ValueError(f"&lattice group is missing from {source_file}")
        lattice_group = lattice_doc[lattice_key]
        if not isinstance(lattice_group, f90nml.Namelist):
            raise ValueError(f"lattice group is invalid in {source_file}")
        raw_a = np.asarray(lattice_group.get("a"), dtype=float)
        raw_crd = np.asarray(lattice_group.get("crd"), dtype=float)
        if raw_a.shape != (3, 3) or raw_crd.ndim != 2 or raw_crd.shape[0] != 3:
            raise ValueError("file orientation requires 3x3 lattice a and 3xN basis crd arrays")
        cell = raw_a.T
        basis = raw_crd.T
        lattice_doc = f90nml.Namelist(lattice_doc)
    else:
        cell, basis, site_types = _builtin_cell(symmetry, c_over_a)
        if len(site_types) != bulk_ntype:
            raise ValueError(
                f"built-in {symmetry} lattice has {len(site_types)} primitive sites but bulk input declares ntype={bulk_ntype}"
            )
        lattice_doc = f90nml.Namelist({
            "lattice": {
                "nbulk_bulk": bulk_ntype,
                "ntot": len(basis),
                "nbas": len(basis),
                "nrec": bulk_ntype,
                "a": cell.T.tolist(),
                "crd": basis.T.tolist(),
                "izp": site_types.tolist(),
                "no": list(range(1, len(basis) + 1)),
                "iu": list(range(1, len(basis) + 1)),
                "ib": list(range(1, bulk_ntype + 1)),
                "irec": list(range(1, bulk_ntype + 1)),
            }
        })
        lattice_group = lattice_doc["lattice"]

    if len(hkl) == 4:
        h, k, i, ell = hkl
        if i != -(h + k):
            raise ValueError("hcp Miller-Bravais indices must satisfy i = -(h + k)")
        # Match build_clusup's conversion to the Cartesian projection vector.
        normal = np.array([2 * h + k, 2 * h + 3 * k, ell], dtype=float)
    else:
        normal = np.asarray(hkl, dtype=float)
    rotation = _rotation_to_global_z(normal)
    rotated_cell = cell @ rotation.T
    rotated_basis = basis @ rotation.T
    lattice_group["a"] = rotated_cell.T.tolist()
    lattice_group["crd"] = rotated_basis.T.tolist()
    return lattice_doc, rotated_basis


def _surface_labels(bulk_labels: list[str], nlay: int, empty_spheres: int) -> list[str]:
    labels: list[str] = []
    if empty_spheres:
        labels.append("ES")
        labels.extend(f"ES-{index}" for index in range(1, empty_spheres))
    host = bulk_labels[0]
    interface_count = nlay - empty_spheres
    labels.extend(f"{host}-S" if index == 0 else f"{host}-S-{index}" for index in range(interface_count))
    return bulk_labels + labels


def _copy_bulk_restart_files(bulk: Path, candidate: Path, host_labels: list[str]) -> list[str]:
    copied: set[str] = set()
    missing: list[str] = []
    for label in host_labels:
        potential = bulk / f"{label}_out.nml"
        if not potential.is_file():
            missing.append(str(potential))
            continue
        shutil.copy2(potential, candidate / potential.name)
        copied.add(potential.name)
        definition = bulk / f"{label}.nml"
        if definition.is_file():
            shutil.copy2(definition, candidate / definition.name)
            copied.add(definition.name)
        for radial in sorted(bulk.glob(f"radial_ground_state_{label}_*.dat")):
            shutil.copy2(radial, candidate / radial.name)
            copied.add(radial.name)
    if missing:
        raise FileNotFoundError(
            "converged bulk potential file(s) are missing: " + ", ".join(missing)
        )
    return sorted(copied)


def _empty_sphere_setup_inputs(
    bulk: Path,
    work: Path,
    source_input: f90nml.Namelist,
    crystal_sym: str,
    c_over_a: float,
    bulk_ntype: int,
) -> tuple[Path, Path]:
    """Stage the bulk cell with explicit empty-sphere labels for atom setup."""
    lattice_file = bulk / "lattice.nml"
    if crystal_sym == "file":
        work.mkdir(parents=True, exist_ok=True)
        staged_lattice = work / "lattice.nml"
        staged_input = work / "input.nml"
        shutil.copy2(lattice_file, staged_lattice)
        input_doc = deepcopy(source_input)
        _set(input_doc[_find_key(input_doc, "lattice")], "crystal_sym", "file")
        _set(input_doc[_find_key(input_doc, "atoms")], "label", ["E1"] * bulk_ntype)
        input_doc.write(str(staged_input), force=True)
        return staged_lattice, staged_input

    cell, basis, site_types = _builtin_cell(crystal_sym, c_over_a)
    if len(site_types) != bulk_ntype:
        raise ValueError(
            f"built-in {crystal_sym} lattice has {len(site_types)} primitive sites, "
            f"but bulk input declares ntype={bulk_ntype}"
        )
    lattice_doc = f90nml.Namelist({
        "lattice": {
            "nbulk_bulk": bulk_ntype,
            "ntot": len(basis),
            "nbas": len(basis),
            "nrec": bulk_ntype,
            "a": cell.T.tolist(),
            "crd": basis.T.tolist(),
            "izp": site_types.tolist(),
            "no": list(range(1, len(basis) + 1)),
            "iu": list(range(1, len(basis) + 1)),
            "ib": list(range(1, bulk_ntype + 1)),
            "irec": list(range(1, bulk_ntype + 1)),
        }
    })
    input_doc = deepcopy(source_input)
    _set(input_doc[_find_key(input_doc, "lattice")], "crystal_sym", "file")
    _set(input_doc[_find_key(input_doc, "atoms")], "label", ["E1"] * bulk_ntype)
    work.mkdir(parents=True, exist_ok=True)
    lattice_path = work / "lattice.nml"
    input_path = work / "input.nml"
    _write_lattice_namelist(lattice_path, lattice_doc, basis)
    input_doc.write(str(input_path), force=True)
    return lattice_path, input_path


def _potential_valence(par: Mapping[str, Any]) -> float:
    ql_key = _find_key(par, "ql")
    ql = par.get(ql_key) if ql_key else None
    try:
        # f90nml exposes the Fortran ql(moment, l, spin) array as
        # ql[spin][l][moment]. Ql(1,l,spin) is the charge in each l channel.
        return float(sum(float(channel[0]) for spin in ql for channel in spin))
    except (TypeError, IndexError):
        return 0.0


def _occupancy_seed(par: Mapping[str, Any], amount: float, *, empty_sphere: bool) -> list[list[list[float]]]:
    ql_key = _find_key(par, "ql")
    ql = par.get(ql_key) if ql_key else None
    if not isinstance(ql, list) or len(ql) != 2 or not all(isinstance(spin, list) and spin for spin in ql):
        raise ValueError("cannot set initial surface occupancy: atom potential has no two-spin ql array")
    seeded = deepcopy(ql)
    if not all(isinstance(channel, list) and len(channel) == 3 for spin in seeded for channel in spin):
        raise ValueError("cannot set initial surface occupancy: atom potential ql array is malformed")
    if empty_sphere:
        # Keep the converged host radial parameters, but start the empty
        # sphere with a small, spin-balanced s charge and no host multipoles.
        for spin in seeded:
            for channel in spin:
                channel[:] = [0.0, 0.0, 0.0]
        seeded[0][0][0] = amount / 2.0
        seeded[1][0][0] = amount / 2.0
    return seeded


def _remove_charge(par: f90nml.Namelist, amount: float) -> list[list[list[float]]]:
    ql_key = _find_key(par, "ql")
    ql = par.get(ql_key) if ql_key else None
    current = _potential_valence(par)
    if not isinstance(ql, list) or amount >= current:
        raise ValueError("initial empty-sphere occupancy exceeds the interface atom's starting valence")
    scale = (current - amount) / current
    seeded = deepcopy(ql)
    # Transfer charge only: preserve the spin and higher-moment starting
    # values while scaling Q0 over all angular-momentum channels.
    for spin in seeded:
        for channel in spin:
            channel[0] = float(channel[0]) * scale
    return seeded


def _write_atom_seed(
    source_path: Path,
    destination: Path,
    label: str,
    *,
    ql: list[list[list[float]]],
    empty_sphere: bool,
    pl: tuple[float, ...] | None = None,
) -> None:
    """Patch an atom file while preserving RS-LMTO's native array indexing.

    f90nml rewrites zero-based orbital arrays using one-based explicit array
    sections. RS-LMTO's reader then indexes the final l channel out of bounds
    (iostat=5010), leaving part of the starting potential unread. Keep the
    original Fortran-generated file and replace only element metadata and QL.
    """
    text = source_path.read_text(encoding="utf-8")
    lines = text.splitlines()
    patched: list[str] = []
    in_element = False
    in_par = False
    ql_inserted = False
    pl_line = re.compile(r"^(\s*)pl\s*\(\s*:\s*,\s*(\d+)\s*\)\s*=", re.IGNORECASE)
    ql_line = re.compile(r"^(\s*)ql\s*\(\s*\d+\s*,\s*:\s*,\s*\d+\s*\)\s*=", re.IGNORECASE)
    assignment = re.compile(r"^(\s*)([a-z_][a-z_0-9]*)\s*=\s*([^!]*?)(\s*!.*)?$", re.IGNORECASE)
    for line in lines:
        stripped = line.strip().lower()
        if stripped.startswith("&element"):
            in_element = True
        if stripped.startswith("&par"):
            in_element = False
            in_par = True
        if in_element:
            match = assignment.match(line)
            if match:
                indent, key, _, comment = match.groups()
                key_lower = key.lower()
                if key_lower == "symbol":
                    line = f"{indent}{key} = '{label}'{comment or ''}"
                elif empty_sphere and key_lower in {"atomic_number", "core", "valence", "f_core"}:
                    line = f"{indent}{key} = 0{comment or ''}"
        if in_par and pl is not None:
            match = pl_line.match(line)
            if match:
                indent, spin_index = match.groups()
                values = ", ".join(f"{float(value):.16E}" for value in pl)
                line = f"{indent}pl(:, {spin_index}) = {values}"
        if in_par and ql_line.match(line):
            if not ql_inserted:
                indent = ql_line.match(line).group(1)
                for spin_index, spin_values in enumerate(ql, start=1):
                    for moment_index in range(3):
                        values = ", ".join(f"{float(channel[moment_index]):.16E}" for channel in spin_values)
                        patched.append(f"{indent}ql({moment_index + 1}, :, {spin_index}) = {values}")
                ql_inserted = True
            continue
        patched.append(line)
    if not ql_inserted:
        raise ValueError(f"atom potential {source_path} has no native ql(moment, :, spin) assignments")
    destination.write_text("\n".join(patched) + "\n", encoding="utf-8")


def _write_surface_atomic_files(
    bulk: Path,
    candidate: Path,
    source_input: f90nml.Namelist,
    crystal_sym: str,
    c_over_a: float,
    bulk_ntype: int,
    bulk_labels: list[str],
    surface_labels: list[str],
    empty_spheres: int,
    charge_per_empty_sphere: float,
    charge_sigma: float,
) -> list[str]:
    """Clone host interface inputs and generate ES radial parameters via setup."""
    written: list[str] = []
    host = bulk_labels[0]
    empty_labels = surface_labels[:empty_spheres]
    interface_labels = surface_labels[empty_spheres:]

    host_atom_path = bulk / f"{host}_out.nml"
    if not host_atom_path.is_file():
        host_atom_path = bulk / f"{host}.nml"
    if not host_atom_path.is_file():
        raise FileNotFoundError(f"no host atom file found for label {host!r}")
    host_atom = f90nml.read(host_atom_path)
    element_key = _find_key(host_atom, "element")
    par_key = _find_key(host_atom, "par")
    if element_key is None or par_key is None:
        raise ValueError(f"host atom file {host_atom_path} must contain &element and &par")
    host_par = host_atom[par_key]
    host_valence = _potential_valence(host_par)

    total_seed = charge_per_empty_sphere * len(empty_labels)
    # A Gaussian is the derivative of an erf step. Integrating it over each
    # unit-width representative layer gives a smooth, discrete transfer profile
    # while conserving the requested ES seed charge over the finite interface.
    if total_seed > 0 and not interface_labels:
        raise ValueError("cannot balance empty-sphere starting charge without an interface layer")
    if interface_labels and total_seed > 0:
        scale = math.sqrt(2.0) * charge_sigma
        weights = [
            0.5 * (math.erf((index + 0.5) / scale) - math.erf((index - 0.5) / scale))
            for index in range(len(interface_labels))
        ]
        norm = sum(weights)
        if norm <= 0:
            raise ValueError("charge_sigma is too large to resolve across the available interface layers")
        removals = [total_seed * weight / norm for weight in weights]
        if any(amount >= host_valence for amount in removals):
            raise ValueError("initial empty-sphere occupancy exceeds the interface starting valence")
    else:
        removals = [0.0] * len(interface_labels)

    for index, label in enumerate(interface_labels):
        ql = (
            _remove_charge(host_par, removals[index])
            if removals[index] > 0
            else _occupancy_seed(host_par, 0.0, empty_sphere=False)
        )
        path = candidate / f"{label}.nml"
        _write_atom_seed(host_atom_path, path, label, ql=ql, empty_sphere=False)
        written.append(path.name)

    if empty_labels:
        generator_work = candidate / ".rsutil-es-generator"
        try:
            lattice_path, input_path = _empty_sphere_setup_inputs(
                bulk, generator_work / "source", source_input, crystal_sym, c_over_a, bulk_ntype
            )
            run_workspace = generator_work / "run"
            result = FortranSetupRunner().run(
                lattice_path,
                input_path,
                SetupSpec(empty_spheres=False),
                workspace=run_workspace,
            )
            if not result.success:
                details = "; ".join(
                    f"{item.message} ({json.dumps(dict(item.details), sort_keys=True)})"
                    for item in result.diagnostics
                )
                worker_output = "\n".join(
                    value.strip() for value in (result.stderr, result.stdout) if value.strip()
                )
                if worker_output:
                    details += "; worker output: " + worker_output[-2000:].replace("\n", " | ")
                raise RuntimeError(f"setup could not generate the empty-sphere atom file: {details}")
            generated = sorted(path for path in run_workspace.glob("E*.nml") if path.stem[1:].isdigit())
            if not generated:
                raise RuntimeError("setup succeeded but produced no E<n>.nml empty-sphere atom file")
            seed_atom = f90nml.read(generated[0])
            seed_par_key = _find_key(seed_atom, "par")
            if seed_par_key is None:
                raise ValueError(f"generated empty-sphere atom file {generated[0]} has no &par group")
            ql = _occupancy_seed(seed_atom[seed_par_key], charge_per_empty_sphere, empty_sphere=True)
            for label in empty_labels:
                path = candidate / f"{label}.nml"
                _write_atom_seed(
                    generated[0],
                    path,
                    label,
                    ql=ql,
                    empty_sphere=True,
                    pl=_EMPTY_SPHERE_FREE_ELECTRON_PL,
                )
                written.append(path.name)
        finally:
            shutil.rmtree(generator_work, ignore_errors=True)
    return sorted(written)


def prepare_surface(
    bulk_dir: str | Path,
    output_dir: str | Path,
    spec: SurfaceSpec | None = None,
) -> dict[str, Any]:
    """Prepare a new surface calculation directory without running RS-LMTO.

    The operation fails if the output directory already exists.  A candidate
    directory is assembled beside the destination and renamed into place only
    after all required bulk restart files have been validated and copied.
    """

    spec = spec or SurfaceSpec()
    bulk = Path(bulk_dir).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    if not bulk.is_dir():
        raise FileNotFoundError(f"bulk directory does not exist: {bulk}")
    if output.exists():
        raise FileExistsError(f"surface output directory already exists: {output}")
    input_path = bulk / "input.nml"
    if not input_path.is_file():
        raise FileNotFoundError(f"bulk input namelist is missing: {input_path}")

    source = f90nml.read(input_path)
    lattice_key = _find_key(source, "lattice")
    atoms_key = _find_key(source, "atoms")
    energy_key = _find_key(source, "energy")
    if lattice_key is None or atoms_key is None or energy_key is None:
        raise ValueError("bulk input.nml must contain &lattice, &atoms, and &energy groups")
    lattice = source[lattice_key]
    atoms = source[atoms_key]
    energy = source[energy_key]

    symmetry_key = _find_key(lattice, "crystal_sym")
    crystal_sym = (spec.crystal_sym or (lattice.get(symmetry_key) if symmetry_key else "file"))
    crystal_sym = str(crystal_sym).strip().lower()
    if crystal_sym not in {"fcc", "bcc", "hcp", "file"}:
        raise ValueError(f"unsupported crystal_sym {crystal_sym!r}; choose fcc, bcc, hcp, or file")
    lattice_file = bulk / "lattice.nml"
    if crystal_sym == "file" and not lattice_file.is_file():
        raise FileNotFoundError("crystal_sym='file' requires lattice.nml in the bulk directory")
    if crystal_sym != "hcp" and len(spec.miller) != 3:
        raise ValueError("fcc, bcc, and file surfaces require three Miller indices")
    if crystal_sym == "hcp" and len(spec.miller) == 4 and spec.miller[2] != -(spec.miller[0] + spec.miller[1]):
        raise ValueError("hcp Miller-Bravais indices must satisfy i = -(h + k)")

    ntype_key = _find_key(lattice, "ntype")
    bulk_ntype = int(lattice.get(ntype_key, 0)) if ntype_key else 0
    if crystal_sym == "file":
        bulk_lattice = f90nml.read(lattice_file)
        bulk_lattice_key = _find_key(bulk_lattice, "lattice")
        if bulk_lattice_key is None:
            raise ValueError(f"&lattice group is missing from {lattice_file}")
        declared_bulk_sites = bulk_lattice[bulk_lattice_key].get("nbulk_bulk")
        if declared_bulk_sites is not None and int(declared_bulk_sites) != bulk_ntype:
            raise ValueError(
                f"bulk input declares ntype={bulk_ntype}, but lattice.nml declares "
                f"nbulk_bulk={int(declared_bulk_sites)}"
            )
    label_key = _find_key(atoms, "label")
    bulk_labels = [str(label).strip() for label in _as_list(atoms.get(label_key) if label_key else None)]
    if bulk_ntype < 1 or len(bulk_labels) != bulk_ntype:
        raise ValueError(
            f"bulk input declares ntype={bulk_ntype} but has {len(bulk_labels)} atom labels"
        )
    if any(not label for label in bulk_labels):
        raise ValueError("bulk labels must be non-empty")
    fermi_key = _find_key(energy, "fermi")
    if fermi_key is None or not math.isfinite(float(energy[fermi_key])):
        raise ValueError("bulk &energy must define a finite fermi value")
    fermi = float(energy[fermi_key])

    output_ntype = bulk_ntype + spec.nlay
    generated_labels = _surface_labels(bulk_labels, spec.nlay, spec.empty_spheres)
    if len(generated_labels) != output_ntype:
        raise ValueError("surface label generation did not match ntype = bulk ntype + nlay")
    if any(len(label) > 8 for label in generated_labels):
        raise ValueError("generated surface labels exceed the RS-LMTO 8-character label limit")
    if len({label[:8].lower() for label in generated_labels}) != len(generated_labels):
        raise ValueError("generated surface labels collide after the RS-LMTO 8-character limit")

    surface_input = f90nml.Namelist(source)
    calc_key = _find_key(surface_input, "calculation")
    if calc_key is None:
        calc_key = "calculation"
        surface_input[calc_key] = f90nml.Namelist()
    _set(surface_input[calc_key], "pre_processing", "buildsurf")

    surface_lattice = surface_input[lattice_key]
    output_symmetry = "file" if spec.orientation == "global-z" else crystal_sym
    hkl = tuple(spec.miller)
    if crystal_sym == "hcp" and len(hkl) == 3:
        hkl = (hkl[0], hkl[1], -hkl[0] - hkl[1], hkl[2])
    hkl = _reduce_miller(hkl)
    surface_miller = (0, 0, 1) if spec.orientation == "global-z" else hkl
    _set(surface_lattice, "crystal_sym", output_symmetry)
    _set(surface_lattice, "surftype", " ".join(str(x) for x in surface_miller))
    _set(surface_lattice, "nlay", spec.nlay)
    _set(surface_lattice, "ntype", output_ntype)
    if spec.rc is not None:
        _set(surface_lattice, "rc", spec.rc)
    if spec.n1 is not None:
        _set(surface_lattice, "n1", spec.n1)
    if spec.n2 is not None:
        _set(surface_lattice, "n2", spec.n2)
    if spec.n3 is not None:
        _set(surface_lattice, "n3", spec.n3)
    ct_key = _find_key(surface_lattice, "ct")
    old_ct = _as_list(surface_lattice.get(ct_key) if ct_key else None)
    default_ct = float(old_ct[0]) if old_ct else 4.0
    _set(surface_lattice, "ct", [float(old_ct[min(i, len(old_ct) - 1)]) if old_ct else default_ct for i in range(output_ntype)])

    surface_atoms = surface_input[atoms_key]
    _set(surface_atoms, "label", generated_labels)

    self_key = _find_key(surface_input, "self")
    if self_key is None:
        self_key = "self"
        surface_input[self_key] = f90nml.Namelist()
    _set(surface_input[self_key], "ws_all", True)
    _set(surface_input[self_key], "nstep", spec.self_steps)

    surface_energy = surface_input[energy_key]
    _set(surface_energy, "fermi", fermi)
    _set(surface_energy, "fix_fermi", True)

    control_key = _find_key(surface_input, "control")
    if control_key is None:
        control_key = "control"
        surface_input[control_key] = f90nml.Namelist()
    _set(surface_input[control_key], "calctype", "S")
    mix_key = _find_key(surface_input, "mix")
    if mix_key is None:
        mix_key = "mix"
        surface_input[mix_key] = f90nml.Namelist()
    _set(surface_input[mix_key], "beta", spec.mixing_beta)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.surface-", dir=output.parent))
    try:
        surface_input.write(str(temporary / "input.nml"), force=True)
        if spec.orientation == "global-z":
            alat_key = _find_key(lattice, "alat")
            if alat_key is None:
                raise ValueError("bulk &lattice must define alat for global-z orientation")
            oriented, oriented_basis = _rotate_file_lattice(
                lattice_file if crystal_sym == "file" else None,
                crystal_sym,
                float(lattice.get(_find_key(lattice, "celldm") or "celldm", 1.633) or 1.633),
                hkl,
                bulk_ntype,
            )
            _write_lattice_namelist(temporary / "lattice.nml", oriented, oriented_basis)
        elif crystal_sym == "file":
            shutil.copy2(lattice_file, temporary / "lattice.nml")
        copied = _copy_bulk_restart_files(bulk, temporary, bulk_labels)
        generated_atomic_files = _write_surface_atomic_files(
            bulk,
            temporary,
            source,
            crystal_sym,
            float(lattice.get(_find_key(lattice, "celldm") or "celldm", 1.633) or 1.633),
            bulk_ntype,
            bulk_labels,
            generated_labels[bulk_ntype:],
            spec.empty_spheres,
            spec.empty_sphere_charge,
            spec.charge_sigma,
        )
        manifest = {
            "schema_version": "rsutil-surface-setup-v1",
            "bulk_dir": str(bulk),
            "crystal_sym": crystal_sym,
            "output_crystal_sym": output_symmetry,
            "miller_requested": list(spec.miller),
            "miller_reduced": list(hkl),
            "orientation": spec.orientation,
            "nlay": spec.nlay,
            "empty_spheres": spec.empty_spheres,
            "empty_sphere_initial_occupancy": spec.empty_sphere_charge,
            "charge_transfer_profile": "erf-derivative layer weights",
            "charge_transfer_sigma_layers": spec.charge_sigma,
            "surface_grid": {"n1": spec.n1, "n2": spec.n2, "n3": spec.n3},
            "rc": spec.rc,
            "self_steps": spec.self_steps,
            "mixing_beta": spec.mixing_beta,
            "ntype_bulk": bulk_ntype,
            "ntype_surface": output_ntype,
            "host_labels": bulk_labels,
            "surface_labels": generated_labels[bulk_ntype:],
            "fermi": fermi,
            "staged_restart_files": copied,
            "generated_atomic_files": generated_atomic_files,
            "geometry_generation": "RS-LMTO buildsurf at calculation runtime",
        }
        (temporary / "surface_setup.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, output)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return manifest


__all__ = ["SurfaceSpec", "prepare_surface"]
