"""Prepare RS-LMTO embedded impurity inputs from a converged host run."""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import f90nml

from .surface import _as_list, _find_key


@dataclass(frozen=True)
class ImpuritySpec:
    """Controls one bulk or surface impurity setup."""

    xyz: str | Path
    scf_radius: float = 0.0
    self_steps: int = 100
    mixing_beta: float = 0.1

    def __post_init__(self) -> None:
        if not math.isfinite(self.scf_radius) or self.scf_radius < 0:
            raise ValueError("scf_radius must be finite and non-negative")
        if self.self_steps < 1:
            raise ValueError("self_steps must be a positive integer")
        if not math.isfinite(self.mixing_beta) or not 0 < self.mixing_beta <= 1:
            raise ValueError("mixing_beta must be in (0, 1]")


def _set(group: f90nml.Namelist, name: str, value: Any) -> None:
    key = _find_key(group, name) or name
    group[key] = value


def _xyz_atoms(path: Path) -> list[tuple[str, tuple[float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"XYZ file is incomplete: {path}")
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError("XYZ file must start with the standard atom-count line") from exc
    if count < 1 or len(lines) < count + 2:
        raise ValueError(f"XYZ declares {count} atoms but does not contain that many coordinate rows")
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for row, line in enumerate(lines[2:2 + count], start=3):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"XYZ row {row} must contain an element and three coordinates")
        element = fields[0].strip()
        if not re.fullmatch(r"[A-Z][a-z]?", element):
            raise ValueError(f"XYZ row {row} must use a chemical element symbol, got {element!r}")
        try:
            xyz = tuple(float(value) for value in fields[1:4])
        except ValueError as exc:
            raise ValueError(f"invalid coordinates on XYZ row {row}") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"non-finite coordinates on XYZ row {row}")
        atoms.append((element, xyz))
    return atoms


def _cluster_sites(path: Path) -> list[tuple[tuple[float, float, float], int]]:
    """Read RS-LMTO's paired five-column ``clust`` records."""
    if not path.is_file():
        raise FileNotFoundError(f"host cluster file is required to select impurity sites: {path}")
    sites: list[tuple[tuple[float, float, float], int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 5:
            continue
        values: list[float] = []
        try:
            values = [float(value.replace("D", "E").replace("d", "e")) for value in fields]
        except ValueError:
            continue
        if len(values) % 5:
            raise ValueError(f"unrecognized clust record in {path}: {line!r}")
        for offset in range(0, len(values), 5):
            x, y, z, typ, _ = values[offset:offset + 5]
            sites.append(((x, y, z), int(typ)))
    if not sites:
        raise ValueError(f"host cluster file contains no sites: {path}")
    return sites


def _match_site(
    xyz_angstrom: tuple[float, float, float],
    sites: list[tuple[tuple[float, float, float], int]],
    alat: float,
) -> int:
    target = tuple(value / alat for value in xyz_angstrom)
    distances = [math.dist(target, site[0]) * alat for site in sites]
    matches = [index for index, distance in enumerate(distances) if distance <= 2.0e-4]
    if not matches:
        nearest = min(distances)
        raise ValueError(
            f"XYZ position ({xyz_angstrom[0]:.6f}, {xyz_angstrom[1]:.6f}, "
            f"{xyz_angstrom[2]:.6f}) Å does not match a site in clust (nearest is {nearest:.6g} Å)"
        )
    if len(matches) != 1:
        raise ValueError(f"XYZ position matches {len(matches)} host sites; the selection is ambiguous")
    return matches[0]


def _element_symbol(path: Path) -> str:
    atom = f90nml.read(path)
    element_key = _find_key(atom, "element")
    if element_key is None:
        raise ValueError(f"atom file has no &element group: {path}")
    symbol = str(atom[element_key].get(_find_key(atom[element_key], "symbol") or "symbol", "")).strip()
    match = re.match(r"([A-Z][a-z]?)", symbol)
    if match is None:
        raise ValueError(f"cannot determine chemical element from atom file {path}")
    return match.group(1)


def _is_empty_sphere(path: Path) -> bool:
    atom = f90nml.read(path)
    element_key = _find_key(atom, "element")
    if element_key is None:
        return False
    number_key = _find_key(atom[element_key], "atomic_number")
    if number_key is not None:
        try:
            return int(atom[element_key][number_key]) == 0
        except (TypeError, ValueError):
            pass
    return str(atom[element_key].get(_find_key(atom[element_key], "symbol") or "symbol", "")).strip().upper() == "ES"


def _unique_label(element: str, kind: str, index: int) -> str:
    label = f"{element}-{kind}{index:03d}"
    if len(label) > 8:
        label = f"{element[:2]}{kind}{index:04d}"
    if len(label) > 8:
        raise ValueError("too many selected sites for the RS-LMTO eight-character label limit")
    return label


def _write_xyz_selection(path: Path, selections: list[tuple[str, tuple[float, float, float]]]) -> None:
    rows = [str(len(selections)), "Generated by rsutil setup impurity; coordinates in host Cartesian frame (Å)"]
    rows.extend(f"{element:<8s} {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}" for element, xyz in selections)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _patch_inclu_assignments(path: Path, coords: list[tuple[float, float, float]]) -> None:
    """Write explicit Fortran inclu(j,:) rows to avoid f90nml dimension ambiguity."""
    lines = path.read_text(encoding="utf-8").splitlines()
    lattice_start = next((i for i, line in enumerate(lines) if re.match(r"\s*&lattice\b", line, re.I)), None)
    if lattice_start is None:
        raise ValueError("input namelist is missing its &lattice group")
    lattice_end = next((i for i in range(lattice_start + 1, len(lines)) if lines[i].strip() == "/"), None)
    if lattice_end is None:
        raise ValueError("input namelist has an unterminated &lattice group")
    lines = [line for line in lines if not re.match(r"\s*inclu(?:\s*\(.*\))?\s*=", line, re.I)]
    lattice_start = next(i for i, line in enumerate(lines) if re.match(r"\s*&lattice\b", line, re.I))
    lattice_end = next(i for i in range(lattice_start + 1, len(lines)) if lines[i].strip() == "/")
    assignments = [
        f"inclu({index}, :) = " + ", ".join(f"{value:.14g}" for value in coord)
        for index, coord in enumerate(coords, start=1)
    ]
    lines[lattice_end:lattice_end] = assignments
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _copy_host_files(host: Path, output: Path, labels: list[str]) -> list[str]:
    copied: set[str] = set()
    missing: list[str] = []
    for label in labels:
        out = host / f"{label}_out.nml"
        if not out.is_file():
            missing.append(str(out))
            continue
        shutil.copy2(out, output / out.name)
        copied.add(out.name)
        seed = host / f"{label}.nml"
        if seed.is_file():
            shutil.copy2(seed, output / seed.name)
            copied.add(seed.name)
        for radial in host.glob(f"radial_ground_state_{label}_*.dat"):
            shutil.copy2(radial, output / radial.name)
            copied.add(radial.name)
    if missing:
        raise FileNotFoundError("converged host potential file(s) are missing: " + ", ".join(missing))
    return sorted(copied)


def _copy_generated_seed(source: Path, destination: Path, label: str) -> None:
    lines = source.read_text(encoding="utf-8").splitlines()
    in_element = False
    changed = False
    patched: list[str] = []
    symbol_assignment = re.compile(r"^(\s*symbol\s*=\s*)('[^']*'|\"[^\"]*\")(.*)$", re.I)
    for line in lines:
        if line.strip().lower().startswith("&element"):
            in_element = True
        elif line.strip().startswith("&"):
            in_element = False
        if in_element:
            match = symbol_assignment.match(line)
            if match:
                line = f"{match.group(1)}'{label}'{match.group(3)}"
                changed = True
        patched.append(line)
    if not changed:
        raise ValueError(f"atom seed has no element symbol assignment: {source}")
    destination.write_text("\n".join(patched) + "\n", encoding="utf-8")


def _generate_element_seed(element: str, work: Path, host_input: f90nml.Namelist, alat: float) -> Path:
    """Generate a neutral atomic seed using the existing Fortran setup worker."""
    from .core.contracts import SetupSpec
    from .setup_runner import FortranSetupRunner

    work.mkdir(parents=True, exist_ok=True)
    lattice_doc = f90nml.Namelist({
        "lattice": {
            "nbulk_bulk": 1, "ntot": 1, "nbas": 1, "nrec": 1,
            "a": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "crd": [[0.0], [0.0], [0.0]], "izp": [1], "no": [1],
            "iu": [1], "ib": [1], "irec": [1],
        }
    })
    lattice_path = work / "lattice.nml"
    lattice_doc.write(str(lattice_path), force=True)
    lines = lattice_path.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if re.match(r"\s*&lattice\b", line, re.I))
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "/")
    dimensions = [line for line in lines[start + 1:end] if re.match(r"\s*(nbulk_bulk|ntot|nbas|nrec)\s*=", line, re.I)]
    lines = [line for line in lines if line not in dimensions]
    start = next(i for i, line in enumerate(lines) if re.match(r"\s*&lattice\b", line, re.I))
    lines[start + 1:start + 1] = dimensions
    lattice_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    input_doc = f90nml.Namelist(host_input)
    lattice_key = _find_key(input_doc, "lattice")
    atoms_key = _find_key(input_doc, "atoms")
    if lattice_key is None or atoms_key is None:
        raise ValueError("host input.nml must contain &lattice and &atoms groups")
    lattice = input_doc[lattice_key]
    _set(lattice, "crystal_sym", "file")
    _set(lattice, "alat", alat)
    _set(lattice, "ntype", 1)
    _set(lattice, "nlay", 0)
    _set(lattice, "nclu", 0)
    _set(lattice, "ct", [3.0])
    _set(lattice, "r2", 9.0)
    _set(input_doc[atoms_key], "label", [element])
    calc_key = _find_key(input_doc, "calculation") or "calculation"
    if calc_key not in input_doc:
        input_doc[calc_key] = f90nml.Namelist()
    _set(input_doc[calc_key], "pre_processing", "bravais")
    control_key = _find_key(input_doc, "control") or "control"
    if control_key not in input_doc:
        input_doc[control_key] = f90nml.Namelist()
    _set(input_doc[control_key], "calctype", "B")
    input_path = work / "input.nml"
    input_doc.write(str(input_path), force=True)
    run = work / "run"
    result = FortranSetupRunner().run(lattice_path, input_path, SetupSpec(empty_spheres=False), workspace=run)
    if not result.success:
        details = "; ".join(item.message for item in result.diagnostics)
        raise RuntimeError(f"could not generate {element} atomic seed with the setup worker: {details}")
    candidates = [run / f"{element}.nml", run / f"{element}_out.nml"]
    seed = next((candidate for candidate in candidates if candidate.is_file()), None)
    if seed is None:
        raise RuntimeError(f"setup worker did not generate an atom file for {element}")
    return seed


def prepare_impurity(
    host_dir: str | Path,
    output_dir: str | Path,
    spec: ImpuritySpec,
) -> dict[str, Any]:
    """Prepare impurity sites and an optional self-consistent embedding radius."""
    host = Path(host_dir).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    xyz_path = Path(spec.xyz).expanduser().resolve()
    if not host.is_dir():
        raise FileNotFoundError(f"host directory does not exist: {host}")
    if not xyz_path.is_file():
        raise FileNotFoundError(f"impurity XYZ file does not exist: {xyz_path}")
    if output.exists():
        raise FileExistsError(f"impurity output directory already exists: {output}")
    input_path = host / "input.nml"
    if not input_path.is_file():
        raise FileNotFoundError(f"host input namelist is missing: {input_path}")

    source = f90nml.read(input_path)
    lattice_key = _find_key(source, "lattice")
    atoms_key = _find_key(source, "atoms")
    energy_key = _find_key(source, "energy")
    control_key = _find_key(source, "control")
    if None in (lattice_key, atoms_key, energy_key, control_key):
        raise ValueError("host input.nml must contain &lattice, &atoms, &energy, and &control groups")
    lattice = source[lattice_key]
    source_atoms = source[atoms_key]
    labels = [str(value).strip() for value in _as_list(source_atoms.get(_find_key(source_atoms, "label") or "label"))]
    ntype_key = _find_key(lattice, "ntype")
    ntype = int(lattice.get(ntype_key, 0)) if ntype_key else 0
    if ntype < 1 or len(labels) != ntype:
        raise ValueError(f"host declares ntype={ntype} but provides {len(labels)} atom labels")
    alat_key = _find_key(lattice, "alat")
    if alat_key is None or float(lattice[alat_key]) <= 0:
        raise ValueError("host &lattice must define a positive alat")
    alat = float(lattice[alat_key])
    fermi_key = _find_key(source[energy_key], "fermi")
    if fermi_key is None or not math.isfinite(float(source[energy_key][fermi_key])):
        raise ValueError("host &energy must define a finite fermi value")
    fermi = float(source[energy_key][fermi_key])
    calctype = str(source[control_key].get(_find_key(source[control_key], "calctype") or "calctype", "")).strip().upper()
    if calctype not in {"B", "S"}:
        raise ValueError("impurity hosts must be converged bulk (calctype='B') or surface (calctype='S') calculations")
    calculation_key = _find_key(source, "calculation")
    preprocessing = "newclubulk" if calctype == "B" else "newclusurf"
    sites = _cluster_sites(host / "clust")
    if any(typ < 1 or typ > len(labels) for _, typ in sites):
        raise ValueError("host clust contains atom type indices outside the input label list")

    requested = _xyz_atoms(xyz_path)
    selected: list[tuple[str, tuple[float, float, float], int, int]] = []
    center_indices: set[int] = set()
    for element, position in requested:
        site_index = _match_site(position, sites, alat)
        if site_index in center_indices:
            raise ValueError("XYZ file selects the same host site more than once")
        center_indices.add(site_index)
        selected.append((element, sites[site_index][0], sites[site_index][1], site_index))

    grown_indices: set[int] = set()
    if spec.scf_radius > 0:
        for index, (coord, typ) in enumerate(sites):
            if index in center_indices:
                continue
            if _is_empty_sphere(host / f"{labels[typ - 1]}_out.nml"):
                continue
            if any(math.dist(coord, center_coord) * alat <= spec.scf_radius + 1.0e-8 for _, center_coord, _, _ in selected):
                grown_indices.add(index)
        if not grown_indices:
            raise ValueError("--scf-radius selected no host sites; increase the radius")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.impurity-", dir=output.parent))
    try:
        copied = _copy_host_files(host, temporary, labels)
        lattice_file = host / "lattice.nml"
        if lattice_file.is_file():
            shutil.copy2(lattice_file, temporary / lattice_file.name)
            copied.append(lattice_file.name)
        selected_types: list[tuple[str, str, tuple[float, float, float], int]] = []
        for index, (element, coord, typ, _site_index) in enumerate(selected, start=1):
            label = _unique_label(element, "I", index)
            selected_types.append((label, element, coord, typ))
        for sequence, site_index in enumerate(sorted(grown_indices), start=1):
            coord, typ = sites[site_index]
            host_label = labels[typ - 1]
            element = _element_symbol(host / f"{host_label}_out.nml")
            label = _unique_label(element, "G", sequence)
            selected_types.append((label, element, coord, typ))

        all_labels = labels + [entry[0] for entry in selected_types]
        if len({label[:8].lower() for label in all_labels}) != len(all_labels):
            raise ValueError("generated labels collide after the RS-LMTO eight-character limit")
        input_doc = f90nml.Namelist(source)
        out_lattice = input_doc[lattice_key]
        _set(out_lattice, "ntype", len(all_labels))
        old_ct = _as_list(out_lattice.get(_find_key(out_lattice, "ct") or "ct"))
        default_ct = float(old_ct[-1]) if old_ct else 3.0
        _set(out_lattice, "ct", old_ct + [default_ct] * len(selected_types))
        _set(out_lattice, "nclu", len(selected_types))
        out_atoms = input_doc[atoms_key]
        _set(out_atoms, "label", all_labels)
        out_calc_key = calculation_key or "calculation"
        if out_calc_key not in input_doc:
            input_doc[out_calc_key] = f90nml.Namelist()
        _set(input_doc[out_calc_key], "pre_processing", preprocessing)
        out_control = input_doc[control_key]
        _set(out_control, "calctype", "I")
        out_self_key = _find_key(input_doc, "self") or "self"
        if out_self_key not in input_doc:
            input_doc[out_self_key] = f90nml.Namelist()
        _set(input_doc[out_self_key], "nstep", spec.self_steps)
        out_mix_key = _find_key(input_doc, "mix") or "mix"
        if out_mix_key not in input_doc:
            input_doc[out_mix_key] = f90nml.Namelist()
        _set(input_doc[out_mix_key], "beta", spec.mixing_beta)
        _set(input_doc[energy_key], "fermi", fermi)
        _set(input_doc[energy_key], "fix_fermi", True)
        coordinates = [entry[2] for entry in selected_types]
        _set(out_lattice, "inclu", [list(coord) for coord in coordinates])
        input_output = temporary / "input.nml"
        input_doc.write(str(input_output), force=True)
        _patch_inclu_assignments(input_output, coordinates)

        generated: list[str] = []
        seed_work = temporary / ".rsutil-atomic-seeds"
        host_symbol_to_seed: dict[str, Path] = {}
        for label in labels:
            atom_path = host / f"{label}_out.nml"
            symbol = _element_symbol(atom_path)
            host_symbol_to_seed.setdefault(symbol.lower(), atom_path)
        generated_element_seeds: dict[str, Path] = {}
        for label, element, _coord, typ in selected_types:
            host_label = labels[typ - 1]
            host_seed = host / f"{host_label}_out.nml"
            if element.lower() == _element_symbol(host_seed).lower():
                seed = host_seed
            else:
                if element.lower() in host_symbol_to_seed:
                    seed = host_symbol_to_seed[element.lower()]
                else:
                    seed = generated_element_seeds.get(element.lower())
                    if seed is None:
                        seed = _generate_element_seed(element, seed_work / element, source, alat)
                        generated_element_seeds[element.lower()] = seed
            destination = temporary / f"{label}.nml"
            _copy_generated_seed(seed, destination, label)
            generated.append(destination.name)
        shutil.rmtree(seed_work, ignore_errors=True)

        xyz_out = temporary / "impurity_sites.xyz"
        selections = [(element, tuple(value * alat for value in coord)) for _label, element, coord, _typ in selected_types]
        _write_xyz_selection(xyz_out, selections)
        manifest = {
            "schema_version": "rsutil-impurity-setup-v1",
            "host_dir": str(host),
            "host_type": "bulk" if calctype == "B" else "surface",
            "xyz_input": str(xyz_path),
            "scf_radius_angstrom": spec.scf_radius,
            "impurity_sites": len(requested),
            "grown_host_sites": len(grown_indices),
            "nclu": len(selected_types),
            "fermi": fermi,
            "pre_processing": preprocessing,
            "ntype_host": ntype,
            "ntype_impurity": len(all_labels),
            "host_labels": labels,
            "added_labels": [entry[0] for entry in selected_types],
            "selected_coordinates_alat": [list(coord) for coord in coordinates],
            "staged_host_files": copied,
            "generated_atom_files": generated,
        }
        (temporary / "impurity_setup.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        shutil.rmtree(seed_work, ignore_errors=True)
        os.replace(temporary, output)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return manifest


__all__ = ["ImpuritySpec", "prepare_impurity"]
