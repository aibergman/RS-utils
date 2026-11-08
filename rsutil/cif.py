"""Convert crystallographic CIF files into RS-LMTO ``lattice.nml``/``input.nml``.

This replaces the external ``cif2cell`` step for ordinary, fully ordered
structures.  CIF parsing and space-group expansion are delegated to pymatgen
(already a core RS-utils dependency); this module only maps the resulting
structure onto the RS-LMTO source-namelist conventions:

* lattice vectors ``a(:, i)`` and site coordinates ``crd(:, j)`` are written
  in units of ``alat``; ``alat`` itself is in angstrom and defaults to the
  length of the first lattice vector, as ``cif2cell`` does;
* symmetry is **not** used by RS-LMTO here, so every atom of the cell is its
  own electronic type: ``izp = no = iu = ib = irec = atom index``, and
  ``ntot = nbas = nrec = nbulk_bulk = ntype = number of atoms``;
* atoms are ordered by descending atomic number, then by CIF site label,
  and labelled ``<Element><running count>`` (``Ga1``, ``Ga2``, ... ``Ni1``).

Partially occupied (disordered) CIFs are rejected rather than silently
collapsed to the majority species.

The cutoff ``ct`` is the neighbour cutoff (angstrom).  Unless given, it is
chosen by :func:`suggest_ct` so every atom has at least ``ct_neighbors`` atoms
(itself included) inside the cutoff, typically 15 for bcc and 19 for fcc.
"""

from __future__ import annotations

import io
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .core import RSStructure, StructureSite
from .core.errors import ValidationError


class CifConversionError(ValueError):
    """Raised when a CIF cannot be converted to a valid RS-LMTO structure."""


@dataclass(frozen=True)
class CifConversion:
    """Result of :func:`read_cif`: the structure plus conversion provenance."""

    structure: RSStructure
    alat: float
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class RSInputFiles:
    """Paths written by :func:`write_rs_input`."""

    directory: Path
    lattice: Path
    input: Path
    ct: float = 0.0
    r2: float = 0.0
    ct_auto: bool = False
    neighbors_min: int = 0
    neighbors_max: int = 0


# Groups written verbatim after &lattice/&atoms when no template is supplied.
# These mirror the cif2cell-generated defaults and are *placeholders* for a
# scientist to review, not derived physics.
DEFAULT_INPUT_GROUPS: dict[str, dict[str, Any]] = {
    "calculation": {"pre_processing": "bravais", "verbose": True},
    "self": {"nstep": 25, "cold": True},
    "energy": {
        "fermi": 0.0,
        "energy_min": -1.5,
        "energy_max": 0.5,
        "channels_ldos": 2000,
    },
    "control": {"calctype": "B", "nsp": 1, "lld": 26, "recur": "block"},
    "mix": {"beta": 0.1, "mixtype": "broyden"},
    "hamiltonian": {"hoh": False},
}


def _natural_key(text: str) -> tuple[Any, ...]:
    return tuple(int(part) if part.isdigit() else part for part in re.split(r"(\d+)", text))


def _parse_structure(path: Path, *, primitive: bool, site_tolerance: float):
    try:
        from pymatgen.io.cif import CifParser
    except ImportError as exc:  # pragma: no cover - pymatgen is a core dependency
        raise CifConversionError("CIF conversion requires pymatgen") from exc

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            parser = CifParser(str(path), site_tolerance=site_tolerance)
            structures = parser.parse_structures(primitive=primitive, on_error="raise")
        except Exception as exc:  # pymatgen raises many exception types
            raise CifConversionError(f"could not parse CIF {path}: {exc}") from exc
    if not structures:
        raise CifConversionError(f"no structure found in CIF {path}")
    messages = tuple(dict.fromkeys(str(item.message) for item in caught))
    return structures[0], messages


def read_cif(
    path: str | Path,
    *,
    alat: float | None = None,
    supercell: Sequence[int] = (1, 1, 1),
    primitive: bool = False,
    site_tolerance: float = 1e-4,
) -> CifConversion:
    """Read ``path`` and return a canonical :class:`RSStructure` in angstrom.

    Parameters
    ----------
    alat:
        Length unit in angstrom for the written namelists.  Defaults to the
        length of the first lattice vector of the (pre-supercell) cell.
    supercell:
        Ordered replication ``(n1, n2, n3)`` applied after parsing.
    primitive:
        Reduce to the primitive cell instead of keeping the CIF cell.
    site_tolerance:
        Fractional distance below which pymatgen merges duplicate sites.
    """

    cif_path = Path(path).expanduser()
    if not cif_path.is_file():
        raise CifConversionError(f"CIF file does not exist: {cif_path}")
    try:
        repeat = tuple(int(n) for n in supercell)
    except (TypeError, ValueError) as exc:
        raise CifConversionError("supercell must be three positive integers") from exc
    if len(repeat) != 3 or any(n < 1 for n in repeat):
        raise CifConversionError("supercell must be three positive integers")

    pmg, parse_warnings = _parse_structure(
        cif_path, primitive=primitive, site_tolerance=site_tolerance
    )

    if not pmg.is_ordered:
        disordered = sorted(
            {
                f"{getattr(site, 'label', None) or site.species_string} "
                f"({site.species_string})"
                for site in pmg
                if not site.is_ordered
            }
        )
        raise CifConversionError(
            "CIF contains partially occupied sites, which are not supported "
            "yet: " + ", ".join(disordered)
        )

    unit_alat = float(np.linalg.norm(pmg.lattice.matrix[0]))
    alat_value = unit_alat if alat is None else float(alat)
    if not np.isfinite(alat_value) or alat_value <= 0:
        raise CifConversionError("alat must be a positive, finite length in angstrom")

    if repeat != (1, 1, 1):
        pmg = pmg.copy()
        pmg.make_supercell(repeat)

    entries = []
    for site in pmg:
        element = site.specie.element if hasattr(site.specie, "element") else site.specie
        frac = np.mod(np.asarray(site.frac_coords, dtype=float), 1.0)
        frac[np.isclose(frac, 1.0, atol=1e-9, rtol=0.0)] = 0.0
        label = getattr(site, "label", None) or ""
        entries.append((element, label, frac))

    entries.sort(
        key=lambda item: (
            -item[0].Z,
            _natural_key(str(item[1])),
            tuple(np.round(item[2], 8)),
        )
    )

    lattice = np.asarray(pmg.lattice.matrix, dtype=float)
    counts: dict[str, int] = {}
    sites = []
    for index, (element, _cif_label, frac) in enumerate(entries):
        symbol = element.symbol
        counts[symbol] = counts.get(symbol, 0) + 1
        cartesian = frac @ lattice
        sites.append(
            StructureSite(
                index=index,
                position=tuple(float(x) for x in cartesian),
                element=symbol,
                atomic_number=int(element.Z),
                rs_label=f"{symbol}{counts[symbol]}",
                rs_type_index=index + 1,
            )
        )

    metadata: dict[str, Any] = {
        "source_cif": cif_path.name,
        "alat_angstrom": alat_value,
        "supercell": list(repeat),
        "primitive": bool(primitive),
    }
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        metadata["spacegroup"] = SpacegroupAnalyzer(pmg, symprec=1e-3).get_space_group_symbol()
    except Exception:  # pragma: no cover - informational only
        pass

    try:
        structure = RSStructure(
            lattice_vectors=tuple(tuple(float(x) for x in row) for row in lattice),
            length_unit="angstrom",
            sites=tuple(sites),
            metadata=metadata,
        )
    except ValidationError as exc:
        raise CifConversionError(f"invalid structure from {cif_path}: {exc}") from exc
    return CifConversion(structure=structure, alat=alat_value, warnings=parse_warnings)


def _pymatgen_structure(structure: RSStructure):
    from pymatgen.core import Lattice, Structure

    return Structure(
        Lattice(np.asarray(structure.lattice_vectors, dtype=float)),
        [site.element for site in structure.sites],
        [site.position for site in structure.sites],
        coords_are_cartesian=True,
    )


def neighbor_counts(structure: RSStructure, ct: float) -> tuple[int, ...]:
    """Atoms within ``ct`` angstrom of each site, counting the site itself."""

    if structure.length_unit.lower() not in {"angstrom", "å"}:
        raise CifConversionError("structure lengths must be in angstrom")
    centers, _, _, _ = _pymatgen_structure(structure).get_neighbor_list(float(ct))
    counts = np.bincount(centers, minlength=len(structure.sites)) + 1
    return tuple(int(c) for c in counts)


def suggest_ct(
    structure: RSStructure,
    *,
    min_atoms: int = 15,
    margin: float = 0.3,
    shell_tol: float = 0.1,
) -> float:
    """Pick a neighbour cutoff ``ct`` (angstrom) for ``structure``.

    The cutoff is the smallest radius at which **every** atom has at least
    ``min_atoms`` atoms within it (itself included), extended to the end of
    that neighbour shell (distances within ``shell_tol`` of each other form
    one shell, so a split shell is never cut through) and then placed
    ``margin`` beyond it, or half the gap to the next shell if that is
    smaller.  This gives 15 atoms for bcc and 19 for fcc.
    """

    if not isinstance(min_atoms, int) or isinstance(min_atoms, bool) or min_atoms < 2:
        raise CifConversionError("ct_neighbors must be an integer >= 2")
    if structure.length_unit.lower() not in {"angstrom", "å"}:
        raise CifConversionError("structure lengths must be in angstrom")
    if len(structure.sites) == 0:
        raise CifConversionError("cannot choose ct for a structure without sites")
    pmg = _pymatgen_structure(structure)
    n_sites = len(structure.sites)
    radius = max(3.0, float(np.cbrt(pmg.volume / n_sites)) * 1.6)
    for _ in range(12):
        centers, _, _, dist = pmg.get_neighbor_list(radius)
        per_site = [np.sort(dist[centers == i]) for i in range(n_sites)]
        if all(len(d) >= min_atoms - 1 for d in per_site):
            reach = max(float(d[min_atoms - 2]) for d in per_site)
            all_d = np.unique(np.round(dist, 6))
            end = reach
            while True:
                later = all_d[all_d > end + 1e-9]
                if later.size and later[0] - end < shell_tol:
                    end = float(later[0])
                    continue
                break
            later = all_d[all_d > end + 1e-9]
            if later.size:  # next shell is inside the search radius: gap is known
                gap = float(later[0]) - end
                return round(end + min(margin, gap / 2.0), 2)
        radius *= 1.5
    raise CifConversionError(
        f"could not choose ct for {min_atoms} neighbours; pass ct explicitly"
    )


def _fmt(value: float) -> str:
    # 10 decimals like cif2cell; normalise negative zero.
    text = f"{value:.10f}"
    return "0.0000000000" if text == "-0.0000000000" else text


def format_lattice_nml(structure: RSStructure, alat: float, *, header: str = "") -> str:
    """Return ``lattice.nml`` text for ``structure`` (one RS type per atom)."""

    if structure.length_unit.lower() not in {"angstrom", "å"}:
        raise CifConversionError("structure lengths must be in angstrom")
    n = len(structure.sites)
    lattice = np.asarray(structure.lattice_vectors, dtype=float) / alat
    lines = [header] if header else []
    lines += [
        "&lattice",
        f"nbulk_bulk  = {n}",
        f"ntot        = {n}",
        f"nbas        = {n}",
        f"nrec        = {n}",
        "!",
    ]
    for i, row in enumerate(lattice, start=1):
        lines.append(f"a(:, {i}) = " + ", ".join(_fmt(x) for x in row))
    lines.append("!")
    for i, site in enumerate(structure.sites, start=1):
        pos = np.asarray(site.position, dtype=float) / alat
        lines.append(f"crd(:, {i}) = " + ", ".join(_fmt(x) for x in pos))
    for name in ("izp", "no", "iu", "ib", "irec"):
        lines.append("!")
        lines.extend(f"{name}({i}) = {i}" for i in range(1, n + 1))
    lines += ["!", "/", ""]
    return "\n".join(lines)


def _load_template(template: str | Path | Mapping[str, Any] | None) -> dict[str, dict[str, Any]]:
    if template is None:
        return {name: dict(group) for name, group in DEFAULT_INPUT_GROUPS.items()}
    if isinstance(template, Mapping):
        source = template
    else:
        import f90nml

        path = Path(template).expanduser()
        if not path.is_file():
            raise CifConversionError(f"input template does not exist: {path}")
        try:
            source = f90nml.read(str(path))
        except Exception as exc:
            raise CifConversionError(f"could not read input template {path}: {exc}") from exc
    return {str(name).lower(): dict(group) for name, group in source.items()}


def _format_other_groups(groups: Mapping[str, Mapping[str, Any]]) -> str:
    import f90nml

    if not groups:
        return ""
    # Namelist's constructor sorts keys; assigning one at a time preserves the
    # template's group and variable order.
    document = f90nml.Namelist()
    for name, group in groups.items():
        block = f90nml.Namelist()
        for key, value in group.items():
            block[key] = value
        document[name] = block
    buffer = io.StringIO()
    document.write(buffer, force=True)
    return buffer.getvalue()


def format_input_nml(
    structure: RSStructure,
    alat: float,
    *,
    ct: float,
    r2: float | None = None,
    rc: float = 50,
    wav: float | None = None,
    template: str | Path | Mapping[str, Any] | None = None,
    header: str = "",
) -> str:  # ``ct`` is resolved by the caller (see write_rs_input)
    """Return ``input.nml`` text.

    ``&lattice`` and ``&atoms`` are generated from the structure; every other
    namelist group comes from ``template`` (or :data:`DEFAULT_INPUT_GROUPS`).
    Template ``&lattice`` entries such as ``wav`` are kept unless overridden.
    """

    if not np.isfinite(ct) or ct <= 0:
        raise CifConversionError("ct must be a positive, finite cutoff")
    r2_value = float(ct) ** 2 if r2 is None else float(r2)
    n = len(structure.sites)
    groups = _load_template(template)
    lattice_extra = {
        key: value
        for key, value in groups.pop("lattice", {}).items()
        if key not in {"alat", "ntype", "ct", "r2", "crystal_sym", "rc"}
    }
    atoms_extra = {
        key: value
        for key, value in groups.pop("atoms", {}).items()
        if key not in {"label"}
    }
    if wav is not None:
        lattice_extra["wav"] = float(wav)

    lines = [header] if header else []
    lines += [
        "&lattice",
        f"rc = {rc:g}",
        f"alat = {alat:.10g}",
        "crystal_sym = 'file'",
        f"ntype       = {n}",
    ]
    for key, value in lattice_extra.items():
        lines.append(f"{key} = {_scalar(value)}")
    lines += [f"ct({i}) = {ct:.2f}" for i in range(1, n + 1)]
    lines += [f"r2 = {r2_value:.10g}", "/", "&atoms"]
    lines.append(f"database = {_scalar(atoms_extra.pop('database', './'))}")
    for key, value in atoms_extra.items():
        lines.append(f"{key} = {_scalar(value)}")
    lines += [
        f"label({i}) = '{site.rs_label}'" for i, site in enumerate(structure.sites, start=1)
    ]
    lines.append("/")
    return "\n".join(lines) + "\n" + _format_other_groups(groups)


def _scalar(value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, (list, tuple)):
        return ", ".join(_scalar(item) for item in value)
    return repr(value)


def write_rs_input(
    conversion: CifConversion,
    directory: str | Path,
    *,
    ct: float | None = None,
    ct_neighbors: int = 15,
    r2: float | None = None,
    rc: float = 50,
    wav: float | None = None,
    template: str | Path | Mapping[str, Any] | None = None,
    overwrite: bool = False,
) -> RSInputFiles:
    """Write ``lattice.nml`` and ``input.nml`` into ``directory``.

    If ``ct`` is omitted it is chosen by :func:`suggest_ct` with
    ``ct_neighbors`` atoms per site.
    """

    destination = Path(directory).expanduser()
    lattice_path = destination / "lattice.nml"
    input_path = destination / "input.nml"
    if not overwrite and (lattice_path.exists() or input_path.exists()):
        raise FileExistsError(
            f"{destination} already contains lattice.nml or input.nml; "
            "pass overwrite=True to replace them"
        )
    structure = conversion.structure
    ct_auto = ct is None
    if ct is None:
        ct = suggest_ct(structure, min_atoms=ct_neighbors)
    elif not np.isfinite(ct) or ct <= 0:
        raise CifConversionError("ct must be a positive, finite cutoff")
    counts = neighbor_counts(structure, float(ct))
    source = structure.metadata.get("source_cif", "CIF")
    header = (
        f"! Generated by rsutil.cif from {source}\n"
        "! One RS electronic type per atom; ct/r2 and non-structural groups are\n"
        "! settings to review, not values derived from the CIF."
    )
    # Build both texts before touching the filesystem so a bad argument does
    # not leave a half-written bundle behind.
    lattice_text = format_lattice_nml(structure, conversion.alat, header=header)
    input_text = format_input_nml(
        structure,
        conversion.alat,
        ct=ct,
        r2=r2,
        rc=rc,
        wav=wav,
        template=template,
        header=header,
    )
    destination.mkdir(parents=True, exist_ok=True)
    lattice_path.write_text(lattice_text, encoding="utf-8")
    input_path.write_text(input_text, encoding="utf-8")
    return RSInputFiles(
        destination,
        lattice_path,
        input_path,
        ct=float(ct),
        r2=float(ct) ** 2 if r2 is None else float(r2),
        ct_auto=ct_auto,
        neighbors_min=min(counts),
        neighbors_max=max(counts),
    )


def convert_cif(
    cif: str | Path,
    directory: str | Path,
    *,
    ct: float | None = None,
    ct_neighbors: int = 15,
    alat: float | None = None,
    supercell: Sequence[int] = (1, 1, 1),
    primitive: bool = False,
    r2: float | None = None,
    rc: float = 50,
    wav: float | None = None,
    template: str | Path | Mapping[str, Any] | None = None,
    overwrite: bool = False,
) -> tuple[RSInputFiles, CifConversion]:
    """Convenience wrapper: :func:`read_cif` then :func:`write_rs_input`."""

    conversion = read_cif(cif, alat=alat, supercell=supercell, primitive=primitive)
    files = write_rs_input(
        conversion,
        directory,
        ct=ct,
        ct_neighbors=ct_neighbors,
        r2=r2,
        rc=rc,
        wav=wav,
        template=template,
        overwrite=overwrite,
    )
    return files, conversion


__all__ = [
    "CifConversion",
    "CifConversionError",
    "DEFAULT_INPUT_GROUPS",
    "RSInputFiles",
    "convert_cif",
    "format_input_nml",
    "format_lattice_nml",
    "neighbor_counts",
    "read_cif",
    "suggest_ct",
    "write_rs_input",
]
