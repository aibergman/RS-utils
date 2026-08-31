"""Optional interoperability between ASE ``Atoms`` and ``RSStructure``.

ASE is intentionally imported only when one of the conversion functions is
called.  The canonical RS-utils structure remains a small dataclass and never
stores an ASE object.

The boundary convention is explicit:

* ASE Cartesian lengths are recorded as ``angstrom`` in ``RSStructure``.
* lattice vectors are stored as rows, matching ASE's ``Cell.array`` and the
  existing RS-utils convention used by periodic pair mapping;
* site positions are Cartesian unless a caller constructs an
  ``RSStructure`` with fractional coordinates before calling :func:`to_ase`;
* PBC flags and atom order are copied exactly.

The adapter preserves JSON-compatible ``Atoms.info`` under
``RSStructure.metadata['ase_info']`` and per-atom arrays (other than the
canonical ``numbers`` and ``positions`` arrays) under
``RSStructure.metadata['ase_arrays']``.  Their dtypes are recorded where
possible.  RS labels/types use the optional ``rsutil_rs_labels`` and
``rsutil_rs_types`` info keys when present.
"""

from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .core import RSStructure, SetupSpec, StructureSite
from .setup_service import SetupService


class ASEDependencyError(ImportError):
    """Raised when an ASE conversion is requested without the optional extra."""


def _ase_atoms_type() -> type:
    try:
        from ase import Atoms
    except ImportError as exc:  # pragma: no cover - exercised with ASE absent
        raise ASEDependencyError(
            "ASE interoperability requires the optional 'ase' package; "
            "install it with 'pip install rsutil[ase]' or 'pip install ase'."
        ) from exc
    return Atoms


def _json_value(value: Any) -> Any:
    """Return a JSON-compatible copy of common ASE metadata values."""

    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        if isinstance(value, float) and not math.isfinite(value):
            raise TypeError("non-finite metadata values are not supported")
        return value
    raise TypeError(f"value of type {type(value).__name__} is not JSON-compatible")


def _json_mapping(values: Mapping[str, Any]) -> tuple[dict[str, Any], tuple[str, ...]]:
    kept: dict[str, Any] = {}
    omitted: list[str] = []
    for key, value in values.items():
        name = str(key)
        try:
            kept[name] = _json_value(value)
        except TypeError:
            omitted.append(name)
    return kept, tuple(sorted(omitted))


def _site_values(info: Mapping[str, Any], *keys: str) -> list[Any] | None:
    for key in keys:
        value = info.get(key)
        if value is not None:
            if isinstance(value, np.ndarray):
                value = value.tolist()
            if isinstance(value, (list, tuple)):
                return list(value)
    return None


def from_ase(atoms: Any) -> RSStructure:
    """Convert an ASE ``Atoms`` object into the canonical RS structure.

    Cartesian positions and cell vectors are copied in Å.  No wrapping or
    reordering is performed, which is important for stable site identity.
    """

    Atoms = _ase_atoms_type()
    if not isinstance(atoms, Atoms):
        raise TypeError(f"expected ase.Atoms, got {type(atoms).__name__}")

    cell = np.asarray(atoms.cell.array, dtype=float)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
        raise ValueError("ASE cell must be a finite 3x3 matrix")
    positions = np.asarray(atoms.get_positions(), dtype=float)
    if positions.shape != (len(atoms), 3) or not np.all(np.isfinite(positions)):
        raise ValueError("ASE Cartesian positions must be finite 3-vectors")

    info, omitted_info = _json_mapping(atoms.info)
    arrays: dict[str, Any] = {}
    dtypes: dict[str, str] = {}
    omitted_arrays: list[str] = []
    for name, value in atoms.arrays.items():
        if name in {"numbers", "positions"}:
            continue
        try:
            arrays[str(name)] = _json_value(value)
            dtypes[str(name)] = str(np.asarray(value).dtype)
        except TypeError:
            omitted_arrays.append(str(name))

    labels = _site_values(atoms.info, "rsutil_rs_labels", "rs_labels")
    types = _site_values(atoms.info, "rsutil_rs_types", "rs_types")
    if labels is not None and len(labels) != len(atoms):
        raise ValueError("ASE RS label metadata must have one entry per atom")
    if types is not None and len(types) != len(atoms):
        raise ValueError("ASE RS type metadata must have one entry per atom")

    symbols = atoms.get_chemical_symbols()
    atomic_numbers = atoms.get_atomic_numbers()
    sites = tuple(
        StructureSite(
            index=index,
            position=tuple(float(value) for value in positions[index]),
            element=str(symbols[index]),
            atomic_number=(int(atomic_numbers[index]) or None),
            rs_label=None if labels is None else str(labels[index]),
            rs_type=None if types is None else str(types[index]),
        )
        for index in range(len(atoms))
    )
    metadata: dict[str, Any] = {
        "ase_metadata_version": 1,
        "ase_info": info,
        "ase_arrays": arrays,
        "ase_array_dtypes": dtypes,
    }
    if omitted_info:
        metadata["ase_info_omitted"] = list(omitted_info)
    if omitted_arrays:
        metadata["ase_arrays_omitted"] = sorted(omitted_arrays)
    return RSStructure(
        lattice_vectors=tuple(tuple(float(value) for value in row) for row in cell),
        length_unit="angstrom",
        sites=sites,
        coordinate_convention="cartesian",
        pbc=tuple(bool(value) for value in atoms.get_pbc()),
        metadata=metadata,
    )


ase_to_rs_structure = from_ase


def _angstrom_factor(length_unit: str) -> float:
    normalized = length_unit.strip().lower().replace("å", "a")
    if normalized in {"angstrom", "ångström", "ang", "a"}:
        return 1.0
    if normalized in {"bohr", "a0", "atomic_unit", "atomic units"}:
        return 0.529177210903
    if normalized in {"nm", "nanometer", "nanometers"}:
        return 10.0
    if normalized in {"pm", "picometer", "picometers"}:
        return 0.01
    raise ValueError(
        f"cannot convert RSStructure length unit {length_unit!r} to ASE Å; "
        "use angstrom, bohr, nm, or pm"
    )


def to_ase(structure: RSStructure) -> Any:
    """Convert an ``RSStructure`` into a new ASE ``Atoms`` object."""

    Atoms = _ase_atoms_type()
    if not isinstance(structure, RSStructure):
        raise TypeError(f"expected RSStructure, got {type(structure).__name__}")
    factor = _angstrom_factor(structure.length_unit)
    cell = np.asarray(structure.lattice_vectors, dtype=float) * factor
    positions = np.asarray(structure.positions, dtype=float)
    if structure.coordinate_convention.value == "fractional":
        positions = positions @ cell
    else:
        positions = positions * factor
    atoms = Atoms(
        symbols=list(structure.elements),
        positions=positions,
        cell=cell,
        pbc=structure.pbc,
    )

    metadata = dict(structure.metadata)
    info = metadata.get("ase_info", {})
    if isinstance(info, Mapping):
        for key, value in info.items():
            try:
                atoms.info[str(key)] = _json_value(value)
            except TypeError:
                continue
    labels = [site.rs_label for site in structure.sites]
    types = [site.rs_type for site in structure.sites]
    if any(value is not None for value in labels):
        atoms.info["rsutil_rs_labels"] = labels
    if any(value is not None for value in types):
        atoms.info["rsutil_rs_types"] = types

    arrays = metadata.get("ase_arrays", {})
    dtypes = metadata.get("ase_array_dtypes", {})
    if isinstance(arrays, Mapping):
        for name, values in arrays.items():
            if str(name) in {"numbers", "positions"}:
                continue
            try:
                dtype = dtypes.get(name) if isinstance(dtypes, Mapping) else None
                array = np.asarray(values, dtype=dtype)
                if array.ndim == 0 or array.shape[0] != len(atoms):
                    continue
                atoms.new_array(str(name), array.copy())
            except (TypeError, ValueError):
                continue
    return atoms


rs_structure_to_ase = to_ase


@dataclass(frozen=True)
class ASESetupSource:
    """Files generated from ASE for consumption by :class:`SetupService`."""

    source_dir: Path
    lattice_path: Path
    input_path: Path
    structure: RSStructure


def _setup_source_text(structure: RSStructure) -> tuple[str, str]:
    """Serialize only the geometry portion of the existing RS input contract."""

    if structure.coordinate_convention.value != "cartesian":
        raise ValueError("ASE setup source serialization requires Cartesian RSStructure positions")
    cell = np.asarray(structure.lattice_vectors, dtype=float)
    determinant = float(np.linalg.det(cell))
    if not math.isfinite(determinant) or abs(determinant) <= 1.0e-12:
        raise ValueError("RS setup source requires a non-singular 3x3 ASE cell")
    positions = np.asarray(structure.positions, dtype=float)
    symbols = list(structure.elements)
    unique_symbols = list(dict.fromkeys(symbols))
    type_numbers = {symbol: index for index, symbol in enumerate(unique_symbols, start=1)}
    pbc = tuple(bool(value) for value in structure.pbc)
    nsites = len(symbols)
    volume_per_site = abs(determinant) / max(nsites, 1)
    wav = volume_per_site ** (1.0 / 3.0)
    ndim = nsites * 1000 if any(pbc) else nsites

    lattice_lines = [
        "&lattice\n",
        f"nbulk_bulk = {nsites}\n",
        f"ntot = {nsites}\n",
        f"ntype = {len(unique_symbols)}\n",
        f"nbas = {nsites}\n",
        f"nrec = {nsites}\n",
    ]
    for index, vector in enumerate(cell, start=1):
        lattice_lines.append(
            f"a(:, {index}) = {vector[0]:.16g}, {vector[1]:.16g}, {vector[2]:.16g}\n"
        )
    for index, position in enumerate(positions, start=1):
        lattice_lines.append(
            f"crd(:, {index}) = {position[0]:.16g}, {position[1]:.16g}, {position[2]:.16g}\n"
        )
    for name, values in (
        ("izp", [type_numbers[symbol] for symbol in symbols]),
        ("no", list(range(1, nsites + 1))),
        ("iu", list(range(1, nsites + 1))),
        ("ib", list(range(1, nsites + 1))),
        ("irec", list(range(1, nsites + 1))),
    ):
        for index, value in enumerate(values, start=1):
            lattice_lines.append(f"{name}({index}) = {value}\n")
    lattice_lines.extend(
        [
            "ct(:) = 0.8\n",
            f"r2 = {(1.005 * 0.8) ** 2:.16g}\n",
            "/\n",
        ]
    )

    input_lines = [
        "&lattice\n",
        "rc = 100\n",
        "alat = 1.0\n",
        "crystal_sym = 'file'\n",
        f"wav = {wav:.16g}\n",
        f"pbc = {'.true.' if any(pbc) else '.false.'}\n",
        f"b1 = {'.true.' if pbc[0] else '.false.'}\n",
        f"b2 = {'.true.' if pbc[1] else '.false.'}\n",
        f"b3 = {'.true.' if pbc[2] else '.false.'}\n",
        "n1 = 10\n",
        "n2 = 10\n",
        "n3 = 10\n",
        f"ndim = {ndim}\n",
        f"ntype = {len(unique_symbols)}\n",
        "ct(:) = 0.8\n",
        f"r2 = {(1.005 * 0.8) ** 2:.16g}\n",
        "/\n",
        "&atoms\n",
        "database = './'\n",
    ]
    input_lines.extend(
        f"label({index}) = '{symbol}'\n" for index, symbol in enumerate(unique_symbols, start=1)
    )
    input_lines.append("/\n")
    return "".join(lattice_lines), "".join(input_lines)


def write_setup_source(
    atoms: Any,
    source_dir: str | Path,
    *,
    overwrite: bool = False,
) -> ASESetupSource:
    """Write an ASE-derived ``lattice.nml``/``input.nml`` source bundle.

    The helper does not invent RS-LMTO element potential files.  Callers that
    invoke the real Fortran setup service must provide those files as explicit
    runner resources, just as they would for a hand-authored source directory.
    """

    structure = from_ase(atoms)
    destination = Path(source_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    lattice_path = destination / "lattice.nml"
    input_path = destination / "input.nml"
    if not overwrite and (lattice_path.exists() or input_path.exists()):
        raise FileExistsError(
            f"ASE setup source already contains lattice.nml or input.nml: {destination}"
        )
    lattice_text, input_text = _setup_source_text(structure)
    lattice_path.write_text(lattice_text, encoding="utf-8")
    input_path.write_text(input_text, encoding="utf-8")
    return ASESetupSource(destination, lattice_path, input_path, structure)


create_setup_source = write_setup_source


def run_setup_from_ase(
    atoms: Any,
    destination_dir: str | Path,
    spec: SetupSpec | Mapping[str, Any] | None = None,
    *,
    source_dir: str | Path | None = None,
    resources: Iterable[str | Path] = (),
    setup_service: SetupService | None = None,
    **kwargs: Any,
) -> Any:
    """Generate a source bundle and delegate setup to the existing service.

    ``resources`` is passed explicitly to the service/runner; no setup logic
    is duplicated here.  A temporary source bundle is used when ``source_dir``
    is omitted and is removed after the service returns.
    """

    service = setup_service or SetupService()
    if source_dir is not None:
        bundle = write_setup_source(atoms, source_dir)
        return service.run(
            bundle.source_dir,
            destination_dir,
            spec,
            resources=tuple(resources),
            **kwargs,
        )
    with tempfile.TemporaryDirectory(prefix="rsutil-ase-source-") as temporary:
        bundle = write_setup_source(atoms, temporary)
        return service.run(
            bundle.source_dir,
            destination_dir,
            spec,
            resources=tuple(resources),
            **kwargs,
        )


__all__ = [
    "ASEDependencyError",
    "ASESetupSource",
    "ase_to_rs_structure",
    "create_setup_source",
    "from_ase",
    "rs_structure_to_ase",
    "run_setup_from_ase",
    "to_ase",
    "write_setup_source",
]
