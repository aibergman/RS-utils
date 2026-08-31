"""Versioned RS-LMTO real-space tight-binding models and ML export.

The versioned text format implemented here is the ``fable_v4``
``rs-tb-meta-v1``/``rs-tb-hoppings-v1`` pair.  Legacy ``rs2txt.dat`` data is
kept in the bands package; conversion is deliberately refused unless the
legacy displacement can be mapped unambiguously to a periodic ``(R, i, j)``
identity.
"""

from __future__ import annotations

import io
import json
import math
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from rsutil.core import RSStructure, StructureSite
from rsutil.core.hashing import hash_file

from .errors import (
    MalformedResultError,
    MissingResultFileError,
    UnsupportedResultFormatError,
)
from .models import ParserProvenance, SourceFile
from .units import ELECTRONVOLT, Quantity, Unit


TB_MODEL_SCHEMA = "rsutil-tb-model-v1"
TB_META_SCHEMA = "rs-tb-meta-v1"
TB_HOPPING_SCHEMA = "rs-tb-hoppings-v1"
ML_SCHEMA = "rsutil-tb-ml-v1"


def _as_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _as_json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_as_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_as_json(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _as_json(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, Unit):
        return value.value
    return value


def _finite_complex(value: complex, name: str) -> complex:
    value = complex(value)
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True, order=True)
class OrbitalMetadata:
    """One local orbital in the exported basis.

    ``orbital_index`` is zero-based within ``site_index``.  The fable_v4 text
    format stores one-based local and global orbital indices; parsing converts
    only the local identity to the canonical zero-based representation.
    """

    site_index: int
    orbital_index: int
    label: str | None = None
    spin_channel: str | None = None
    scalar_orbital_index: int | None = None

    def __post_init__(self) -> None:
        if self.site_index < 0 or self.orbital_index < 0:
            raise ValueError("orbital indices must be non-negative")
        if self.scalar_orbital_index is not None and self.scalar_orbital_index < 0:
            raise ValueError("scalar_orbital_index must be non-negative")

    def to_dict(self) -> dict[str, Any]:
        return {
            "site_index": self.site_index,
            "orbital_index": self.orbital_index,
            "label": self.label,
            "spin_channel": self.spin_channel,
            "scalar_orbital_index": self.scalar_orbital_index,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "OrbitalMetadata":
        return cls(
            site_index=int(data["site_index"]),
            orbital_index=int(data["orbital_index"]),
            label=data.get("label"),
            spin_channel=data.get("spin_channel"),
            scalar_orbital_index=(
                None
                if data.get("scalar_orbital_index") is None
                else int(data["scalar_orbital_index"])
            ),
        )


@dataclass(frozen=True, order=True)
class HoppingTerm:
    """One complex matrix element of ``H(R)``.

    The target site ``site_j`` is in the cell translated by ``R`` relative to
    the reference cell.  Thus the identity is independent of text line order:
    ``(R, site_i, site_j, orbital_i, orbital_j)``.
    """

    R: tuple[int, int, int]
    site_i: int
    site_j: int
    orbital_i: int
    orbital_j: int
    value: complex = field(compare=False)
    global_orbital_i: int | None = field(default=None, compare=False)
    global_orbital_j: int | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        if len(self.R) != 3 or any(int(value) != value for value in self.R):
            raise ValueError("R must contain three integer lattice translations")
        object.__setattr__(self, "R", tuple(int(value) for value in self.R))
        for name in ("site_i", "site_j", "orbital_i", "orbital_j"):
            if int(getattr(self, name)) < 0:
                raise ValueError(f"{name} must be non-negative")
        object.__setattr__(self, "value", _finite_complex(self.value, "hopping value"))

    def identity(self) -> tuple[int, int, int, int, int, int, int]:
        return (*self.R, self.site_i, self.site_j, self.orbital_i, self.orbital_j)

    def to_dict(self) -> dict[str, Any]:
        return {
            "R": list(self.R),
            "site_i": self.site_i,
            "site_j": self.site_j,
            "orbital_i": self.orbital_i,
            "orbital_j": self.orbital_j,
            "value": [self.value.real, self.value.imag],
            "global_orbital_i": self.global_orbital_i,
            "global_orbital_j": self.global_orbital_j,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "HoppingTerm":
        value = data["value"]
        if isinstance(value, (list, tuple)) and len(value) == 2:
            value = complex(float(value[0]), float(value[1]))
        return cls(
            R=tuple(int(item) for item in data["R"]),
            site_i=int(data["site_i"]),
            site_j=int(data["site_j"]),
            orbital_i=int(data["orbital_i"]),
            orbital_j=int(data["orbital_j"]),
            value=value,
            global_orbital_i=data.get("global_orbital_i"),
            global_orbital_j=data.get("global_orbital_j"),
        )


def _structure_from_lattice_reader(lattice_reader: Any) -> RSStructure:
    """Convert the existing lightweight ``LatticeReader`` contract."""

    # LatticeReader stores a(:, 1:3) as columns; RSStructure stores vectors
    # as rows, matching the fable metadata and public spglib convention.
    vectors = np.asarray(lattice_reader.lattice_vectors, dtype=float).T
    positions = np.asarray(lattice_reader.atomic_positions, dtype=float)
    # ``LatticeReader.atomic_numbers`` is the legacy ``izp``/type array, not
    # guaranteed chemical atomic numbers (e.g. Fe/Co is encoded as 1/2).
    numbers = list(np.asarray(lattice_reader.atomic_numbers).astype(int))
    if vectors.shape != (3, 3) or positions.ndim != 2 or positions.shape[1] != 3:
        raise MalformedResultError("legacy lattice reader lacks complete geometry")
    if len(positions) != int(lattice_reader.nbas) or len(numbers) != len(positions):
        raise MalformedResultError("legacy lattice reader geometry counts are inconsistent")
    sites = []
    for index, (position, atomic_number) in enumerate(zip(positions, numbers)):
        # Legacy lattice.nml does not provide a reliable element symbol here;
        # preserve the type number as metadata and keep ML identity honest.
        sites.append(StructureSite(index, tuple(position), "X", None, metadata={"legacy_type": atomic_number}))
    return RSStructure(
        lattice_vectors=tuple(tuple(vector) for vector in vectors),
        length_unit="unknown",
        coordinate_convention="cartesian",
        sites=tuple(sites),
        metadata={"source": str(getattr(lattice_reader, "filename", ""))},
    )


@dataclass(frozen=True)
class TightBindingModel:
    """Typed, sparse/list representation of a real-space TB Hamiltonian."""

    structure: RSStructure
    orbitals: tuple[OrbitalMetadata, ...]
    hoppings: tuple[HoppingTerm, ...]
    energy_unit: Unit | str = ELECTRONVOLT
    spin_representation: str = "unknown"
    soc_representation: str = "not_encoded"
    control_nsp: int | None = None
    soc_enabled: bool | None = None
    gbt_compatible: bool | None = None
    magnetic_representation: str = "unknown"
    fermi_energy: Quantity | None = None
    provenance: ParserProvenance | Mapping[str, Any] | None = None
    overlap: tuple[HoppingTerm, ...] | None = None
    schema_version: str = TB_MODEL_SCHEMA

    def __post_init__(self) -> None:
        if self.schema_version != TB_MODEL_SCHEMA:
            raise ValueError(f"unsupported TB model schema: {self.schema_version}")
        unit = self.energy_unit if isinstance(self.energy_unit, Unit) else Unit(str(self.energy_unit))
        object.__setattr__(self, "energy_unit", unit)
        orbitals = tuple(sorted(self.orbitals))
        hoppings = tuple(sorted(self.hoppings, key=lambda term: term.identity()))
        if not orbitals:
            raise ValueError("a tight-binding model needs at least one orbital")
        if len({(item.site_index, item.orbital_index) for item in orbitals}) != len(orbitals):
            raise ValueError("orbital identities must be unique")
        site_count = len(self.structure.sites)
        valid_orbitals = {(item.site_index, item.orbital_index) for item in orbitals}
        for term in hoppings:
            if term.site_i >= site_count or term.site_j >= site_count:
                raise ValueError("hopping site index is outside the structure")
            if (term.site_i, term.orbital_i) not in valid_orbitals or (term.site_j, term.orbital_j) not in valid_orbitals:
                raise ValueError("hopping orbital index is outside the basis")
        if self.overlap is not None:
            object.__setattr__(self, "overlap", tuple(sorted(self.overlap, key=lambda term: term.identity())))
        object.__setattr__(self, "orbitals", orbitals)
        object.__setattr__(self, "hoppings", hoppings)
        if isinstance(self.provenance, Mapping):
            object.__setattr__(self, "provenance", dict(self.provenance))
        if self.control_nsp is not None and self.control_nsp not in {1, 2, 3, 4}:
            raise ValueError("control%nsp must be one of 1, 2, 3, or 4")
        if self.control_nsp is not None:
            expected_soc = self.control_nsp in {2, 4}
            if self.soc_enabled is not None and self.soc_enabled != expected_soc:
                raise ValueError("soc_enabled disagrees with control%nsp")
            object.__setattr__(self, "soc_enabled", expected_soc)
            object.__setattr__(self, "gbt_compatible", self.control_nsp == 1)

    @property
    def orbitals_per_site(self) -> int:
        counts = [sum(item.site_index == site for item in self.orbitals) for site in range(len(self.structure.sites))]
        if not counts or len(set(counts)) != 1:
            raise ValueError("H(k) reconstruction requires an equal orbital count per site")
        return counts[0]

    @property
    def dimension(self) -> int:
        return len(self.orbitals)

    def hamiltonian(self, k_point: Sequence[float]) -> np.ndarray:
        """Reconstruct ``H(k)`` for fractional reciprocal coordinates."""

        k = np.asarray(tuple(k_point), dtype=float)
        if k.shape != (3,) or not np.all(np.isfinite(k)):
            raise ValueError("k_point must contain three finite fractional coordinates")
        orbital_count = self.orbitals_per_site
        matrix = np.zeros((self.dimension, self.dimension), dtype=np.complex128)
        for term in self.hoppings:
            row = term.site_i * orbital_count + term.orbital_i
            column = term.site_j * orbital_count + term.orbital_j
            phase = np.exp(2j * np.pi * float(np.dot(k, term.R)))
            matrix[row, column] += term.value * phase
        return matrix

    def eigenvalues(self, k_point: Sequence[float]) -> np.ndarray:
        """Return sorted band energies reconstructed from ``H(k)``."""

        matrix = self.hamiltonian(k_point)
        if np.allclose(matrix, matrix.conj().T):
            return np.linalg.eigvalsh(matrix)
        return np.sort(np.linalg.eigvals(matrix).real)

    @classmethod
    def from_fable_v4(cls, metadata_path: str | Path, hopping_path: str | Path | None = None, **kwargs: Any) -> "TightBindingModel":
        return read_fable_v4_tb(metadata_path, hopping_path, **kwargs)

    @classmethod
    def from_ml_dataset(cls, npz_path: str | Path, metadata_path: str | Path | None = None) -> "TightBindingModel":
        return tight_binding_from_ml_dataset(npz_path, metadata_path)

    def to_dict(self, *, include_hoppings: bool = True) -> dict[str, Any]:
        provenance = self.provenance
        if isinstance(provenance, ParserProvenance):
            provenance = provenance.to_dict()
        result = {
            "schema_version": self.schema_version,
            "structure": self.structure.to_dict(),
            "orbitals": [item.to_dict() for item in self.orbitals],
            "energy_unit": self.energy_unit.value,
            "spin_representation": self.spin_representation,
            "soc_representation": self.soc_representation,
            "control_nsp": self.control_nsp,
            "soc_enabled": self.soc_enabled,
            "gbt_compatible": self.gbt_compatible,
            "magnetic_representation": self.magnetic_representation,
            "fermi_energy": None if self.fermi_energy is None else self.fermi_energy.to_dict(),
            "provenance": None if provenance is None else _as_json(provenance),
            "overlap_present": self.overlap is not None,
        }
        if include_hoppings:
            result["hoppings"] = [item.to_dict() for item in self.hoppings]
        return result


def _read_metadata(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != f"# {TB_META_SCHEMA}":
        raise UnsupportedResultFormatError(f"{path} is not {TB_META_SCHEMA}")
    metadata: dict[str, Any] = {"lattice_vectors": [], "positions": []}
    mode = None
    for line_number, raw in enumerate(lines[1:], start=2):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        key = parts[0].lower()
        values = parts[1:]
        if key == "lattice_vectors_cart":
            mode = "lattice"
            continue
        if key == "positions_cart":
            mode = "positions"
            continue
        if mode == "lattice" and len(parts) == 3:
            try:
                metadata["lattice_vectors"].append(tuple(float(value) for value in parts))
                continue
            except ValueError:
                pass
        if mode == "positions" and len(values) == 3:
            try:
                index = int(key)
                metadata["positions"].append((index, tuple(float(value) for value in values)))
                continue
            except ValueError:
                pass
        mode = None
        if key in {"n_atoms", "norb_scalar", "block_size", "index_base"} and len(values) == 1:
            metadata[key] = int(values[0])
        elif key in {"energy_unit", "length_unit", "spin_representation", "soc_representation"} and len(values) == 1:
            metadata[key] = values[0]
        elif key in {"nsp", "control_nsp"} and len(values) == 1:
            metadata["control_nsp"] = int(values[0])
        elif key in {"fermi_energy", "fermi_energy_ev"} and values:
            if len(values) > 2:
                raise MalformedResultError(f"malformed Fermi energy on metadata line {line_number}")
            metadata["fermi_energy"] = Quantity(float(values[0]), values[1] if len(values) == 2 else ELECTRONVOLT)
        else:
            raise MalformedResultError(f"unsupported TB metadata field {key!r} on line {line_number}")
    required = {"n_atoms", "norb_scalar", "block_size", "index_base", "energy_unit", "length_unit"}
    missing = sorted(required.difference(metadata))
    if missing or len(metadata["lattice_vectors"]) != 3 or len(metadata["positions"]) != metadata.get("n_atoms", -1):
        raise MalformedResultError(f"incomplete TB metadata: {', '.join(missing) or 'geometry/counts'}")
    if metadata["index_base"] != 1:
        raise UnsupportedResultFormatError("only one-based fable_v4 TB metadata is supported")
    return metadata


def _structure_from_metadata(metadata: Mapping[str, Any], lattice: Any = None) -> RSStructure:
    if lattice is not None:
        if isinstance(lattice, RSStructure):
            return lattice
        return _structure_from_lattice_reader(lattice)
    positions = sorted(metadata["positions"], key=lambda item: item[0])
    expected_indices = list(range(1, len(positions) + 1))
    if [index for index, _position in positions] != expected_indices:
        raise MalformedResultError("TB position indices must be one-based and contiguous")
    sites = tuple(
        StructureSite(index, position, "X", None, metadata={"identity_only": True})
        for index, (_source_index, position) in enumerate(positions)
    )
    return RSStructure(
        lattice_vectors=tuple(metadata["lattice_vectors"]),
        length_unit=str(metadata["length_unit"]),
        coordinate_convention="cartesian",
        sites=sites,
        metadata={"species_metadata_available": False, "format": TB_META_SCHEMA},
    )


def _orbital_metadata(metadata: Mapping[str, Any]) -> tuple[OrbitalMetadata, ...]:
    n_atoms = int(metadata["n_atoms"])
    scalar_count = int(metadata["norb_scalar"])
    block_size = int(metadata["block_size"])
    if scalar_count <= 0 or block_size <= 0:
        raise MalformedResultError("TB orbital counts must be positive")
    spinor = block_size == 2 * scalar_count
    result = []
    for site in range(n_atoms):
        for orbital in range(block_size):
            channel = None
            scalar_orbital = orbital
            if spinor:
                channel = "up" if orbital < scalar_count else "down"
                scalar_orbital = orbital % scalar_count
            result.append(OrbitalMetadata(site, orbital, f"orbital_{orbital + 1}", channel, scalar_orbital))
    return tuple(result)


def control_nsp_from_input(input_path: str | Path) -> int:
    """Read the RS-LMTO ``control%nsp`` selector from ``input.nml``."""

    path = Path(input_path)
    if not path.is_file():
        raise MissingResultFileError(f"RS-LMTO input file not found: {path}")
    text = re.sub(r"!.*$", "", path.read_text(encoding="utf-8"), flags=re.MULTILINE)
    control = re.search(
        r"&control\b(?P<body>.*?)(?:^\s*/|\Z)",
        text,
        flags=re.IGNORECASE | re.DOTALL | re.MULTILINE,
    )
    if control is None:
        raise MalformedResultError(f"{path} does not contain a control namelist")
    match = re.search(r"\bnsp\s*=\s*([1-4])\b", control.group("body"), flags=re.IGNORECASE)
    if match is None:
        raise MalformedResultError(f"{path} does not define control%nsp")
    return int(match.group(1))


def _spin_metadata(nsp: int | None, block_size: int, scalar_count: int) -> tuple[str, str, bool | None, str]:
    if nsp == 1:
        return "collinear_magnetic", "disabled", False, "collinear_magnetic"
    if nsp == 2:
        return "collinear_magnetic", "enabled", True, "collinear_magnetic"
    if nsp == 4:
        return "noncollinear_spinor", "enabled", True, "noncollinear"
    if nsp == 3:
        return "unknown", "not_encoded", None, "unknown"
    return (
        "spinor" if block_size == 2 * scalar_count else "scalar",
        "not_encoded",
        None,
        "unknown",
    )


def read_fable_v4_tb(
    metadata_path: str | Path,
    hopping_path: str | Path | None = None,
    *,
    lattice: Any = None,
    input_path: str | Path | None = None,
    nsp: int | None = None,
) -> TightBindingModel:
    """Read the exact ``fable_v4`` versioned metadata/hopping pair."""

    metadata_file = Path(metadata_path)
    if not metadata_file.is_file():
        raise MissingResultFileError(f"TB metadata file not found: {metadata_file}")
    hopping_file = Path(hopping_path) if hopping_path is not None else metadata_file.with_suffix(".tb")
    if not hopping_file.is_file():
        raise MissingResultFileError(f"TB hopping file not found: {hopping_file}")
    metadata = _read_metadata(metadata_file)
    input_nsp = control_nsp_from_input(input_path) if input_path is not None else None
    selected_nsp = nsp if nsp is not None else (input_nsp if input_nsp is not None else metadata.get("control_nsp"))
    if selected_nsp is not None:
        selected_nsp = int(selected_nsp)
        if selected_nsp not in {1, 2, 3, 4}:
            raise MalformedResultError("control%nsp must be one of 1, 2, 3, or 4")
    structure = _structure_from_metadata(metadata, lattice)
    if len(structure.sites) != metadata["n_atoms"]:
        raise MalformedResultError("TB metadata/site structure count mismatch")
    hopping_lines = hopping_file.read_text(encoding="utf-8").splitlines()
    first_content = next((line.strip() for line in hopping_lines if line.strip()), "")
    if first_content != f"# {TB_HOPPING_SCHEMA}":
        raise UnsupportedResultFormatError(f"{hopping_file} is not {TB_HOPPING_SCHEMA}")
    terms: list[HoppingTerm] = []
    for line_number, raw in enumerate(hopping_lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            if line == f"# {TB_HOPPING_SCHEMA}":
                continue
            continue
        parts = line.split()
        if len(parts) != 12:
            raise MalformedResultError(f"expected 12 TB hopping fields on line {line_number}")
        try:
            R = tuple(int(item) for item in parts[0:3])
            site_i, site_j, _k_neigh, local_i, local_j, global_i, global_j = (int(item) for item in parts[3:10])
            value = complex(float(parts[10]), float(parts[11]))
        except ValueError as exc:
            raise MalformedResultError(f"malformed TB hopping on line {line_number}") from exc
        if not (1 <= site_i <= metadata["n_atoms"] and 1 <= site_j <= metadata["n_atoms"]):
            raise MalformedResultError(f"TB site index out of range on line {line_number}")
        if not (1 <= local_i <= metadata["block_size"] and 1 <= local_j <= metadata["block_size"]):
            raise MalformedResultError(f"TB local orbital index out of range on line {line_number}")
        terms.append(HoppingTerm(R, site_i - 1, site_j - 1, local_i - 1, local_j - 1, value, global_i, global_j))
    if len({term.identity() for term in terms}) != len(terms):
        raise MalformedResultError("duplicate TB hopping identity")
    block_size = metadata["block_size"]
    spin_representation, soc_representation, soc_enabled, magnetic_representation = _spin_metadata(
        selected_nsp, block_size, metadata["norb_scalar"]
    )
    spin_representation = metadata.get("spin_representation", spin_representation)
    soc_representation = metadata.get("soc_representation", soc_representation)
    source_files = [
        SourceFile(str(metadata_file), hash_file(metadata_file), "tb_metadata"),
        SourceFile(str(hopping_file), hash_file(hopping_file), "tb_hoppings"),
    ]
    if input_path is not None:
        input_file = Path(input_path)
        source_files.append(SourceFile(str(input_file), hash_file(input_file), "input_namelist"))
    provenance = ParserProvenance(
        source_path=str(hopping_file),
        content_hash=hash_file(hopping_file),
        format="rs-lmto-fable-v4-tb",
        format_version=f"{TB_META_SCHEMA}+{TB_HOPPING_SCHEMA}",
        source_files=tuple(source_files),
    )
    return TightBindingModel(
        structure=structure,
        orbitals=_orbital_metadata(metadata),
        hoppings=tuple(terms),
        energy_unit=metadata["energy_unit"],
        spin_representation=spin_representation,
        soc_representation=soc_representation,
        control_nsp=selected_nsp,
        soc_enabled=soc_enabled,
        gbt_compatible=(selected_nsp == 1) if selected_nsp is not None else None,
        magnetic_representation=magnetic_representation,
        fermi_energy=metadata.get("fermi_energy"),
        provenance=provenance,
    )


def _match_legacy_translation(structure: RSStructure, site_i: int, displacement: Sequence[float], tolerance: float = 1e-6) -> tuple[int, int, int, int] | None:
    """Return the unique target site and R matching a legacy Cartesian vector."""

    if structure.coordinate_convention.value != "cartesian":
        raise MalformedResultError("legacy TB conversion requires Cartesian site positions")
    vectors = np.asarray(structure.lattice_vectors, dtype=float)
    origin = np.asarray(structure.sites[site_i].position, dtype=float)
    displacement = np.asarray(displacement, dtype=float)
    matches = []
    for site_j, site in enumerate(structure.sites):
        position = np.asarray(site.position, dtype=float)
        for r1 in range(-5, 6):
            for r2 in range(-5, 6):
                for r3 in range(-5, 6):
                    R = np.array((r1, r2, r3), dtype=float)
                    expected = origin - (position + R @ vectors)
                    if np.linalg.norm(expected - displacement) <= tolerance:
                        matches.append((site_j, (r1, r2, r3)))
    if len(matches) != 1:
        return None
    return (*matches[0][1], matches[0][0])


def tight_binding_from_legacy_reader(hamiltonian_reader: Any, lattice_reader: Any) -> TightBindingModel:
    """Convert a legacy ``rs2txt.dat`` reader only with provable periodic identity."""

    structure = _structure_from_lattice_reader(lattice_reader)
    orbital_count = int(hamiltonian_reader.orbital_size)
    orbitals = tuple(
        OrbitalMetadata(site, orbital, f"orbital_{orbital + 1}")
        for site in range(len(structure.sites))
        for orbital in range(orbital_count)
    )
    terms = []
    for source_type, neighbors in hamiltonian_reader.hamiltonian_data.items():
        site_i = int(source_type) - 1
        if not 0 <= site_i < len(structure.sites):
            raise MalformedResultError("legacy TB source site cannot be established")
        for _neighbor, data in neighbors.items():
            match = _match_legacy_translation(structure, site_i, data["displacement"])
            if match is None:
                raise UnsupportedResultFormatError(
                    "legacy rs2txt displacement does not establish one periodic site/R identity"
                )
            R = match[:3]
            site_j = match[3]
            matrix = np.asarray(data["matrix"], dtype=np.complex128)
            if matrix.shape != (orbital_count, orbital_count):
                raise MalformedResultError("legacy Hamiltonian matrix shape is inconsistent")
            for orbital_i in range(orbital_count):
                for orbital_j in range(orbital_count):
                    if matrix[orbital_i, orbital_j] != 0:
                        terms.append(HoppingTerm(R, site_i, site_j, orbital_i, orbital_j, matrix[orbital_i, orbital_j]))
    provenance = ParserProvenance(
        source_path=str(hamiltonian_reader.filename),
        content_hash=hash_file(hamiltonian_reader.filename),
        format="rs-lmto-legacy-rs2txt",
        format_version=None,
        source_files=(SourceFile(str(hamiltonian_reader.filename), hash_file(hamiltonian_reader.filename), "legacy_hamiltonian"),),
    )
    return TightBindingModel(
        structure=structure,
        orbitals=orbitals,
        hoppings=tuple(terms),
        energy_unit="unknown",
        spin_representation="spinor" if orbital_count == 18 else "unknown",
        soc_representation="not_encoded",
        control_nsp=None,
        soc_enabled=None,
        gbt_compatible=None,
        magnetic_representation="unknown",
        provenance=provenance,
    )


def _output_paths(output_prefix: str | Path) -> tuple[Path, Path]:
    prefix = Path(output_prefix)
    if prefix.suffix == ".npz":
        return prefix, prefix.with_suffix(".json")
    return prefix.with_suffix(".npz"), prefix.with_suffix(".json")


def _npy_bytes(array: np.ndarray) -> bytes:
    stream = io.BytesIO()
    np.lib.format.write_array(stream, np.asarray(array), allow_pickle=False)
    return stream.getvalue()


def export_ml_dataset(model: TightBindingModel, output_prefix: str | Path) -> tuple[Path, Path]:
    """Write deterministic compressed NPZ arrays and JSON metadata.

    Native ``complex128`` values are stored in the NPZ; JSON contains only
    schema, structure, basis, units and provenance metadata.
    """

    npz_path, json_path = _output_paths(output_prefix)
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    terms = model.hoppings
    arrays = {
        "R": np.asarray([term.R for term in terms], dtype=np.int64).reshape((-1, 3)),
        "site_i": np.asarray([term.site_i for term in terms], dtype=np.int64),
        "site_j": np.asarray([term.site_j for term in terms], dtype=np.int64),
        "orbital_i": np.asarray([term.orbital_i for term in terms], dtype=np.int64),
        "orbital_j": np.asarray([term.orbital_j for term in terms], dtype=np.int64),
        "value": np.asarray([term.value for term in terms], dtype=np.complex128),
    }
    with zipfile.ZipFile(npz_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in ("R", "site_i", "site_j", "orbital_i", "orbital_j", "value"):
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, _npy_bytes(arrays[name]))
    metadata = {
        "schema_version": ML_SCHEMA,
        "arrays": {
            "R": {"dtype": "int64", "shape": [len(terms), 3], "meaning": "target-cell translation"},
            "site_i": {"dtype": "int64", "meaning": "zero-based reference-cell row site"},
            "site_j": {"dtype": "int64", "meaning": "zero-based reference-cell column site"},
            "orbital_i": {"dtype": "int64", "meaning": "zero-based local row orbital"},
            "orbital_j": {"dtype": "int64", "meaning": "zero-based local column orbital"},
            "value": {"dtype": "complex128", "meaning": "H(R) matrix element"},
        },
        "complex_encoding": "native_numpy_complex128",
        "term_count": len(terms),
        "model": model.to_dict(include_hoppings=False),
    }
    json_path.write_text(json.dumps(_as_json(metadata), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return npz_path, json_path


def tight_binding_from_ml_dataset(npz_path: str | Path, metadata_path: str | Path | None = None) -> TightBindingModel:
    """Reconstruct a model from the dependency-light ML dataset pair."""

    npz_file = Path(npz_path)
    json_file = Path(metadata_path) if metadata_path is not None else npz_file.with_suffix(".json")
    if not npz_file.is_file() or not json_file.is_file():
        raise MissingResultFileError("ML TB dataset requires both .npz and .json files")
    metadata = json.loads(json_file.read_text(encoding="utf-8"))
    if metadata.get("schema_version") != ML_SCHEMA:
        raise UnsupportedResultFormatError("unsupported ML TB dataset schema")
    with np.load(npz_file, allow_pickle=False) as arrays:
        required = ("R", "site_i", "site_j", "orbital_i", "orbital_j", "value")
        if any(name not in arrays for name in required):
            raise MalformedResultError("ML TB dataset is missing an identity/value array")
        count = int(metadata.get("term_count", -1))
        if count < 0 or any(len(arrays[name]) != count for name in required if name != "R") or arrays["R"].shape != (count, 3):
            raise MalformedResultError("ML TB array lengths do not agree")
        terms = tuple(
            HoppingTerm(
                tuple(int(item) for item in arrays["R"][index]),
                int(arrays["site_i"][index]),
                int(arrays["site_j"][index]),
                int(arrays["orbital_i"][index]),
                int(arrays["orbital_j"][index]),
                complex(arrays["value"][index]),
            )
            for index in range(count)
        )
    model_data = metadata.get("model")
    if not isinstance(model_data, Mapping):
        raise MalformedResultError("ML TB metadata lacks model metadata")
    provenance = model_data.get("provenance")
    if isinstance(provenance, Mapping) and "source_path" in provenance:
        provenance = ParserProvenance.from_dict(provenance)
    fermi = model_data.get("fermi_energy")
    return TightBindingModel(
        structure=RSStructure.from_dict(model_data["structure"]),
        orbitals=tuple(OrbitalMetadata.from_dict(item) for item in model_data["orbitals"]),
        hoppings=terms,
        energy_unit=model_data["energy_unit"],
        spin_representation=model_data.get("spin_representation", "unknown"),
        soc_representation=model_data.get("soc_representation", "not_encoded"),
        control_nsp=model_data.get("control_nsp"),
        soc_enabled=model_data.get("soc_enabled"),
        gbt_compatible=model_data.get("gbt_compatible"),
        magnetic_representation=model_data.get("magnetic_representation", "unknown"),
        fermi_energy=None if fermi is None else Quantity.from_dict(fermi),
        provenance=provenance,
    )


def load_tight_binding(source: str | Path, *, lattice: Any = None) -> TightBindingModel:
    """Load a versioned fable_v4 pair or a previously exported ML dataset."""

    path = Path(source)
    if path.suffix == ".npz":
        return tight_binding_from_ml_dataset(path)
    if path.suffix == ".json":
        raise UnsupportedResultFormatError("pass the ML TB .npz file as the source")
    if path.suffix == ".tb":
        return read_fable_v4_tb(path.with_suffix(".meta"), path, lattice=lattice)
    if path.suffix == ".meta":
        return read_fable_v4_tb(path, path.with_suffix(".tb"), lattice=lattice)
    raise UnsupportedResultFormatError(f"cannot detect a TB format from {path}")


read_tight_binding = load_tight_binding
export_tight_binding_ml = export_ml_dataset


__all__ = [
    "HoppingTerm",
    "ML_SCHEMA",
    "OrbitalMetadata",
    "TB_HOPPING_SCHEMA",
    "TB_META_SCHEMA",
    "TB_MODEL_SCHEMA",
    "TightBindingModel",
    "export_ml_dataset",
    "export_tight_binding_ml",
    "load_tight_binding",
    "read_fable_v4_tb",
    "read_tight_binding",
    "tight_binding_from_legacy_reader",
    "tight_binding_from_ml_dataset",
    "control_nsp_from_input",
]
