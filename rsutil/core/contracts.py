"""Small, dependency-light data contracts used by HT services."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Sequence

from .errors import ValidationError


class CoordinateConvention(str, Enum):
    CARTESIAN = "cartesian"
    FRACTIONAL = "fractional"


class DiagnosticSeverity(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


class OperationStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"


class ErrorCode(str, Enum):
    INVALID_INPUT = "invalid_input"
    VALIDATION_FAILED = "validation_failed"
    MISSING_INPUT = "missing_input"
    WORKSPACE_ERROR = "workspace_error"
    WORKER_FAILED = "worker_failed"
    OUTPUT_VALIDATION_FAILED = "output_validation_failed"
    PARSE_ERROR = "parse_error"
    EXCHANGE_ERROR = "exchange_error"
    UNKNOWN = "unknown"


def _vector(value: Sequence[float], name: str) -> tuple[float, float, float]:
    if isinstance(value, (str, bytes)) or len(value) != 3:
        raise ValidationError(f"{name} must contain exactly three coordinates")
    result = tuple(float(item) for item in value)
    if not all(math.isfinite(item) for item in result):
        raise ValidationError(f"{name} must contain only finite numbers")
    return result  # type: ignore[return-value]


def _mapping(value: Mapping[str, Any] | None) -> dict[str, Any]:
    return dict(value or {})


_ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh "
    "Fl Mc Lv Ts Og"
).split()
_ATOMIC_NUMBERS = {symbol: number for number, symbol in enumerate(_ELEMENTS, start=1)}


@dataclass(frozen=True)
class StructureSite:
    """One site in the canonical, stable ordering of an ``RSStructure``.

    Site indices are zero-based and must be contiguous in ``RSStructure``.
    ``element`` is the chemical identity; RS-specific labels/types are kept
    separately so labels such as ``Fe1`` do not become chemical symbols.
    """

    index: int
    position: tuple[float, float, float]
    element: str
    atomic_number: int | None = None
    rs_label: str | None = None
    rs_type: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.index, int) or isinstance(self.index, bool) or self.index < 0:
            raise ValidationError("site index must be a non-negative integer")
        object.__setattr__(self, "position", _vector(self.position, "site position"))
        if not isinstance(self.element, str) or not self.element.strip():
            raise ValidationError("site element must be a non-empty string")
        element = self.element.strip()
        if not re.fullmatch(r"[A-Za-z][A-Za-z]{0,2}", element):
            raise ValidationError(f"invalid chemical element symbol: {self.element!r}")
        canonical_element = element[0].upper() + element[1:].lower()
        object.__setattr__(self, "element", canonical_element)
        if self.atomic_number is not None:
            if (
                not isinstance(self.atomic_number, int)
                or isinstance(self.atomic_number, bool)
                or not 1 <= self.atomic_number <= 118
            ):
                raise ValidationError("atomic number must be an integer in [1, 118]")
            expected = _ATOMIC_NUMBERS.get(canonical_element)
            if expected is not None and expected != self.atomic_number:
                raise ValidationError(
                    f"atomic number {self.atomic_number} does not match {canonical_element}"
                )
        object.__setattr__(self, "metadata", _mapping(self.metadata))

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "position": list(self.position),
            "element": self.element,
            "atomic_number": self.atomic_number,
            "rs_label": self.rs_label,
            "rs_type": self.rs_type,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "StructureSite":
        return cls(
            index=int(data["index"]),
            position=data["position"],
            element=str(data["element"]),
            atomic_number=data.get("atomic_number"),
            rs_label=data.get("rs_label"),
            rs_type=data.get("rs_type"),
            metadata=data.get("metadata"),
        )


Site = StructureSite


@dataclass(frozen=True)
class RSStructure:
    """Canonical RS structure independent of ASE and pymatgen."""

    lattice_vectors: tuple[tuple[float, float, float], ...]
    length_unit: str
    sites: tuple[StructureSite, ...]
    coordinate_convention: CoordinateConvention | str = CoordinateConvention.CARTESIAN
    pbc: tuple[bool, bool, bool] = (True, True, True)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.length_unit, str) or not self.length_unit.strip():
            raise ValidationError("length_unit must be a non-empty string")
        vectors = tuple(_vector(vector, "lattice vector") for vector in self.lattice_vectors)
        if len(vectors) != 3:
            raise ValidationError("lattice_vectors must contain exactly three vectors")
        object.__setattr__(self, "lattice_vectors", vectors)
        try:
            convention = CoordinateConvention(self.coordinate_convention)
        except ValueError as exc:
            raise ValidationError("coordinate_convention must be cartesian or fractional") from exc
        object.__setattr__(self, "coordinate_convention", convention)
        if len(self.pbc) != 3:
            raise ValidationError("pbc must contain exactly three booleans")
        if any(not isinstance(value, bool) for value in self.pbc):
            raise ValidationError("pbc must contain booleans")
        pbc = tuple(self.pbc)
        object.__setattr__(self, "pbc", pbc)
        sites = tuple(
            site if isinstance(site, StructureSite) else StructureSite.from_dict(site)
            for site in self.sites
        )
        if tuple(site.index for site in sites) != tuple(range(len(sites))):
            raise ValidationError("site indices must be stable, zero-based, and contiguous")
        object.__setattr__(self, "sites", sites)
        object.__setattr__(self, "metadata", _mapping(self.metadata))

    @property
    def positions(self) -> tuple[tuple[float, float, float], ...]:
        return tuple(site.position for site in self.sites)

    @property
    def elements(self) -> tuple[str, ...]:
        return tuple(site.element for site in self.sites)

    @property
    def atomic_numbers(self) -> tuple[int | None, ...]:
        return tuple(site.atomic_number for site in self.sites)

    def to_dict(self) -> dict[str, Any]:
        return {
            "lattice_vectors": [list(vector) for vector in self.lattice_vectors],
            "length_unit": self.length_unit,
            "coordinate_convention": self.coordinate_convention.value,
            "pbc": list(self.pbc),
            "sites": [site.to_dict() for site in self.sites],
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "RSStructure":
        return cls(
            lattice_vectors=data["lattice_vectors"],
            length_unit=str(data["length_unit"]),
            coordinate_convention=data.get("coordinate_convention", "cartesian"),
            pbc=data.get("pbc", (True, True, True)),
            sites=tuple(StructureSite.from_dict(site) for site in data["sites"]),
            metadata=data.get("metadata"),
        )


@dataclass(frozen=True)
class SetupSpec:
    """Scientific setup parameters exposed by the legacy ``-S`` workflow."""

    cluster_type: str = "periodic"
    empty_spheres: bool = False
    min_empty_radius: float = 1.2
    max_empty_radius: float = 2.5
    lmx: int = 2
    rec_method: str = "block"
    rec_steps: int = 20
    lattice_dist: float = 0.0
    basis_dist: float = 0.0
    alat_scale: float = 1.0
    seed: int | None = None
    mode_2d: bool = False
    z_iso_factor: float | None = None
    supercell: tuple[int, int, int] = (1, 1, 1)
    sqs: str | None = None
    sqs_shells: int = 3
    sqs_steps: int = 5000
    sqs_seed: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.cluster_type, str) or self.cluster_type not in {"periodic", "free"}:
            raise ValidationError("cluster_type must be 'periodic' or 'free'")
        if not isinstance(self.rec_method, str) or self.rec_method.lower() not in {
            "block",
            "chebyshev",
            "lanczos",
        }:
            raise ValidationError("rec_method must be block, chebyshev, or lanczos")
        if not 0 <= self.lmx <= 3:
            raise ValidationError("lmx must be in [0, 3]")
        for field_name in ("min_empty_radius", "max_empty_radius", "lattice_dist", "basis_dist", "alat_scale"):
            value = float(getattr(self, field_name))
            if not math.isfinite(value):
                raise ValidationError(f"{field_name} must be finite")
        if self.min_empty_radius <= 0 or self.max_empty_radius <= 0:
            raise ValidationError("empty-sphere radii must be positive")
        if self.min_empty_radius > self.max_empty_radius:
            raise ValidationError("min_empty_radius cannot exceed max_empty_radius")
        if self.lattice_dist < 0 or self.basis_dist < 0 or self.alat_scale <= 0:
            raise ValidationError("distortions must be non-negative and alat_scale positive")
        if self.rec_steps <= 0 or self.sqs_shells < 1 or self.sqs_steps < 0:
            raise ValidationError("rec_steps and sqs_shells must be positive; sqs_steps non-negative")
        if len(self.supercell) != 3 or any(int(value) != value or int(value) < 1 for value in self.supercell):
            raise ValidationError("supercell must contain three positive integers")
        object.__setattr__(self, "supercell", tuple(int(value) for value in self.supercell))
        if self.z_iso_factor is not None and (not math.isfinite(self.z_iso_factor) or self.z_iso_factor <= 0):
            raise ValidationError("z_iso_factor must be positive when supplied")
        if self.z_iso_factor is not None and not self.mode_2d:
            raise ValidationError("z_iso_factor requires mode_2d")
        if self.seed is not None and not isinstance(self.seed, int):
            raise ValidationError("seed must be an integer or None")
        if self.sqs_seed is not None and not isinstance(self.sqs_seed, int):
            raise ValidationError("sqs_seed must be an integer or None")
        object.__setattr__(self, "rec_method", self.rec_method.lower())

    @property
    def sqs_spec(self) -> str | None:
        """Compatibility spelling for the SQS composition specification."""

        return self.sqs

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SetupSpec":
        return cls(**dict(data))

    def to_dict(self) -> dict[str, Any]:
        return {
            "cluster_type": self.cluster_type,
            "empty_spheres": self.empty_spheres,
            "min_empty_radius": self.min_empty_radius,
            "max_empty_radius": self.max_empty_radius,
            "lmx": self.lmx,
            "rec_method": self.rec_method,
            "rec_steps": self.rec_steps,
            "lattice_dist": self.lattice_dist,
            "basis_dist": self.basis_dist,
            "alat_scale": self.alat_scale,
            "seed": self.seed,
            "mode_2d": self.mode_2d,
            "z_iso_factor": self.z_iso_factor,
            "supercell": list(self.supercell),
            "sqs": self.sqs,
            "sqs_shells": self.sqs_shells,
            "sqs_steps": self.sqs_steps,
            "sqs_seed": self.sqs_seed,
        }


@dataclass(frozen=True)
class Diagnostic:
    severity: DiagnosticSeverity | str
    code: ErrorCode | str
    message: str
    path: str | None = None
    details: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        try:
            object.__setattr__(self, "severity", DiagnosticSeverity(self.severity))
        except ValueError as exc:
            raise ValidationError("diagnostic severity must be info, warning, or error") from exc
        object.__setattr__(self, "code", self.code.value if isinstance(self.code, ErrorCode) else str(self.code))
        if not self.message.strip():
            raise ValidationError("diagnostic message must not be empty")
        object.__setattr__(self, "details", _mapping(self.details))

    def to_dict(self) -> dict[str, Any]:
        return {
            "severity": self.severity.value,
            "code": self.code,
            "message": self.message,
            "path": self.path,
            "details": dict(self.details),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Diagnostic":
        return cls(
            severity=data["severity"],
            code=data["code"],
            message=str(data["message"]),
            path=data.get("path"),
            details=data.get("details"),
        )


@dataclass(frozen=True)
class ValidationResult:
    valid: bool
    diagnostics: tuple[Diagnostic, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {"valid": self.valid, "diagnostics": [item.to_dict() for item in self.diagnostics]}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ValidationResult":
        return cls(
            valid=bool(data["valid"]),
            diagnostics=tuple(Diagnostic.from_dict(item) for item in data.get("diagnostics", ())),
        )


@dataclass(frozen=True)
class OperationResult:
    status: OperationStatus = OperationStatus.SUCCESS
    diagnostics: tuple[Diagnostic, ...] = ()
    stdout: str = ""
    stderr: str = ""
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        try:
            object.__setattr__(self, "status", OperationStatus(self.status))
        except ValueError as exc:
            raise ValidationError("status must be a known operation status") from exc
        object.__setattr__(self, "provenance", _mapping(self.provenance))

    @property
    def success(self) -> bool:
        return self.status is OperationStatus.SUCCESS

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status.value,
            "success": self.success,
            "diagnostics": [item.to_dict() for item in self.diagnostics],
            "stdout": self.stdout,
            "stderr": self.stderr,
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class SetupResult(OperationResult):
    candidate_files: tuple[str, ...] = ()
    semantic_validation: ValidationResult | None = None
    sphere_outputs: Mapping[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        result = super().to_dict()
        result.update(
            {
                "candidate_files": list(self.candidate_files),
                "semantic_validation": (
                    None if self.semantic_validation is None else self.semantic_validation.to_dict()
                ),
                "sphere_outputs": dict(self.sphere_outputs),
            }
        )
        return result


@dataclass(frozen=True)
class PairKey:
    """Canonical periodic pair identity ``(site_i, site_j, lattice_R)``."""

    i: int
    j: int
    R: tuple[int, int, int]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.i, int)
            or isinstance(self.i, bool)
            or not isinstance(self.j, int)
            or isinstance(self.j, bool)
            or self.i < 0
            or self.j < 0
        ):
            raise ValidationError("pair site indices must be non-negative integers")
        if len(self.R) != 3 or any(
            not isinstance(value, int) or isinstance(value, bool) for value in self.R
        ):
            raise ValidationError("pair translation R must contain three integers")
        object.__setattr__(self, "R", tuple(self.R))

    def reverse(self) -> "PairKey":
        return PairKey(self.j, self.i, tuple(-value for value in self.R))

    def to_dict(self) -> dict[str, Any]:
        return {"i": self.i, "j": self.j, "R": list(self.R)}
