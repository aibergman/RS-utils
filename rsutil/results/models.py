"""Versioned, JSON-compatible models for parsed RS-LMTO runs."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping

from .units import Quantity, Unit


class ConvergenceState(str, Enum):
    CONVERGED = "converged"
    NOT_CONVERGED = "not_converged"
    UNKNOWN = "unknown"


class ParseSeverity(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True)
class ParseDiagnostic:
    """A parser observation; it is deliberately separate from validation."""

    severity: ParseSeverity | str
    code: str
    message: str
    source: str | None = None
    line: int | None = None
    details: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "severity", ParseSeverity(self.severity))
        if not self.code.strip() or not self.message.strip():
            raise ValueError("parse diagnostic code and message must not be empty")
        object.__setattr__(self, "details", dict(self.details))

    def to_dict(self) -> dict[str, Any]:
        return {
            "severity": self.severity.value,
            "code": self.code,
            "message": self.message,
            "source": self.source,
            "line": self.line,
            "details": _json_value(self.details),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ParseDiagnostic":
        return cls(
            severity=data["severity"],
            code=str(data["code"]),
            message=str(data["message"]),
            source=data.get("source"),
            line=data.get("line"),
            details=data.get("details", {}),
        )


@dataclass(frozen=True)
class ParseValidation:
    """Validation evidence for a parsed model, independent of parse notices."""

    valid: bool
    diagnostics: tuple[ParseDiagnostic, ...] = ()
    schema_version: str = "rsutil-parse-validation-v1"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "valid": self.valid,
            "diagnostics": [item.to_dict() for item in self.diagnostics],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ParseValidation":
        return cls(
            valid=bool(data["valid"]),
            diagnostics=tuple(ParseDiagnostic.from_dict(item) for item in data.get("diagnostics", ())),
            schema_version=str(data.get("schema_version", "rsutil-parse-validation-v1")),
        )


@dataclass(frozen=True)
class SourceFile:
    path: str
    content_hash: str
    role: str

    def to_dict(self) -> dict[str, str]:
        return {"path": self.path, "content_hash": self.content_hash, "role": self.role}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SourceFile":
        return cls(str(data["path"]), str(data["content_hash"]), str(data["role"]))


@dataclass(frozen=True)
class ParserProvenance:
    """Raw-source identity retained alongside derived result metadata."""

    source_path: str
    content_hash: str
    format: str = "rs-lmto-report"
    format_version: str | None = None
    parser_version: str = "rsutil-result-v1"
    source_files: tuple[SourceFile, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_path": self.source_path,
            "content_hash": self.content_hash,
            "format": self.format,
            "format_version": self.format_version,
            "parser_version": self.parser_version,
            "source_files": [item.to_dict() for item in self.source_files],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ParserProvenance":
        return cls(
            source_path=str(data["source_path"]),
            content_hash=str(data["content_hash"]),
            format=str(data.get("format", "rs-lmto-report")),
            format_version=data.get("format_version"),
            parser_version=str(data.get("parser_version", "rsutil-result-v1")),
            source_files=tuple(SourceFile.from_dict(item) for item in data.get("source_files", ())),
        )


@dataclass(frozen=True)
class SiteMoment:
    """One atom's parsed moment, using zero-based canonical site indexing."""

    site_index: int
    magnitude: Quantity | None = None
    vector: Quantity | None = None
    orientation: tuple[float, float, float] | None = None
    orientation_representation: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "site_index": self.site_index,
            "magnitude": None if self.magnitude is None else self.magnitude.to_dict(),
            "vector": None if self.vector is None else self.vector.to_dict(),
            "orientation": None if self.orientation is None else list(self.orientation),
            "orientation_representation": self.orientation_representation,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SiteMoment":
        orientation = data.get("orientation")
        return cls(
            site_index=int(data["site_index"]),
            magnitude=None if data.get("magnitude") is None else Quantity.from_dict(data["magnitude"]),
            vector=None if data.get("vector") is None else Quantity.from_dict(data["vector"]),
            orientation=None if orientation is None else tuple(float(value) for value in orientation),
            orientation_representation=data.get("orientation_representation"),
        )


@dataclass(frozen=True)
class SiteQuantity:
    site_index: int
    quantity: Quantity

    def to_dict(self) -> dict[str, Any]:
        return {"site_index": self.site_index, "quantity": self.quantity.to_dict()}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SiteQuantity":
        return cls(int(data["site_index"]), Quantity.from_dict(data["quantity"]))


@dataclass(frozen=True)
class ElectronicResult:
    """Electronic observables confidently parsed from an RS-LMTO run."""

    schema_version: str = "rsutil-electronic-result-v1"
    fermi_energy: Quantity | None = None
    total_energy: Quantity | None = None
    band_energy: Quantity | None = None
    site_occupations: tuple[SiteQuantity, ...] = ()
    site_up_orbital_occupations: tuple[SiteQuantity, ...] = ()
    site_down_orbital_occupations: tuple[SiteQuantity, ...] = ()
    site_charges: tuple[SiteQuantity, ...] = ()
    charge_transfers: tuple[SiteQuantity, ...] = ()

    @property
    def fermi(self) -> float | None:
        return None if self.fermi_energy is None else float(self.fermi_energy.value)

    def _electronic_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "fermi_energy": _quantity_dict(self.fermi_energy),
            "total_energy": _quantity_dict(self.total_energy),
            "band_energy": _quantity_dict(self.band_energy),
            "site_occupations": [item.to_dict() for item in self.site_occupations],
            "site_up_orbital_occupations": [item.to_dict() for item in self.site_up_orbital_occupations],
            "site_down_orbital_occupations": [item.to_dict() for item in self.site_down_orbital_occupations],
            "site_charges": [item.to_dict() for item in self.site_charges],
            "charge_transfers": [item.to_dict() for item in self.charge_transfers],
        }

    def to_dict(self) -> dict[str, Any]:
        return self._electronic_dict()

    def to_json(self, *, indent: int | None = None) -> str:
        return json.dumps(self.to_dict(), sort_keys=True, indent=indent, separators=None if indent else (",", ":"))


@dataclass(frozen=True)
class SCFResult(ElectronicResult):
    """Versioned SCF/converged-run result model."""

    schema_version: str = "rsutil-scf-result-v1"
    convergence: ConvergenceState | str = ConvergenceState.UNKNOWN
    total_spin_moment: Quantity | None = None
    site_spin_moments: tuple[SiteMoment, ...] = ()
    total_orbital_moment: Quantity | None = None
    site_orbital_moments: tuple[SiteMoment, ...] = ()
    total_rms_difference: Quantity | None = None
    diagnostics: tuple[ParseDiagnostic, ...] = ()
    validation: ParseValidation = field(default_factory=lambda: ParseValidation(True))
    provenance: ParserProvenance = field(
        default_factory=lambda: ParserProvenance(source_path="", content_hash="")
    )
    input_path: str | None = None
    run_manifest: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "convergence", ConvergenceState(self.convergence))
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))
        if self.run_manifest is not None:
            object.__setattr__(self, "run_manifest", _json_value(dict(self.run_manifest)))

    @property
    def converged(self) -> bool | None:
        if self.convergence is ConvergenceState.UNKNOWN:
            return None
        return self.convergence is ConvergenceState.CONVERGED

    @property
    def site_magnetic_moments(self) -> tuple[SiteMoment, ...]:
        """Descriptive alias for the spin-moment records."""

        return self.site_spin_moments

    @property
    def convergence_status(self) -> ConvergenceState:
        """Compatibility spelling for consumers using a status noun."""

        return self.convergence

    @property
    def total_magnetic_moment(self) -> Quantity | None:
        return self.total_spin_moment

    @property
    def site_moments(self) -> tuple[SiteMoment, ...]:
        return self.site_spin_moments

    def to_dict(self) -> dict[str, Any]:
        result = self._electronic_dict()
        result.update(
            {
                "convergence": self.convergence.value,
                "total_spin_moment": _quantity_dict(self.total_spin_moment),
                "site_spin_moments": [item.to_dict() for item in self.site_spin_moments],
                "total_orbital_moment": _quantity_dict(self.total_orbital_moment),
                "site_orbital_moments": [item.to_dict() for item in self.site_orbital_moments],
                "total_rms_difference": _quantity_dict(self.total_rms_difference),
                "diagnostics": [item.to_dict() for item in self.diagnostics],
                "validation": self.validation.to_dict(),
                "provenance": self.provenance.to_dict(),
                "input_path": self.input_path,
                "run_manifest": _json_value(self.run_manifest),
            }
        )
        return result

    def to_json(self, *, indent: int | None = None) -> str:
        return json.dumps(self.to_dict(), sort_keys=True, indent=indent, separators=None if indent else (",", ":"))

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SCFResult":
        return cls(
            schema_version=str(data.get("schema_version", "rsutil-scf-result-v1")),
            fermi_energy=_quantity_from(data.get("fermi_energy")),
            total_energy=_quantity_from(data.get("total_energy")),
            band_energy=_quantity_from(data.get("band_energy")),
            site_occupations=tuple(SiteQuantity.from_dict(item) for item in data.get("site_occupations", ())),
            site_up_orbital_occupations=tuple(SiteQuantity.from_dict(item) for item in data.get("site_up_orbital_occupations", ())),
            site_down_orbital_occupations=tuple(SiteQuantity.from_dict(item) for item in data.get("site_down_orbital_occupations", ())),
            site_charges=tuple(SiteQuantity.from_dict(item) for item in data.get("site_charges", ())),
            charge_transfers=tuple(SiteQuantity.from_dict(item) for item in data.get("charge_transfers", ())),
            convergence=data.get("convergence", ConvergenceState.UNKNOWN),
            total_spin_moment=_quantity_from(data.get("total_spin_moment")),
            site_spin_moments=tuple(SiteMoment.from_dict(item) for item in data.get("site_spin_moments", ())),
            total_orbital_moment=_quantity_from(data.get("total_orbital_moment")),
            site_orbital_moments=tuple(SiteMoment.from_dict(item) for item in data.get("site_orbital_moments", ())),
            total_rms_difference=_quantity_from(data.get("total_rms_difference")),
            diagnostics=tuple(ParseDiagnostic.from_dict(item) for item in data.get("diagnostics", ())),
            validation=ParseValidation.from_dict(data.get("validation", {"valid": True})),
            provenance=ParserProvenance.from_dict(data.get("provenance", {"source_path": "", "content_hash": ""})),
            input_path=data.get("input_path"),
            run_manifest=data.get("run_manifest"),
        )

    @classmethod
    def from_json(cls, value: str) -> "SCFResult":
        return cls.from_dict(json.loads(value))


def _quantity_dict(value: Quantity | None) -> dict[str, Any] | None:
    return None if value is None else value.to_dict()


def _quantity_from(value: Any) -> Quantity | None:
    return None if value is None else Quantity.from_dict(value)


def _json_value(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Quantity):
        return value.to_dict()
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    return value


__all__ = [
    "ConvergenceState",
    "ElectronicResult",
    "ParseDiagnostic",
    "ParseSeverity",
    "ParseValidation",
    "ParserProvenance",
    "SCFResult",
    "SiteMoment",
    "SiteQuantity",
    "SourceFile",
]
