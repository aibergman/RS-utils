"""Evidence-first contracts for RS-LMTO transport and damping results.

The public fable_v4 source contains legacy conductivity and damping writers,
but their text outputs do not carry enough machine-readable units and context
for safe HT parsing.  This module therefore provides the versioned result
shape and an intentionally small registry, but does not guess at legacy text
columns.  Files that cannot be identified as a documented stable format
produce an explicit :class:`UnsupportedTransportResult`.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from .errors import MissingResultFileError
from .models import (
    ParseDiagnostic,
    ParseSeverity,
    ParseValidation,
    ParserProvenance,
    SourceFile,
)
from .units import DIMENSIONLESS, Quantity, Unit


class TransportSupport(str, Enum):
    """Evidence status of a transport/damping output format."""

    SUPPORTED_STABLE = "SUPPORTED_STABLE"
    SUPPORTED_EXPERIMENTAL = "SUPPORTED_EXPERIMENTAL"
    UNSUPPORTED_AMBIGUOUS = "UNSUPPORTED_AMBIGUOUS"


class TransportObservable(str, Enum):
    CONDUCTIVITY = "conductivity"
    DAMPING = "damping"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ResultAxis:
    """One explicit result axis, including its role and unit."""

    name: str
    values: tuple[Any, ...]
    unit: Unit | str = DIMENSIONLESS
    role: str = "index"

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("result axis name must not be empty")
        if not self.values:
            raise ValueError(f"result axis {self.name!r} must not be empty")
        object.__setattr__(self, "values", tuple(self.values))
        object.__setattr__(self, "unit", Unit(self.unit))
        object.__setattr__(self, "role", str(self.role))

    @property
    def length(self) -> int:
        return len(self.values)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "values": list(self.values),
            "unit": self.unit.value,
            "role": self.role,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ResultAxis":
        return cls(
            name=str(data["name"]),
            values=tuple(data["values"]),
            unit=data.get("unit", DIMENSIONLESS),
            role=str(data.get("role", "index")),
        )


def _shape(value: Any) -> tuple[int, ...]:
    if not isinstance(value, (tuple, list)):
        return ()
    if not value:
        return (0,)
    child_shapes = {_shape(item) for item in value}
    if len(child_shapes) != 1:
        raise ValueError("result values must form a rectangular array")
    return (len(value),) + next(iter(child_shapes))


def _require_text(name: str, value: str) -> str:
    value = str(value).strip()
    if not value:
        raise ValueError(f"{name} must not be empty")
    return value


@dataclass(frozen=True)
class ConductivityResult:
    """Typed conductivity contract for a future evidenced parser.

    Values are ordered as ``axes + component_order``.  The component order is
    explicit because no transpose or Cartesian convention may be inferred
    from a file.  ``method``, ``channel`` and ``approximation`` are required
    interpretation metadata, not presentation labels.
    """

    values: Quantity
    axes: tuple[ResultAxis, ...]
    component_order: tuple[str, ...]
    method: str
    channel: str
    approximation: str
    provenance: ParserProvenance
    options: Mapping[str, Any] = field(default_factory=dict)
    diagnostics: tuple[ParseDiagnostic, ...] = ()
    validation: ParseValidation = field(default_factory=lambda: ParseValidation(True))
    schema_version: str = "rsutil-conductivity-result-v1"

    def __post_init__(self) -> None:
        if not isinstance(self.values, Quantity):
            raise ValueError("conductivity values must be a Quantity")
        object.__setattr__(self, "axes", tuple(self.axes))
        object.__setattr__(self, "component_order", tuple(str(item) for item in self.component_order))
        if not self.component_order or any(not item.strip() for item in self.component_order):
            raise ValueError("conductivity component_order must be non-empty")
        expected = tuple(axis.length for axis in self.axes) + (len(self.component_order),)
        if _shape(self.values.value) != expected:
            raise ValueError(f"conductivity values shape {_shape(self.values.value)} != expected {expected}")
        for name, value in (("method", self.method), ("channel", self.channel), ("approximation", self.approximation)):
            _require_text(name, value)
        object.__setattr__(self, "options", dict(self.options))
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "observable": TransportObservable.CONDUCTIVITY.value,
            "values": self.values.to_dict(),
            "axes": [axis.to_dict() for axis in self.axes],
            "component_order": list(self.component_order),
            "method": self.method,
            "channel": self.channel,
            "approximation": self.approximation,
            "options": self.options,
            "diagnostics": [item.to_dict() for item in self.diagnostics],
            "validation": self.validation.to_dict(),
            "provenance": self.provenance.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ConductivityResult":
        return cls(
            values=Quantity.from_dict(data["values"]),
            axes=tuple(ResultAxis.from_dict(item) for item in data.get("axes", ())),
            component_order=tuple(data["component_order"]),
            method=str(data["method"]),
            channel=str(data["channel"]),
            approximation=str(data["approximation"]),
            options=data.get("options", {}),
            diagnostics=tuple(ParseDiagnostic.from_dict(item) for item in data.get("diagnostics", ())),
            validation=ParseValidation.from_dict(data.get("validation", {"valid": True})),
            provenance=ParserProvenance.from_dict(data["provenance"]),
            schema_version=str(data.get("schema_version", "rsutil-conductivity-result-v1")),
        )

    def to_json(self, *, indent: int | None = None) -> str:
        return json.dumps(self.to_dict(), sort_keys=True, indent=indent)


@dataclass(frozen=True)
class DampingResult:
    """Typed site/pair damping contract for a future evidenced parser.

    ``index_kind`` must be ``site`` or ``pair``.  ``index_labels`` identify
    the corresponding last-but-one axis; ``component_order`` identifies the
    scalar or tensor components on the last axis.  This prevents a parser
    from calling an unlabeled scalar "Gilbert damping".
    """

    values: Quantity
    axes: tuple[ResultAxis, ...]
    index_kind: str
    index_labels: tuple[str, ...]
    component_order: tuple[str, ...]
    method: str
    channel: str
    approximation: str
    provenance: ParserProvenance
    options: Mapping[str, Any] = field(default_factory=dict)
    diagnostics: tuple[ParseDiagnostic, ...] = ()
    validation: ParseValidation = field(default_factory=lambda: ParseValidation(True))
    schema_version: str = "rsutil-damping-result-v1"

    def __post_init__(self) -> None:
        if not isinstance(self.values, Quantity):
            raise ValueError("damping values must be a Quantity")
        object.__setattr__(self, "axes", tuple(self.axes))
        object.__setattr__(self, "index_kind", str(self.index_kind))
        object.__setattr__(self, "index_labels", tuple(str(item) for item in self.index_labels))
        object.__setattr__(self, "component_order", tuple(str(item) for item in self.component_order))
        if self.index_kind not in {"site", "pair"}:
            raise ValueError("damping index_kind must be 'site' or 'pair'")
        if not self.index_labels or not self.component_order:
            raise ValueError("damping index labels and component_order must be non-empty")
        expected = tuple(axis.length for axis in self.axes) + (len(self.index_labels), len(self.component_order))
        if _shape(self.values.value) != expected:
            raise ValueError(f"damping values shape {_shape(self.values.value)} != expected {expected}")
        for name, value in (("method", self.method), ("channel", self.channel), ("approximation", self.approximation)):
            _require_text(name, value)
        object.__setattr__(self, "options", dict(self.options))
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "observable": TransportObservable.DAMPING.value,
            "values": self.values.to_dict(),
            "axes": [axis.to_dict() for axis in self.axes],
            "index_kind": self.index_kind,
            "index_labels": list(self.index_labels),
            "component_order": list(self.component_order),
            "method": self.method,
            "channel": self.channel,
            "approximation": self.approximation,
            "options": self.options,
            "diagnostics": [item.to_dict() for item in self.diagnostics],
            "validation": self.validation.to_dict(),
            "provenance": self.provenance.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "DampingResult":
        return cls(
            values=Quantity.from_dict(data["values"]),
            axes=tuple(ResultAxis.from_dict(item) for item in data.get("axes", ())),
            index_kind=str(data["index_kind"]),
            index_labels=tuple(data["index_labels"]),
            component_order=tuple(data["component_order"]),
            method=str(data["method"]),
            channel=str(data["channel"]),
            approximation=str(data["approximation"]),
            options=data.get("options", {}),
            diagnostics=tuple(ParseDiagnostic.from_dict(item) for item in data.get("diagnostics", ())),
            validation=ParseValidation.from_dict(data.get("validation", {"valid": True})),
            provenance=ParserProvenance.from_dict(data["provenance"]),
            schema_version=str(data.get("schema_version", "rsutil-damping-result-v1")),
        )

    def to_json(self, *, indent: int | None = None) -> str:
        return json.dumps(self.to_dict(), sort_keys=True, indent=indent)


@dataclass(frozen=True)
class UnsupportedTransportResult:
    """A deliberate non-result for an absent, ambiguous, or unstable format."""

    observable: TransportObservable | str
    status: TransportSupport
    reason: str
    provenance: ParserProvenance
    diagnostics: tuple[ParseDiagnostic, ...] = ()
    input_path: str | None = None
    run_manifest: Mapping[str, Any] | None = None
    schema_version: str = "rsutil-transport-damping-result-v1"

    def __post_init__(self) -> None:
        object.__setattr__(self, "observable", TransportObservable(self.observable))
        object.__setattr__(self, "status", TransportSupport(self.status))
        _require_text("reason", self.reason)
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))
        if self.run_manifest is not None:
            object.__setattr__(self, "run_manifest", dict(self.run_manifest))

    @property
    def supported(self) -> bool:
        return self.status in {TransportSupport.SUPPORTED_STABLE, TransportSupport.SUPPORTED_EXPERIMENTAL}

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "observable": self.observable.value,
            "status": self.status.value,
            "reason": self.reason,
            "diagnostics": [item.to_dict() for item in self.diagnostics],
            "provenance": self.provenance.to_dict(),
            "input_path": self.input_path,
            "run_manifest": self.run_manifest,
        }

    def to_json(self, *, indent: int | None = None) -> str:
        return json.dumps(self.to_dict(), sort_keys=True, indent=indent)


@dataclass(frozen=True)
class ParserSpec:
    """One explicit registry entry; this is intentionally not a plugin API."""

    name: str
    observable: TransportObservable
    status: TransportSupport
    evidence: str
    parser: Callable[..., Any] | None = None


PARSER_REGISTRY: tuple[ParserSpec, ...] = (
    ParserSpec(
        name="rs-lmto-conductivity-text",
        observable=TransportObservable.CONDUCTIVITY,
        status=TransportSupport.UNSUPPORTED_AMBIGUOUS,
        evidence=(
            "fable_v4 source/conductivity.f90::calculate_conductivity_tensor writes "
            "energy/real/imag text, but the files do not declare units, operator "
            "channel, tensor convention, or calculation options."
        ),
    ),
    ParserSpec(
        name="rs-lmto-damping-text",
        observable=TransportObservable.DAMPING,
        status=TransportSupport.UNSUPPORTED_AMBIGUOUS,
        evidence=(
            "fable_v4 source/exchange_dynamics.f90::calculate_gilbert_damping writes "
            "pair and energy damping tables, but the files do not declare units, "
            "normalization, SOC/input identity, or a complete run context."
        ),
    ),
)


_CANDIDATE_NAMES = {
    "conductivity.out",
    "transport.out",
    "kubo.out",
    "kpm.out",
    "damping.out",
    "gilbert.out",
    "alpha.out",
    "report.out",
}


def _is_candidate_name(name: str) -> bool:
    """Return whether *name* is a known legacy or fable_v4 result filename."""

    lower = name.lower()
    if lower in _CANDIDATE_NAMES:
        return True
    return lower in {
        "cond_total.out",
        "cond_total_orb_real.out",
        "cond_total_orb_im.out",
        "damping-energy.out",
        "alldampings.out",
    } or lower.endswith(("_cond.out", "_cond_orb_real.out", "_cond_orb_im.out"))


def parse_transport_result(path: str | Path, *, observable: TransportObservable | str | None = None) -> UnsupportedTransportResult:
    """Inspect a transport/damping result without guessing an unstable schema.

    No parser is currently production-supported.  Existing candidate files
    return a structured unsupported result; a missing explicit file remains a
    normal missing-file error.
    """

    supplied = Path(path).expanduser()
    requested = TransportObservable.UNKNOWN if observable is None else TransportObservable(observable)
    if supplied.is_dir():
        root = supplied.resolve()
        candidates = tuple(sorted(item for item in root.iterdir() if item.is_file() and _is_candidate_name(item.name)))
        if not candidates:
            raise MissingResultFileError(f"no transport/damping result candidate found in {root}")
        source = candidates[0]
        ambiguous = len(candidates) > 1
    elif supplied.is_file():
        source = supplied.resolve()
        root = source.parent
        candidates = (source,)
        ambiguous = False
    else:
        raise MissingResultFileError(f"result path does not exist: {supplied}")

    detected = _detect_observable(source.name)
    if requested is not TransportObservable.UNKNOWN and detected is not TransportObservable.UNKNOWN and requested is not detected:
        status = TransportSupport.UNSUPPORTED_AMBIGUOUS
        reason = f"file identifies {detected.value}, not requested {requested.value}; no values were parsed"
    else:
        status = TransportSupport.UNSUPPORTED_AMBIGUOUS
        reason = "no evidenced stable RS-LMTO transport/damping text schema is available; no values were parsed"
    if ambiguous:
        reason += "; multiple candidate outputs make the run selection ambiguous"

    source_files = [SourceFile(str(item), _sha256(item), "transport_candidate") for item in candidates]
    for name, role in (("input.nml", "input"), ("lattice.nml", "lattice")):
        candidate = root / name
        if candidate.is_file():
            source_files.append(SourceFile(str(candidate), _sha256(candidate), role))
    manifest = _read_manifest(root)
    if manifest is not None:
        manifest_path, manifest_data = manifest
        source_files.append(SourceFile(str(manifest_path), _sha256(manifest_path), "run_manifest"))
    else:
        manifest_data = None

    diagnostics = (
        ParseDiagnostic(
            ParseSeverity.WARNING,
            "transport_format_unsupported",
            reason,
            source=str(source),
            details={"candidate_files": [str(item) for item in candidates]},
        ),
    )
    provenance = ParserProvenance(
        source_path=str(source),
        content_hash=_sha256(source),
        format="rs-lmto-transport-damping",
        format_version=None,
        parser_version="rsutil-transport-result-v1",
        source_files=tuple(source_files),
    )
    return UnsupportedTransportResult(
        observable=requested if requested is not TransportObservable.UNKNOWN else detected,
        status=status,
        reason=reason,
        diagnostics=diagnostics,
        provenance=provenance,
        input_path=str(root / "input.nml") if (root / "input.nml").is_file() else None,
        run_manifest=manifest_data,
    )


def parse_conductivity_result(path: str | Path) -> UnsupportedTransportResult:
    return parse_transport_result(path, observable=TransportObservable.CONDUCTIVITY)


def parse_damping_result(path: str | Path) -> UnsupportedTransportResult:
    return parse_transport_result(path, observable=TransportObservable.DAMPING)


def _detect_observable(name: str) -> TransportObservable:
    lower = name.lower()
    if any(token in lower for token in ("damp", "gilbert", "alpha")):
        return TransportObservable.DAMPING
    if any(token in lower for token in ("conduct", "transport", "kubo", "kpm", "sigma")):
        return TransportObservable.CONDUCTIVITY
    return TransportObservable.UNKNOWN


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_manifest(root: Path) -> tuple[Path, Mapping[str, Any]] | None:
    for name in ("run_manifest.json", "rsutil-run.json", "manifest.json", "rsutil-exchange-manifest.json"):
        candidate = root / name
        if not candidate.is_file():
            continue
        try:
            data = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            return None
        if isinstance(data, Mapping):
            return candidate, data
    return None


__all__ = [
    "ConductivityResult",
    "DampingResult",
    "PARSER_REGISTRY",
    "ParserSpec",
    "ResultAxis",
    "TransportObservable",
    "TransportSupport",
    "UnsupportedTransportResult",
    "parse_conductivity_result",
    "parse_damping_result",
    "parse_transport_result",
]
