"""Conservative parser for the legacy RS-LMTO text result files."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .errors import MalformedResultError, MissingResultFileError
from .models import (
    ConvergenceState,
    ParseDiagnostic,
    ParseSeverity,
    ParseValidation,
    ParserProvenance,
    SCFResult,
    SiteMoment,
    SiteQuantity,
    SourceFile,
)
from .units import ELECTRONS, Quantity, RYDBERG, UNKNOWN_UNIT, scalar, vector


_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_NUMBER_RE = re.compile(rf"^{_NUMBER}$")
_FERMI_RE = re.compile(rf"^\s*Fermi\s+energy\s*:\s*(?P<value>{_NUMBER})\s*$", re.I)
_ENERGY_RE = re.compile(rf"^\s*(?P<label>Total|Band) energy of system\s*:\s*(?P<value>{_NUMBER})\s*$", re.I)
_TOTAL_VECTOR_RE = re.compile(
    rf"^\s*Total\s+(?P<label>spin|orbital)\s+moment\s*:\s*(?P<values>.*)$", re.I
)
_SITE_SCALAR_RE = re.compile(
    rf"^\s*(?P<label>Spin moment|Orbital moment)\s+of\s+atom\s+(?P<index>\d+)\s*:\s*(?P<value>{_NUMBER})\s*$",
    re.I,
)
_SITE_VECTOR_RE = re.compile(
    r"^\s*(?P<label>Spin moment projections|Orbital moment projections)\s+of\s+atom\s+(?P<index>\d+)\s*:\s*(?P<values>.*)$",
    re.I,
)
_SITE_DIRECTION_RE = re.compile(
    r"^\s*(?P<label>Spin moment direction|Orbital moment direction)\s+of\s+atom\s+(?P<index>\d+)\s*:\s*(?P<values>.*)$",
    re.I,
)
_SITE_ARRAY_RE = re.compile(
    r"^\s*(?P<label>Up orbital occupation|Down orbital occupation)\s+at\s+atom\s+(?P<index>\d+)\s*:\s*(?P<values>.*)$",
    re.I,
)
_SITE_VALUE_RE = re.compile(
    rf"^\s*(?P<label>Occupation|Charge|Charge transfer)\s+at\s+atom\s+(?P<index>\d+)\s*:\s*(?P<value>{_NUMBER})\s*$",
    re.I,
)
_RMS_RE = re.compile(rf"^\s*Total RMS Diff\s*:\s*(?P<value>{_NUMBER})\s*$", re.I)
_CONVERGENCE_RE = re.compile(r"(?P<state>Not\s+converged!|Converged!)", re.I)

_REPORT_NAMES = ("report.out",)
_LOG_NAMES = ("out.scf", "out.cold")
_MANIFEST_NAMES = (
    "run_manifest.json",
    "rsutil-run.json",
    "manifest.json",
    "rsutil-exchange-manifest.json",
)


def parse_run(path: str | Path) -> SCFResult:
    """Parse one explicit RS-LMTO run directory or ``report.out`` path.

    The report is the required source.  A neighboring ``out.scf`` or
    ``out.cold`` is used only as explicit convergence evidence; a report's
    presence or its ``Calculation finished`` log line is not evidence of SCF
    convergence.
    """

    supplied = Path(path).expanduser()
    if supplied.is_dir():
        root = supplied.resolve()
        report_path = root / _REPORT_NAMES[0]
    elif supplied.is_file():
        report_path = supplied.resolve()
        root = report_path.parent
    else:
        raise MissingResultFileError(f"result path does not exist: {supplied}")

    if not report_path.is_file():
        raise MissingResultFileError(f"required report file is missing: {report_path}")

    report_text = _read_text(report_path)
    parsed = _parse_report(report_text, report_path)

    auxiliary: list[tuple[Path, str]] = []
    for name in _LOG_NAMES:
        candidate = root / name
        if candidate.is_file():
            auxiliary.append((candidate, "run_log"))
    orientation_path = root / "minfo.out"
    if orientation_path.is_file():
        auxiliary.append((orientation_path, "spin_orientation"))
        _merge_orientations(parsed["spin_moments"], _parse_directions(_read_text(orientation_path), orientation_path), "spin")
    orbital_orientation_path = root / "linfo.out"
    if orbital_orientation_path.is_file():
        auxiliary.append((orbital_orientation_path, "orbital_orientation"))
        _merge_orientations(parsed["orbital_moments"], _parse_directions(_read_text(orbital_orientation_path), orbital_orientation_path), "orbital")

    convergence, convergence_diagnostics = _convergence_from_logs(auxiliary)
    diagnostics: list[ParseDiagnostic] = parsed["diagnostics"] + convergence_diagnostics
    if convergence is ConvergenceState.UNKNOWN:
        diagnostics.append(
            ParseDiagnostic(
                ParseSeverity.WARNING,
                "convergence_unknown",
                "no explicit Converged!/Not converged! evidence was found",
                source=str(report_path),
            )
        )

    source_files = [SourceFile(str(report_path), _sha256(report_path), "report")]
    for auxiliary_path, role in auxiliary:
        source_files.append(SourceFile(str(auxiliary_path), _sha256(auxiliary_path), role))
    run_manifest, manifest_diagnostic = _read_manifest(root)
    if manifest_diagnostic is not None:
        diagnostics.append(manifest_diagnostic)
    for name in _MANIFEST_NAMES:
        manifest_path = root / name
        if manifest_path.is_file():
            source_files.append(SourceFile(str(manifest_path), _sha256(manifest_path), "run_manifest"))
            break

    provenance = ParserProvenance(
        source_path=str(report_path),
        content_hash=_sha256(report_path),
        format="rs-lmto-report",
        format_version=None,
        source_files=tuple(source_files),
    )

    validation = _validate(parsed, provenance)
    return SCFResult(
        convergence=convergence,
        fermi_energy=parsed["fermi_energy"],
        total_energy=parsed["total_energy"],
        band_energy=parsed["band_energy"],
        site_occupations=tuple(sorted(parsed["occupations"].values(), key=lambda item: item.site_index)),
        site_up_orbital_occupations=tuple(sorted(parsed["up_occupations"].values(), key=lambda item: item.site_index)),
        site_down_orbital_occupations=tuple(sorted(parsed["down_occupations"].values(), key=lambda item: item.site_index)),
        site_charges=tuple(sorted(parsed["charges"].values(), key=lambda item: item.site_index)),
        charge_transfers=tuple(sorted(parsed["charge_transfers"].values(), key=lambda item: item.site_index)),
        total_spin_moment=parsed["total_spin_moment"],
        site_spin_moments=tuple(sorted(parsed["spin_moments"].values(), key=lambda item: item.site_index)),
        total_orbital_moment=parsed["total_orbital_moment"],
        site_orbital_moments=tuple(sorted(parsed["orbital_moments"].values(), key=lambda item: item.site_index)),
        total_rms_difference=parsed["total_rms_difference"],
        diagnostics=tuple(diagnostics),
        validation=validation,
        provenance=provenance,
        input_path=str(root / "input.nml") if (root / "input.nml").is_file() else None,
        run_manifest=run_manifest,
    )


def parse_report(path: str | Path) -> SCFResult:
    """Descriptive alias for :func:`parse_run`."""

    return parse_run(path)


parse_scf = parse_run


def _parse_report(text: str, source: Path) -> dict:
    values: dict = {
        "fermi_energy": None,
        "total_energy": None,
        "band_energy": None,
        "total_spin_moment": None,
        "total_orbital_moment": None,
        "total_rms_difference": None,
        "spin_moments": {},
        "orbital_moments": {},
        "occupations": {},
        "up_occupations": {},
        "down_occupations": {},
        "charges": {},
        "charge_transfers": {},
        "diagnostics": [],
    }
    recognized = 0
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        match = _FERMI_RE.match(line)
        if match:
            if values["fermi_energy"] is not None:
                raise MalformedResultError(f"multiple Fermi energy records in {source}")
            values["fermi_energy"] = scalar(_number(match.group("value")), RYDBERG)
            recognized += 1
            continue
        if "fermi energy" in stripped.lower() and ":" in stripped:
            raise MalformedResultError(f"malformed Fermi energy at {source}:{line_number}")

        match = _ENERGY_RE.match(line)
        if match:
            quantity = scalar(_number(match.group("value")), UNKNOWN_UNIT)
            key = "total_energy" if match.group("label").lower() == "total" else "band_energy"
            if values[key] is not None:
                raise MalformedResultError(f"multiple {key.replace('_', ' ')} records in {source}")
            values[key] = quantity
            recognized += 1
            continue
        if re.match(r"^\s*(Total|Band) energy of system\s*:", line, re.I):
            raise MalformedResultError(f"malformed energy record at {source}:{line_number}")

        match = _TOTAL_VECTOR_RE.match(line)
        if match:
            parsed_vector = _numbers(match.group("values"), expected=3, source=source, line_number=line_number)
            key = "total_spin_moment" if match.group("label").lower() == "spin" else "total_orbital_moment"
            if values[key] is not None:
                raise MalformedResultError(f"multiple {key.replace('_', ' ')} records in {source}")
            values[key] = vector(parsed_vector, UNKNOWN_UNIT)
            recognized += 1
            continue

        match = _SITE_VECTOR_RE.match(line)
        if match:
            site = _site_index(match.group("index"), source, line_number)
            parsed_vector = _numbers(match.group("values"), expected=3, source=source, line_number=line_number)
            target = values["spin_moments"] if match.group("label").lower().startswith("spin") else values["orbital_moments"]
            record = target.get(site, SiteMoment(site_index=site))
            if record.vector is not None:
                raise MalformedResultError(f"multiple moment-vector records for atom {site + 1} in {source}")
            target[site] = SiteMoment(
                site_index=site,
                magnitude=record.magnitude,
                vector=vector(parsed_vector, UNKNOWN_UNIT),
                orientation=record.orientation,
                orientation_representation=record.orientation_representation,
            )
            recognized += 1
            continue

        match = _SITE_SCALAR_RE.match(line)
        if match:
            site = _site_index(match.group("index"), source, line_number)
            target = values["spin_moments"] if match.group("label").lower().startswith("spin") else values["orbital_moments"]
            record = target.get(site, SiteMoment(site_index=site))
            if record.magnitude is not None:
                raise MalformedResultError(f"multiple moment-magnitude records for atom {site + 1} in {source}")
            target[site] = SiteMoment(
                site_index=site,
                magnitude=scalar(_number(match.group("value")), UNKNOWN_UNIT),
                vector=record.vector,
                orientation=record.orientation,
                orientation_representation=record.orientation_representation,
            )
            recognized += 1
            continue
        if re.match(r"^\s*(Spin|Orbital) moment(?: projections)?\s+of\s+atom", line, re.I):
            raise MalformedResultError(f"malformed site moment record at {source}:{line_number}")

        match = _SITE_VALUE_RE.match(line)
        if match:
            site = _site_index(match.group("index"), source, line_number)
            label = match.group("label").lower()
            quantity = scalar(_number(match.group("value")), ELECTRONS)
            if label == "occupation":
                if site in values["occupations"]:
                    raise MalformedResultError(f"multiple occupation records for atom {site + 1} in {source}")
                values["occupations"][site] = SiteQuantity(site, quantity)
            elif label == "charge":
                if site in values["charges"]:
                    raise MalformedResultError(f"multiple charge records for atom {site + 1} in {source}")
                values["charges"][site] = SiteQuantity(site, quantity)
            else:
                if site in values["charge_transfers"]:
                    raise MalformedResultError(f"multiple charge-transfer records for atom {site + 1} in {source}")
                values["charge_transfers"][site] = SiteQuantity(site, quantity)
            recognized += 1
            continue
        if re.match(r"^\s*(Occupation|Charge|Charge transfer)\s+at\s+atom", line, re.I):
            raise MalformedResultError(f"malformed site charge record at {source}:{line_number}")

        match = _SITE_ARRAY_RE.match(line)
        if match:
            site = _site_index(match.group("index"), source, line_number)
            numbers = _numbers(match.group("values"), expected=None, source=source, line_number=line_number)
            quantity = SiteQuantity(site, Quantity(tuple(numbers), ELECTRONS))
            if match.group("label").lower().startswith("up"):
                if site in values["up_occupations"]:
                    raise MalformedResultError(f"multiple up-occupation records for atom {site + 1} in {source}")
                values["up_occupations"][site] = quantity
            else:
                if site in values["down_occupations"]:
                    raise MalformedResultError(f"multiple down-occupation records for atom {site + 1} in {source}")
                values["down_occupations"][site] = quantity
            recognized += 1
            continue

        match = _RMS_RE.match(line)
        if match:
            if values["total_rms_difference"] is not None:
                raise MalformedResultError(f"multiple RMS-difference records in {source}")
            values["total_rms_difference"] = scalar(_number(match.group("value")), UNKNOWN_UNIT)
            recognized += 1
            continue
        if re.match(r"^\s*Total RMS Diff\s*:", line, re.I):
            raise MalformedResultError(f"malformed RMS-difference record at {source}:{line_number}")

    if values["fermi_energy"] is None:
        raise MissingResultFileError(f"Fermi energy record is missing from {source}")
    if recognized == 0:
        raise MalformedResultError(f"no supported records found in {source}")
    return values


def _parse_directions(text: str, source: Path) -> dict[int, tuple[float, float, float]]:
    result: dict[int, tuple[float, float, float]] = {}
    for line_number, line in enumerate(text.splitlines(), start=1):
        match = _SITE_DIRECTION_RE.match(line)
        if not match:
            if "moment direction" in line.lower() and ":" in line:
                raise MalformedResultError(f"malformed direction record at {source}:{line_number}")
            continue
        site = _site_index(match.group("index"), source, line_number)
        if site in result:
            raise MalformedResultError(f"multiple moment directions for atom {site + 1} in {source}")
        values = _numbers(match.group("values"), expected=3, source=source, line_number=line_number)
        norm = math.sqrt(sum(value * value for value in values))
        if not math.isfinite(norm) or norm == 0:
            raise MalformedResultError(f"zero/non-finite moment direction at {source}:{line_number}")
        result[site] = tuple(values)  # type: ignore[assignment]
    return result


def _merge_orientations(records: dict[int, SiteMoment], directions: Mapping[int, Sequence[float]], kind: str) -> None:
    for site, direction in directions.items():
        record = records.get(site, SiteMoment(site_index=site))
        records[site] = SiteMoment(
            site_index=site,
            magnitude=record.magnitude,
            vector=record.vector,
            orientation=tuple(float(value) for value in direction),
            orientation_representation="cartesian_unit_vector",
        )


def _convergence_from_logs(auxiliary: Iterable[tuple[Path, str]]) -> tuple[ConvergenceState, list[ParseDiagnostic]]:
    evidence: list[tuple[ConvergenceState, Path]] = []
    diagnostics: list[ParseDiagnostic] = []
    for path, role in auxiliary:
        if role != "run_log":
            continue
        text = _read_text(path)
        matches = list(_CONVERGENCE_RE.finditer(text))
        if matches:
            match = matches[-1]
            state = (
                ConvergenceState.NOT_CONVERGED
                if match.group("state").lower().startswith("not")
                else ConvergenceState.CONVERGED
            )
            evidence.append((state, path))
    if not evidence:
        return ConvergenceState.UNKNOWN, diagnostics
    states = {item[0] for item in evidence}
    if len(states) > 1:
        diagnostics.append(
            ParseDiagnostic(
                ParseSeverity.WARNING,
                "conflicting_convergence_evidence",
                "convergence logs contain both converged and not-converged markers",
                source=str(evidence[-1][1]),
            )
        )
        return ConvergenceState.UNKNOWN, diagnostics
    return evidence[-1][0], diagnostics


def _validate(parsed: Mapping, provenance: ParserProvenance) -> ParseValidation:
    diagnostics: list[ParseDiagnostic] = []
    fermi = parsed["fermi_energy"]
    if fermi is None or fermi.unit is not RYDBERG:
        diagnostics.append(ParseDiagnostic(ParseSeverity.ERROR, "fermi_unit_invalid", "Fermi energy does not carry the established Ry unit", source=provenance.source_path))
    for name in ("spin_moments", "orbital_moments"):
        for site, record in parsed[name].items():
            if site < 0:
                diagnostics.append(ParseDiagnostic(ParseSeverity.ERROR, "site_index_invalid", "site index is negative", source=provenance.source_path))
            if record.vector is not None and len(record.vector.value) != 3:
                diagnostics.append(ParseDiagnostic(ParseSeverity.ERROR, "vector_shape_invalid", "moment vector does not have three components", source=provenance.source_path))
    return ParseValidation(valid=not any(item.severity is ParseSeverity.ERROR for item in diagnostics), diagnostics=tuple(diagnostics))


def _read_manifest(root: Path) -> tuple[dict | None, ParseDiagnostic | None]:
    for name in _MANIFEST_NAMES:
        path = root / name
        if not path.is_file():
            continue
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            return None, ParseDiagnostic(ParseSeverity.WARNING, "manifest_unreadable", f"could not read run manifest: {exc}", source=str(path))
        if not isinstance(value, dict):
            return None, ParseDiagnostic(ParseSeverity.WARNING, "manifest_invalid", "run manifest is not a JSON object", source=str(path))
        return value, None
    return None, None


def _number(value: str) -> float:
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError(f"non-finite numeric value: {value}")
    return result


def _numbers(value: str, *, expected: int | None, source: Path, line_number: int) -> list[float]:
    tokens = value.split()
    if expected is not None and len(tokens) != expected:
        raise MalformedResultError(f"expected {expected} numeric values at {source}:{line_number}")
    if expected is None and not tokens:
        raise MalformedResultError(f"expected numeric values at {source}:{line_number}")
    if any(not _NUMBER_RE.match(token) for token in tokens):
        raise MalformedResultError(f"malformed numeric value at {source}:{line_number}")
    return [_number(token) for token in tokens]


def _site_index(value: str, source: Path, line_number: int) -> int:
    parsed = int(value)
    if parsed < 1:
        raise MalformedResultError(f"atom index must be one-based and positive at {source}:{line_number}")
    return parsed - 1


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise MissingResultFileError(f"could not read result file {path}: {exc}") from exc


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["parse_report", "parse_run", "parse_scf"]
