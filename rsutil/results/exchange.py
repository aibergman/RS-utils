"""Parsers for RS-LMTO magnetic-interaction result tables.

The supported text tables are deliberately small and strict.  Their ``r_ij``
columns are Cartesian displacements in the same coordinate scale as the
reference structure; they are never used as the canonical pair identity.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from rsutil.core import ExchangePair, PairKey, RSStructure
from rsutil.core.hashing import hash_json

from .errors import (
    InteractionMappingError,
    MalformedResultError,
    MissingResultFileError,
    SourceUnitMismatchError,
    UnsupportedResultFormatError,
)
from .models import (
    ExchangeSet,
    InteractionRecord,
    ParseDiagnostic,
    ParseSeverity,
    ParseValidation,
    ParserProvenance,
    SourceFile,
)
from .units import MILLI_RYDBERG, Quantity, UNKNOWN_UNIT, Unit, scalar, vector


_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_NUMBER_RE = re.compile(rf"^{_NUMBER}$")
_MANIFEST_NAME = "rsutil_exchange_manifest.json"
_SCHEMAS = frozenset({"rsutil-exchange-v1", "rsutil-exchange-v2", "rsutil-exchange-v3"})
_FABLE_V4_MRY_STEMS = frozenset(
    {
        "jij",
        "jijso",
        "jijfo",
        "dij",
        "dijso",
        "dijfo",
        "aij",
        "aijso",
        "aijfo",
    }
)
_TENSOR_STEMS = frozenset({"aij", "aijso", "aijfo"})


class _RawRow:
    def __init__(self, values: list[float], line: int, source: Path) -> None:
        self.values = values
        self.line = line
        self.source = source


def parse_exchange_results(
    path: str | Path | Sequence[str | Path] | None = None,
    manifest: str | Path | Mapping[str, Any] | None = None,
    structure: RSStructure | None = None,
    *,
    jij_path: str | Path | None = None,
    dmi_path: str | Path | None = None,
    tensor_path: str | Path | None = None,
    tolerance: float = 1e-6,
    strict_physics: bool = False,
) -> ExchangeSet:
    """Parse one or more supported interaction files into an ``ExchangeSet``.

    A directory is searched only for the explicitly supported conventional
    names.  When a job manifest is present beside a file it is used
    automatically.  Mapping is by the one-based ``itype/jtype`` columns and
    the Cartesian displacement against the manifest's canonical pair set;
    row order is never used as identity.
    """

    paths = _input_paths(path, jij_path=jij_path, dmi_path=dmi_path, tensor_path=tensor_path)
    if not paths:
        raise MissingResultFileError("no exchange result file was supplied")
    loaded_manifest = _load_manifest(manifest)
    if loaded_manifest is None:
        for candidate in paths:
            sibling = candidate.parent / _MANIFEST_NAME
            if sibling.is_file():
                loaded_manifest = _load_manifest(sibling)
                break

    records: dict[PairKey, InteractionRecord] = {}
    parse_diagnostics: list[ParseDiagnostic] = []
    source_files: list[SourceFile] = []
    for result_path in paths:
        kind, metadata, headers = _classify(result_path)
        raw_rows = _read_rows(result_path, kind=kind, headers=headers)
        if not raw_rows:
            raise UnsupportedResultFormatError(f"no supported interaction rows found in {result_path}")
        source_files.append(SourceFile(str(result_path), _sha256(result_path), f"interaction_{kind}"))
        for row in raw_rows:
            key = _map_row(
                row,
                loaded_manifest,
                structure,
                kind=kind,
                tolerance=tolerance,
            )
            record = _record_from_row(
                row,
                key,
                kind=kind,
                metadata=metadata,
                structure=structure,
            )
            previous = records.get(key)
            if previous is None:
                records[key] = record
            else:
                records[key] = _merge_record(previous, record)

    validation = validate_exchange_set(
        ExchangeSet(records=records),
        strict=strict_physics,
        atol=tolerance,
        rtol=1e-5,
    )
    combined_hash = _combined_hash(source_files)
    provenance = ParserProvenance(
        source_path=str(paths[0]),
        content_hash=combined_hash,
        format="rs-lmto-exchange",
        format_version="rsutil-exchange-result-v1",
        source_files=tuple(source_files),
    )
    return ExchangeSet(
        records=records,
        diagnostics=tuple(parse_diagnostics),
        validation=validation,
        provenance=provenance,
        run_manifest=loaded_manifest,
    )


def parse_exchange_file(
    path: str | Path,
    manifest: str | Path | Mapping[str, Any] | None = None,
    structure: RSStructure | None = None,
    *,
    tolerance: float = 1e-6,
    strict_physics: bool = False,
) -> ExchangeSet:
    """Descriptive single-file alias for :func:`parse_exchange_results`."""

    return parse_exchange_results(
        path,
        manifest=manifest,
        structure=structure,
        tolerance=tolerance,
        strict_physics=strict_physics,
    )


parse_exchange = parse_exchange_results
parse_interactions = parse_exchange_results
parse_exchange_set = parse_exchange_results


def write_exchange_json(result: ExchangeSet, path: str | Path, *, indent: int = 2) -> Path:
    """Write the versioned interaction set for database ingestion."""

    destination = Path(path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(result.to_json(indent=indent) + "\n", encoding="utf-8")
    return destination


def validate_exchange_set(
    result: ExchangeSet,
    *,
    strict: bool = False,
    atol: float = 1e-8,
    rtol: float = 1e-5,
) -> ParseValidation:
    """Validate only relationships justified by pairs present in the set.

    Reverse checks are applied only when the reverse key is actually present.
    By default numerical disagreement is a physics warning; ``strict=True``
    promotes it to an error.  No inversion or crystal-symmetry constraint is
    inferred from the interaction table alone.
    """

    diagnostics: list[ParseDiagnostic] = []
    seen: set[PairKey] = set()
    for key in sorted(result.records, key=lambda item: item.sort_key):
        if key in seen or key.reverse() not in result.records:
            continue
        reverse = key.reverse()
        seen.add(key)
        seen.add(reverse)
        first = result.records[key]
        second = result.records[reverse]
        severity = ParseSeverity.ERROR if strict else ParseSeverity.WARNING
        if first.J is not None and second.J is not None:
            if not _close(first.J.value, second.J.value, atol=atol, rtol=rtol):
                diagnostics.append(_reverse_diagnostic("isotropic_exchange", key, reverse, severity))
        if first.D is not None and second.D is not None:
            expected = tuple(-float(value) for value in first.D.value)
            if not _close(first.D.value, tuple(-float(value) for value in second.D.value), atol=atol, rtol=rtol):
                diagnostics.append(_reverse_diagnostic("dmi_vector", key, reverse, severity, expected=expected))
        if first.tensor is not None and second.tensor is not None:
            expected = tuple(tuple(float(value) for value in row) for row in np.asarray(first.tensor.value).T)
            actual = tuple(tuple(float(value) for value in row) for row in np.asarray(second.tensor.value))
            if not _close(expected, actual, atol=atol, rtol=rtol):
                diagnostics.append(_reverse_diagnostic("anisotropic_exchange_tensor_transpose", key, reverse, severity))
    return ParseValidation(valid=not any(item.severity is ParseSeverity.ERROR for item in diagnostics), diagnostics=tuple(diagnostics))


def _input_paths(
    path: str | Path | Sequence[str | Path] | None,
    *,
    jij_path: str | Path | None,
    dmi_path: str | Path | None,
    tensor_path: str | Path | None,
) -> tuple[Path, ...]:
    values: list[str | Path] = []
    if path is not None:
        if isinstance(path, (str, Path)):
            values.append(path)
        else:
            values.extend(path)
    values.extend(value for value in (jij_path, dmi_path, tensor_path) if value is not None)
    if len(values) == 1 and Path(values[0]).expanduser().is_dir():
        root = Path(values[0]).expanduser().resolve()
        names = (
            "jijso.out",
            "jijfo.out",
            "dijso.out",
            "dijfo.out",
            "aij.out",
            "aijso.out",
            "aijfo.out",
            "dmi.out",
            "dmi.dat",
            "exchange_tensor.out",
        )
        return tuple(root / name for name in names if (root / name).is_file())
    result = tuple(Path(value).expanduser().resolve() for value in values)
    for item in result:
        if not item.is_file():
            raise MissingResultFileError(f"exchange result file does not exist: {item}")
    return result


def _load_manifest(value: str | Path | Mapping[str, Any] | None) -> dict[str, Any] | None:
    if value is None:
        return None
    if isinstance(value, Mapping):
        data = dict(value)
    else:
        path = Path(value).expanduser()
        if path.is_dir():
            path = path / _MANIFEST_NAME
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise MalformedResultError(f"could not read exchange manifest {path}: {exc}") from exc
    if not isinstance(data, dict) or data.get("schema_version") not in _SCHEMAS:
        raise UnsupportedResultFormatError("unsupported or malformed exchange manifest schema")
    spec = data.get("spec")
    if spec is not None:
        if not isinstance(spec, Mapping) or not isinstance(data.get("spec_hash"), str):
            raise MalformedResultError("exchange manifest has an incomplete specification")
        if hash_json(spec) != data["spec_hash"]:
            raise MalformedResultError("exchange manifest specification hash is invalid")
        keys = spec.get("pair_keys")
        if not isinstance(keys, list):
            raise MalformedResultError("exchange manifest spec has no pair_keys")
        for name in ("canonical_pair_list", "pair_keys"):
            view = data.get(name)
            if view is not None and view != keys:
                raise MalformedResultError(
                    f"exchange manifest {name} does not match spec.pair_keys"
                )
        if data.get("legacy_ijpair") is not None and data["legacy_ijpair"] != spec.get("legacy_ijpair"):
            raise MalformedResultError("exchange manifest legacy_ijpair does not match spec")
        pair_records = spec.get("pairs")
        if pair_records is not None:
            if not isinstance(pair_records, list) or len(pair_records) != len(keys):
                raise MalformedResultError("exchange manifest spec has invalid pair evidence")
            if data.get("pairs") is not None and data["pairs"] != pair_records:
                raise MalformedResultError("exchange manifest pairs do not match spec")
            for key, record in zip(keys, pair_records):
                if not isinstance(record, Mapping) or record.get("key") != key:
                    raise MalformedResultError("exchange manifest pair evidence does not match pair_keys")
                for name in ("rs_type_i", "rs_type_j"):
                    type_index = record.get(name)
                    if type_index is not None and (
                        isinstance(type_index, bool)
                        or not isinstance(type_index, int)
                        or type_index < 1
                    ):
                        raise MalformedResultError(f"exchange manifest contains invalid {name}")
    else:
        # Older hand-authored parser manifests may contain only one pair-list
        # view.  If compatibility views are duplicated, they must still agree.
        keys = data.get("canonical_pair_list", data.get("pair_keys"))
        if data.get("canonical_pair_list") is not None and data.get("pair_keys") is not None:
            if data["canonical_pair_list"] != data["pair_keys"]:
                raise MalformedResultError("exchange manifest pair-list views disagree")
    if not isinstance(keys, list):
        raise MalformedResultError("exchange manifest has no canonical_pair_list")
    try:
        data["_canonical_keys"] = tuple(PairKey.from_dict(item) for item in keys)
    except (KeyError, TypeError, ValueError) as exc:
        raise MalformedResultError("exchange manifest contains invalid canonical pair keys") from exc
    pair_records = spec.get("pairs") if isinstance(spec, Mapping) else None
    if pair_records is None:
        pair_records = [{"key": key} for key in keys]
    data["_pair_records"] = tuple(dict(record) for record in pair_records)
    return data


def _classify(path: Path) -> tuple[str, dict[str, Any], tuple[str, ...]]:
    stem = path.stem.lower()
    headers = tuple(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not _is_numeric_line(line))
    header_text = " ".join(headers).lower()
    if "dmi" in stem or stem.startswith("dij") or "d_ij" in header_text:
        kind = "dmi"
    elif stem in _TENSOR_STEMS or "tensor" in stem or "tensor" in header_text or "full exchange" in header_text:
        kind = "tensor"
    elif "j_ij" in header_text or stem.startswith("jij"):
        kind = "exchange"
    else:
        raise UnsupportedResultFormatError(
            f"cannot identify supported interaction format from filename/header: {path.name}"
        )
    suffix = None
    match = re.search(r"(?:jij|dij|aij)(so|fo)$", stem)
    if match:
        suffix = match.group(1)
    declared_units = re.findall(r"\b(mry|ry|rydberg|ev|mev|hartree|ha)\b", header_text, re.I)
    declared_unit_names = {unit.lower() for unit in declared_units}
    if declared_unit_names and not declared_unit_names.issubset({"ry", "rydberg", "mry"}):
        raise UnsupportedResultFormatError(
            f"interaction file declares an energy unit without a supported conversion: {path}"
        )
    normalized_declared_units = {"mry" if unit == "mry" else "ry" for unit in declared_unit_names}
    if len(normalized_declared_units) > 1:
        raise UnsupportedResultFormatError(f"interaction file declares conflicting energy units: {path}")
    if stem in _FABLE_V4_MRY_STEMS:
        if normalized_declared_units and normalized_declared_units != {"mry"}:
            raise SourceUnitMismatchError(
                f"SOURCE_UNIT_MISMATCH: {path.name} is a known fable_v4 interaction format "
                f"written in mRy, but declares {', '.join(sorted(normalized_declared_units))}"
            )
        energy_unit = MILLI_RYDBERG
        energy_unit_evidence = "fable_v4_exchange_writer_scaling"
    elif normalized_declared_units:
        energy_unit = MILLI_RYDBERG if "mry" in normalized_declared_units else Unit.RY
        energy_unit_evidence = "explicit_file_header"
    else:
        energy_unit = UNKNOWN_UNIT
        energy_unit_evidence = "unverified_source_unit"
    metadata: dict[str, Any] = {
        "format": kind,
        "channel": suffix,
        "method": None,
        "soc": None,
        "energy_unit": energy_unit.value,
        "energy_unit_evidence": energy_unit_evidence,
    }
    for name in ("method", "channel"):
        match = re.search(rf"\b{name}\s*[:=]\s*([A-Za-z0-9_.+-]+)", header_text, re.I)
        if match:
            metadata[name] = match.group(1)
    match = re.search(r"\bsoc\s*[:=]\s*(true|false|yes|no|on|off|1|0)", header_text, re.I)
    if match:
        metadata["soc"] = match.group(1).lower() in {"true", "yes", "on", "1"}
    return kind, metadata, headers


def _read_rows(path: Path, *, kind: str, headers: Sequence[str]) -> list[_RawRow]:
    rows: list[_RawRow] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("!"):
            continue
        if not _is_numeric_line(stripped):
            continue
        tokens = stripped.split()
        if any(not _NUMBER_RE.match(token) for token in tokens):
            raise MalformedResultError(f"malformed numeric row at {path}:{line_number}")
        values = [_number(token) for token in tokens]
        expected = {"exchange": (7, 14, 15), "dmi": (9,), "tensor": (14, 15)}[kind]
        if len(values) not in expected:
            raise MalformedResultError(
                f"unsupported {kind} row with {len(values)} columns at {path}:{line_number}; expected {expected}"
            )
        if kind == "exchange" and len(values) in (14, 15):
            # A 14/15-column J table is a tensor table only when the format is
            # identified by a tensor filename/header; otherwise fail closed.
            raise MalformedResultError(
                f"tensor-shaped row requires an explicit tensor filename/header at {path}:{line_number}"
            )
        rows.append(_RawRow(values, line_number, path))
    if kind == "exchange" and any(len(row.values) in (14, 15) for row in rows):
        raise MalformedResultError(f"tensor-shaped rows are not valid in scalar exchange file {path}")
    return rows


def _map_row(
    row: _RawRow,
    manifest: Mapping[str, Any] | None,
    structure: RSStructure | None,
    *,
    kind: str,
    tolerance: float,
) -> PairKey:
    values = row.values
    i_one, j_one = _positive_int(values[0], row, "itype"), _positive_int(values[1], row, "jtype")
    observed_distance = values[-1] if kind != "tensor" or len(values) == 15 else None
    if manifest is not None:
        records = list(manifest["_pair_records"])
        typed = any(record.get("rs_type_i") is not None for record in records)
        if typed:
            records = [
                record
                for record in records
                if record.get("rs_type_i") == i_one and record.get("rs_type_j") == j_one
            ]
        else:
            # v1/v2 manifests predate the source-faithful type contract.  Their
            # site-number mapping remains readable only as a compatibility path.
            records = [
                record for record in records
                if PairKey.from_dict(record["key"]).i + 1 == i_one
                and PairKey.from_dict(record["key"]).j + 1 == j_one
            ]
        if structure is not None:
            if typed and any(
                type_index is not None
                and type_index != structure.sites[index].rs_type_index
                for record in records
                for index, type_index in (
                    (PairKey.from_dict(record["key"]).i, record.get("rs_type_i")),
                    (PairKey.from_dict(record["key"]).j, record.get("rs_type_j")),
                )
            ):
                raise InteractionMappingError(
                    f"row at {row.source}:{row.line} conflicts with structure RS type evidence"
                )
            records = [
                record for record in records
                if _geometry_matches(
                    PairKey.from_dict(record["key"]),
                    values[2:5],
                    observed_distance,
                    structure,
                    tolerance,
                )
            ]
        elif typed:
            records = [
                record for record in records
                if _manifest_geometry_matches(record, values[2:5], observed_distance, tolerance)
            ]
        candidates = [PairKey.from_dict(record["key"]) for record in records]
    elif structure is not None:
        candidates = _structure_candidates(
            structure, i_one, j_one, values[2:5], observed_distance, tolerance
        )
    else:
        candidates = []
    if len(candidates) != 1:
        if not candidates:
            raise InteractionMappingError(
                f"row at {row.source}:{row.line} cannot be mapped uniquely to the exchange manifest"
            )
        raise InteractionMappingError(
            f"row at {row.source}:{row.line} maps ambiguously to canonical pairs: {candidates}"
        )
    return candidates[0]


def _structure_candidates(
    structure: RSStructure,
    i_one: int,
    j_one: int,
    displacement: Sequence[float],
    distance: float | None,
    tolerance: float,
) -> list[PairKey]:
    """Find site pairs by source type and geometry without assuming site=type."""

    all_types = [site.rs_type_index for site in structure.sites]
    if any(type_index is not None for type_index in all_types):
        i_sites = [index for index, type_index in enumerate(all_types) if type_index == i_one]
        j_sites = [index for index, type_index in enumerate(all_types) if type_index == j_one]
    else:
        # Non-RS structures retain the historical site-index fallback.
        i_sites = [i_one - 1] if 1 <= i_one <= len(structure.sites) else []
        j_sites = [j_one - 1] if 1 <= j_one <= len(structure.sites) else []
    candidates: list[PairKey] = []
    for i_index in i_sites:
        for j_index in j_sites:
            try:
                key = _infer_key_from_geometry(
                    i_index, j_index, displacement, structure, tolerance
                )
            except InteractionMappingError:
                continue
            if _geometry_matches(key, displacement, distance, structure, tolerance):
                candidates.append(key)
    return list(dict.fromkeys(candidates))


def _manifest_geometry_matches(
    record: Mapping[str, Any],
    displacement: Sequence[float],
    distance: float | None,
    tolerance: float,
) -> bool:
    expected = record.get("cartesian_displacement")
    if expected is None:
        return True
    if not np.allclose(np.asarray(expected, dtype=float), np.asarray(displacement, dtype=float), atol=tolerance, rtol=0):
        return False
    recorded_distance = record.get("distance")
    if distance is not None and recorded_distance is not None:
        if not math.isclose(float(distance), float(recorded_distance), abs_tol=tolerance, rel_tol=1e-5):
            return False
    return True


def _record_from_row(
    row: _RawRow,
    key: PairKey,
    *,
    kind: str,
    metadata: Mapping[str, Any],
    structure: RSStructure | None,
) -> InteractionRecord:
    values = row.values
    pair = None
    if structure is not None:
        displacement = _pair_displacement(structure, key)
        pair = ExchangePair(
            key=key,
            cartesian_displacement=displacement,
            distance=float(np.linalg.norm(displacement)),
            species_i=structure.sites[key.i].element,
            species_j=structure.sites[key.j].element,
            length_unit=structure.length_unit,
            label_i=structure.sites[key.i].rs_label,
            label_j=structure.sites[key.j].rs_label,
            rs_type_i=structure.sites[key.i].rs_type_index,
            rs_type_j=structure.sites[key.j].rs_type_index,
            provenance={"geometry_source": "reference_structure"},
        )
    quantities: dict[str, Quantity | None] = {
        "isotropic_exchange": None,
        "dmi_vector": None,
        "anisotropic_exchange_tensor": None,
    }
    if kind == "exchange":
        if len(values) != 7:
            raise MalformedResultError(f"scalar exchange row has unexpected shape at {row.source}:{row.line}")
        quantities["isotropic_exchange"] = scalar(values[5], metadata["energy_unit"])
    elif kind == "dmi":
        quantities["dmi_vector"] = vector(values[5:8], metadata["energy_unit"])
    else:
        tensor_values = values[5:14]
        # A whole Fortran array is emitted in column-major storage order.
        # Thus flat values 1..9 mean (1,1),(2,1),(3,1),(1,2),... .
        quantities["anisotropic_exchange_tensor"] = Quantity(
            tuple(
                tuple(float(tensor_values[row + 3 * column]) for column in range(3))
                for row in range(3)
            ),
            metadata["energy_unit"],
        )
    return InteractionRecord(
        key=key,
        pair=pair,
        method=metadata.get("method"),
        channel=metadata.get("channel"),
        soc=metadata.get("soc"),
        provenance={
            "source_path": str(row.source),
            "line": row.line,
            "format": metadata["format"],
            "coordinate_columns": "cartesian_displacement",
            "energy_unit": metadata["energy_unit"],
            "energy_unit_evidence": metadata["energy_unit_evidence"],
        },
        **quantities,
    )


def _merge_record(first: InteractionRecord, second: InteractionRecord) -> InteractionRecord:
    if first.channel != second.channel or first.method != second.method or first.soc != second.soc:
        raise MalformedResultError(
            f"multiple interaction channels/methods map to the same PairKey {first.key}; parse them as separate datasets"
        )
    fields = ("isotropic_exchange", "dmi_vector", "anisotropic_exchange_tensor")
    for field_name in fields:
        if getattr(first, field_name) is not None and getattr(second, field_name) is not None:
            raise MalformedResultError(f"duplicate {field_name} for canonical pair {first.key}")
    return replace(
        first,
        isotropic_exchange=first.isotropic_exchange or second.isotropic_exchange,
        dmi_vector=first.dmi_vector or second.dmi_vector,
        anisotropic_exchange_tensor=(
            first.anisotropic_exchange_tensor or second.anisotropic_exchange_tensor
        ),
        provenance={"sources": [first.provenance, second.provenance]},
    )


def _infer_key_from_geometry(i: int, j: int, displacement: Sequence[float], structure: RSStructure, tolerance: float) -> PairKey:
    lattice = np.asarray(structure.lattice_vectors, dtype=float)
    base = _pair_displacement(structure, PairKey(i, j, (0, 0, 0)))
    translation = (np.asarray(displacement, dtype=float) - base) @ np.linalg.inv(lattice)
    rounded = np.rint(translation).astype(int)
    if not np.allclose(translation, rounded, atol=tolerance, rtol=0):
        raise InteractionMappingError("Cartesian exchange displacement does not imply an integer lattice translation")
    if any(not periodic and value != 0 for periodic, value in zip(structure.pbc, rounded)):
        raise InteractionMappingError("exchange displacement uses a non-periodic lattice translation")
    return PairKey(i, j, tuple(int(value) for value in rounded))


def _geometry_matches(
    key: PairKey,
    displacement: Sequence[float],
    distance: float | None,
    structure: RSStructure,
    tolerance: float,
) -> bool:
    expected = _pair_displacement(structure, key)
    observed = np.asarray(displacement, dtype=float)
    if not np.allclose(expected, observed, atol=tolerance, rtol=0):
        return False
    return distance is None or math.isclose(
        float(np.linalg.norm(observed)), float(distance), abs_tol=tolerance, rel_tol=1e-5
    )


def _pair_displacement(structure: RSStructure, key: PairKey) -> tuple[float, float, float]:
    lattice = np.asarray(structure.lattice_vectors, dtype=float)
    positions = np.asarray(structure.positions, dtype=float)
    if structure.coordinate_convention.value == "fractional":
        positions = positions @ lattice
    displacement = positions[key.j] + np.asarray(key.R, dtype=float) @ lattice - positions[key.i]
    return tuple(float(value) for value in displacement)


def _positive_int(value: float, row: _RawRow, label: str) -> int:
    integer = int(value)
    if value != integer or integer < 1:
        raise MalformedResultError(f"{label} must be a positive one-based integer at {row.source}:{row.line}")
    return integer


def _is_numeric_line(line: str) -> bool:
    first = line.split()[0] if line.split() else ""
    return bool(_NUMBER_RE.match(first))


def _number(value: str) -> float:
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError(f"non-finite numeric value: {value}")
    return result


def _close(first: Any, second: Any, *, atol: float, rtol: float) -> bool:
    return bool(np.allclose(np.asarray(first, dtype=float), np.asarray(second, dtype=float), atol=atol, rtol=rtol))


def _reverse_diagnostic(field: str, key: PairKey, reverse: PairKey, severity: ParseSeverity, *, expected: Any = None) -> ParseDiagnostic:
    details = {"pair": key.to_dict(), "reverse_pair": reverse.to_dict()}
    if expected is not None:
        details["expected_reverse"] = expected
    return ParseDiagnostic(
        severity=severity,
        code="reverse_pair_inconsistent",
        message=f"{field} is inconsistent for reverse pairs {key} and {reverse}",
        details=details,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _combined_hash(source_files: Sequence[SourceFile]) -> str:
    digest = hashlib.sha256()
    for item in source_files:
        digest.update(item.path.encode("utf-8"))
        digest.update(item.content_hash.encode("ascii"))
    return digest.hexdigest()


__all__ = [
    "parse_exchange",
    "parse_exchange_file",
    "parse_exchange_results",
    "parse_exchange_set",
    "parse_interactions",
    "validate_exchange_set",
    "write_exchange_json",
]
