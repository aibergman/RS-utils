"""Lightweight Slurm, status, and result-collection helpers for campaigns.

The functions in this module deliberately operate on the relocatable
``campaign_index.json`` produced by HT-16.  They render ordinary shell/Slurm
artifacts and inspect files; they do not submit jobs, query Slurm, retry work,
or model dependencies.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shlex
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..core.hashing import hash_file
from ..results import parse_exchange_results, parse_run
from .generation import CAMPAIGN_INDEX_FILENAME
from .models import CampaignIndex, CampaignIndexEntry


SLURM_SCHEMA_VERSION = "rsutil-slurm-array-v1"
CAMPAIGN_RESULTS_SCHEMA_VERSION = "rsutil-campaign-results-v1"
RUN_LIST_FILENAME = "slurm_run_list.tsv"
SLURM_SCRIPT_FILENAME = "submit_array.sbatch"
CAMPAIGN_RESULTS_FILENAME = "campaign_results.json"

GENERATED = "GENERATED"
LOCKED = "LOCKED"
COMPLETED_VALID = "COMPLETED_VALID"
COMPLETED_INVALID = "COMPLETED_INVALID"
FAILED = "FAILED"
UNKNOWN = "UNKNOWN"

_LOCK_NAMES = (".rsutil-running", ".rsutil-running.lock", "run.lock", "RUNNING")
_FAILURE_NAMES = ("FAILED", "FAILED.json", "run.failed", "rsutil_failed.json")
_EXCHANGE_MANIFEST = "rsutil_exchange_manifest.json"
_EXCHANGE_RESULT_NAMES = {
    "jijso.out",
    "jijfo.out",
    "dijso.out",
    "dijfo.out",
    "dmi.out",
    "dmi.dat",
    "exchange_tensor.out",
}


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent), text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            stream.write(content)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    _atomic_write(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _load_index(value: CampaignIndex | str | Path, *, campaign_root: str | Path | None = None) -> tuple[CampaignIndex, Path, Path]:
    if isinstance(value, CampaignIndex):
        if campaign_root is None:
            raise ValueError("campaign_root is required when passing a CampaignIndex object")
        root = Path(campaign_root).expanduser().resolve()
        return value, root, root / CAMPAIGN_INDEX_FILENAME

    supplied = Path(value).expanduser().resolve()
    index_path = supplied / CAMPAIGN_INDEX_FILENAME if supplied.is_dir() else supplied
    try:
        index = CampaignIndex.from_dict(json.loads(index_path.read_text(encoding="utf-8")))
    except Exception as exc:
        raise ValueError(f"could not read campaign index {index_path}: {exc}") from exc
    return index, index_path.parent, index_path


def load_campaign_index(value: str | Path) -> CampaignIndex:
    """Load a relocatable ``campaign_index.json`` from a file or directory."""

    return _load_index(value)[0]


def _index_path(root: Path, value: str, label: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"campaign index contains unsafe {label}: {value!r}")
    result = (root / relative).resolve()
    if result != root and root not in result.parents:
        raise ValueError(f"campaign index {label} escapes campaign root: {value!r}")
    return result


def _run_dir(root: Path, entry: CampaignIndexEntry) -> Path:
    return _index_path(root, entry.run_path, "run path")


def _manifest_path(root: Path, entry: CampaignIndexEntry) -> Path:
    return _index_path(root, entry.manifest_path, "manifest path")


def _command_tokens(value: str | Path | Sequence[str] | None, *, required: bool, name: str) -> list[str]:
    if value is None:
        if required:
            raise ValueError(f"{name} is required")
        return []
    if isinstance(value, Path):
        tokens = [str(value)]
    elif isinstance(value, str):
        try:
            tokens = shlex.split(value)
        except ValueError as exc:
            raise ValueError(f"invalid {name}: {exc}") from exc
    else:
        tokens = [str(item) for item in value]
    if not tokens or any("\x00" in token or "\n" in token or "\r" in token for token in tokens):
        raise ValueError(f"{name} must contain non-empty shell-safe tokens")
    return tokens


def _executable_tokens(value: str | Path | Sequence[str] | None) -> list[str]:
    """Treat an existing/path-like string as one executable token."""

    if isinstance(value, str) and ("/" in value or "\\" in value or Path(value).expanduser().is_file()):
        return _command_tokens([value], required=True, name="executable")
    return _command_tokens(value, required=True, name="executable")


def _preamble_lines(value: str | Sequence[str] | None) -> list[str]:
    if value is None:
        return []
    lines = value.splitlines() if isinstance(value, str) else [str(item) for item in value]
    if any("\x00" in line for line in lines):
        raise ValueError("environment preamble cannot contain NUL bytes")
    return lines


def _safe_directive(value: Any, name: str) -> str:
    text = str(value)
    if not text or any(char in text for char in "\r\n"):
        raise ValueError(f"Slurm {name} must be a non-empty single-line value")
    return text


def _job_name(value: Any) -> str:
    text = _safe_directive(value, "job-name")
    safe = re.sub(r"[^A-Za-z0-9._-]+", "-", text).strip("-.")
    return safe or "rsutil-campaign"


@dataclass(frozen=True)
class SlurmRenderResult:
    campaign_root: Path
    run_list_path: Path
    script_path: Path
    array_size: int
    schema_version: str = SLURM_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "campaign_root": str(self.campaign_root),
            "run_list_path": str(self.run_list_path),
            "script_path": str(self.script_path),
            "array_size": self.array_size,
        }


def render_slurm_array(
    campaign_index: CampaignIndex | str | Path,
    *,
    executable: str | Path | Sequence[str] | None = None,
    launcher: str | Sequence[str] | None = None,
    arguments: Sequence[str] = (),
    campaign_root: str | Path | None = None,
    output_dir: str | Path | None = None,
    run_list_name: str = RUN_LIST_FILENAME,
    script_name: str = SLURM_SCRIPT_FILENAME,
    job_name: str | None = None,
    partition: str | None = None,
    account: str | None = None,
    time_limit: str | None = None,
    resources: Mapping[str, Any] | None = None,
    environment_preamble: str | Sequence[str] | None = None,
) -> SlurmRenderResult:
    """Render a stable run list and readable ``sbatch`` array script.

    ``executable`` and ``launcher`` are render parameters.  No site module,
    MPI launcher, or RS-LMTO installation path is selected by this function.
    The output script executes one command in the selected run directory and
    exports the array mapping as ``RSUTIL_CAMPAIGN_*`` variables.
    """

    index, root, _ = _load_index(campaign_index, campaign_root=campaign_root)
    destination = root if output_dir is None else Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    run_list_path = destination / run_list_name
    script_path = destination / script_name
    if Path(run_list_name).name != run_list_name or Path(script_name).name != script_name:
        raise ValueError("run-list and script names must be simple file names")
    executable_tokens = _executable_tokens(executable)
    launcher_tokens = _command_tokens(launcher, required=False, name="launcher")
    argument_tokens = [str(item) for item in arguments]
    if any("\x00" in token or "\n" in token or "\r" in token for token in argument_tokens):
        raise ValueError("command arguments cannot contain NUL or newline characters")
    preamble = _preamble_lines(environment_preamble)
    resource_values = dict(resources or {})

    rows: list[str] = ["array_index\tentry_id\trun_path\tmanifest_path"]
    for array_index, entry in enumerate(index.entries):
        _run_dir(root, entry)
        _manifest_path(root, entry)
        for value in (entry.entry_id, entry.run_path, entry.manifest_path):
            if any(character in value for character in "\t\r\n"):
                raise ValueError("campaign index paths and IDs cannot contain tabs or newlines")
        rows.append(f"{array_index}\t{entry.entry_id}\t{entry.run_path}\t{entry.manifest_path}")
    _atomic_write(run_list_path, "\n".join(rows) + "\n")

    job = _job_name(job_name or index.campaign_name)
    directives = [f"#SBATCH --job-name={job}", f"#SBATCH --array=0-{len(index.entries) - 1}"]
    for option, value in (("partition", partition), ("account", account), ("time", time_limit)):
        if value is not None:
            directives.append(f"#SBATCH --{option}={_safe_directive(value, option)}")
    known_resources = {
        "cpus_per_task": "cpus-per-task",
        "mem": "mem",
        "gpus": "gpus",
        "nodes": "nodes",
        "ntasks": "ntasks",
        "ntasks_per_node": "ntasks-per-node",
        "constraint": "constraint",
        "qos": "qos",
    }
    for key, option in known_resources.items():
        if key in resource_values and resource_values[key] is not None:
            directives.append(f"#SBATCH --{option}={_safe_directive(resource_values[key], option)}")
    for key, value in resource_values.items():
        if key not in known_resources:
            option = key.replace("_", "-")
            directives.append(f"#SBATCH --{option}={_safe_directive(value, option)}")

    command = launcher_tokens + executable_tokens + argument_tokens
    quoted_command = " ".join(shlex.quote(token) for token in command)
    root_reference = shlex.quote(root.resolve().as_posix())
    lines = [
        "#!/usr/bin/env bash",
        *directives,
        "set -euo pipefail",
        "",
    ]
    if preamble:
        lines.extend(["# User-supplied environment preamble"] + preamble + [""])
    lines.extend(
        [
            f"CAMPAIGN_ROOT={root_reference}",
            f"RUN_LIST={shlex.quote(run_list_path.resolve().as_posix())}",
            'TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"',
            'if ! [[ "$TASK_ID" =~ ^[0-9]+$ ]]; then echo "invalid array task ID: $TASK_ID" >&2; exit 2; fi',
            'LINE_NUMBER=$((TASK_ID + 2))',
            'IFS=$\'\\t\' read -r ARRAY_INDEX ENTRY_ID RUN_PATH MANIFEST_PATH < <(sed -n "${LINE_NUMBER}p" "$RUN_LIST")',
            'if [[ -z "${RUN_PATH:-}" || "$ARRAY_INDEX" != "$TASK_ID" ]]; then echo "array task is not present in run list: $TASK_ID" >&2; exit 2; fi',
            'RUN_DIR="$CAMPAIGN_ROOT/$RUN_PATH"',
            'if [[ ! -d "$RUN_DIR" ]]; then echo "run directory does not exist: $RUN_DIR" >&2; exit 3; fi',
            'export RSUTIL_CAMPAIGN_ARRAY_INDEX="$ARRAY_INDEX"',
            'export RSUTIL_CAMPAIGN_ENTRY_ID="$ENTRY_ID"',
            'export RSUTIL_CAMPAIGN_RUN_DIR="$RUN_DIR"',
            'export RSUTIL_CAMPAIGN_MANIFEST_PATH="$CAMPAIGN_ROOT/$MANIFEST_PATH"',
            'cd -- "$RUN_DIR"',
            f"echo \"rsutil campaign entry $ENTRY_ID (array $ARRAY_INDEX) in $RUN_DIR\"",
            f"exec {quoted_command}",
            "",
        ]
    )
    _atomic_write(script_path, "\n".join(lines))
    script_path.chmod(script_path.stat().st_mode | 0o111)
    return SlurmRenderResult(root, run_list_path, script_path, len(index.entries))


render_slurm = render_slurm_array


@dataclass(frozen=True)
class CampaignStatusEntry:
    array_index: int
    entry_id: str
    name: str
    run_path: str
    manifest_path: str
    manifest_hash: str
    status: str
    evidence: str
    parser_error: Mapping[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "array_index": self.array_index,
            "entry_id": self.entry_id,
            "name": self.name,
            "run_path": self.run_path,
            "manifest_path": self.manifest_path,
            "manifest_hash": self.manifest_hash,
            "status": self.status,
            "evidence": self.evidence,
            "parser_error": None if self.parser_error is None else dict(self.parser_error),
        }


@dataclass(frozen=True)
class CampaignStatusReport:
    campaign_id: str
    campaign_name: str
    entries: tuple[CampaignStatusEntry, ...]
    scheduler_state: str = "not_queried"
    schema_version: str = "rsutil-campaign-status-v1"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "campaign_id": self.campaign_id,
            "campaign_name": self.campaign_name,
            "scheduler_state": self.scheduler_state,
            "entries": [entry.to_dict() for entry in self.entries],
        }


def _error_payload(exc: BaseException) -> dict[str, Any]:
    code = getattr(exc, "error_code", None) or getattr(exc, "code", None) or exc.__class__.__name__.lower()
    return {"code": str(code), "type": exc.__class__.__name__, "message": str(exc)}


def _read_json(path: Path) -> Mapping[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError(f"JSON document is not an object: {path}")
    return data


def _manifest_matches(manifest: Mapping[str, Any], entry: CampaignIndexEntry) -> bool:
    campaign = manifest.get("campaign")
    return isinstance(campaign, Mapping) and campaign.get("entry_id") == entry.entry_id and campaign.get("manifest_hash") == entry.manifest_hash


def _first_marker(run_dir: Path, names: Sequence[str]) -> Path | None:
    for name in names:
        path = run_dir / name
        if path.exists():
            return path
    return None


def _status_entry(root: Path, array_index: int, entry: CampaignIndexEntry) -> CampaignStatusEntry:
    run_dir = _run_dir(root, entry)
    manifest_path = _manifest_path(root, entry)
    manifest: Mapping[str, Any] | None = None
    manifest_error: dict[str, Any] | None = None
    if manifest_path.is_file():
        try:
            manifest = _read_json(manifest_path)
            if not _manifest_matches(manifest, entry):
                manifest_error = {"code": "campaign_manifest_mismatch", "message": "manifest does not match the campaign index"}
        except Exception as exc:
            manifest_error = _error_payload(exc)
    else:
        manifest_error = {"code": "campaign_manifest_missing", "message": f"manifest is missing: {manifest_path}"}

    failure_marker = _first_marker(run_dir, _FAILURE_NAMES) if run_dir.is_dir() else None
    if failure_marker is not None:
        return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, FAILED, f"failure marker {failure_marker.name}", manifest_error)
    lock_marker = _first_marker(run_dir, _LOCK_NAMES) if run_dir.is_dir() else None
    if lock_marker is not None:
        return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, LOCKED, f"lock marker {lock_marker.name}", manifest_error)
    report_path = run_dir / "report.out"
    if report_path.is_file():
        if manifest_error is not None:
            return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, COMPLETED_INVALID, "result exists but campaign manifest is invalid", manifest_error)
        try:
            result = parse_run(run_dir)
            status = COMPLETED_VALID if result.validation.valid else COMPLETED_INVALID
            evidence = "report.out parsed and validation passed" if status == COMPLETED_VALID else "report.out parsed with validation errors"
            return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, status, evidence)
        except Exception as exc:
            return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, COMPLETED_INVALID, "report.out is malformed or unsupported", _error_payload(exc))
    if manifest_error is not None:
        return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, UNKNOWN, "campaign manifest missing or invalid", manifest_error)
    return CampaignStatusEntry(array_index, entry.entry_id, entry.name, entry.run_path, entry.manifest_path, entry.manifest_hash, GENERATED, "owned campaign run generated; no result evidence")


def status_campaign(campaign_index: CampaignIndex | str | Path, *, campaign_root: str | Path | None = None) -> CampaignStatusReport:
    """Classify runs using local evidence only; no scheduler state is queried."""

    index, root, _ = _load_index(campaign_index, campaign_root=campaign_root)
    return CampaignStatusReport(
        campaign_id=index.campaign_id,
        campaign_name=index.campaign_name,
        entries=tuple(_status_entry(root, position, entry) for position, entry in enumerate(index.entries)),
    )


campaign_status = status_campaign


def _relative(root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _raw_files(root: Path, run_dir: Path, excluded: set[Path]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    if not run_dir.is_dir():
        return result
    for path in sorted(item for item in run_dir.rglob("*") if item.is_file()):
        resolved = path.resolve()
        if resolved in excluded or any(directory in resolved.parents for directory in excluded if directory.is_dir()):
            continue
        try:
            result.append({"path": _relative(root, path), "sha256": hash_file(path), "size": path.stat().st_size})
        except OSError as exc:
            result.append({"path": _relative(root, path), "error": _error_payload(exc)})
    return result


def _logicalize(value: Any, root: Path) -> Any:
    if isinstance(value, Mapping):
        return {key: _logicalize(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_logicalize(item, root) for item in value]
    if isinstance(value, str) and value.startswith("/"):
        try:
            return Path(value).resolve().relative_to(root.resolve()).as_posix()
        except ValueError:
            return value
    return value


def _scalar_summary(result: Any) -> dict[str, Any]:
    summary: dict[str, Any] = {"convergence": getattr(result.convergence, "value", result.convergence)}
    for field_name in ("fermi_energy", "total_energy", "band_energy", "total_rms_difference"):
        quantity = getattr(result, field_name, None)
        if quantity is not None and not isinstance(quantity.value, (list, tuple, dict)):
            summary[field_name] = quantity.to_dict()
        else:
            summary[field_name] = None
    summary["validation_valid"] = bool(result.validation.valid)
    return summary


def _exchange_subjobs(root: Path, run_dir: Path, entry_id: str, sidecar_root: Path) -> list[dict[str, Any]]:
    subjobs: list[dict[str, Any]] = []
    for manifest_path in sorted(run_dir.rglob(_EXCHANGE_MANIFEST)) if run_dir.is_dir() else ():
        item: dict[str, Any] = {
            "parent_run_id": entry_id,
            "manifest_path": _relative(root, manifest_path),
            "manifest_hash": hash_file(manifest_path),
        }
        try:
            manifest = _read_json(manifest_path)
            pairs = manifest.get("canonical_pair_list", manifest.get("pair_keys", []))
            if not isinstance(pairs, list):
                raise ValueError("exchange manifest canonical pair list is not a list")
            item["canonical_pair_list"] = pairs
            result_paths = sorted(path for path in manifest_path.parent.iterdir() if path.is_file() and path.name.lower() in _EXCHANGE_RESULT_NAMES)
            if result_paths:
                parsed = parse_exchange_results(result_paths, manifest=manifest_path)
                sidecar = sidecar_root / "exchange" / entry_id / f"{manifest_path.parent.name}-{hash_file(manifest_path)[:16]}.json"
                _atomic_json(sidecar, _logicalize(parsed.to_dict(), root))
                item["result_sidecar"] = _relative(root, sidecar)
                item["result_files"] = [_relative(root, path) for path in result_paths]
            else:
                item["result_files"] = []
        except Exception as exc:
            item["parse_error"] = _error_payload(exc)
        subjobs.append(item)
    return subjobs


@dataclass(frozen=True)
class CampaignCollection:
    aggregate_path: Path
    entries: tuple[Mapping[str, Any], ...]
    csv_path: Path | None = None
    jsonl_path: Path | None = None
    schema_version: str = CAMPAIGN_RESULTS_SCHEMA_VERSION

    @property
    def failures(self) -> tuple[Mapping[str, Any], ...]:
        return tuple(entry for entry in self.entries if entry.get("parse_error") or entry.get("validation", {}).get("valid") is False)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "aggregate_path": str(self.aggregate_path),
            "entries": [dict(entry) for entry in self.entries],
            "csv_path": None if self.csv_path is None else str(self.csv_path),
            "jsonl_path": None if self.jsonl_path is None else str(self.jsonl_path),
        }


def collect_campaign(
    campaign_index: CampaignIndex | str | Path,
    *,
    campaign_root: str | Path | None = None,
    output: str | Path | None = None,
    csv_output: str | Path | None = None,
    jsonl_output: str | Path | None = None,
    sidecar_dir: str | Path | None = None,
) -> CampaignCollection:
    """Parse every run independently and write provenance-preserving results.

    A malformed or missing result is recorded on its entry and does not stop
    collection of later entries.  Full parsed models are written to JSON
    sidecars, leaving arrays and tensors structured rather than flattened into
    a scalar table.
    """

    index, root, index_path = _load_index(campaign_index, campaign_root=campaign_root)
    aggregate_path = (root / CAMPAIGN_RESULTS_FILENAME) if output is None else Path(output).expanduser().resolve()
    sidecars = (root / "collected") if sidecar_dir is None else Path(sidecar_dir).expanduser().resolve()
    sidecars.mkdir(parents=True, exist_ok=True)
    status_report = status_campaign(index, campaign_root=root)
    status_by_id = {entry.entry_id: entry for entry in status_report.entries}
    aggregate_exclusions = {aggregate_path.resolve(), sidecars.resolve()}
    entries: list[dict[str, Any]] = []
    for array_index, index_entry in enumerate(index.entries):
        run_dir = _run_dir(root, index_entry)
        status = status_by_id[index_entry.entry_id]
        manifest_path = _manifest_path(root, index_entry)
        manifest_content_hash = hash_file(manifest_path) if manifest_path.is_file() else None
        collected: dict[str, Any] = {
            "array_index": array_index,
            "entry_id": index_entry.entry_id,
            "name": index_entry.name,
            "run_path": index_entry.run_path,
            "manifest_path": index_entry.manifest_path,
            "manifest_hash": index_entry.manifest_hash,
            "manifest_content_hash": manifest_content_hash,
            "status": status.status,
            "status_evidence": status.evidence,
            "raw_files": _raw_files(root, run_dir, aggregate_exclusions),
            "exchange_subjobs": _exchange_subjobs(root, run_dir, index_entry.entry_id, sidecars),
        }
        report_path = run_dir / "report.out"
        if report_path.is_file():
            try:
                parsed = parse_run(run_dir)
                sidecar_path = sidecars / "results" / f"{index_entry.entry_id}.json"
                _atomic_json(sidecar_path, _logicalize(parsed.to_dict(), root))
                collected["result_sidecar"] = _relative(root, sidecar_path)
                collected["summary"] = _scalar_summary(parsed)
                collected["validation"] = parsed.validation.to_dict()
                collected["parser_provenance"] = _logicalize(parsed.provenance.to_dict(), root)
            except Exception as exc:
                collected["parse_error"] = _error_payload(exc)
        else:
            collected["parse_error"] = {"code": "missing_result_file", "type": "MissingResultFileError", "message": "report.out is not present"}
        entries.append(collected)

    aggregate: dict[str, Any] = {
        "schema_version": CAMPAIGN_RESULTS_SCHEMA_VERSION,
        "campaign_id": index.campaign_id,
        "campaign_name": index.campaign_name,
        "campaign_spec_hash": index.campaign_spec_hash,
        "collected_at": datetime.now(timezone.utc).isoformat(),
        "source_index": {"path": _relative(root, index_path), "sha256": hash_file(index_path)},
        "entries": entries,
    }
    _atomic_json(aggregate_path, aggregate)

    csv_path = None if csv_output is None else Path(csv_output).expanduser().resolve()
    if csv_path is not None:
        fields = ["array_index", "entry_id", "name", "run_path", "status", "convergence", "fermi_energy", "total_energy", "band_energy", "validation_valid", "result_sidecar", "parse_error"]
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", prefix=f".{csv_path.name}.", dir=str(csv_path.parent), delete=False) as stream:
            temporary_csv = Path(stream.name)
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for entry in entries:
                summary = entry.get("summary", {})
                writer.writerow({
                    "array_index": entry["array_index"], "entry_id": entry["entry_id"], "name": entry["name"], "run_path": entry["run_path"], "status": entry["status"],
                    "convergence": summary.get("convergence"),
                    "fermi_energy": json.dumps(summary.get("fermi_energy"), sort_keys=True),
                    "total_energy": json.dumps(summary.get("total_energy"), sort_keys=True),
                    "band_energy": json.dumps(summary.get("band_energy"), sort_keys=True),
                    "validation_valid": entry.get("validation", {}).get("valid"), "result_sidecar": entry.get("result_sidecar"), "parse_error": json.dumps(entry.get("parse_error"), sort_keys=True),
                })
        os.replace(temporary_csv, csv_path)

    jsonl_path = None if jsonl_output is None else Path(jsonl_output).expanduser().resolve()
    if jsonl_path is not None:
        _atomic_write(jsonl_path, "".join(json.dumps(entry, sort_keys=True) + "\n" for entry in entries))
    return CampaignCollection(aggregate_path, tuple(entries), csv_path, jsonl_path)


campaign_collect = collect_campaign


__all__ = [
    "CAMPAIGN_RESULTS_FILENAME",
    "CAMPAIGN_RESULTS_SCHEMA_VERSION",
    "COMPLETED_INVALID",
    "COMPLETED_VALID",
    "CampaignCollection",
    "CampaignStatusEntry",
    "CampaignStatusReport",
    "FAILED",
    "GENERATED",
    "LOCKED",
    "RUN_LIST_FILENAME",
    "SLURM_SCHEMA_VERSION",
    "SLURM_SCRIPT_FILENAME",
    "SlurmRenderResult",
    "UNKNOWN",
    "campaign_collect",
    "campaign_status",
    "collect_campaign",
    "load_campaign_index",
    "render_slurm",
    "render_slurm_array",
    "status_campaign",
]
