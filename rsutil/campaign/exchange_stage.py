"""A small, explicit exchange stage for generated campaigns.

This module is intentionally narrower than a workflow engine.  It selects
eligible parent runs, delegates exchange planning and job serialization to
``ExchangeService``, and records enough orchestration metadata to render and
inspect a Slurm array later.  Exchange-job manifests and result parsers remain
the scientific source of truth.
"""

from __future__ import annotations

import json
import math
import os
import shlex
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..core.hashing import hash_file, hash_json
from ..exchange.jobs import (
    EXCHANGE_MANIFEST_FILENAME,
    ExchangeJobSpec,
    read_exchange_manifest,
    write_exchange_jobs,
)
from ..exchange.reference import read_reference_calculation
from ..exchange.service import ExchangeService
from ..exchange.symmetry import (
    MAGNETIC_REFERENCE_COLLINEAR_FM,
    MagneticSymmetryContext,
)
from ..results import (
    UnsupportedResultFormatError,
    parse_exchange_results,
)
from .ht17 import (
    COMPLETED_INVALID,
    COMPLETED_VALID,
    EXECUTION_EXITED_NONZERO,
    EXECUTION_EXITED_ZERO,
    EXECUTION_NOT_RUN,
    EXECUTION_RUNNING_OR_UNKNOWN,
    EXECUTION_SCHEMA_VERSION,
    EXECUTION_SIGNALED,
    FAILED,
    GENERATED,
    OUTPUT_NO_OUTPUT,
    OUTPUT_PARSE_INVALID,
    OUTPUT_PARSE_UNSUPPORTED,
    OUTPUT_PARSE_VALID,
    UNKNOWN,
    _atomic_json,
    _atomic_write,
    _command_tokens,
    _error_payload,
    _executable_tokens,
    _job_name,
    _load_index,
    _preamble_lines,
    _safe_directive,
    status_campaign,
)
from .models import _validate_resource_values


EXCHANGE_STAGE_MANIFEST_FILENAME = "exchange_stage_manifest.json"
EXCHANGE_STAGE_MANIFEST_SCHEMA = "rsutil-exchange-stage-v1"
EXCHANGE_STAGE_RESULTS_FILENAME = "exchange_stage_results.json"
EXCHANGE_STAGE_RESULTS_SCHEMA = "rsutil-exchange-stage-results-v1"
EXCHANGE_STAGE_RUN_LIST_FILENAME = "exchange_slurm_run_list.tsv"
EXCHANGE_STAGE_SLURM_SCRIPT_FILENAME = "exchange_submit_array.sbatch"

EXCHANGE_OUTPUT_PRESENT = "OUTPUT_PRESENT"
EXCHANGE_PARSER_VALID = OUTPUT_PARSE_VALID
EXCHANGE_PARSER_INVALID = OUTPUT_PARSE_INVALID
EXCHANGE_PARSER_UNSUPPORTED = OUTPUT_PARSE_UNSUPPORTED


class ExchangeStageError(RuntimeError):
    """Base error for invalid or unsafe exchange-stage operations."""


class UnsupportedMagneticContextError(ExchangeStageError, ValueError):
    """A requested symmetry reduction lacks an established magnetic context."""


@dataclass(frozen=True)
class ExchangeStagePreparation:
    stage_root: Path
    manifest_path: Path
    stage_id: str
    jobs: tuple[Mapping[str, Any], ...]
    skipped_parents: tuple[Mapping[str, Any], ...]

    @property
    def array_size(self) -> int:
        return len(self.jobs)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": EXCHANGE_STAGE_MANIFEST_SCHEMA,
            "stage_root": str(self.stage_root),
            "manifest_path": str(self.manifest_path),
            "stage_id": self.stage_id,
            "array_size": self.array_size,
            "jobs": [dict(item) for item in self.jobs],
            "skipped_parents": [dict(item) for item in self.skipped_parents],
        }


@dataclass(frozen=True)
class ExchangeSlurmRenderResult:
    stage_root: Path
    run_list_path: Path
    script_path: Path
    array_size: int
    schema_version: str = "rsutil-exchange-slurm-v1"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "stage_root": str(self.stage_root),
            "run_list_path": str(self.run_list_path),
            "script_path": str(self.script_path),
            "array_size": self.array_size,
        }


@dataclass(frozen=True)
class ExchangeStageStatusEntry:
    array_index: int
    parent_entry_id: str
    parent_name: str
    parent_run_path: str
    exchange_manifest_path: str
    atom_path: str
    job_identity: str
    status: str
    evidence: str
    exchange_spec_hash: str = ""
    execution_status: str = EXECUTION_NOT_RUN
    output_status: str = OUTPUT_NO_OUTPUT
    parser_status: str = EXCHANGE_PARSER_INVALID
    parser_error: Mapping[str, Any] | None = None
    execution_evidence: Mapping[str, Any] | None = None
    result_files: tuple[str, ...] = ()
    pair_keys: tuple[Mapping[str, Any], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "array_index": self.array_index,
            "parent_entry_id": self.parent_entry_id,
            "parent_name": self.parent_name,
            "parent_run_path": self.parent_run_path,
            "exchange_manifest_path": self.exchange_manifest_path,
            "exchange_job_manifest": self.exchange_manifest_path,
            "atom_path": self.atom_path,
            "job_identity": self.job_identity,
            "canonical_exchange_job_identity": self.exchange_spec_hash or self.job_identity,
            "exchange_spec_hash": self.exchange_spec_hash,
            "status": self.status,
            "evidence": self.evidence,
            "execution_status": self.execution_status,
            "output_status": self.output_status,
            "exchange_output_status": self.output_status,
            "output_present": self.output_status == EXCHANGE_OUTPUT_PRESENT,
            "parser_status": self.parser_status,
            "parser_valid": self.parser_status == EXCHANGE_PARSER_VALID,
            "parser_error": None if self.parser_error is None else dict(self.parser_error),
            "execution_evidence": (
                None if self.execution_evidence is None else dict(self.execution_evidence)
            ),
            "result_files": list(self.result_files),
            "pair_keys": [dict(item) for item in self.pair_keys],
        }


@dataclass(frozen=True)
class ExchangeStageStatusReport:
    stage_id: str
    campaign_id: str
    campaign_name: str
    entries: tuple[ExchangeStageStatusEntry, ...]
    scheduler_state: str = "not_queried"
    schema_version: str = "rsutil-exchange-stage-status-v1"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "stage_id": self.stage_id,
            "campaign_id": self.campaign_id,
            "campaign_name": self.campaign_name,
            "scheduler_state": self.scheduler_state,
            "entries": [entry.to_dict() for entry in self.entries],
        }


@dataclass(frozen=True)
class ExchangeStageCollection:
    aggregate_path: Path
    entries: tuple[Mapping[str, Any], ...]
    schema_version: str = EXCHANGE_STAGE_RESULTS_SCHEMA

    @property
    def successful(self) -> tuple[Mapping[str, Any], ...]:
        return tuple(item for item in self.entries if item.get("status") == COMPLETED_VALID)

    @property
    def failures(self) -> tuple[Mapping[str, Any], ...]:
        return tuple(item for item in self.entries if item.get("status") != COMPLETED_VALID)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "aggregate_path": str(self.aggregate_path),
            "entries": [dict(item) for item in self.entries],
        }


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _version() -> str | None:
    try:
        import rsutil

        return getattr(rsutil, "__version__", None)
    except Exception:
        return None


def _stage_root(value: str | Path, parent_root: Path) -> Path:
    candidate = Path(value).expanduser()
    root = candidate.resolve() if candidate.is_absolute() else (parent_root / candidate).resolve()
    if root == parent_root:
        raise ExchangeStageError("exchange stage destination must be distinct from the parent campaign")
    if root.exists() and not root.is_dir():
        raise ExchangeStageError(f"exchange stage destination is not a directory: {root}")
    return root


def _stage_directory(destination: str | Path, parent_root: Path, stage_id: str) -> tuple[Path, Path]:
    """Resolve a persistent directory for one stage without orphaning others.

    New stages live below ``destination/<stage_id>``.  A legacy stage whose
    manifest is directly in ``destination`` remains readable, and the first
    new nested stage gets a compatibility manifest at the destination root.
    The compatibility file is never replaced once it exists.
    """

    base = _stage_root(destination, parent_root)
    nested = base / stage_id
    nested_manifest = nested / EXCHANGE_STAGE_MANIFEST_FILENAME

    if nested_manifest.is_file():
        try:
            existing = json.loads(nested_manifest.read_text(encoding="utf-8"))
        except Exception as exc:
            raise ExchangeStageError(
                f"could not read existing exchange stage manifest {nested_manifest}: {exc}"
            ) from exc
        if not isinstance(existing, Mapping) or existing.get("stage_id") != stage_id:
            raise ExchangeStageError(
                f"exchange stage directory {nested} belongs to a different stage_id"
            )
        return nested, base
    if nested.exists():
        if not nested.is_dir():
            raise ExchangeStageError(f"exchange stage path is not a directory: {nested}")
        if any(nested.iterdir()):
            raise ExchangeStageError(
                f"exchange stage directory has no valid manifest; refusing to use it: {nested}"
            )
        return nested, base

    root_manifest = base / EXCHANGE_STAGE_MANIFEST_FILENAME
    if root_manifest.is_file():
        try:
            existing = json.loads(root_manifest.read_text(encoding="utf-8"))
        except Exception:
            existing = None
        if isinstance(existing, Mapping) and existing.get("stage_id") == stage_id:
            return base, base
    return nested, base


def _write_stage_manifest_alias(base: Path, stage_root: Path) -> None:
    """Keep the historical destination-root entry point without overwriting it."""

    if stage_root == base:
        return
    alias = base / EXCHANGE_STAGE_MANIFEST_FILENAME
    if alias.exists() or alias.is_symlink():
        return
    alias.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(os.path.relpath(stage_root / EXCHANGE_STAGE_MANIFEST_FILENAME, alias.parent), alias)


def _write_stage_job_alias(base: Path, stage_root: Path, stage_id: str) -> None:
    """Expose the historical job path while retaining one actual job tree."""

    if stage_root == base:
        return
    actual = stage_root / "jobs" / stage_id
    if not actual.is_dir():
        return
    alias = base / "jobs" / stage_id
    if alias.exists() or alias.is_symlink():
        if alias.is_symlink() and alias.resolve() == actual.resolve():
            return
        raise ExchangeStageError(
            f"compatibility exchange job path already exists and is not this stage: {alias}"
        )
    alias.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(os.path.relpath(actual, alias.parent), alias, target_is_directory=True)


def _stage_path(root: Path, value: str, label: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ExchangeStageError(f"exchange stage contains unsafe {label}: {value!r}")
    result = (root / candidate).resolve()
    if result != root and root not in result.parents:
        raise ExchangeStageError(f"exchange stage {label} escapes its root: {value!r}")
    return result


def _relative(root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _load_stage(value: str | Path) -> tuple[dict[str, Any], Path, Path]:
    supplied = Path(value).expanduser().resolve()
    manifest_path = supplied / EXCHANGE_STAGE_MANIFEST_FILENAME if supplied.is_dir() else supplied
    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ExchangeStageError(f"could not read exchange stage manifest {manifest_path}: {exc}") from exc
    if not isinstance(data, dict) or data.get("schema_version") != EXCHANGE_STAGE_MANIFEST_SCHEMA:
        raise ExchangeStageError(f"unsupported exchange stage manifest: {manifest_path}")
    if not isinstance(data.get("stage_id"), str) or not isinstance(data.get("jobs"), list):
        raise ExchangeStageError(f"exchange stage manifest is incomplete: {manifest_path}")
    return data, manifest_path.parent, manifest_path


def _candidate_file(run_dir: Path, names: Sequence[str]) -> Path | None:
    for name in names:
        candidate = run_dir / name
        if candidate.is_file():
            return candidate
    return None


def _reference_paths(run_dir: Path) -> dict[str, Path | None]:
    return {
        "lattice": _candidate_file(run_dir, ("lattice.nml", "lattice_new.nml")),
        "input": _candidate_file(run_dir, ("input.nml", "input_new.nml")),
        "cluster": _candidate_file(run_dir, ("clust",)),
        "report": _candidate_file(run_dir, ("report.out",)),
    }


def _parent_diagnostics(status: Any) -> list[dict[str, Any]]:
    diagnostics: list[dict[str, Any]] = []
    if status.execution_status != EXECUTION_EXITED_ZERO:
        diagnostics.append({
            "code": "parent_execution_not_successful",
            "message": (
                "parent execution is not a recorded successful zero exit "
                f"({status.execution_status})"
            ),
            "execution_status": status.execution_status,
        })
    if status.output_status != OUTPUT_PARSE_VALID:
        diagnostics.append({
            "code": "parent_parser_invalid",
            "message": f"parent result parser status is {status.output_status}",
            "output_status": status.output_status,
        })
    if status.scientific_status != "CONVERGED":
        diagnostics.append({
            "code": "parent_scf_not_converged",
            "message": f"parent SCF status is {status.scientific_status}",
            "scientific_status": status.scientific_status,
        })
    if status.parser_error is not None:
        diagnostics.append({
            "code": "parent_status_diagnostic",
            "message": str(status.parser_error.get("message", "parent status has a diagnostic")),
            "details": dict(status.parser_error),
        })
    return diagnostics


def _canonical_exchange_options(
    *,
    cutoff: float,
    unique: bool,
    full_pairs: bool | None,
    shell_width: float,
    soc: bool | None,
    magnetic_axis: Sequence[float] | None,
    magnetic_reference: str | None,
    magnetic_context: MagneticSymmetryContext | None,
    max_pairs_per_job: int = 0,
) -> dict[str, Any]:
    if isinstance(max_pairs_per_job, bool) or int(max_pairs_per_job) != max_pairs_per_job or max_pairs_per_job < 0:
        raise ValueError("max_pairs_per_job must be zero or a positive integer")
    if not math.isfinite(float(cutoff)) or float(cutoff) <= 0:
        raise ValueError("exchange cutoff must be a finite positive number")
    if not math.isfinite(float(shell_width)) or float(shell_width) < 0:
        raise ValueError("exchange shell_width must be a finite non-negative number")
    effective_full_pairs = (not unique) if full_pairs is None else bool(full_pairs)
    effective_unique = bool(unique and not effective_full_pairs)
    axis = None if magnetic_axis is None else [float(value) for value in magnetic_axis]
    if axis is not None and len(axis) != 3:
        raise ValueError("magnetic_axis must contain exactly three values")
    if axis is not None and not all(math.isfinite(value) for value in axis):
        raise ValueError("magnetic_axis must contain finite values")
    reference = (
        MAGNETIC_REFERENCE_COLLINEAR_FM
        if magnetic_reference is None
        else str(magnetic_reference).strip().lower()
    )
    options: dict[str, Any] = {
        "cutoff": float(cutoff),
        "shell_width": float(shell_width),
        "full_pairs": effective_full_pairs,
        "unique": effective_unique,
        "pair_policy": "unique" if effective_unique else "full",
        "soc": None if soc is None else bool(soc),
        "magnetic_axis": axis,
        "magnetic_reference": reference,
        "magnetic_context": (
            None
            if magnetic_context is None
            else {
                "reference_class": magnetic_context.reference_class,
                "soc": magnetic_context.soc,
                "magnetic_axis": (
                    None
                    if magnetic_context.magnetic_axis is None
                    else list(magnetic_context.magnetic_axis)
                ),
                "nsp": magnetic_context.nsp,
                "nsp_source": magnetic_context.nsp_source,
            }
        ),
    }
    if max_pairs_per_job:
        # Recorded only when slicing so unsliced stage identities are unchanged.
        options["max_pairs_per_job"] = int(max_pairs_per_job)
    return options


def _check_reduction_policy(options: Mapping[str, Any]) -> None:
    if not options["unique"]:
        return
    reference = str(options["magnetic_reference"])
    context = options.get("magnetic_context")
    context_reference = (
        None if not isinstance(context, Mapping) else context.get("reference_class")
    )
    if context_reference not in (None, MAGNETIC_REFERENCE_COLLINEAR_FM):
        raise UnsupportedMagneticContextError(
            "unique exchange reduction requires an explicit proven "
            "magnetic context with reference_class='collinear_fm'"
        )
    if context_reference is None and reference != MAGNETIC_REFERENCE_COLLINEAR_FM:
        raise UnsupportedMagneticContextError(
            "unique exchange reduction requires an explicit proven "
            "magnetic_reference='collinear_fm'"
        )


def _stage_identity(parent_records: Sequence[Mapping[str, Any]], options: Mapping[str, Any]) -> str:
    return hash_json({
        "schema_version": EXCHANGE_STAGE_MANIFEST_SCHEMA,
        "parents": [
            {
                "entry_id": item["entry_id"],
                "exchange_source_identity": item["exchange_source_identity"],
            }
            for item in parent_records
            if item.get("exchange_source_identity")
        ],
        "exchange_options": dict(options),
    })


def _job_record(
    *,
    stage_root: Path,
    parent: Mapping[str, Any],
    job: ExchangeJobSpec,
    manifest_path: Path,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    relative_manifest = _relative(stage_root, manifest_path)
    atom_path = manifest_path.parent
    pair_keys = list(manifest["spec"]["pair_keys"])
    exchange_identity = str(manifest["spec_hash"])
    # The exchange manifest identity is scientific.  Scope the array identity
    # by parent and Atom index as well, so two identical parent runs remain
    # distinct executable rows while still reusing the canonical manifest.
    job_identity = hash_json({
        "parent_entry_id": parent["entry_id"],
        "atom_index": int(job.atom_index),
        "exchange_spec_hash": exchange_identity,
    })
    return {
        "parent_entry_id": parent["entry_id"],
        "parent_name": parent["name"],
        "parent_run_path": parent["run_path"],
        "parent_manifest_path": parent["manifest_path"],
        "parent_manifest_hash": parent["manifest_hash"],
        "atom_index": int(job.atom_index),
        **(
            {"chunk_index": int(job.chunk_index), "chunk_count": int(job.chunk_count), "chunk_label": job.chunk_label}
            if job.chunk_label
            else {}
        ),
        "exchange_manifest_path": relative_manifest,
        "exchange_job_manifest": relative_manifest,
        "exchange_manifest_hash": hash_file(manifest_path),
        "atom_path": _relative(stage_root, atom_path),
        "atom_directory": _relative(stage_root, atom_path),
        "canonical_exchange_job_identity": exchange_identity,
        "exchange_job_identity": job_identity,
        "job_identity": job_identity,
        "exchange_spec_hash": exchange_identity,
        "canonical_pair_list": pair_keys,
        "pair_key_inventory": pair_keys,
    }


def prepare_exchange_stage(
    campaign_index: str | Path | Any,
    destination: str | Path = "exchange-stage",
    *,
    campaign_root: str | Path | None = None,
    cutoff: float = 1.0,
    unique: bool = False,
    full_pairs: bool | None = None,
    shell_width: float = 0.1,
    soc: bool | None = None,
    magnetic_axis: Sequence[float] | None = None,
    magnetic_reference: str | None = None,
    magnetic_context: MagneticSymmetryContext | None = None,
    max_pairs_per_job: int = 0,
) -> ExchangeStagePreparation:
    """Prepare exchange jobs from only successful, parsed, converged parents."""

    if not isinstance(campaign_index, (str, Path)) and campaign_root is None:
        destination_path = Path(destination).expanduser()
        campaign_root = (
            destination_path.resolve().parent
            if destination_path.is_absolute()
            else Path.cwd()
        )
    index, parent_root, index_path = _load_index(
        campaign_index, campaign_root=campaign_root
    )
    options = _canonical_exchange_options(
        cutoff=cutoff,
        unique=unique,
        full_pairs=full_pairs,
        shell_width=shell_width,
        soc=soc,
        magnetic_axis=magnetic_axis,
        magnetic_reference=magnetic_reference,
        magnetic_context=magnetic_context,
        max_pairs_per_job=max_pairs_per_job,
    )
    _check_reduction_policy(options)

    parent_status = status_campaign(index, campaign_root=parent_root)
    status_by_id = {item.entry_id: item for item in parent_status.entries}
    parent_records: list[dict[str, Any]] = []
    for entry in index.entries:
        status = status_by_id[entry.entry_id]
        diagnostics = _parent_diagnostics(status)
        parent_records.append({
            "entry_id": entry.entry_id,
            "name": entry.name,
            "run_path": entry.run_path,
            "manifest_path": entry.manifest_path,
            "manifest_hash": entry.manifest_hash,
            "status": status.status,
            "execution_status": status.execution_status,
            "output_status": status.output_status,
            "scientific_status": status.scientific_status,
            "parent_status_eligible": not diagnostics,
            "eligible": False,
            "diagnostics": diagnostics,
        })

    # Plan first so the stage identity can use the exchange service's
    # authoritative source-run identity.  This identity includes the current
    # parent report/Fermi source and required converged atom files.
    planned: dict[int, tuple[dict[str, Any], dict[str, Path | None], Any]] = {}
    skipped: list[dict[str, Any]] = []
    for parent_position, (entry, parent) in enumerate(zip(index.entries, parent_records)):
        if not parent["parent_status_eligible"]:
            skipped.append(dict(parent))
            continue
        run_dir = (parent_root / entry.run_path).resolve()
        paths = _reference_paths(run_dir)
        missing = [name for name, path in paths.items() if path is None]
        if missing:
            parent["diagnostics"].extend({
                "code": "exchange_source_missing",
                "message": f"eligible parent is missing exchange source file: {name}",
                "path": str(run_dir / name),
            } for name in missing)
            skipped.append(dict(parent))
            continue
        try:
            plan = ExchangeService(run_dir).plan(
                cutoff=options["cutoff"],
                unique=options["unique"],
                full_pairs=options["full_pairs"],
                shell_width=options["shell_width"],
                soc=options["soc"],
                magnetic_axis=(
                    None
                    if options["magnetic_axis"] is None
                    else tuple(options["magnetic_axis"])
                ),
                # Leave an unspecified reference to the service so each
                # parent resolves its own default (noncollinear -> unknown).
                magnetic_reference=(
                    None
                    if magnetic_context is not None or magnetic_reference is None
                    else options["magnetic_reference"]
                ),
                magnetic_context=magnetic_context,
                max_pairs_per_job=options.get("max_pairs_per_job", 0),
                input_path=paths["input"],
                lattice_path=paths["lattice"],
                cluster_path=paths["cluster"],
                report_path=paths["report"],
            )
            if options["unique"] and not all(
                job.symmetry_reduction_applied for job in plan.jobs
            ):
                reasons = sorted({
                    job.symmetry_reduction_reason
                    for job in plan.jobs
                    if not job.symmetry_reduction_applied
                })
                reason = "; ".join(reasons) or "the exchange service did not prove reduction support"
                raise UnsupportedMagneticContextError(
                    f"unique exchange reduction is unsupported for parent {entry.entry_id}: {reason}"
                )
            if not plan.jobs:
                parent["diagnostics"].append({
                    "code": "exchange_jobs_empty",
                    "message": "exchange service generated no Atom_* jobs",
                })
                skipped.append(dict(parent))
                continue
            source_identities = {
                str(job.source_run_identity)
                for job in plan.jobs
                if str(job.source_run_identity)
            }
            if len(source_identities) != 1:
                raise ExchangeStageError(
                    "exchange service did not provide one authoritative parent source-run identity"
                )
            parent["exchange_source_identity"] = next(iter(source_identities))
            parent["eligible"] = True
            parent["reference_inputs"] = {
                name: path.name for name, path in paths.items() if path is not None
            }
            planned[parent_position] = (parent, paths, plan)
        except UnsupportedMagneticContextError:
            raise
        except Exception as exc:
            parent["diagnostics"].append({
                "code": "exchange_preparation_failed",
                "message": str(exc),
                "type": exc.__class__.__name__,
            })
            skipped.append(dict(parent))

    stage_id = _stage_identity(parent_records, options)
    stage_root, stage_base = _stage_directory(destination, parent_root, stage_id)
    jobs_by_parent = stage_root / "jobs" / stage_id
    generated: list[tuple[int, tuple[int, int], dict[str, Any]]] = []
    for parent_position, (parent, paths, plan) in planned.items():
        entry = index.entries[parent_position]
        try:
            parent_job_root = jobs_by_parent / entry.entry_id
            run_dir = (parent_root / entry.run_path).resolve()
            rebased_jobs = tuple(
                replace(job, output_dir=parent_job_root / job.directory_name)
                for job in plan.jobs
            )
            write_exchange_jobs(
                rebased_jobs,
                lattice_path=paths["lattice"],  # type: ignore[arg-type]
                source_dir=run_dir,
                destination_root=parent_job_root,
                collision_policy="reuse_if_identical",
            )
            for job in rebased_jobs:
                atom_dir = job.output_dir.resolve()
                manifest_path = atom_dir / EXCHANGE_MANIFEST_FILENAME
                manifest = read_exchange_manifest(manifest_path)
                generated.append((
                    parent_position,
                    (int(job.atom_index), int(job.chunk_index)),
                    _job_record(
                        stage_root=stage_root,
                        parent=parent,
                        job=job,
                        manifest_path=manifest_path,
                        manifest=manifest,
                    ),
                ))
            parent["eligible"] = True
        except UnsupportedMagneticContextError:
            raise
        except Exception as exc:
            generated[:] = [
                item for item in generated if item[0] != parent_position
            ]
            parent["diagnostics"].append({
                "code": "exchange_preparation_failed",
                "message": str(exc),
                "type": exc.__class__.__name__,
            })
            skipped.append(dict(parent))

    generated.sort(key=lambda item: (item[0], item[1], item[2]["atom_path"]))
    job_records: list[dict[str, Any]] = []
    for array_index, (_parent_position, _atom_index, record) in enumerate(generated):
        record = dict(record)
        record["array_index"] = array_index
        job_records.append(record)

    selected_parents = [item for item in parent_records if item.get("eligible")]
    generation = {
        "version": EXCHANGE_STAGE_MANIFEST_SCHEMA,
        "timestamp": _now(),
        "rsutil_version": _version(),
    }
    manifest: dict[str, Any] = {
        "schema_version": EXCHANGE_STAGE_MANIFEST_SCHEMA,
        "owner": "rsutil",
        "operation": "campaign_exchange_stage",
        "stage_id": stage_id,
        "campaign_id": index.campaign_id,
        "campaign_name": index.campaign_name,
        "parent_campaign_id": index.campaign_id,
        "parent_campaign_index": _relative(stage_root, index_path),
        # Keep the stage relocatable when it is moved together with its
        # parent campaign.  This is provenance, not a stage-owned path, so it
        # may legitimately contain ``..`` components.
        "parent_campaign_root": _relative(stage_root, parent_root),
        "parent_entry_ids": [item["entry_id"] for item in selected_parents],
        "parent_run_paths": [item["run_path"] for item in selected_parents],
        "exchange_options": dict(options),
        "parent_entries": parent_records,
        "skipped_parents": skipped,
        "jobs": job_records,
        "generated_exchange_job_manifests": [item["exchange_manifest_path"] for item in job_records],
        "atom_paths": [item["atom_path"] for item in job_records],
        "pair_key_inventory": [
            {
                "parent_entry_id": item["parent_entry_id"],
                "job_identity": item["job_identity"],
                "exchange_spec_hash": item["exchange_spec_hash"],
                "atom_path": item["atom_path"],
                "pair_keys": item["canonical_pair_list"],
            }
            for item in job_records
        ],
        "generation": generation,
        "generated_at": generation["timestamp"],
        "generation_version": generation["version"],
    }
    _write_stage_job_alias(stage_base, stage_root, stage_id)
    manifest_path = stage_root / EXCHANGE_STAGE_MANIFEST_FILENAME
    _atomic_json(manifest_path, manifest)
    _write_stage_manifest_alias(stage_base, stage_root)
    return ExchangeStagePreparation(
        stage_root=stage_root,
        manifest_path=manifest_path,
        stage_id=stage_id,
        jobs=tuple(job_records),
        skipped_parents=tuple(skipped),
    )


prepare_campaign_exchange = prepare_exchange_stage
exchange_stage_prepare = prepare_exchange_stage
prepare_exchange = prepare_exchange_stage


def _exchange_result_names() -> tuple[str, ...]:
    return (
        "jij.out", "jijso.out", "jijfo.out",
        "dij.out", "dijso.out", "dijfo.out", "dmi.out", "dmi.dat",
        "aij.out", "aijso.out", "aijfo.out", "jtens.out", "exchange_tensor.out",
    )


def _result_paths(atom_dir: Path) -> tuple[Path, ...]:
    return tuple(atom_dir / name for name in _exchange_result_names() if (atom_dir / name).is_file())


def _read_exchange_execution(
    atom_dir: Path,
    *,
    stage: Mapping[str, Any],
    job: Mapping[str, Any],
) -> tuple[str, Mapping[str, Any] | None, dict[str, Any] | None]:
    path = atom_dir / "rsutil_execution.json"
    if not path.is_file():
        return EXECUTION_NOT_RUN, None, None
    try:
        evidence = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(evidence, Mapping):
            raise ValueError("execution evidence is not an object")
        if evidence.get("schema_version") != EXECUTION_SCHEMA_VERSION:
            raise ValueError(f"unsupported execution evidence schema: {evidence.get('schema_version')!r}")
        identity = evidence.get("identity")
        if not isinstance(identity, Mapping):
            raise ValueError("execution evidence identity is missing")
        expected = {
            "exchange_stage_id": stage["stage_id"],
            "parent_entry_id": job["parent_entry_id"],
            "exchange_manifest_hash": job["exchange_manifest_hash"],
            "job_identity": job["job_identity"],
            "array_index": str(job["array_index"]),
        }
        if any(identity.get(key) != value for key, value in expected.items()):
            return EXECUTION_RUNNING_OR_UNKNOWN, evidence, {
                "code": "execution_identity_mismatch",
                "message": "exchange execution evidence does not match the stage job",
            }
        status = evidence.get("status")
        if status == "running":
            return EXECUTION_RUNNING_OR_UNKNOWN, evidence, None
        return_code = evidence.get("return_code")
        if not isinstance(return_code, int):
            raise ValueError("final execution evidence has no integer return_code")
        if status == "success" and return_code == 0:
            return EXECUTION_EXITED_ZERO, evidence, None
        if status == "failure" and return_code != 0:
            if isinstance(evidence.get("signal"), int) and evidence["signal"] > 0:
                return EXECUTION_SIGNALED, evidence, None
            return EXECUTION_EXITED_NONZERO, evidence, None
        raise ValueError("execution evidence status and return_code disagree")
    except Exception as exc:
        return EXECUTION_RUNNING_OR_UNKNOWN, None, _error_payload(exc)


def _exchange_summary(execution_status: str, output_status: str, parser_status: str) -> str:
    if execution_status in (EXECUTION_EXITED_NONZERO, EXECUTION_SIGNALED):
        return FAILED
    if parser_status in (EXCHANGE_PARSER_INVALID, EXCHANGE_PARSER_UNSUPPORTED):
        return COMPLETED_INVALID
    if output_status == OUTPUT_NO_OUTPUT:
        return COMPLETED_INVALID if execution_status == EXECUTION_EXITED_ZERO else GENERATED
    if execution_status != EXECUTION_EXITED_ZERO:
        return UNKNOWN
    return COMPLETED_VALID


def _stage_status_entry(stage_root: Path, stage: Mapping[str, Any], job: Mapping[str, Any]) -> ExchangeStageStatusEntry:
    atom_dir = _stage_path(stage_root, str(job["atom_path"]), "atom path")
    manifest_path = _stage_path(stage_root, str(job["exchange_manifest_path"]), "exchange manifest path")
    result_paths = _result_paths(atom_dir) if atom_dir.is_dir() else ()
    execution_status, evidence, execution_error = _read_exchange_execution(
        atom_dir, stage=stage, job=job
    )
    parser_status = EXCHANGE_PARSER_INVALID
    parser_error: Mapping[str, Any] | None = execution_error
    pair_keys: tuple[Mapping[str, Any], ...] = ()
    manifest_error: dict[str, Any] | None = None
    if not manifest_path.is_file():
        manifest_error = {"code": "exchange_manifest_missing", "message": f"manifest is missing: {manifest_path}"}
    else:
        try:
            if hash_file(manifest_path) != job["exchange_manifest_hash"]:
                raise ValueError("exchange manifest content hash does not match the stage manifest")
            manifest = read_exchange_manifest(manifest_path)
            pair_keys = tuple(manifest["spec"]["pair_keys"])
        except Exception as exc:
            manifest_error = _error_payload(exc)
    output_status = EXCHANGE_OUTPUT_PRESENT if result_paths else OUTPUT_NO_OUTPUT
    if result_paths and manifest_error is None:
        try:
            parsed = parse_exchange_results(atom_dir, manifest=manifest_path)
            parser_status = EXCHANGE_PARSER_VALID if parsed.validation.valid else EXCHANGE_PARSER_INVALID
            if parser_status != EXCHANGE_PARSER_VALID:
                parser_error = {
                    "code": "exchange_result_validation_failed",
                    "message": "canonical exchange result validation is invalid",
                    "validation": parsed.validation.to_dict(),
                }
        except Exception as exc:
            parser_status = EXCHANGE_PARSER_UNSUPPORTED if isinstance(exc, UnsupportedResultFormatError) else EXCHANGE_PARSER_INVALID
            parser_error = _error_payload(exc)
    elif manifest_error is not None:
        parser_error = manifest_error
    status = _exchange_summary(execution_status, output_status, parser_status)
    if manifest_error is not None:
        status = COMPLETED_INVALID if result_paths or execution_status == EXECUTION_EXITED_ZERO else UNKNOWN
        evidence_text = "exchange manifest is missing, stale, or invalid"
    elif execution_error is not None:
        evidence_text = "execution evidence is missing, unfinished, or has the wrong identity"
    elif execution_status in (EXECUTION_EXITED_NONZERO, EXECUTION_SIGNALED):
        evidence_text = "exchange execution evidence records failure"
    elif output_status == OUTPUT_NO_OUTPUT and execution_status == EXECUTION_NOT_RUN:
        evidence_text = "exchange job generated; no execution or result evidence"
    elif output_status == OUTPUT_NO_OUTPUT:
        evidence_text = "exchange execution completed without supported output"
    elif parser_status == EXCHANGE_PARSER_VALID and execution_status == EXECUTION_EXITED_ZERO:
        evidence_text = "execution exited zero and canonical exchange results parsed"
    elif parser_status == EXCHANGE_PARSER_VALID:
        evidence_text = "canonical exchange results parsed; execution evidence is incomplete"
    else:
        evidence_text = "exchange output is malformed or unsupported"
    return ExchangeStageStatusEntry(
        array_index=int(job["array_index"]),
        parent_entry_id=str(job["parent_entry_id"]),
        parent_name=str(job["parent_name"]),
        parent_run_path=str(job["parent_run_path"]),
        exchange_manifest_path=str(job["exchange_manifest_path"]),
        atom_path=str(job["atom_path"]),
        job_identity=str(job["job_identity"]),
        exchange_spec_hash=str(job["exchange_spec_hash"]),
        status=status,
        evidence=evidence_text,
        execution_status=execution_status,
        output_status=output_status,
        parser_status=parser_status,
        parser_error=parser_error,
        execution_evidence=evidence,
        result_files=tuple(_relative(stage_root, path) for path in result_paths),
        pair_keys=pair_keys,
    )


def status_exchange_stage(stage: str | Path) -> ExchangeStageStatusReport:
    """Inspect exchange execution, output presence, and parser validity separately."""

    data, root, _manifest_path = _load_stage(stage)
    campaign = data.get("campaign_id") or data.get("parent_campaign_id")
    if not isinstance(campaign, str):
        raise ExchangeStageError("exchange stage manifest has no campaign_id")
    entries = tuple(
        _stage_status_entry(root, data, job)
        for job in sorted(data["jobs"], key=lambda item: int(item.get("array_index", 0)))
    )
    return ExchangeStageStatusReport(
        stage_id=str(data["stage_id"]),
        campaign_id=campaign,
        campaign_name=str(data.get("campaign_name", "campaign")),
        entries=entries,
    )


exchange_stage_status = status_exchange_stage
status_exchange = status_exchange_stage


def _slurm_resources(resources: Mapping[str, Any] | None) -> dict[str, Any]:
    values = dict(resources or {})
    aliases = {"memory": "mem", "time_limit": "time"}
    result: dict[str, Any] = {}
    for key, value in values.items():
        if value is None:
            continue
        normalized = aliases.get(str(key).lower().replace("-", "_"), str(key).lower().replace("-", "_"))
        result[normalized] = value
    _validate_resource_values(result)
    return result


def render_exchange_slurm(
    stage: str | Path,
    *,
    executable: str | Path | Sequence[str] | None = None,
    launcher: str | Sequence[str] | None = None,
    arguments: Sequence[str] = (),
    output_dir: str | Path | None = None,
    run_list_name: str = EXCHANGE_STAGE_RUN_LIST_FILENAME,
    script_name: str = EXCHANGE_STAGE_SLURM_SCRIPT_FILENAME,
    job_name: str | None = None,
    partition: str | None = None,
    account: str | None = None,
    time_limit: str | None = None,
    resources: Mapping[str, Any] | None = None,
    environment_preamble: str | Sequence[str] | None = None,
) -> ExchangeSlurmRenderResult:
    """Render a stable Slurm array for the prepared Atom_* directories."""

    data, root, _manifest_path = _load_stage(stage)
    jobs = sorted(data["jobs"], key=lambda item: int(item.get("array_index", 0)))
    if not jobs:
        raise ExchangeStageError("exchange stage contains no generated jobs")
    if Path(run_list_name).name != run_list_name or Path(script_name).name != script_name:
        raise ValueError("run-list and script names must be simple file names")
    executable_tokens = _executable_tokens(executable)
    launcher_tokens = _command_tokens(launcher, required=False, name="launcher")
    argument_tokens = [str(item) for item in arguments]
    if any("\x00" in token or "\n" in token or "\r" in token for token in argument_tokens):
        raise ValueError("command arguments cannot contain NUL or newline characters")
    preamble = _preamble_lines(environment_preamble)
    destination = root if output_dir is None else Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    run_list_path = destination / run_list_name
    script_path = destination / script_name

    headers = (
        "array_index\tparent_entry_id\tparent_name\tparent_run_path\t"
        "exchange_manifest_path\texchange_manifest_hash\tatom_path\tjob_identity\texchange_spec_hash"
    )
    rows = [headers]
    for position, job in enumerate(jobs):
        if int(job.get("array_index", -1)) != position:
            raise ExchangeStageError("exchange stage array indices are not contiguous")
        fields = [
            str(position),
            str(job["parent_entry_id"]),
            str(job["parent_name"]),
            str(job["parent_run_path"]),
            str(job["exchange_manifest_path"]),
            str(job["exchange_manifest_hash"]),
            str(job["atom_path"]),
            str(job["job_identity"]),
            str(job["exchange_spec_hash"]),
        ]
        if any(any(character in field for character in "\t\r\n") for field in fields):
            raise ExchangeStageError("exchange stage run-list fields cannot contain tabs or newlines")
        _stage_path(root, fields[6], "atom path")
        _stage_path(root, fields[4], "exchange manifest path")
        rows.append("\t".join(fields))
    _atomic_write(run_list_path, "\n".join(rows) + "\n")

    job = _job_name(job_name or f"{data.get('campaign_name', 'campaign')}-exchange")
    directives = [f"#SBATCH --job-name={job}", f"#SBATCH --array=0-{len(jobs) - 1}"]
    for option, value in (("partition", partition), ("account", account), ("time", time_limit)):
        if value is not None:
            directives.append(f"#SBATCH --{option}={_safe_directive(value, option)}")
    known = {
        "cpus_per_task": "cpus-per-task",
        "mem": "mem",
        "gpus": "gpus",
        "nodes": "nodes",
        "ntasks": "ntasks",
        "ntasks_per_node": "ntasks-per-node",
        "constraint": "constraint",
        "qos": "qos",
    }
    resource_values = _slurm_resources(resources)
    for key, option in known.items():
        if key in resource_values:
            directives.append(f"#SBATCH --{option}={_safe_directive(resource_values[key], option)}")
    for key, value in resource_values.items():
        if key not in known and key != "time":
            directives.append(f"#SBATCH --{key.replace('_', '-')}={_safe_directive(value, key)}")

    command = launcher_tokens + executable_tokens + argument_tokens
    root_reference = shlex.quote(root.resolve().as_posix())
    lines = ["#!/usr/bin/env bash", *directives, "set -euo pipefail", ""]
    if preamble:
        lines.extend(["# User-supplied environment preamble", *preamble, ""])
    lines.extend([
        f"STAGE_ROOT={root_reference}",
        f"RUN_LIST={shlex.quote(run_list_path.resolve().as_posix())}",
        'TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"',
        'if ! [[ "$TASK_ID" =~ ^[0-9]+$ ]]; then echo "invalid array task ID: $TASK_ID" >&2; exit 2; fi',
        'LINE_NUMBER=$((TASK_ID + 2))',
        'IFS=$\'\t\' read -r ARRAY_INDEX PARENT_ENTRY_ID PARENT_NAME PARENT_RUN_PATH MANIFEST_PATH MANIFEST_HASH ATOM_PATH JOB_IDENTITY EXCHANGE_SPEC_HASH < <(sed -n "${LINE_NUMBER}p" "$RUN_LIST")',
        'if [[ -z "${ATOM_PATH:-}" || "$ARRAY_INDEX" != "$TASK_ID" ]]; then echo "array task is not present in exchange run list: $TASK_ID" >&2; exit 2; fi',
        'RUN_DIR="$STAGE_ROOT/$ATOM_PATH"',
        'if [[ ! -d "$RUN_DIR" ]]; then echo "Atom directory does not exist: $RUN_DIR" >&2; exit 3; fi',
        'export RSUTIL_EXCHANGE_STAGE_ID=' + shlex.quote(str(data["stage_id"])),
        'export RSUTIL_EXCHANGE_CAMPAIGN_ID=' + shlex.quote(str(data.get("campaign_id", ""))),
        'export RSUTIL_EXCHANGE_ARRAY_INDEX="$ARRAY_INDEX"',
        'export RSUTIL_EXCHANGE_PARENT_ENTRY_ID="$PARENT_ENTRY_ID"',
        'export RSUTIL_EXCHANGE_PARENT_RUN_PATH="$PARENT_RUN_PATH"',
        'export RSUTIL_EXCHANGE_MANIFEST_PATH="$STAGE_ROOT/$MANIFEST_PATH"',
        'export RSUTIL_EXCHANGE_MANIFEST_HASH="$MANIFEST_HASH"',
        'export RSUTIL_EXCHANGE_JOB_IDENTITY="$JOB_IDENTITY"',
        'export RSUTIL_EXECUTION_PATH="$RUN_DIR/rsutil_execution.json"',
        'cd -- "$RUN_DIR"',
        f"echo \"rsutil exchange job $JOB_IDENTITY (array $ARRAY_INDEX) in $RUN_DIR\"",
        f"COMMAND_JSON={shlex.quote(json.dumps(command, separators=(',', ':')))}",
        'python - "$COMMAND_JSON" <<\'PY\'',
        "import json",
        "import os",
        "import socket",
        "import subprocess",
        "import sys",
        "import tempfile",
        "from datetime import datetime, timezone",
        "from pathlib import Path",
        "",
        'execution_path = Path(os.environ["RSUTIL_EXECUTION_PATH"])',
        'command = json.loads(sys.argv[1])',
        'started_at = datetime.now(timezone.utc).isoformat()',
        'identity = {',
        '    "exchange_stage_id": os.environ.get("RSUTIL_EXCHANGE_STAGE_ID"),',
        '    "campaign_id": os.environ.get("RSUTIL_EXCHANGE_CAMPAIGN_ID"),',
        '    "parent_entry_id": os.environ.get("RSUTIL_EXCHANGE_PARENT_ENTRY_ID"),',
        '    "parent_run_path": os.environ.get("RSUTIL_EXCHANGE_PARENT_RUN_PATH"),',
        '    "exchange_manifest_hash": os.environ.get("RSUTIL_EXCHANGE_MANIFEST_HASH"),',
        '    "job_identity": os.environ.get("RSUTIL_EXCHANGE_JOB_IDENTITY"),',
        '    "array_index": os.environ.get("RSUTIL_EXCHANGE_ARRAY_INDEX"),',
        '}',
        'slurm = {key: value for key, value in {',
        '    "job_id": os.environ.get("SLURM_JOB_ID"),',
        '    "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),',
        '    "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),',
        '}.items() if value is not None}',
        'def atomic_record(payload):',
        '    execution_path.parent.mkdir(parents=True, exist_ok=True)',
        '    fd, temporary_name = tempfile.mkstemp(prefix="." + execution_path.name + ".", dir=str(execution_path.parent), text=True)',
        '    temporary = Path(temporary_name)',
        '    try:',
        '        with os.fdopen(fd, "w", encoding="utf-8") as stream:',
        '            json.dump(payload, stream, indent=2, sort_keys=True)',
        '            stream.write("\\n")',
        '        os.replace(temporary, execution_path)',
        '    finally:',
        '        if temporary.exists():',
        '            temporary.unlink()',
        'base = {',
        '    "schema_version": "' + EXECUTION_SCHEMA_VERSION + '",',
        '    "identity": identity,',
        '    "command": command,',
        '    "cwd": os.getcwd(),',
        '    "hostname": socket.gethostname(),',
        '    "started_at": started_at,',
        '}',
        'atomic_record({**base, "status": "running", "finished_at": None, "return_code": None, "signal": None})',
        'try:',
        '    completed = subprocess.run(command, cwd=os.getcwd(), check=False)',
        '    return_code = completed.returncode',
        '    error = None',
        'except OSError as exc:',
        '    return_code = 127',
        '    error = {"type": exc.__class__.__name__, "message": str(exc)}',
        'finished_at = datetime.now(timezone.utc).isoformat()',
        'signal = -return_code if return_code < 0 else None',
        'payload = {**base, "status": "success" if return_code == 0 else "failure", "finished_at": finished_at, "return_code": return_code, "signal": signal}',
        'if error is not None:',
        '    payload["error"] = error',
        'atomic_record(payload)',
        'raise SystemExit(return_code)',
        'PY',
        "",
    ])
    _atomic_write(script_path, "\n".join(lines))
    script_path.chmod(script_path.stat().st_mode | 0o111)
    return ExchangeSlurmRenderResult(root, run_list_path, script_path, len(jobs))


exchange_stage_render_slurm = render_exchange_slurm
render_slurm_exchange = render_exchange_slurm


def _reference_structure(stage_root: Path, stage: Mapping[str, Any], job: Mapping[str, Any]) -> Any:
    parent_root_value = stage.get("parent_campaign_root")
    roots: list[Path] = []
    if isinstance(parent_root_value, str):
        candidate = Path(parent_root_value).expanduser()
        roots.append(
            candidate.resolve()
            if candidate.is_absolute()
            else (stage_root / candidate).resolve()
        )
    inferred = stage_root.parent.resolve()
    if inferred not in roots:
        roots.append(inferred)
    for parent_root in roots:
        run_dir = (parent_root / str(job["parent_run_path"])).resolve()
        if not run_dir.is_dir():
            continue
        paths = _reference_paths(run_dir)
        if all(paths.get(name) is not None for name in ("lattice", "input", "cluster", "report")):
            try:
                reference = read_reference_calculation(
                    run_dir,
                    input_path=paths["input"],
                    lattice_path=paths["lattice"],
                    cluster_path=paths["cluster"],
                    report_path=paths["report"],
                )
                return reference.canonical_structure
            except Exception:
                continue
    return None


def _result_payload(parsed: Any) -> dict[str, Any]:
    # The parser attaches private lookup caches (including PairKey objects) to
    # the loaded manifest.  They are useful while mapping rows but are not a
    # JSON result contract, so remove those caches before serializing the
    # canonical ExchangeSet.
    run_manifest = getattr(parsed, "run_manifest", None)
    if isinstance(run_manifest, Mapping):
        safe_manifest = {
            key: value for key, value in run_manifest.items() if not str(key).startswith("_")
        }
        from dataclasses import replace as dataclass_replace

        parsed = dataclass_replace(parsed, run_manifest=safe_manifest)
    payload = parsed.to_dict()
    enriched_records: list[dict[str, Any]] = []
    for record in parsed.interactions:
        item = record.to_dict()
        full = record.full_exchange_tensor()
        if full is not None:
            item["derived_full_tensor"] = full.to_dict()
        enriched_records.append(item)
    # Keep the exact canonical result alongside the focused records.  The
    # latter only adds the derived tensor serialization supplied by the model.
    payload["records"] = enriched_records
    return payload


def collect_exchange_stage(
    stage: str | Path,
    *,
    output: str | Path | None = None,
    sidecar_dir: str | Path | None = None,
) -> ExchangeStageCollection:
    """Collect canonical exchange results without introducing another parser."""

    data, root, manifest_path = _load_stage(stage)
    aggregate_path = root / EXCHANGE_STAGE_RESULTS_FILENAME if output is None else Path(output).expanduser().resolve()
    sidecars = root / "collected" if sidecar_dir is None else Path(sidecar_dir).expanduser().resolve()
    sidecars.mkdir(parents=True, exist_ok=True)
    status_report = status_exchange_stage(root)
    status_by_id = {entry.job_identity: entry for entry in status_report.entries}
    entries: list[dict[str, Any]] = []
    structure_cache: dict[str, Any] = {}
    for job in sorted(data["jobs"], key=lambda item: int(item.get("array_index", 0))):
        identity = str(job["job_identity"])
        status = status_by_id[identity]
        atom_dir = _stage_path(root, str(job["atom_path"]), "atom path")
        result_files = _result_paths(atom_dir)
        collected: dict[str, Any] = {
            **dict(job),
            "status": status.status,
            "status_evidence": status.evidence,
            "execution_status": status.execution_status,
            "output_status": status.output_status,
            "exchange_output_status": status.output_status,
            "output_present": status.output_status == EXCHANGE_OUTPUT_PRESENT,
            "parser_status": status.parser_status,
            "parser_valid": status.parser_status == EXCHANGE_PARSER_VALID,
            "parser_error": None if status.parser_error is None else dict(status.parser_error),
            "execution_evidence": status.execution_evidence,
            "result_files": [_relative(root, path) for path in result_files],
        }
        if result_files and status.parser_status == EXCHANGE_PARSER_VALID:
            try:
                if job["parent_entry_id"] not in structure_cache:
                    structure_cache[job["parent_entry_id"]] = _reference_structure(root, data, job)
                parsed = parse_exchange_results(
                    atom_dir,
                    manifest=_stage_path(root, str(job["exchange_manifest_path"]), "exchange manifest path"),
                    structure=structure_cache[job["parent_entry_id"]],
                )
                parsed_payload = _result_payload(parsed)
                sidecar = sidecars / "exchange" / f"{identity}.json"
                _atomic_json(sidecar, parsed_payload)
                collected["result_sidecar"] = _relative(root, sidecar)
                collected["result"] = parsed_payload
                collected["records"] = parsed_payload["records"]
                collected["validation"] = parsed.validation.to_dict()
                collected["parser_provenance"] = parsed.provenance.to_dict()
            except Exception as exc:
                collected["parser_status"] = EXCHANGE_PARSER_INVALID
                collected["parser_error"] = _error_payload(exc)
                if collected["status"] == COMPLETED_VALID:
                    collected["status"] = COMPLETED_INVALID
        elif not result_files:
            collected["parse_error"] = {
                "code": "missing_exchange_result_file",
                "type": "MissingResultFileError",
                "message": "no supported exchange result file is present",
            }
        entries.append(collected)

    aggregate = {
        "schema_version": EXCHANGE_STAGE_RESULTS_SCHEMA,
        "stage_id": data["stage_id"],
        "campaign_id": data.get("campaign_id"),
        "campaign_name": data.get("campaign_name"),
        "exchange_options": data.get("exchange_options", {}),
        "collected_at": _now(),
        "source_stage_manifest": {
            "path": _relative(root, manifest_path),
            "sha256": hash_file(manifest_path),
        },
        "pair_key_inventory": data.get("pair_key_inventory", []),
        "skipped_parents": data.get("skipped_parents", []),
        "entries": entries,
    }
    _atomic_json(aggregate_path, aggregate)
    return ExchangeStageCollection(aggregate_path, tuple(entries))


exchange_stage_collect = collect_exchange_stage
collect_exchange = collect_exchange_stage


__all__ = [
    "EXCHANGE_OUTPUT_PRESENT",
    "EXCHANGE_PARSER_INVALID",
    "EXCHANGE_PARSER_UNSUPPORTED",
    "EXCHANGE_PARSER_VALID",
    "EXCHANGE_STAGE_MANIFEST_FILENAME",
    "EXCHANGE_STAGE_MANIFEST_SCHEMA",
    "EXCHANGE_STAGE_RESULTS_FILENAME",
    "EXCHANGE_STAGE_RESULTS_SCHEMA",
    "EXCHANGE_STAGE_RUN_LIST_FILENAME",
    "EXCHANGE_STAGE_SLURM_SCRIPT_FILENAME",
    "ExchangeStageCollection",
    "ExchangeStageError",
    "ExchangeStagePreparation",
    "ExchangeStageStatusEntry",
    "ExchangeStageStatusReport",
    "ExchangeSlurmRenderResult",
    "UnsupportedMagneticContextError",
    "collect_exchange_stage",
    "collect_exchange",
    "exchange_stage_collect",
    "exchange_stage_prepare",
    "exchange_stage_render_slurm",
    "exchange_stage_status",
    "prepare_exchange",
    "prepare_campaign_exchange",
    "prepare_exchange_stage",
    "render_exchange_slurm",
    "render_slurm_exchange",
    "status_exchange",
    "status_exchange_stage",
]
