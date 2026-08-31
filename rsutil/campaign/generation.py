"""Specification loading, validation, and deterministic campaign generation."""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from ..core.contracts import (
    Diagnostic,
    DiagnosticSeverity,
    ErrorCode,
    OperationStatus,
)
from ..core.hashing import hash_file, scientific_hash
from ..setup_service import SetupService
from ..setup_validation import validate_setup_namelists
from .models import Campaign, CampaignEntry, CampaignIndex, CampaignIndexEntry


CAMPAIGN_INDEX_FILENAME = "campaign_index.json"
RUN_MANIFEST_FILENAME = "run_manifest.json"


class CampaignError(RuntimeError):
    """Base error raised for invalid or unpublishable campaigns."""


class CampaignCollisionError(CampaignError):
    """An existing path is not the identical run requested by the campaign."""


def _diagnostic(code: ErrorCode | str, message: str, *, path: str | None = None, **details: Any) -> Diagnostic:
    return Diagnostic(DiagnosticSeverity.ERROR, code, message, path=path, details=details)


@dataclass(frozen=True)
class CampaignValidation:
    """Machine-readable result of campaign-wide preflight validation."""

    valid: bool
    diagnostics: tuple[Diagnostic, ...] = ()
    entry_ids: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "entry_ids", dict(self.entry_ids or {}))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": "rsutil-campaign-validation-v1",
            "valid": self.valid,
            "entry_ids": dict(self.entry_ids),
            "diagnostics": [item.to_dict() for item in self.diagnostics],
        }


def _load_json_or_yaml(path: Path) -> Mapping[str, Any]:
    try:
        if path.suffix.lower() in {".yaml", ".yml"}:
            try:
                import yaml  # type: ignore[import-not-found]
            except ImportError as exc:
                raise CampaignError("YAML campaign files require PyYAML; use the primary JSON format") from exc
            data = yaml.safe_load(path.read_text(encoding="utf-8"))
        else:
            data = json.loads(path.read_text(encoding="utf-8"))
    except CampaignError:
        raise
    except Exception as exc:
        raise CampaignError(f"could not parse campaign specification {path}: {exc}") from exc
    if not isinstance(data, Mapping):
        raise CampaignError("campaign specification must be a JSON/YAML object")
    return data


def load_campaign(path: str | Path) -> Campaign:
    """Load the primary JSON campaign format (and YAML when PyYAML is present)."""

    source = Path(path).resolve()
    if not source.is_file():
        raise CampaignError(f"campaign specification does not exist: {source}")
    try:
        return Campaign.from_dict(_load_json_or_yaml(source))
    except CampaignError:
        raise
    except Exception as exc:
        raise CampaignError(f"invalid campaign specification {source}: {exc}") from exc


def _resolve_path(value: str | Path, base_dir: Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base_dir / path).resolve()


def _source_value(campaign: Campaign, entry: CampaignEntry) -> Any:
    if entry.structure is not None:
        if entry.structure not in campaign.structures:
            raise CampaignError(f"entry {entry.name!r} references unknown structure {entry.structure!r}")
        return campaign.structures[entry.structure]
    return entry.source


def _source_files(campaign: Campaign, entry: CampaignEntry, base_dir: Path) -> tuple[Path, Path]:
    value = _source_value(campaign, entry)
    if isinstance(value, Mapping):
        lattice = value.get("lattice", value.get("lattice_nml"))
        input_nml = value.get("input", value.get("input_nml"))
        if lattice is None or input_nml is None:
            raise CampaignError(f"entry {entry.name!r} source mapping needs lattice and input files")
        return _resolve_path(str(lattice), base_dir), _resolve_path(str(input_nml), base_dir)
    if not isinstance(value, (str, Path)):
        raise CampaignError(f"entry {entry.name!r} source must be a directory or file mapping")
    directory = _resolve_path(value, base_dir)
    if not directory.is_dir():
        raise CampaignError(f"entry {entry.name!r} source directory does not exist: {directory}")
    return directory / "lattice.nml", directory / "input.nml"


def _source_hashes(lattice: Path, input_nml: Path) -> dict[str, Any]:
    return {
        "lattice.nml": {"sha256": hash_file(lattice), "size": lattice.stat().st_size},
        "input.nml": {"sha256": hash_file(input_nml), "size": input_nml.stat().st_size},
    }


def _effective_setup(entry: CampaignEntry) -> Any:
    return entry.setup_spec


def _validate_determinism(entry: CampaignEntry) -> list[Diagnostic]:
    spec = _effective_setup(entry)
    diagnostics: list[Diagnostic] = []
    if (spec.lattice_dist or spec.basis_dist) and spec.seed is None:
        diagnostics.append(_diagnostic(ErrorCode.UNSEEDED_SETUP, "stochastic lattice/basis distortion requires setup.seed", path=entry.name))
    if spec.sqs is not None and spec.sqs_seed is None:
        diagnostics.append(_diagnostic(ErrorCode.UNSEEDED_SETUP, "stochastic SQS setup requires setup.sqs_seed", path=entry.name))
    return diagnostics


def _validate_campaign(campaign: Campaign, base_dir: Path) -> CampaignValidation:
    diagnostics: list[Diagnostic] = []
    entry_ids: dict[str, str] = {}
    for entry in campaign.entries:
        try:
            lattice, input_nml = _source_files(campaign, entry, base_dir)
            if not lattice.is_file() or not input_nml.is_file():
                missing = [str(path) for path in (lattice, input_nml) if not path.is_file()]
                diagnostics.append(_diagnostic(ErrorCode.CAMPAIGN_SOURCE_MISSING, "campaign source file is missing", path=entry.name, missing=missing))
                continue
            source_hashes = _source_hashes(lattice, input_nml)
            entry_ids[entry.name] = entry.identity_hash(
                source_hashes=source_hashes,
                setup_template=campaign.setup_template,
                calculation_template=campaign.calculation_template,
            )
            diagnostics.extend(_validate_determinism(entry))
            validation = validate_setup_namelists(lattice, input_nml, spec=entry.setup_spec)
            diagnostics.extend(validation.diagnostics)
        except CampaignError as exc:
            diagnostics.append(_diagnostic(ErrorCode.CAMPAIGN_SPEC_INVALID, str(exc), path=entry.name))
        except Exception as exc:
            diagnostics.append(_diagnostic(ErrorCode.CAMPAIGN_SPEC_INVALID, f"could not validate entry: {exc}", path=entry.name))
    return CampaignValidation(
        valid=not any(item.severity is DiagnosticSeverity.ERROR for item in diagnostics),
        diagnostics=tuple(diagnostics),
        entry_ids=entry_ids,
    )


def validate_campaign(campaign: Campaign, *, base_dir: str | Path = ".") -> CampaignValidation:
    """Validate all entries without creating or changing a destination."""

    return _validate_campaign(campaign, Path(base_dir).resolve())


def _slug(value: str) -> str:
    result = re.sub(r"[^A-Za-z0-9._-]+", "-", value.strip()).strip(".-")
    return result or "entry"


def _relative_path(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent), text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


@contextmanager
def _materialize_source(lattice: Path, input_nml: Path):
    """Yield a setup source directory, copying explicitly named source files."""

    if lattice.parent == input_nml.parent and lattice.name == "lattice.nml" and input_nml.name == "input.nml":
        yield lattice.parent
        return
    with tempfile.TemporaryDirectory(prefix="rsutil-campaign-source-") as temporary:
        root = Path(temporary)
        shutil.copy2(lattice, root / "lattice.nml")
        shutil.copy2(input_nml, root / "input.nml")
        yield root


def _read_run_manifest(path: Path) -> Mapping[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise CampaignCollisionError(f"could not read existing run manifest {path}: {exc}") from exc
    if not isinstance(data, Mapping) or not isinstance(data.get("campaign"), Mapping):
        raise CampaignCollisionError(f"existing directory is not an rsutil campaign run: {path.parent}")
    return data


def _logicalize_manifest_paths(value: Any, *, lattice: Path, input_nml: Path) -> Any:
    """Remove machine-local paths from diagnostic records in a run manifest."""

    if isinstance(value, Mapping):
        return {
            key: _logicalize_manifest_paths(item, lattice=lattice, input_nml=input_nml)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_logicalize_manifest_paths(item, lattice=lattice, input_nml=input_nml) for item in value]
    if isinstance(value, str):
        for source, logical in ((str(lattice), "source/lattice.nml"), (str(input_nml), "source/input.nml")):
            if value == source or value.startswith(source + ":"):
                return logical + value[len(source):]
        # Diagnostic paths are evidence labels, not scientific identity. Keep
        # their basename when a worker/temp path is the only available value.
        if value.startswith("/"):
            path_part, separator, suffix = value.partition(":")
            return "generated/" + Path(path_part).name + (separator + suffix if separator else "")
    return value


def _logical_source(campaign: Campaign, entry: CampaignEntry, base_dir: Path) -> Any:
    if entry.structure is not None:
        return entry.structure
    if isinstance(entry.source, Mapping):
        result = {}
        for key, value in entry.source.items():
            if isinstance(value, (str, Path)):
                path = _resolve_path(value, base_dir)
                try:
                    result[key] = path.relative_to(base_dir).as_posix()
                except ValueError:
                    result[key] = path.name
            else:
                result[key] = value
        return result
    if isinstance(entry.source, (str, Path)):
        path = _resolve_path(entry.source, base_dir)
        try:
            return path.relative_to(base_dir).as_posix()
        except ValueError:
            return path.name
    return entry.source


class CampaignGenerator:
    """Generate setup-owned run directories from a validated Campaign."""

    def __init__(self, setup_service: SetupService | None = None) -> None:
        self.setup_service = setup_service or SetupService()

    def validate(self, campaign: Campaign, *, base_dir: str | Path = ".") -> CampaignValidation:
        return validate_campaign(campaign, base_dir=base_dir)

    def generate(
        self,
        campaign: Campaign,
        destination: str | Path | None = None,
        *,
        base_dir: str | Path = ".",
        collision_policy: str = "reuse_if_identical",
    ) -> CampaignIndex:
        if collision_policy not in {"error", "reuse_if_identical"}:
            raise CampaignError("collision_policy must be 'error' or 'reuse_if_identical'")
        base = Path(base_dir).resolve()
        root = _resolve_path(destination or campaign.destination or f"{_slug(campaign.name)}-campaign", base)
        validation = self.validate(campaign, base_dir=base)
        if not validation.valid:
            details = "; ".join(f"{item.path or 'campaign'}: {item.message}" for item in validation.diagnostics if item.severity is DiagnosticSeverity.ERROR)
            raise CampaignError(f"campaign preflight validation failed: {details}")
        if root.exists() and not root.is_dir():
            raise CampaignCollisionError(f"campaign destination is not a directory: {root}")

        index_path = root / CAMPAIGN_INDEX_FILENAME
        if index_path.exists():
            try:
                old_index = CampaignIndex.from_dict(json.loads(index_path.read_text(encoding="utf-8")))
            except Exception as exc:
                raise CampaignCollisionError(f"existing campaign index is invalid: {index_path}") from exc
            if old_index.campaign_name != campaign.name:
                raise CampaignCollisionError(
                    f"destination already contains campaign {old_index.campaign_name!r}; refusing to reuse {root}"
                )

        entry_ids = dict(validation.entry_ids)
        campaign_id = campaign.campaign_hash(list(entry_ids.values()))
        spec_hash = scientific_hash({"campaign": campaign.to_dict(), "entry_ids": entry_ids})
        index_entries: list[CampaignIndexEntry] = []
        for entry in campaign.entries:
            entry_id = entry_ids[entry.name]
            relative_run = entry.destination or campaign.destination_naming.format(
                name=_slug(entry.name), entry_id=entry_id, campaign_id=campaign_id
            )
            run_path = Path(relative_run)
            if run_path.is_absolute() or ".." in run_path.parts:
                raise CampaignError(f"unsafe destination naming for entry {entry.name!r}: {relative_run!r}")
            run_dir = (root / run_path).resolve()
            if root not in run_dir.parents:
                raise CampaignError(f"destination naming escapes campaign root: {relative_run!r}")
            manifest_hash = scientific_hash({"schema_version": "rsutil-run-v1", "entry_id": entry_id})
            manifest_path = run_dir / RUN_MANIFEST_FILENAME
            if run_dir.exists():
                if collision_policy == "error":
                    raise CampaignCollisionError(f"run directory already exists: {run_dir}")
                if not manifest_path.is_file():
                    raise CampaignCollisionError(f"existing run directory has no campaign manifest: {run_dir}")
                existing = _read_run_manifest(manifest_path)
                existing_campaign = existing["campaign"]
                if existing_campaign.get("entry_id") != entry_id or existing_campaign.get("manifest_hash") != manifest_hash:
                    raise CampaignCollisionError(f"existing run directory is incompatible with entry {entry.name!r}: {run_dir}")
                index_entries.append(CampaignIndexEntry(
                    entry_id=entry_id,
                    name=entry.name,
                    run_path=_relative_path(run_dir, root),
                    manifest_path=_relative_path(manifest_path, root),
                    manifest_hash=manifest_hash,
                ))
                continue
            lattice, input_nml = _source_files(campaign, entry, base)
            with _materialize_source(lattice, input_nml) as source_dir:
                result = self.setup_service.run(
                    source_dir,
                    run_dir,
                    entry.setup_spec,
                    collision_policy="error",
                )
            if not result.success:
                message = "; ".join(item.message for item in result.diagnostics) or "setup service failed"
                raise CampaignError(f"entry {entry.name!r} generation failed: {message}")
            setup_manifest = result.provenance.get("manifest", {})
            if not isinstance(setup_manifest, Mapping):
                setup_manifest = {}
            run_manifest = dict(setup_manifest)
            run_manifest.setdefault("schema_version", "rsutil-run-v1")
            run_manifest.setdefault("operation", "setup")
            run_manifest.setdefault("status", OperationStatus.SUCCESS.value)
            run_manifest["campaign"] = {
                "campaign_name": campaign.name,
                "entry_name": entry.name,
                "entry_id": entry_id,
                "manifest_hash": manifest_hash,
                "source": _logical_source(campaign, entry, base),
                "calculation": dict(campaign.calculation_template) | dict(entry.calculation),
                "exchange": None if entry.exchange is None else dict(entry.exchange),
            }
            run_manifest = _logicalize_manifest_paths(
                run_manifest, lattice=lattice, input_nml=input_nml
            )
            _atomic_json(manifest_path, run_manifest)
            index_entries.append(CampaignIndexEntry(
                entry_id=entry_id,
                name=entry.name,
                run_path=_relative_path(run_dir, root),
                manifest_path=_relative_path(manifest_path, root),
                manifest_hash=manifest_hash,
            ))

        index = CampaignIndex(
            campaign_id=campaign_id,
            campaign_name=campaign.name,
            campaign_spec_hash=spec_hash,
            entries=tuple(index_entries),
        )
        _atomic_json(index_path, index.to_dict())
        return index


def generate_campaign(
    campaign: Campaign,
    destination: str | Path | None = None,
    *,
    base_dir: str | Path = ".",
    setup_service: SetupService | None = None,
    collision_policy: str = "reuse_if_identical",
) -> CampaignIndex:
    return CampaignGenerator(setup_service=setup_service).generate(
        campaign, destination, base_dir=base_dir, collision_policy=collision_policy
    )


__all__ = [
    "CAMPAIGN_INDEX_FILENAME",
    "RUN_MANIFEST_FILENAME",
    "CampaignCollisionError",
    "CampaignError",
    "CampaignGenerator",
    "CampaignValidation",
    "generate_campaign",
    "load_campaign",
    "validate_campaign",
]
