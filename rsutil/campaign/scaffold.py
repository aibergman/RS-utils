"""Discover RS-LMTO source bundles and scaffold campaign specifications."""

from __future__ import annotations

import fnmatch
import hashlib
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..core.hashing import hash_file
from .generation import CampaignError, CampaignValidation, _atomic_json, validate_campaign
from .models import CAMPAIGN_SCHEMA_VERSION, Campaign


class CampaignScaffoldError(CampaignError):
    """Raised when a source tree cannot be converted into a campaign."""


@dataclass(frozen=True)
class SourceBundle:
    """A complete or incomplete source bundle candidate."""

    directory: Path
    relative_path: str
    has_input: bool
    has_lattice: bool
    has_source_cif: bool
    cif_files: tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        return self.has_input and self.has_lattice

    @property
    def missing(self) -> tuple[str, ...]:
        return tuple(
            filename
            for filename, present in (
                ("input.nml", self.has_input),
                ("lattice.nml", self.has_lattice),
            )
            if not present
        )


@dataclass(frozen=True)
class CampaignScaffoldResult:
    """Details of a successfully discovered and validated scaffold."""

    campaign: Campaign
    output: Path
    source_root: Path
    discovered: int
    complete: int
    with_source_cif: int
    validation: CampaignValidation
    ambiguous_cif: tuple[tuple[str, tuple[str, ...]], ...] = ()


def _safe_component(value: str) -> str:
    """Convert one path component into a portable campaign identifier."""

    result = re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip(".-")
    return result or "entry"


def _bundle_alias(bundle: SourceBundle) -> str:
    """Return a stable alias, retaining simple directory names verbatim."""

    parts = Path(bundle.relative_path).parts
    if not parts or parts == (".",):
        return "root"
    components = [_safe_component(part) for part in parts]
    if len(components) == 1:
        return components[0]

    # Keep nested aliases readable while making path-derived identity robust
    # against ambiguous separators and sanitized component names.
    relative = Path(bundle.relative_path).as_posix()
    digest = hashlib.sha256(relative.encode("utf-8")).hexdigest()[:8]
    return "--".join((*components, digest))


def _matches_filter(bundle: SourceBundle, patterns: tuple[str, ...]) -> bool:
    """Match a bundle against aliases and POSIX relative paths."""

    if not patterns:
        return False
    alias = _bundle_alias(bundle)
    return any(
        fnmatch.fnmatchcase(candidate, pattern)
        for pattern in patterns
        for candidate in (bundle.relative_path, alias)
    )


def discover_source_bundles(
    source_root: str | Path,
    *,
    include: tuple[str, ...] = (),
    exclude: tuple[str, ...] = (),
) -> tuple[SourceBundle, ...]:
    """Recursively discover complete and incomplete RS-LMTO candidates.

    A candidate is any directory containing at least one of ``input.nml`` and
    ``lattice.nml``.  Directories containing neither file are ignored.
    ``include`` and ``exclude`` use case-sensitive shell-style glob matching
    against both the POSIX relative path and the generated source alias.
    Results are sorted by relative path for reproducible campaign ordering.
    """

    root = Path(source_root).resolve()
    if not root.exists():
        raise CampaignScaffoldError(f"source root does not exist: {root}")
    if not root.is_dir():
        raise CampaignScaffoldError(f"source root is not a directory: {root}")

    directories = [root]
    directories.extend(
        path
        for path in root.rglob("*")
        if path.is_dir() and not path.is_symlink()
    )
    candidates: list[SourceBundle] = []
    for directory in sorted(directories, key=lambda path: path.relative_to(root).as_posix()):
        has_input = (directory / "input.nml").is_file()
        has_lattice = (directory / "lattice.nml").is_file()
        if not has_input and not has_lattice:
            continue
        relative = directory.relative_to(root).as_posix()
        candidate = SourceBundle(
            directory=directory,
            relative_path=relative,
            has_input=has_input,
            has_lattice=has_lattice,
            has_source_cif=(directory / "source.cif").is_file(),
            cif_files=tuple(sorted(path.name for path in directory.glob("*.cif") if path.is_file())),
        )
        if include and not _matches_filter(candidate, include):
            continue
        if exclude and _matches_filter(candidate, exclude):
            continue
        candidates.append(candidate)
    return tuple(candidates)


def _relative_source_path(directory: Path, output: Path) -> str:
    """Represent a source directory relative to the campaign JSON."""

    return Path(os.path.relpath(directory, output.parent)).as_posix()


def _incomplete_message(candidates: tuple[SourceBundle, ...]) -> str:
    lines = ["incomplete source bundles found:"]
    for candidate in candidates:
        relative = candidate.relative_path or "."
        lines.append(f"  {relative}/ missing {', '.join(candidate.missing)}")
    return "\n".join(lines)


def scaffold_campaign(
    source_root: str | Path,
    *,
    name: str,
    output: str | Path,
    include: tuple[str, ...] = (),
    exclude: tuple[str, ...] = (),
    dry_run: bool = False,
) -> CampaignScaffoldResult:
    """Create and validate a campaign specification from source bundles.

    The campaign JSON is written only after the existing campaign validator
    accepts the in-memory specification, unless ``dry_run`` is true. Source
    paths are relative to the generated JSON, keeping a moved campaign and
    source tree relocatable.
    """

    if not isinstance(name, str) or not name.strip():
        raise CampaignScaffoldError("campaign name must be non-empty")

    root = Path(source_root).resolve()
    output_path = Path(output).resolve()
    candidates = discover_source_bundles(root, include=include, exclude=exclude)
    incomplete = tuple(candidate for candidate in candidates if not candidate.complete)
    if incomplete:
        raise CampaignScaffoldError(_incomplete_message(incomplete))
    if not candidates:
        raise CampaignScaffoldError(f"no RS-LMTO source bundles found under: {root}")

    structures: dict[str, str] = {}
    aliases: dict[str, str] = {}
    ambiguous_cif: list[tuple[str, tuple[str, ...]]] = []
    structure_sources: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        alias = _bundle_alias(candidate)
        previous = aliases.get(alias)
        if previous is not None:
            raise CampaignScaffoldError(
                f"source alias collision for {alias!r}: {previous} and {candidate.relative_path}"
            )
        aliases[alias] = candidate.relative_path
        structures[alias] = _relative_source_path(candidate.directory, output_path)
        if candidate.has_source_cif:
            source_cif = candidate.directory / "source.cif"
            try:
                structure_sources[alias] = {
                    "kind": "cif",
                    "file": "source.cif",
                    "sha256": hash_file(source_cif),
                }
            except OSError as exc:
                raise CampaignScaffoldError(f"could not hash source.cif for {candidate.relative_path}: {exc}") from exc
        elif len(candidate.cif_files) > 1:
            ambiguous_cif.append((candidate.relative_path or ".", candidate.cif_files))

    raw_spec: dict[str, Any] = {
        "schema_version": CAMPAIGN_SCHEMA_VERSION,
        "name": name,
        "structures": structures,
        "metadata": {"structure_sources": structure_sources},
        "entries": [
            {
                "name": alias,
                "structure": alias,
                "calculation": {"mode": "scf"},
            }
            for alias in structures
        ],
    }
    try:
        campaign = Campaign.from_dict(raw_spec)
    except Exception as exc:
        raise CampaignScaffoldError(f"could not build campaign specification: {exc}") from exc

    validation = validate_campaign(campaign, base_dir=output_path.parent)
    if not validation.valid:
        details = "; ".join(
            f"{item.path or 'campaign'}: {item.message}"
            for item in validation.diagnostics
            if getattr(item.severity, "value", item.severity) == "error"
        )
        raise CampaignScaffoldError(f"campaign validation failed: {details}")

    if output_path.exists() and output_path.is_dir():
        raise CampaignScaffoldError(f"output path is a directory: {output_path}")
    if not dry_run:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        _atomic_json(output_path, campaign.to_dict())
    with_source_cif = sum(1 for candidate in candidates if candidate.has_source_cif)
    return CampaignScaffoldResult(
        campaign=campaign,
        output=output_path,
        source_root=root,
        discovered=len(candidates),
        complete=len(candidates),
        with_source_cif=with_source_cif,
        validation=validation,
        ambiguous_cif=tuple(ambiguous_cif),
    )


__all__ = [
    "CampaignScaffoldError",
    "CampaignScaffoldResult",
    "SourceBundle",
    "discover_source_bundles",
    "scaffold_campaign",
]
