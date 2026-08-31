"""Exchange job specifications and legacy-compatible serialization."""

from __future__ import annotations

import copy
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import f90nml

from rsutil.core import ExchangePair, PairKey, RSStructure
from rsutil.core.hashing import hash_file, hash_json

from .pairs import map_cluster_neighbor


EXCHANGE_MANIFEST_FILENAME = "rsutil_exchange_manifest.json"
EXCHANGE_MANIFEST_SCHEMA = "rsutil-exchange-v2"
_SUPPORTED_EXCHANGE_MANIFEST_SCHEMAS = frozenset({"rsutil-exchange-v1", EXCHANGE_MANIFEST_SCHEMA})
_COLLISION_POLICIES = frozenset({"error", "replace_generated", "reuse_if_identical"})
_RUNTIME_OUTPUT_NAMES = frozenset(
    {
        "report.out",
        "out.scf",
        "out.cold",
        "rsutil_execution.json",
        "jij.out",
        "jijso.out",
        "jijfo.out",
        "dij.out",
        "dijso.out",
        "dijfo.out",
        "aij.out",
        "aijso.out",
        "aijfo.out",
        "jtens.out",
    }
)
_RUNTIME_OUTPUT_SUFFIXES = frozenset({".log", ".out"})


def _is_runtime_output_name(name: str) -> bool:
    """Return whether an unowned directory entry is a post-generation artifact."""

    return name in _RUNTIME_OUTPUT_NAMES or Path(name).suffix.lower() in _RUNTIME_OUTPUT_SUFFIXES


class ExchangeJobError(RuntimeError):
    """Base class for transactional exchange-job publication failures."""


class ExchangeCollisionError(ExchangeJobError, FileExistsError):
    """A destination is not safely replaceable under the selected policy."""


class ExchangeManifestError(ExchangeJobError, ValueError):
    """A generated or existing exchange manifest is invalid."""


@dataclass(frozen=True)
class ExchangeWriteResult:
    """Publication evidence for one transactional exchange-write operation."""

    written: tuple[Path, ...] = ()
    reused: tuple[Path, ...] = ()


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def required_atom_files_from_input(input_nml: Mapping[str, Any]) -> tuple[str, ...]:
    """Return the explicit converged atom files required by an input namelist.

    The RS-LMTO convention is one ``<atoms.label>_out.nml`` per atom type.
    This list is derived from the supplied input object, never from the
    process working directory or an unrestricted glob.
    """

    try:
        labels = _as_list(input_nml["atoms"]["label"])
    except (KeyError, TypeError) as exc:
        raise ValueError("exchange input must define atoms.label") from exc
    result: list[str] = []
    for label in labels:
        filename = f"{str(label).strip()}_out.nml"
        if not str(label).strip() or Path(filename).name != filename:
            raise ValueError(f"invalid atom output label: {label!r}")
        if filename not in result:
            result.append(filename)
    if not result:
        raise ValueError("exchange input must define at least one atom output file")
    return tuple(result)


def build_exchange_input(
    input_nml: Mapping[str, Any],
    *,
    atom_index: int,
    matching_indices: Sequence[int],
    fermi_energy: float,
) -> dict[str, Any]:
    """Create one exchange input namelist without reading or writing files."""

    result = copy.deepcopy(input_nml)
    result.setdefault("energy", f90nml.Namelist())["fermi"] = fermi_energy
    calculation = result.setdefault("calculation", f90nml.Namelist())
    calculation["post_processing"] = "exchange"
    calculation.pop("pre_processing", None)

    lattice = result["lattice"]
    original_ct = lattice.get("ct")
    original_r2 = lattice.get("r2")
    rebuilt = f90nml.Namelist()
    for key, value in lattice.items():
        if key not in {"njij", "ijpair", "ct", "r2"}:
            rebuilt[key] = value
    rebuilt["njij"] = len(matching_indices)
    rebuilt["ijpair"] = [[atom_index] * len(matching_indices), list(matching_indices)]
    if original_ct is not None:
        rebuilt["ct"] = original_ct
    if original_r2 is not None:
        rebuilt["r2"] = original_r2
    result["lattice"] = rebuilt
    return result


@dataclass(frozen=True)
class ExchangeJobSpec:
    """A planned exchange calculation, independent of its disk layout."""

    atom_index: int
    center_cluster_index: int
    center_site_index: int
    matching_indices: tuple[int, ...]
    pairs: tuple[ExchangePair, ...]
    input_nml: Mapping[str, Any]
    output_dir: Path
    cutoff: float = 1.0
    reduction: str = "full"
    soc: bool = False
    magnetic_axis: tuple[float, float, float] | None = None
    magnetic_reference: str = "unknown"
    fermi_energy_source: str = "report.out"
    required_atom_files: tuple[str, ...] = ()
    source_run_identity: str = ""
    source_run_hash: str = ""
    nsp: int | None = None
    nsp_source: str = "unknown"
    is_collinear: bool | None = None
    is_noncollinear: bool | None = None
    soc_source: str = "unknown"
    symmetry_reduction_applied: bool = False
    symmetry_reduction_reason: str = ""

    @property
    def pair_keys(self) -> tuple[PairKey, ...]:
        return tuple(pair.key for pair in self.pairs)

    @property
    def canonical_pair_list(self) -> tuple[PairKey, ...]:
        """The complete pair set in manifest/canonical ordering."""

        return tuple(sorted(self.pair_keys, key=lambda key: key.sort_key))

    @property
    def spec_hash(self) -> str:
        """Deterministic identity of this job's scientific/job specification."""

        return hash_json(_spec_payload(self))

    @property
    def fermi_source(self) -> str:
        """Compatibility alias for the explicit Fermi-energy source path."""

        return self.fermi_energy_source


def construct_exchange_jobs(
    input_nml: Mapping[str, Any],
    cluster: Sequence[Any],
    neighbors: Mapping[int, Sequence[Any]],
    recursion_sites: Sequence[int],
    canonical_structure: RSStructure,
    *,
    fermi_energy: float,
    output_root: str | Path,
    cutoff: float = 1.0,
    reduction: str = "full",
    soc: bool = False,
    magnetic_axis: tuple[float, float, float] | None = None,
    magnetic_reference: str = "unknown",
    fermi_energy_source: str = "report.out",
    required_atom_files: Sequence[str] | None = None,
    source_run_identity: str = "",
    source_run_hash: str = "",
    nsp: int | None = None,
    nsp_source: str = "unknown",
    is_collinear: bool | None = None,
    is_noncollinear: bool | None = None,
    soc_source: str = "unknown",
    symmetry_reduction_applied: bool = False,
    symmetry_reduction_reason: str = "",
) -> tuple[ExchangeJobSpec, ...]:
    """Construct canonical pair-backed jobs while preserving legacy indices."""

    root = Path(output_root).expanduser().resolve()
    output_files = tuple(required_atom_files or required_atom_files_from_input(input_nml))
    jobs: list[ExchangeJobSpec] = []
    for center_cluster_index, atom_index in enumerate(recursion_sites):
        if atom_index <= 0 or atom_index > len(canonical_structure.sites):
            raise ValueError(f"recursion site {atom_index} is outside the reference cell")
        entries = tuple(neighbors.get(center_cluster_index, ()))
        pairs = tuple(
            map_cluster_neighbor(
                canonical_structure,
                cluster,
                center_cluster_index,
                neighbor,
                center_site_index=atom_index - 1,
            )
            for neighbor in entries
        )
        matching_indices = tuple(int(neighbor.index) + 1 for neighbor in entries)
        job_input = build_exchange_input(
            input_nml,
            atom_index=atom_index,
            matching_indices=matching_indices,
            fermi_energy=fermi_energy,
        )
        jobs.append(
            ExchangeJobSpec(
                atom_index=atom_index,
                center_cluster_index=center_cluster_index,
                center_site_index=atom_index - 1,
                matching_indices=matching_indices,
                pairs=pairs,
                input_nml=job_input,
                output_dir=root / f"Atom_{atom_index}",
                cutoff=float(cutoff),
                reduction=str(reduction),
                soc=bool(soc),
                magnetic_axis=(
                    None if magnetic_axis is None else tuple(float(value) for value in magnetic_axis)
                ),
                magnetic_reference=str(magnetic_reference),
                fermi_energy_source=str(fermi_energy_source),
                required_atom_files=output_files,
                source_run_identity=str(source_run_identity),
                source_run_hash=str(source_run_hash),
                nsp=None if nsp is None else int(nsp),
                nsp_source=str(nsp_source),
                is_collinear=is_collinear,
                is_noncollinear=is_noncollinear,
                soc_source=str(soc_source),
                symmetry_reduction_applied=bool(symmetry_reduction_applied),
                symmetry_reduction_reason=str(symmetry_reduction_reason),
            )
        )
    return tuple(jobs)


def inventory_pairs(jobs: Sequence[ExchangeJobSpec]) -> tuple[ExchangePair, ...]:
    """Return the complete canonical pair inventory in deterministic order."""

    pairs = [pair for job in jobs for pair in job.pairs]
    pairs.sort(key=lambda pair: pair.key.sort_key)
    if len({pair.key for pair in pairs}) != len(pairs):
        raise ValueError("exchange jobs contain duplicate canonical pair identities")
    return tuple(pairs)


def _logical_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _spec_payload(job: ExchangeJobSpec) -> dict[str, Any]:
    """Return the deterministic identity payload for a planned job."""

    pair_keys = [key.to_dict() for key in job.canonical_pair_list]
    return {
        "atom_index": job.atom_index,
        "center_cluster_index": job.center_cluster_index,
        "center_site_index": job.center_site_index,
        "cutoff": job.cutoff,
        "reduction": job.reduction,
        "soc": job.soc,
        "magnetic_axis": None if job.magnetic_axis is None else list(job.magnetic_axis),
        "magnetic_reference": job.magnetic_reference,
        "fermi_energy": float(job.input_nml["energy"]["fermi"]),
        "fermi_energy_source": job.fermi_energy_source,
        "required_atom_files": list(job.required_atom_files),
        "source_run_identity": job.source_run_identity,
        "source_run_hash": job.source_run_hash,
        "nsp": job.nsp,
        "nsp_source": job.nsp_source,
        "is_collinear": job.is_collinear,
        "is_noncollinear": job.is_noncollinear,
        "soc_source": job.soc_source,
        "symmetry_reduction_applied": job.symmetry_reduction_applied,
        "symmetry_reduction_reason": job.symmetry_reduction_reason,
        "pair_keys": pair_keys,
        "legacy_ijpair": [
            [job.atom_index] * len(job.matching_indices),
            list(job.matching_indices),
        ],
    }


def build_exchange_manifest(
    job: ExchangeJobSpec,
    *,
    root: str | Path,
    file_hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Build a versioned manifest for one staged exchange job directory."""

    root_path = Path(root).expanduser().resolve()
    spec = _spec_payload(job)
    return {
        "schema_version": EXCHANGE_MANIFEST_SCHEMA,
        "owner": "rsutil",
        "operation": "exchange_job",
        "spec_hash": hash_json(spec),
        "spec": spec,
        "canonical_pair_list": spec["pair_keys"],
        "pair_keys": spec["pair_keys"],
        "legacy_ijpair": spec["legacy_ijpair"],
        "source_run": {
            "identity": job.source_run_identity,
            "hash": job.source_run_hash,
            "fermi_energy_source": job.fermi_energy_source,
        },
        "files": dict(sorted((file_hashes or {}).items())),
        "paths": {
            "job": _logical_path(job.output_dir, root_path),
            "required_atom_files": list(job.required_atom_files),
        },
    }


def read_exchange_manifest(path: str | Path) -> dict[str, Any]:
    """Read and minimally validate an rsutil-owned exchange manifest."""

    manifest_path = Path(path)
    try:
        with manifest_path.open("r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ExchangeManifestError(f"could not read exchange manifest {manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise ExchangeManifestError(f"exchange manifest is not an object: {manifest_path}")
    if manifest.get("schema_version") not in _SUPPORTED_EXCHANGE_MANIFEST_SCHEMAS:
        raise ExchangeManifestError(
            f"unsupported exchange manifest schema in {manifest_path}: "
            f"{manifest.get('schema_version')!r}"
        )
    if manifest.get("owner") != "rsutil" or manifest.get("operation") != "exchange_job":
        raise ExchangeManifestError(f"manifest is not an rsutil exchange job: {manifest_path}")
    if not isinstance(manifest.get("spec"), dict) or not isinstance(manifest.get("spec_hash"), str):
        raise ExchangeManifestError(f"exchange manifest has no valid specification: {manifest_path}")
    if hash_json(manifest["spec"]) != manifest["spec_hash"]:
        raise ExchangeManifestError(f"exchange manifest specification hash is invalid: {manifest_path}")
    _validate_manifest_pair_identity(manifest, manifest_path)
    return manifest


def _validate_manifest_pair_identity(
    manifest: Mapping[str, Any], source: str | Path = "exchange manifest"
) -> list[dict[str, Any]]:
    """Validate pair compatibility views against the hashed specification.

    ``spec.pair_keys`` is the sole scientific identity.  The top-level views
    are retained for compatibility with older consumers, but are never
    authoritative.
    """

    spec = manifest.get("spec")
    if not isinstance(spec, Mapping) or not isinstance(spec.get("pair_keys"), list):
        raise ExchangeManifestError(f"{source} has no authoritative spec.pair_keys")
    authoritative = spec["pair_keys"]
    for name in ("canonical_pair_list", "pair_keys"):
        view = manifest.get(name)
        if view is not None and view != authoritative:
            raise ExchangeManifestError(
                f"{source} {name} does not match authoritative spec.pair_keys"
            )
    legacy = manifest.get("legacy_ijpair")
    if legacy is not None and legacy != spec.get("legacy_ijpair"):
        raise ExchangeManifestError(
            f"{source} legacy_ijpair does not match the hashed specification"
        )
    try:
        for item in authoritative:
            PairKey.from_dict(item)
    except (KeyError, TypeError, ValueError) as exc:
        raise ExchangeManifestError(f"{source} contains invalid canonical pair keys") from exc
    return authoritative


def _validate_exchange_input(path: Path, job: ExchangeJobSpec) -> None:
    try:
        namelist = f90nml.read(str(path))
        lattice = namelist["lattice"]
        njij = int(lattice["njij"])
        ijpair = lattice["ijpair"]
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise ExchangeManifestError(f"invalid generated exchange input: {path}") from exc
    if njij != len(job.pairs):
        raise ExchangeManifestError(
            f"{path}: njij={njij} does not match manifest pair count {len(job.pairs)}"
        )
    try:
        rows = [_as_list(row) for row in ijpair]
    except TypeError as exc:
        raise ExchangeManifestError(f"{path}: ijpair is not a two-row array") from exc
    if njij == 0 and rows == [[None], [None]]:
        # f90nml represents an intentionally empty indexed array as one
        # ``None`` placeholder per row when it is written back out.
        rows = [[], []]
    expected = [
        [job.atom_index] * len(job.matching_indices),
        list(job.matching_indices),
    ]
    if rows != expected:
        raise ExchangeManifestError(
            f"{path}: ijpair does not match the one-based RS/Fortran job convention"
        )
    if any(pair.key.i != job.center_site_index for pair in job.pairs):
        raise ExchangeManifestError(
            f"{path}: canonical pairs do not all use the declared reference center"
        )


def _validate_staged_job(stage_dir: Path, job: ExchangeJobSpec, manifest: Mapping[str, Any]) -> None:
    required = ("lattice.nml", "input.nml", *job.required_atom_files, EXCHANGE_MANIFEST_FILENAME)
    for filename in required:
        path = stage_dir / filename
        if not path.is_file():
            raise ExchangeManifestError(f"staged exchange job is missing {filename}")
    _validate_exchange_input(stage_dir / "input.nml", job)
    try:
        f90nml.read(str(stage_dir / "lattice.nml"))
        for filename in job.required_atom_files:
            f90nml.read(str(stage_dir / filename))
    except (OSError, ValueError, TypeError, KeyError) as exc:
        raise ExchangeManifestError(f"staged exchange job contains unreadable namelists") from exc
    expected_pairs = [key.to_dict() for key in job.canonical_pair_list]
    authoritative_pairs = _validate_manifest_pair_identity(manifest)
    if authoritative_pairs != expected_pairs:
        raise ExchangeManifestError("exchange manifest canonical pair list is inconsistent")
    if manifest.get("legacy_ijpair") != [
        [job.atom_index] * len(job.matching_indices),
        list(job.matching_indices),
    ]:
        raise ExchangeManifestError("exchange manifest ijpair is inconsistent")


def _owned_manifest(destination: Path) -> dict[str, Any] | None:
    if not destination.is_dir() or destination.is_symlink():
        return None
    marker = destination / EXCHANGE_MANIFEST_FILENAME
    if not marker.is_file():
        return None
    return read_exchange_manifest(marker)


def _destination_is_safe(destination: Path, root: Path) -> None:
    try:
        destination.resolve(strict=False).relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError(f"exchange job destination escapes source directory: {destination}") from exc


def _write_staged_job(
    stage_dir: Path,
    job: ExchangeJobSpec,
    *,
    lattice: Path,
    source_root: Path,
) -> dict[str, Any]:
    stage_dir.mkdir(parents=True, exist_ok=False)
    shutil.copy2(lattice, stage_dir / "lattice.nml")
    for filename in job.required_atom_files:
        source = source_root / filename
        if not source.is_file():
            raise FileNotFoundError(f"required converged atom file is missing: {source}")
        shutil.copy2(source, stage_dir / filename)
    f90nml.write(job.input_nml, str(stage_dir / "input.nml"), force=True)
    hashes = {
        filename: hash_file(stage_dir / filename)
        for filename in ("lattice.nml", "input.nml", *job.required_atom_files)
    }
    manifest = build_exchange_manifest(job, root=source_root, file_hashes=hashes)
    with (stage_dir / EXCHANGE_MANIFEST_FILENAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    _validate_staged_job(stage_dir, job, manifest)
    return manifest


def _existing_files_match_manifest(
    destination: Path, manifest: Mapping[str, Any], expected_names: set[str]
) -> bool:
    files = manifest.get("files")
    if not isinstance(files, Mapping):
        return False
    if set(files) != expected_names - {EXCHANGE_MANIFEST_FILENAME}:
        return False
    try:
        return all(hash_file(destination / name) == str(digest) for name, digest in files.items())
    except (OSError, TypeError):
        return False


def write_exchange_jobs(
    jobs: Sequence[ExchangeJobSpec],
    *,
    lattice_path: str | Path,
    source_dir: str | Path,
    collision_policy: str = "error",
    force_legacy_collision: bool = False,
) -> ExchangeWriteResult:
    """Stage, validate, and transactionally publish familiar ``Atom_*`` jobs.

    ``error`` is the safe default.  ``replace_generated`` can replace only a
    directory carrying a valid rsutil ownership marker.  ``reuse_if_identical``
    leaves an identical owned directory untouched and rejects a stale one.
    An unmarked directory is never removed unless ``force_legacy_collision``
    is explicitly supplied.
    """

    return _write_exchange_jobs(
        jobs,
        lattice_path=lattice_path,
        source_dir=source_dir,
        collision_policy=collision_policy,
        force_legacy_collision=force_legacy_collision,
    )


def _write_exchange_jobs(
    jobs: Sequence[ExchangeJobSpec],
    *,
    lattice_path: str | Path,
    source_dir: str | Path,
    collision_policy: str,
    force_legacy_collision: bool,
) -> ExchangeWriteResult:
    """Implementation for the public writer and service collision policies."""

    if collision_policy not in _COLLISION_POLICIES:
        raise ValueError(f"collision_policy must be one of {sorted(_COLLISION_POLICIES)}")
    if not jobs:
        return ExchangeWriteResult()

    root = Path(source_dir).expanduser().resolve()
    lattice = Path(lattice_path).expanduser()
    if not lattice.is_absolute():
        lattice = root / lattice
    if not lattice.is_file():
        raise FileNotFoundError(lattice)
    root.mkdir(parents=True, exist_ok=True)

    stage_root = Path(tempfile.mkdtemp(prefix=".rsutil-exchange-stage-", dir=str(root)))
    prepared: list[tuple[ExchangeJobSpec, Path, dict[str, Any], Path]] = []
    reused: list[Path] = []
    try:
        for index, job in enumerate(jobs):
            destination = Path(job.output_dir).expanduser()
            if not destination.is_absolute():
                destination = root / destination
            destination = destination.resolve(strict=False)
            _destination_is_safe(destination, root)
            stage_dir = stage_root / f"job-{index}"
            manifest = _write_staged_job(
                stage_dir,
                job,
                lattice=lattice,
                source_root=root,
            )
            existing = destination.exists() or destination.is_symlink()
            if existing:
                owned = _owned_manifest(destination)
                if collision_policy == "error":
                    raise ExchangeCollisionError(
                        f"exchange job destination exists; use collision_policy='replace_generated' "
                        f"or 'reuse_if_identical': {destination}"
                    )
                if owned is None and (
                    not force_legacy_collision or not destination.is_dir() or destination.is_symlink()
                ):
                    raise ExchangeCollisionError(
                        f"refusing to replace unmarked legacy/user directory: {destination}"
                    )
                if collision_policy == "reuse_if_identical":
                    if owned is None or owned.get("spec_hash") != manifest.get("spec_hash"):
                        raise ExchangeCollisionError(
                            f"existing exchange job is not identical; refusing stale reuse: {destination}"
                        )
                    expected_names = {
                        "lattice.nml",
                        "input.nml",
                        *job.required_atom_files,
                        EXCHANGE_MANIFEST_FILENAME,
                    }
                    actual_names = {item.name for item in destination.iterdir()}
                    unexpected = {
                        name for name in actual_names - expected_names
                        if not _is_runtime_output_name(name)
                    }
                    if unexpected or not _existing_files_match_manifest(
                        destination, manifest, expected_names
                    ):
                        raise ExchangeCollisionError(
                            f"existing exchange job contains unexpected/stale files or hash mismatches against "
                            f"freshly generated inputs; "
                            f"use collision_policy='replace_generated': {destination}"
                        )
                    reused.append(destination)
                    continue
            prepared.append((job, destination, manifest, stage_dir))

        backups: list[tuple[Path, Path]] = []
        committed: list[tuple[Path, Path | None]] = []
        try:
            for commit_index, (_job, destination, _manifest, stage_dir) in enumerate(prepared):
                backup: Path | None = None
                if destination.exists() or destination.is_symlink():
                    backup = root / f".rsutil-exchange-backup-{commit_index}-{os.getpid()}"
                    suffix = 0
                    while backup.exists():
                        suffix += 1
                        backup = root / f".rsutil-exchange-backup-{commit_index}-{os.getpid()}-{suffix}"
                    os.replace(destination, backup)
                    backups.append((destination, backup))
                os.replace(stage_dir, destination)
                committed.append((destination, backup))
        except Exception:
            for destination, backup in reversed(committed):
                if destination.exists() or destination.is_symlink():
                    shutil.rmtree(destination)
                if backup is not None and backup.exists():
                    os.replace(backup, destination)
            for destination, backup in reversed(backups):
                if backup.exists() and not destination.exists():
                    os.replace(backup, destination)
            raise
        else:
            for _destination, backup in backups:
                if backup.exists():
                    shutil.rmtree(backup)
        written = tuple(destination for _job, destination, _manifest, _stage_dir in prepared)
        for destination in written:
            print(f"Exchange job committed to {destination}")
        return ExchangeWriteResult(written=written, reused=tuple(reused))
    finally:
        shutil.rmtree(stage_root, ignore_errors=True)


__all__ = [
    "ExchangeJobSpec",
    "ExchangeJobError",
    "ExchangeCollisionError",
    "ExchangeManifestError",
    "ExchangeWriteResult",
    "EXCHANGE_MANIFEST_FILENAME",
    "EXCHANGE_MANIFEST_SCHEMA",
    "build_exchange_input",
    "build_exchange_manifest",
    "construct_exchange_jobs",
    "inventory_pairs",
    "read_exchange_manifest",
    "required_atom_files_from_input",
    "write_exchange_jobs",
    "_write_exchange_jobs",
]
