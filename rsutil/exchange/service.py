"""Reusable, path-explicit exchange planning and generation service."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from rsutil.core import ExchangePair, PairKey
from rsutil.core.hashing import hash_file, hash_json

from .jobs import (
    ExchangeJobSpec,
    ExchangeWriteResult,
    construct_exchange_jobs,
    inventory_pairs,
    required_atom_files_from_input,
    write_exchange_jobs,
)
from .neighbors import discover_neighbors
from .pairs import legacy_displacement_mapping
from .reference import ReferenceCalculation, read_reference_calculation
from .writer import write_Rij_mappings


@dataclass(frozen=True)
class ExchangePlan:
    """Planned exchange jobs and their canonical periodic pair inventory."""

    reference: ReferenceCalculation
    neighbors: dict[int, list]
    jobs: tuple[ExchangeJobSpec, ...]
    pairs: tuple[ExchangePair, ...]

    @property
    def pair_keys(self) -> tuple[PairKey, ...]:
        return tuple(pair.key for pair in self.pairs)


@dataclass(frozen=True)
class ExchangeResult:
    """Result of exchange-job generation."""

    plan: ExchangePlan
    mapping_path: Path | None = None
    write_result: ExchangeWriteResult | None = None

    @property
    def pairs(self) -> tuple[ExchangePair, ...]:
        return self.plan.pairs

    @property
    def jobs(self) -> tuple[ExchangeJobSpec, ...]:
        return self.plan.jobs

    @property
    def pair_keys(self) -> tuple[PairKey, ...]:
        return self.plan.pair_keys


def _resolve_output(root: Path, path: str | Path) -> Path:
    candidate = Path(path).expanduser()
    return candidate if candidate.is_absolute() else root / candidate


def _legacy_mapping_records(plan: ExchangePlan) -> list[dict]:
    """Keep legacy cluster indices at the CLI serialization boundary."""

    records: list[dict] = []
    for job in plan.jobs:
        for pair in job.pairs:
            record = legacy_displacement_mapping(pair)
            record["i"] = job.atom_index
            if pair.cluster_index is not None:
                record["j"] = pair.cluster_index + 1
            records.append(record)
    return records


def _source_run_identity(reference: ReferenceCalculation, required_files: tuple[str, ...]) -> str:
    """Hash the explicit source evidence used to construct every job."""

    paths = (
        reference.input_path,
        reference.lattice_path,
        reference.cluster_path,
        reference.report_path,
        *(reference.source_dir / name for name in required_files),
    )
    records = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"required exchange source file is missing: {path}")
        records.append(
            {
                "path": _logical_source_path(path, reference.source_dir),
                "sha256": hash_file(path),
            }
        )
    return hash_json(records)


def _logical_source_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


class ExchangeService:
    """Plan and serialize exchange jobs without changing process CWD."""

    def __init__(self, source_dir: str | Path = ".") -> None:
        self.source_dir = Path(source_dir).expanduser().resolve()

    def plan(
        self,
        *,
        cutoff: float = 1.0,
        unique: bool = False,
        shell_width: float = 0.1,
        soc: bool = False,
        full_pairs: bool = False,
        magnetic_axis: tuple[float, float, float] | None = None,
        verbose: bool = False,
        input_path: str | Path | None = None,
        lattice_path: str | Path | None = None,
        cluster_path: str | Path | None = None,
        report_path: str | Path | None = None,
    ) -> ExchangePlan:
        """Read inputs, discover neighbors, and construct canonical jobs."""

        reference = read_reference_calculation(
            self.source_dir,
            input_path=input_path,
            lattice_path=lattice_path,
            cluster_path=cluster_path,
            report_path=report_path,
            verbose=verbose,
        )
        required_files = required_atom_files_from_input(reference.input_nml)
        source_identity = _source_run_identity(reference, required_files)
        neighbors = discover_neighbors(
            reference.cluster,
            reference.recursion_sites,
            cutoff=cutoff,
            unique=unique,
            shell_width=shell_width,
            analyzer=reference.analyzer,
            structure=reference.structure,
            canonical_structure=reference.canonical_structure,
            soc=soc,
            magnetic_axis=magnetic_axis,
            full_pairs=full_pairs,
            verbose=verbose,
        )
        jobs = construct_exchange_jobs(
            reference.input_nml,
            reference.cluster,
            neighbors,
            reference.recursion_sites,
            reference.canonical_structure,
            fermi_energy=reference.fermi_energy,
            output_root=self.source_dir,
            cutoff=cutoff,
            reduction="unique" if unique and not full_pairs else "full",
            soc=soc,
            magnetic_axis=magnetic_axis,
            fermi_energy_source=_logical_source_path(reference.report_path, self.source_dir),
            required_atom_files=required_files,
            source_run_identity=source_identity,
            source_run_hash=source_identity,
        )
        return ExchangePlan(
            reference=reference,
            neighbors=neighbors,
            jobs=jobs,
            pairs=inventory_pairs(jobs),
        )

    def generate(
        self,
        *,
        cutoff: float = 1.0,
        unique: bool = False,
        shell_width: float = 0.1,
        soc: bool = False,
        full_pairs: bool = False,
        magnetic_axis: tuple[float, float, float] | None = None,
        rij_out: str | Path | None = None,
        rij_fmt: str = "json",
        verbose: bool = False,
        input_path: str | Path | None = None,
        lattice_path: str | Path | None = None,
        cluster_path: str | Path | None = None,
        report_path: str | Path | None = None,
        collision_policy: str = "error",
        force_legacy_collision: bool = False,
    ) -> ExchangeResult:
        """Generate legacy ``Atom_*`` jobs and an optional mapping artifact."""

        plan = self.plan(
            cutoff=cutoff,
            unique=unique,
            shell_width=shell_width,
            soc=soc,
            full_pairs=full_pairs,
            magnetic_axis=magnetic_axis,
            verbose=verbose,
            input_path=input_path,
            lattice_path=lattice_path,
            cluster_path=cluster_path,
            report_path=report_path,
        )
        mapping_path = None
        if rij_out is not None:
            mapping_path = _resolve_output(self.source_dir, rij_out)
            # Publish an explicitly requested mapping before job directories.
            # A write failure therefore cannot look like a successful job run.
            write_Rij_mappings(
                _legacy_mapping_records(plan),
                mapping_path,
                fmt=rij_fmt,
            )
            if verbose:
                print(f"Wrote R_ij mappings to {mapping_path} (format={rij_fmt})")
        write_result = write_exchange_jobs(
            plan.jobs,
            lattice_path=plan.reference.lattice_path,
            source_dir=self.source_dir,
            collision_policy=collision_policy,
            force_legacy_collision=force_legacy_collision,
        )
        return ExchangeResult(
            plan=plan,
            mapping_path=mapping_path,
            write_result=write_result,
        )


def plan_exchange_jobs(source_dir: str | Path = ".", **options) -> ExchangePlan:
    """Convenience API for pair/job inventory generation."""

    return ExchangeService(source_dir).plan(**options)


def generate_exchange_jobs(source_dir: str | Path = ".", **options) -> ExchangeResult:
    """Convenience API for path-explicit exchange generation."""

    return ExchangeService(source_dir).generate(**options)


run_exchange = generate_exchange_jobs


__all__ = [
    "ExchangePlan",
    "ExchangeResult",
    "ExchangeService",
    "generate_exchange_jobs",
    "plan_exchange_jobs",
    "run_exchange",
]
