"""Exchange job specifications and legacy-compatible serialization."""

from __future__ import annotations

import copy
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import f90nml

from rsutil.core import ExchangePair, PairKey, RSStructure

from .pairs import map_cluster_neighbor


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

    @property
    def pair_keys(self) -> tuple[PairKey, ...]:
        return tuple(pair.key for pair in self.pairs)


def construct_exchange_jobs(
    input_nml: Mapping[str, Any],
    cluster: Sequence[Any],
    neighbors: Mapping[int, Sequence[Any]],
    recursion_sites: Sequence[int],
    canonical_structure: RSStructure,
    *,
    fermi_energy: float,
    output_root: str | Path,
) -> tuple[ExchangeJobSpec, ...]:
    """Construct canonical pair-backed jobs while preserving legacy indices."""

    root = Path(output_root).expanduser().resolve()
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


def write_exchange_jobs(
    jobs: Sequence[ExchangeJobSpec],
    *,
    lattice_path: str | Path,
    source_dir: str | Path,
) -> None:
    """Write familiar ``Atom_*`` directories from explicit source paths."""

    root = Path(source_dir).expanduser().resolve()
    lattice = Path(lattice_path).expanduser()
    if not lattice.is_absolute():
        lattice = root / lattice
    output_files = sorted(root.glob("*_out.nml"))
    for job in jobs:
        job.output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(lattice, job.output_dir / "lattice.nml")
        for output_file in output_files:
            shutil.copy2(output_file, job.output_dir / output_file.name)
        f90nml.write(job.input_nml, str(job.output_dir / "input.nml"), force=True)
        print(f"Updated input_nml saved to {job.output_dir / 'input.nml'}")
        print(f"lattice.nml copied to {job.output_dir / 'lattice.nml'}")


__all__ = [
    "ExchangeJobSpec",
    "build_exchange_input",
    "construct_exchange_jobs",
    "inventory_pairs",
    "write_exchange_jobs",
]
