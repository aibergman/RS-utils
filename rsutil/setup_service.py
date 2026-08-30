"""Reusable, path-explicit orchestration for the legacy ``-S`` workflow.

The numerical setup kernel remains in :mod:`rsutil` and is invoked through
``FortranSetupRunner``.  This module owns the filesystem and transform policy
around that kernel; it never changes the caller's working directory.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .core.contracts import (
    Diagnostic,
    DiagnosticSeverity,
    ErrorCode,
    OperationStatus,
    SetupResult,
    SetupSpec,
    ValidationResult,
)
from .core.paths import Workspace
from .setup_runner import FortranSetupRunner
from .setup_validation import build_run_manifest, validate_setup_namelists, validate_setup_outputs


_COLLISION_POLICIES = frozenset({"error", "replace"})


def _diagnostic(code: ErrorCode, message: str, **details: Any) -> Diagnostic:
    return Diagnostic(DiagnosticSeverity.ERROR, code, message, details=details)


def parse_sqs_spec(sqs_spec: str) -> tuple[list[str], np.ndarray]:
    """Parse an SQS composition such as ``Fe=0.5,Ni=0.5``."""
    parts = [part.strip() for part in sqs_spec.split(",") if part.strip()]
    if not parts:
        raise ValueError("Empty --sqs specification")
    species: list[str] = []
    concentrations: list[float] = []
    for part in parts:
        if "=" not in part:
            raise ValueError(f"Invalid --sqs token '{part}', expected Label=fraction")
        label, value = (item.strip() for item in part.split("=", 1))
        if not label:
            raise ValueError(f"Invalid --sqs token '{part}', missing label")
        fraction = float(value)
        if fraction < 0:
            raise ValueError(f"Invalid --sqs fraction for {label}: {fraction} (must be >= 0)")
        species.append(label)
        concentrations.append(fraction)
    values = np.asarray(concentrations, dtype=float)
    total = values.sum()
    if total <= 0:
        raise ValueError("Sum of --sqs fractions must be > 0")
    return species, values / total


def _allocate_species_counts(concentrations: np.ndarray, n_sites: int) -> np.ndarray:
    raw = concentrations * n_sites
    counts = np.floor(raw).astype(int)
    remainder = n_sites - counts.sum()
    if remainder > 0:
        for index in np.argsort(raw - counts)[::-1][:remainder]:
            counts[index] += 1
    elif remainder < 0:
        for index in np.argsort(raw - counts)[: -remainder]:
            if counts[index] > 0:
                counts[index] -= 1
    return counts


def _build_shell_pairs(positions: np.ndarray, n_shells: int, distance_tol: float = 1e-5):
    pairs: list[tuple[int, int]] = []
    distances: list[float] = []
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            pairs.append((i, j))
            distances.append(float(np.linalg.norm(positions[i] - positions[j])))
    if not pairs:
        return []
    order = np.argsort(distances)
    shells = [[pairs[int(order[0])]]]
    reference = distances[int(order[0])]
    for raw_index in order[1:]:
        index = int(raw_index)
        distance = distances[index]
        if abs(distance - reference) <= distance_tol:
            shells[-1].append(pairs[index])
        else:
            if len(shells) >= n_shells:
                break
            shells.append([pairs[index]])
            reference = distance
    return shells[:n_shells]


def _sqs_objective(assignment: np.ndarray, n_species: int, shells, concentrations: np.ndarray) -> float:
    target = np.zeros((n_species, n_species), dtype=float)
    for a in range(n_species):
        for b in range(a, n_species):
            target[a, b] = concentrations[a] * concentrations[b] * (1.0 if a == b else 2.0)
    total = 0.0
    for shell in shells:
        observed = np.zeros_like(target)
        for i, j in shell:
            a, b = sorted((int(assignment[i]), int(assignment[j])))
            observed[a, b] += 1.0
        difference = observed / max(len(shell), 1) - target
        total += float(np.sum(difference * difference))
    return total


def build_sqs_assignment(
    positions: np.ndarray,
    species: list[str],
    concentrations: np.ndarray,
    n_shells: int = 3,
    n_steps: int = 5000,
    seed: int | None = None,
):
    """Build the historical seeded MC SQS-like assignment."""
    n_sites = positions.shape[0]
    counts = _allocate_species_counts(concentrations, n_sites)
    if counts.sum() != n_sites:
        raise RuntimeError("Internal error: SQS counts do not sum to n_sites")
    rng = np.random.default_rng(seed)
    assignment = np.concatenate(
        [np.full(count, index, dtype=int) for index, count in enumerate(counts)]
    )
    rng.shuffle(assignment)
    shells = _build_shell_pairs(positions, n_shells)
    if not shells:
        raise ValueError("Could not construct SQS pair shells from structure")
    best = assignment.copy()
    best_objective = _sqs_objective(best, len(species), shells, concentrations)
    current_objective = best_objective
    for _ in range(max(n_steps, 0)):
        i, j = (int(value) for value in rng.integers(0, n_sites, size=2))
        if i == j or assignment[i] == assignment[j]:
            continue
        assignment[i], assignment[j] = assignment[j], assignment[i]
        objective = _sqs_objective(assignment, len(species), shells, concentrations)
        if objective <= current_objective:
            current_objective = objective
            if objective < best_objective:
                best_objective = objective
                best = assignment.copy()
        else:
            assignment[i], assignment[j] = assignment[j], assignment[i]
    return best, counts, best_objective


def validate_setup_inputs(source_dir: str | Path = ".") -> tuple[Path, Path]:
    source = Path(source_dir).resolve()
    paths = (source / "lattice.nml", source / "input.nml")
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    return paths


def validate_setup_parameters(spec: SetupSpec) -> None:
    """Apply the legacy ``-S`` parameter constraints to a typed spec."""
    if spec.min_empty_radius >= spec.max_empty_radius:
        raise ValueError(
            f"--min_empty_radius ({spec.min_empty_radius}) must be less than "
            f"--max_empty_radius ({spec.max_empty_radius})"
        )
    if spec.sqs is not None:
        parse_sqs_spec(spec.sqs)
        if spec.supercell == (1, 1, 1):
            raise ValueError("--sqs currently requires --supercell")


def _set_lattice_scalar(text: str, key: str, value: str) -> str:
    block = re.search(r"&lattice\b(.*?)/", text, flags=re.IGNORECASE | re.DOTALL)
    if not block:
        return text
    body = block.group(0)
    body, count = re.subn(
        rf"(^\s*{key}\s*=\s*).*$",
        rf"\g<1>{value}",
        body,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    if count == 0:
        body = re.sub(r"/\s*$", f"{key} = {value}\n/", body, flags=re.MULTILINE)
    return text[: block.start()] + body + text[block.end() :]


def _insert_lattice_lines(text: str, lines: list[str]) -> str:
    block = re.search(r"&lattice\b(.*?)/", text, flags=re.IGNORECASE | re.DOTALL)
    if not block:
        return text
    body = re.sub(r"/\s*$", "".join(lines) + "/", block.group(0), flags=re.MULTILINE)
    return text[: block.start()] + body + text[block.end() :]


def _replace_lines(text: str, first_pattern: str, line_pattern: str, lines: list[str]) -> str:
    first = re.search(first_pattern, text, flags=re.IGNORECASE)
    if not first:
        raise ValueError(f"could not find section matching {first_pattern!r}")
    start = first.start()
    section_lines = text[start:].splitlines(keepends=True)
    count = sum(1 for line in section_lines if re.search(line_pattern, line, flags=re.IGNORECASE))
    position = start
    for _ in range(count):
        newline = text.find("\n", position)
        position = len(text) if newline < 0 else newline + 1
    return text[:start] + "".join(lines) + text[position:]


def create_supercell(
    nx: int,
    ny: int,
    nz: int,
    lattice_path: str | Path = "lattice.nml",
    input_path: str | Path = "input.nml",
    output_lattice_path: str | Path = "lattice.nml",
    output_input_path: str | Path = "input.nml",
    in_place: bool = True,
    sqs_spec: str | None = None,
    sqs_shells: int = 3,
    sqs_steps: int = 5000,
    sqs_seed: int | None = None,
) -> bool:
    """Create the legacy supercell files without changing process CWD."""
    if min(nx, ny, nz) < 1:
        print("ERROR: Supercell dimensions must be positive integers", file=sys.stderr)
        return False
    lattice_file = Path(lattice_path).resolve()
    input_file = Path(input_path).resolve()
    output_lattice = Path(output_lattice_path).resolve()
    output_input = Path(output_input_path).resolve()
    try:
        content = lattice_file.read_text(encoding="utf-8")
        vectors = []
        for index in (1, 2, 3):
            match = re.search(
                rf"a\s*\(\s*:\s*,\s*{index}\s*\)\s*=\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)",
                content,
            )
            if not match:
                raise ValueError("Could not parse all 3 lattice vectors")
            vectors.append(np.asarray([float(value) for value in match.groups()]))
        positions = [
            np.asarray([float(value) for value in match.groups()[1:]])
            for match in re.finditer(
                r"crd\s*\(\s*:\s*,\s*(\d+)\s*\)\s*=\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)",
                content,
            )
        ]
        if not positions:
            raise ValueError("Could not parse atomic positions from lattice.nml")
        print(f"\nCreating {nx}×{ny}×{nz} supercell...")
        print(f"  Original unit cell: {len(positions)} atoms")
        expanded = [
            position + ix * vectors[0] + iy * vectors[1] + iz * vectors[2]
            for ix in range(nx)
            for iy in range(ny)
            for iz in range(nz)
            for position in positions
        ]
        n_atoms = len(expanded)
        print(f"  Supercell: {n_atoms} atoms ({nx}×{ny}×{nz} = {nx * ny * nz} unit cells)")
        for key in ("nbulk_bulk", "ntot", "nbas", "nrec"):
            content = re.sub(rf"({key}\s*=\s*)\d+", rf"\g<1>{n_atoms}", content, flags=re.IGNORECASE)
        for index, scale in enumerate((nx, ny, nz), start=1):
            content = re.sub(
                rf"(a\s*\(\s*:\s*,\s*{index}\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)",
                lambda match, scale=scale: match.group(1)
                + ", ".join(f"{float(value) * scale:.10f}" for value in match.groups()[1:]),
                content,
                flags=re.IGNORECASE,
            )
        content = _replace_lines(
            content,
            r"crd\s*\(\s*:\s*,\s*1\s*\)",
            r"crd\s*\(",
            [f"crd(:, {index}) = {position[0]:.10f}, {position[1]:.10f}, {position[2]:.10f}\n" for index, position in enumerate(expanded, 1)],
        )
        for name in ("izp", "no", "iu", "ib", "irec"):
            if re.search(rf"{name}\s*\(\s*1\s*\)", content, flags=re.IGNORECASE):
                content = _replace_lines(
                    content,
                    rf"{name}\s*\(\s*1\s*\)",
                    rf"{name}\s*\(",
                    [f"{name}({index}) = {index}\n" for index in range(1, n_atoms + 1)],
                )
        output_lattice.parent.mkdir(parents=True, exist_ok=True)
        if in_place:
            shutil.copy2(lattice_file, lattice_file.with_name("lattice_original.nml.bak"))
        output_lattice.write_text(content, encoding="utf-8")

        if input_file.is_file():
            input_content = input_file.read_text(encoding="utf-8")
            labels = {int(match.group(1)): match.group(2) for match in re.finditer(r"label\s*\(\s*(\d+)\s*\)\s*=\s*['\"]([^'\"]+)['\"]", input_content)}
            ct_values = {int(match.group(1)): match.group(2) for match in re.finditer(r"ct\s*\(\s*(\d+)\s*\)\s*=\s*([-\d.eE+]+)", input_content)}
            if labels:
                ordered = [labels[index] for index in sorted(labels) if index <= len(positions)]
                ordered_ct = [ct_values[index] for index in sorted(labels) if index in ct_values]
                source_labels: list[str] = []
                super_labels: list[str] = []
                super_ct: list[str] = []
                counts: dict[str, int] = {}
                if sqs_spec:
                    species, concentrations = parse_sqs_spec(sqs_spec)
                    assignment, _, objective = build_sqs_assignment(
                        np.asarray(expanded), species, concentrations, sqs_shells, sqs_steps, sqs_seed
                    )
                    print(f"  SQS enabled: {sqs_spec}")
                    print(f"  SQS shells={sqs_shells}, steps={sqs_steps}, seed={sqs_seed}")
                    print(f"  SQS objective={objective:.6e}")
                    for assignment_index in assignment:
                        base = species[int(assignment_index)]
                        counts[base] = counts.get(base, 0) + 1
                        super_labels.append(f"{base}{counts[base]}")
                        source_labels.append(base)
                        if ordered_ct:
                            super_ct.append(ordered_ct[0])
                else:
                    metadata = []
                    for label in ordered:
                        match = re.match(r"^(.*?)(\d+)$", label)
                        metadata.append((match.group(1), int(match.group(2)), label) if match else (label, 0, label))
                    for _ in range(nx * ny * nz):
                        for index, (base, start, source_label) in enumerate(metadata):
                            counts.setdefault(base, start - 1 if start else 0)
                            counts[base] += 1
                            super_labels.append(f"{base}{counts[base]}")
                            source_labels.append(source_label)
                            if ordered_ct:
                                super_ct.append(ordered_ct[index])
                input_content = _replace_lines(
                    input_content,
                    r"label\s*\(\s*1\s*\)",
                    r"label\s*\(",
                    [f"label({index}) = '{label}'\n" for index, label in enumerate(super_labels, 1)],
                )
                input_content = _set_lattice_scalar(input_content, "ntype", str(len(super_labels)))
                if super_ct:
                    ct_match = re.search(r"ct\s*\(\s*(?::|\d+)\s*\)", input_content, flags=re.IGNORECASE)
                    ct_lines = [f"ct({index}) = {value}\n" for index, value in enumerate(super_ct, 1)]
                    if ct_match:
                        input_content = _replace_lines(input_content, r"ct\s*\(\s*(?::|\d+)\s*\)", r"ct\s*\(", ct_lines)
                    else:
                        input_content = _insert_lattice_lines(input_content, ct_lines)
                output_input.parent.mkdir(parents=True, exist_ok=True)
                if in_place:
                    shutil.copy2(input_file, input_file.with_name("input_original.nml.bak"))
                output_input.write_text(input_content, encoding="utf-8")
                if in_place:
                    for new_label, source_label in zip(super_labels, source_labels):
                        source_file = input_file.parent / f"{source_label}.nml"
                        target_file = input_file.parent / f"{new_label}.nml"
                        if source_file.is_file() and not target_file.exists():
                            shutil.copy2(source_file, target_file)
        print("  ✓ Supercell creation complete")
        return True
    except Exception as exc:
        print(f"ERROR creating supercell: {exc}", file=sys.stderr)
        return False


class SetupService:
    """Execute setup from explicit source/destination directories.

    Candidate files are never published until the isolated worker has exited
    successfully and validated its outputs.  The default ``error`` collision
    policy protects existing destination files; callers that own the output
    directory may explicitly request ``replace``.
    """

    def __init__(self, runner: FortranSetupRunner | None = None) -> None:
        self.runner = runner or FortranSetupRunner()

    def run(
        self,
        source_dir: str | Path,
        destination_dir: str | Path,
        spec: SetupSpec | Mapping[str, Any] | None = None,
        *,
        in_place: bool = False,
        backup: bool = False,
        workspace: str | Path | None = None,
        collision_policy: str = "error",
    ) -> SetupResult:
        source = Path(source_dir).resolve()
        destination = Path(destination_dir).resolve()
        try:
            if collision_policy not in _COLLISION_POLICIES:
                raise ValueError(
                    f"collision_policy must be one of {sorted(_COLLISION_POLICIES)}"
                )
            if destination.exists() and not destination.is_dir():
                raise ValueError(f"destination is not a directory: {destination}")
            setup_spec = spec if isinstance(spec, SetupSpec) else SetupSpec.from_dict(spec or {})
            validate_setup_parameters(setup_spec)
            lattice, input_nml = validate_setup_inputs(source)
        except FileNotFoundError as exc:
            return SetupResult(status=OperationStatus.FAILED, diagnostics=(_diagnostic(ErrorCode.MISSING_INPUT, "required setup input is missing", paths=str(exc)),))
        except Exception as exc:
            code = ErrorCode.SETUP_SPEC_INVALID if isinstance(spec, Mapping) else ErrorCode.INVALID_INPUT
            return SetupResult(status=OperationStatus.FAILED, diagnostics=(_diagnostic(code, f"invalid setup input: {exc}"),))

        preflight_validation = validate_setup_namelists(
            lattice,
            input_nml,
            spec=setup_spec,
        )
        if not preflight_validation.valid:
            failed = SetupResult(
                status=OperationStatus.FAILED,
                diagnostics=tuple(
                    diagnostic
                    for diagnostic in preflight_validation.diagnostics
                    if diagnostic.severity is DiagnosticSeverity.ERROR
                ),
                semantic_validation=preflight_validation,
            )
            return self._with_manifest(failed, source, input_nml, setup_spec)

        transform_workspace: Workspace | None = None
        runner_workspace: Workspace | None = None
        try:
            run_lattice, run_input = lattice, input_nml
            if setup_spec.supercell != (1, 1, 1):
                transform_workspace = Workspace.create(prefix="rsutil-setup-transform-")
                run_lattice = transform_workspace.path / "lattice.nml"
                run_input = transform_workspace.path / "input.nml"
                if not create_supercell(*setup_spec.supercell, lattice, input_nml, run_lattice, run_input, False, setup_spec.sqs, setup_spec.sqs_shells, setup_spec.sqs_steps, setup_spec.sqs_seed):
                    raise ValueError("failed to create supercell")
            if workspace is None:
                runner_workspace = Workspace.create(prefix="rsutil-setup-run-")
                run_workspace = runner_workspace.path
            else:
                run_workspace = Path(workspace).resolve()
                run_workspace.mkdir(parents=True, exist_ok=True)
            result = self.runner.run(run_lattice, run_input, setup_spec, workspace=run_workspace)
            if not result.success:
                return self._with_manifest(result, source, input_nml, setup_spec)

            if result.semantic_validation is not None:
                semantic_validation = ValidationResult(
                    valid=preflight_validation.valid and result.semantic_validation.valid,
                    diagnostics=preflight_validation.diagnostics + result.semantic_validation.diagnostics,
                )
                result = replace(result, semantic_validation=semantic_validation)
                if not semantic_validation.valid:
                    return self._with_manifest(
                        replace(
                            result,
                            status=OperationStatus.FAILED,
                            diagnostics=tuple(
                                diagnostic
                                for diagnostic in semantic_validation.diagnostics
                                if diagnostic.severity is DiagnosticSeverity.ERROR
                            ),
                        ),
                        source,
                        input_nml,
                        setup_spec,
                    )

            if isinstance(self.runner, FortranSetupRunner) and setup_spec.supercell != (1, 1, 1):
                expansion_validation = validate_setup_outputs(
                    run_workspace,
                    lattice_path=lattice,
                    input_path=input_nml,
                    spec=setup_spec,
                    check_supercell_expansion=True,
                )
                expansion_errors = tuple(
                    diagnostic
                    for diagnostic in expansion_validation.diagnostics
                    if diagnostic.code == ErrorCode.SUPERCELL_INVALID.value
                    and diagnostic.severity is DiagnosticSeverity.ERROR
                )
                if expansion_errors:
                    return self._with_manifest(
                        replace(result, status=OperationStatus.FAILED, diagnostics=expansion_errors),
                        source,
                        input_nml,
                        setup_spec,
                    )

            committed_files = _commit_setup_outputs(
                run_workspace,
                destination,
                result.candidate_files,
                in_place=in_place,
                backup=backup or in_place,
                collision_policy=collision_policy,
            )
            provenance = dict(result.provenance)
            reproducible = setup_spec.seed is not None and (
                setup_spec.sqs is None or setup_spec.sqs_seed is not None
            )
            output_paths = {
                name: (
                    destination / "lattice.nml"
                    if in_place and name == "lattice_new.nml"
                    else destination / "input.nml"
                    if in_place and name == "input_new.nml"
                    else destination / name
                )
                for name in committed_files
            }
            manifest = build_run_manifest(
                lattice,
                input_nml,
                setup_spec,
                status=result.status,
                validation=result.semantic_validation,
                output_paths=output_paths,
                base_provenance=result.provenance,
            )
            provenance.update(
                {
                    "source_dir": str(source),
                    "destination_dir": str(destination),
                    "setup_spec": setup_spec.to_dict(),
                    "collision_policy": collision_policy,
                    "in_place": in_place,
                    "backup": backup or in_place,
                    "committed_files": list(committed_files),
                    "reproducible": reproducible,
                    "randomness": {
                        "seed": setup_spec.seed,
                        "sqs_seed": setup_spec.sqs_seed,
                        "reproducible": reproducible,
                    },
                    "manifest": manifest.to_dict(),
                }
            )
            result = replace(result, provenance=provenance)
            return result
        except _OutputCollisionError as exc:
            return SetupResult(
                status=OperationStatus.FAILED,
                diagnostics=(_diagnostic(ErrorCode.OUTPUT_COLLISION, str(exc)),),
            )
        except _CommitError as exc:
            return SetupResult(
                status=OperationStatus.FAILED,
                diagnostics=(_diagnostic(ErrorCode.COMMIT_FAILED, str(exc)),),
            )
        except Exception as exc:
            return SetupResult(status=OperationStatus.FAILED, diagnostics=(_diagnostic(ErrorCode.WORKSPACE_ERROR, f"setup service failed: {exc}"),))
        finally:
            if transform_workspace is not None:
                transform_workspace.cleanup()
            if runner_workspace is not None:
                runner_workspace.cleanup()

    @staticmethod
    def _with_manifest(
        result: SetupResult,
        source: Path,
        input_path: Path,
        spec: SetupSpec,
    ) -> SetupResult:
        """Attach the same machine-readable provenance to success and failure."""

        try:
            manifest = build_run_manifest(
                source / "lattice.nml",
                input_path,
                spec,
                status=result.status,
                validation=result.semantic_validation,
                base_provenance=result.provenance,
            )
            provenance = dict(result.provenance)
            provenance["manifest"] = manifest.to_dict()
            return replace(result, provenance=provenance)
        except Exception:
            # Provenance must never hide the original scientific/setup result.
            return result

    setup = run


class _OutputCollisionError(RuntimeError):
    """A destination contains a file that the selected policy will not replace."""


class _CommitError(RuntimeError):
    """Candidate publication failed after worker validation."""


def _safe_candidate_name(name: str) -> Path:
    candidate = Path(name)
    if (
        not name
        or candidate.is_absolute()
        or candidate == Path(".")
        or ".." in candidate.parts
    ):
        raise _CommitError(f"unsafe candidate output path: {name!r}")
    return candidate


def _commit_setup_outputs(
    runner_workspace: Path,
    destination: Path,
    candidate_files: tuple[str, ...],
    *,
    in_place: bool,
    backup: bool,
    collision_policy: str,
) -> tuple[str, ...]:
    """Publish validated worker outputs with all-or-nothing file replacement."""

    mappings: list[tuple[str, Path]] = []
    for name in sorted(candidate_files):
        relative = _safe_candidate_name(name)
        candidate = runner_workspace / relative
        if not candidate.is_file():
            raise _CommitError(f"validated candidate output is missing: {name}")
        if in_place and name == "lattice_new.nml":
            target = destination / "lattice.nml"
        elif in_place and name == "input_new.nml":
            target = destination / "input.nml"
        else:
            target = destination / relative
        mappings.append((name, target))

    targets = [target for _, target in mappings]
    if len(set(targets)) != len(targets):
        raise _CommitError("candidate outputs map to duplicate destination paths")

    # In-place input replacement is the operation explicitly requested by the
    # caller.  All other existing outputs still require an explicit replace
    # policy, which protects unrelated atom files and reports.
    allowed_existing = {destination / "lattice.nml", destination / "input.nml"} if in_place else set()
    collisions = [
        target.relative_to(destination).as_posix()
        for target in targets
        if target.exists() and target not in allowed_existing
    ]
    if collision_policy == "error" and collisions:
        raise _OutputCollisionError(
            "destination output collision; existing files require collision_policy='replace': "
            + ", ".join(sorted(collisions))
        )

    backup_mappings: list[tuple[str, Path]] = []
    if backup:
        for original_name in ("lattice.nml", "input.nml"):
            original = destination / original_name
            if original.is_file():
                backup_mappings.append((f"{original_name}.bak", original.with_name(original.name + ".bak")))

    all_mappings = mappings + backup_mappings
    if len({target for _, target in all_mappings}) != len(all_mappings):
        raise _CommitError("candidate outputs and backups map to duplicate destination paths")

    commit_workspace = Workspace.create(
        same_filesystem_as=destination,
        prefix=".rsutil-setup-commit-",
    )
    created_dirs: list[Path] = []
    snapshots: list[tuple[Path, Path | None]] = []
    staged: list[tuple[Path, Path]] = []
    try:
        payload_root = commit_workspace.path / "payload"
        rollback_root = commit_workspace.path / "rollback"
        for index, (name, target) in enumerate(all_mappings):
            if target.exists() and target.is_dir():
                raise _CommitError(f"destination target is a directory: {target}")
            if name in dict(mappings):
                source = runner_workspace / _safe_candidate_name(name)
            else:
                source = destination / name.removesuffix(".bak")
            payload = payload_root / str(index)
            payload.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, payload)
            staged.append((payload, target))
            snapshot = None
            if target.exists():
                snapshot = rollback_root / str(index)
                snapshot.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(target, snapshot)
            snapshots.append((target, snapshot))

        for payload, target in staged:
            missing: list[Path] = []
            parent = target.parent
            while not parent.exists():
                missing.append(parent)
                parent = parent.parent
            if missing:
                target.parent.mkdir(parents=True, exist_ok=True)
                created_dirs.extend(missing)
            os.replace(payload, target)
    except Exception as exc:
        rollback_error: Exception | None = None
        for target, snapshot in reversed(snapshots):
            try:
                if snapshot is not None and snapshot.exists():
                    os.replace(snapshot, target)
                elif target.exists():
                    target.unlink()
            except Exception as restore_exc:  # pragma: no cover - filesystem catastrophe
                rollback_error = restore_exc
        for directory in sorted(created_dirs, key=lambda path: len(path.parts), reverse=True):
            try:
                directory.rmdir()
            except OSError:
                pass
        if rollback_error is not None:
            raise _CommitError(f"commit failed ({exc}); rollback also failed ({rollback_error})") from exc
        if isinstance(exc, _CommitError):
            raise
        raise _CommitError(f"commit failed: {exc}") from exc
    finally:
        commit_workspace.cleanup()

    return tuple(name for name, _ in mappings)


def run_setup(source_dir: str | Path, destination_dir: str | Path, spec: SetupSpec | Mapping[str, Any] | None = None, **kwargs: Any) -> SetupResult:
    return SetupService().run(source_dir, destination_dir, spec, **kwargs)


__all__ = [
    "SetupService",
    "run_setup",
    "create_supercell",
    "parse_sqs_spec",
    "build_sqs_assignment",
    "validate_setup_inputs",
    "validate_setup_parameters",
]
