"""Process-isolated adapter for the legacy Fortran setup kernel."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Mapping

from .core.contracts import (
    Diagnostic,
    DiagnosticSeverity,
    ErrorCode,
    OperationStatus,
    SetupResult,
    SetupSpec,
    ValidationResult,
)
from .core.paths import Workspace, stage_file
from .setup_validation import validate_setup_outputs


_WORKER_PARAMETERS = (
    "cluster_type",
    "empty_spheres",
    "min_empty_radius",
    "max_empty_radius",
    "lmx",
    "lattice_dist",
    "basis_dist",
    "alat_scale",
    "seed",
    "mode_2d",
    "z_iso_factor",
    "rec_steps",
    "rec_method",
)
_REQUIRED_OUTPUTS = ("lattice_new.nml", "input_new.nml", "radii.out")
_SPHERE_OUTPUTS = ("radii.out", "spheres.log")


def _parameters(spec: SetupSpec) -> dict[str, Any]:
    values = spec.to_dict()
    return {name: values[name] for name in _WORKER_PARAMETERS}


def _diagnostic(code: ErrorCode, message: str, **details: Any) -> Diagnostic:
    return Diagnostic(
        severity=DiagnosticSeverity.ERROR,
        code=code,
        message=message,
        details=details,
    )


class FortranSetupRunner:
    """Run the existing setup call in a child process and isolated CWD.

    ``lattice_path`` and ``input_path`` are copied into the workspace under
    the fixed names expected by Fortran.  By default the workspace is removed
    after the call; pass ``workspace=...`` when a caller needs to inspect or
    commit candidate files itself.  A caller-provided workspace is never
    removed by this class.
    """

    def __init__(
        self,
        *,
        python_executable: str | os.PathLike[str] | None = None,
        timeout: float | None = None,
        workspace_parent: str | os.PathLike[str] | None = None,
        keep_workspace: bool = False,
        environment: Mapping[str, str] | None = None,
    ) -> None:
        if timeout is not None and timeout <= 0:
            raise ValueError("timeout must be positive when supplied")
        self.python_executable = str(python_executable or sys.executable)
        self.timeout = timeout
        self.workspace_parent = Path(workspace_parent) if workspace_parent is not None else None
        self.keep_workspace = keep_workspace
        self.environment = dict(environment or {})

    def run(
        self,
        lattice_path: str | os.PathLike[str] = "lattice.nml",
        input_path: str | os.PathLike[str] = "input.nml",
        spec: SetupSpec | Mapping[str, Any] | None = None,
        *,
        workspace: str | os.PathLike[str] | None = None,
        resources: Iterable[str | os.PathLike[str]] = (),
        timeout: float | None = None,
    ) -> SetupResult:
        """Execute setup and return structured worker/output diagnostics.

        ``resources`` is an explicit extension point for future legacy
        resources.  The inspected setup call graph currently needs only the
        two namelists; resource files are staged by basename and never found
        implicitly by scanning the source directory.
        """

        started = time.time()
        try:
            setup_spec = spec if isinstance(spec, SetupSpec) else SetupSpec.from_dict(spec or {})
        except Exception as exc:
            return self._failed(
                _diagnostic(ErrorCode.INVALID_INPUT, f"invalid setup specification: {exc}"),
                started=started,
            )

        source_paths = (Path(lattice_path), Path(input_path))
        missing = [str(path) for path in source_paths if not path.is_file()]
        if missing:
            return self._failed(
                _diagnostic(ErrorCode.MISSING_INPUT, "required setup input is missing", paths=missing),
                started=started,
            )

        owned_workspace: Workspace | None = None
        if workspace is None:
            owned_workspace = Workspace.create(
                parent=self.workspace_parent,
                prefix="rsutil-fortran-",
                keep=self.keep_workspace,
            )
            workspace_path = owned_workspace.path
        else:
            workspace_path = Path(workspace).resolve()
            workspace_path.mkdir(parents=True, exist_ok=True)

        try:
            return self._run_in_workspace(
                workspace_path,
                source_paths,
                setup_spec,
                resources=resources,
                timeout=self.timeout if timeout is None else timeout,
                started=started,
            )
        except Exception as exc:
            return self._failed(
                _diagnostic(ErrorCode.WORKSPACE_ERROR, f"could not prepare setup workspace: {exc}"),
                started=started,
            )
        finally:
            if owned_workspace is not None:
                owned_workspace.cleanup()

    def _run_in_workspace(
        self,
        workspace: Path,
        source_paths: tuple[Path, Path],
        spec: SetupSpec,
        *,
        resources: Iterable[str | os.PathLike[str]],
        timeout: float | None,
        started: float,
    ) -> SetupResult:
        # A supplied workspace is an explicit isolation boundary.  Reject
        # pre-existing candidate names so stale outputs cannot be mistaken for
        # this run's results.
        stale = [name for name in (*_REQUIRED_OUTPUTS, "spheres.log") if (workspace / name).exists()]
        if stale:
            return self._failed(
                _diagnostic(ErrorCode.WORKSPACE_ERROR, "workspace contains stale candidate outputs", paths=stale),
                started=started,
            )

        staged_names = {"lattice.nml", "input.nml"}
        stage_file(source_paths[0], workspace, name="lattice.nml")
        stage_file(source_paths[1], workspace, name="input.nml")
        for resource in resources:
            staged_names.add(stage_file(resource, workspace).name)

        request_path = workspace / "rsutil_setup_request.json"
        request_path.write_text(
            json.dumps({"version": 1, "parameters": _parameters(spec)}, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        environment = os.environ.copy()
        environment.update(self.environment)
        package_root = str(Path(__file__).resolve().parent.parent)
        python_path = environment.get("PYTHONPATH")
        environment["PYTHONPATH"] = package_root if not python_path else os.pathsep.join((package_root, python_path))
        command = [
            self.python_executable,
            "-m",
            "rsutil.worker",
            "--request",
            str(request_path),
            "--workdir",
            str(workspace),
        ]
        try:
            completed = subprocess.run(
                command,
                cwd=workspace,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout,
            )
            stdout, stderr, returncode = completed.stdout, completed.stderr, completed.returncode
        except subprocess.TimeoutExpired as exc:
            stdout = self._text(exc.stdout)
            stderr = self._text(exc.stderr)
            return self._result(
                status=OperationStatus.FAILED,
                diagnostics=(_diagnostic(ErrorCode.WORKER_FAILED, "Fortran setup worker timed out", timeout=timeout),),
                workspace=workspace,
                stdout=stdout,
                stderr=stderr,
                started=started,
            )

        candidate_files = self._inventory(workspace, excluded={*staged_names, request_path.name})
        sphere_outputs = {
            name: name for name in _SPHERE_OUTPUTS if (workspace / name).is_file()
        }
        if returncode != 0:
            diagnostic = _diagnostic(
                ErrorCode.WORKER_FAILED,
                "Fortran setup worker failed",
                returncode=returncode,
                candidate_files=list(candidate_files),
            )
            return self._result(
                status=OperationStatus.FAILED,
                diagnostics=(diagnostic,),
                workspace=workspace,
                stdout=stdout,
                stderr=stderr,
                candidate_files=candidate_files,
                sphere_outputs=sphere_outputs,
                started=started,
            )

        validation = self._validate_outputs(workspace, source_paths, spec)
        diagnostics = () if validation.valid else tuple(
            Diagnostic(
                severity=diagnostic.severity,
                code=(
                    ErrorCode.OUTPUT_VALIDATION_FAILED.value
                    if diagnostic.code
                    in {
                        ErrorCode.OUTPUT_MISSING.value,
                        ErrorCode.OUTPUT_UNREADABLE.value,
                        ErrorCode.OUTPUT_NAMELIST_INVALID.value,
                        ErrorCode.SPHERE_OUTPUT_INVALID.value,
                    }
                    else diagnostic.code
                ),
                message=diagnostic.message,
                path=diagnostic.path,
                details=diagnostic.details,
            )
            for diagnostic in validation.diagnostics
        )
        return self._result(
            status=OperationStatus.SUCCESS if validation.valid else OperationStatus.FAILED,
            diagnostics=diagnostics,
            workspace=workspace,
            stdout=stdout,
            stderr=stderr,
            candidate_files=candidate_files,
            semantic_validation=validation,
            sphere_outputs=sphere_outputs,
            started=started,
        )

    @staticmethod
    def _inventory(workspace: Path, *, excluded: set[str]) -> tuple[str, ...]:
        return tuple(
            sorted(
                path.relative_to(workspace).as_posix()
                for path in workspace.rglob("*")
                if path.is_file() and path.relative_to(workspace).as_posix() not in excluded
            )
        )

    @staticmethod
    def _validate_outputs(
        workspace: Path,
        source_paths: tuple[Path, Path],
        spec: SetupSpec,
    ) -> ValidationResult:
        """Apply semantic validation to the complete candidate output set."""

        return validate_setup_outputs(
            workspace,
            lattice_path=source_paths[0],
            input_path=source_paths[1],
            spec=spec,
        )

    @staticmethod
    def _text(value: str | bytes | None) -> str:
        if value is None:
            return ""
        return value.decode(errors="replace") if isinstance(value, bytes) else value

    @staticmethod
    def _result(
        *,
        status: OperationStatus,
        diagnostics: tuple[Diagnostic, ...],
        workspace: Path,
        stdout: str = "",
        stderr: str = "",
        candidate_files: tuple[str, ...] = (),
        semantic_validation: ValidationResult | None = None,
        sphere_outputs: Mapping[str, str] | None = None,
        started: float,
    ) -> SetupResult:
        return SetupResult(
            status=status,
            diagnostics=diagnostics,
            stdout=stdout,
            stderr=stderr,
            candidate_files=candidate_files,
            semantic_validation=semantic_validation,
            sphere_outputs=dict(sphere_outputs or {}),
            provenance={
                "worker": "python -m rsutil.worker",
                "workspace": str(workspace),
                "elapsed_seconds": time.time() - started,
            },
        )

    def _failed(self, diagnostic: Diagnostic, *, started: float) -> SetupResult:
        return SetupResult(
            status=OperationStatus.FAILED,
            diagnostics=(diagnostic,),
            provenance={"worker": "python -m rsutil.worker", "elapsed_seconds": time.time() - started},
        )


__all__ = ["FortranSetupRunner"]
