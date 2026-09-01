"""Private subprocess entry point for unsafe legacy setup execution."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the isolated RS-utils Fortran setup kernel")
    parser.add_argument("--request", required=True, type=Path, help="JSON request file")
    parser.add_argument("--workdir", required=True, type=Path, help="Isolated kernel working directory")
    return parser


def _read_parameters(request_path: Path) -> dict[str, Any]:
    with request_path.open("r", encoding="utf-8") as handle:
        request = json.load(handle)
    if not isinstance(request, dict):
        raise ValueError("worker request must be a JSON object")
    parameters = request.get("parameters")
    if not isinstance(parameters, dict):
        raise ValueError("worker request must contain a parameters object")
    return parameters


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    workdir = args.workdir.resolve()
    request = args.request.resolve()
    if not workdir.is_dir():
        raise FileNotFoundError(f"worker working directory does not exist: {workdir}")

    os.chdir(workdir)
    parameters = _read_parameters(request)

    # Import after changing CWD so all legacy fixed-name I/O occurs in the
    # supplied workspace.  The import itself loads the same shared library as
    # the parent; no numerical implementation is duplicated here.
    from . import _process_lattice_unsafe

    _process_lattice_unsafe(**parameters)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
