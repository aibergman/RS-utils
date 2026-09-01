"""Public setup APIs for legacy and high-throughput workflows."""

from .setup_runner import FortranSetupRunner
from .setup_service import SetupService, run_setup
from .setup_validation import (
    build_run_manifest,
    read_setup_report,
    setup_report,
    validate_setup,
    validate_setup_namelists,
    validate_setup_outputs,
    write_setup_report,
)

__all__ = [
    "FortranSetupRunner",
    "SetupService",
    "build_run_manifest",
    "read_setup_report",
    "run_setup",
    "setup_report",
    "validate_setup",
    "validate_setup_namelists",
    "validate_setup_outputs",
    "write_setup_report",
]
