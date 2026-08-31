"""Versioned, unit-aware result parsing for RS-LMTO runs."""

from .errors import (
    MalformedResultError,
    MissingResultFileError,
    ResultParserError,
    UnsupportedResultFormatError,
)
from .models import (
    ConvergenceState,
    ElectronicResult,
    ParseDiagnostic,
    ParseSeverity,
    ParseValidation,
    ParserProvenance,
    SCFResult,
    SiteMoment,
    SiteQuantity,
    SourceFile,
)
from .parser import parse_report, parse_run, parse_scf
from .units import (
    ELECTRONS,
    MU_B,
    Quantity,
    RY,
    RYDBERG,
    UNKNOWN_UNIT,
    Unit,
    scalar,
    vector,
)

__all__ = [
    "ConvergenceState",
    "ELECTRONS",
    "ElectronicResult",
    "MU_B",
    "MalformedResultError",
    "MissingResultFileError",
    "ParseDiagnostic",
    "ParseSeverity",
    "ParseValidation",
    "ParserProvenance",
    "Quantity",
    "RY",
    "RYDBERG",
    "ResultParserError",
    "SCFResult",
    "SiteMoment",
    "SiteQuantity",
    "SourceFile",
    "UNKNOWN_UNIT",
    "Unit",
    "UnsupportedResultFormatError",
    "parse_report",
    "parse_run",
    "parse_scf",
    "scalar",
    "vector",
]
