"""Structured parsing of the Fermi energy from an RS-LMTO report."""

from __future__ import annotations

from pathlib import Path

from rsutil.core.errors import ParserError
from rsutil.results import ResultParserError, parse_run


class FermiEnergyParseError(ResultParserError):
    """The report does not contain one valid Fermi-energy record."""

    code = "parse_error"
    error_code = code


def parse_fermi_energy(report_path: str | Path) -> float:
    """Read one finite Fermi energy from an explicit report path.

    The parser intentionally accepts the historical ``Fermi energy: value``
    record only. Missing, malformed, or non-finite values are reported as a
    structured parser exception instead of being silently ignored.
    """

    try:
        result = parse_run(report_path)
    except (ParserError, OSError, ValueError) as exc:
        raise FermiEnergyParseError(str(exc)) from exc
    if result.fermi_energy is None:
        raise FermiEnergyParseError(f"Fermi energy not found in report file: {report_path}")
    return float(result.fermi_energy.value)


# Descriptive alias for callers that use a noun-based parser name.
read_fermi_energy = parse_fermi_energy


__all__ = ["FermiEnergyParseError", "parse_fermi_energy", "read_fermi_energy"]
