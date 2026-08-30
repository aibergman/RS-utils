"""Structured parsing of the Fermi energy from an RS-LMTO report."""

from __future__ import annotations

import re
from pathlib import Path

from rsutil.core.errors import ParserError


class FermiEnergyParseError(ParserError):
    """The report does not contain one valid Fermi-energy record."""

    code = "parse_error"
    error_code = code


_FERMI_PATTERN = re.compile(r"^\s*Fermi\s+energy\s*:\s*(?P<value>[^\s]+)", re.IGNORECASE)


def parse_fermi_energy(report_path: str | Path) -> float:
    """Read one finite Fermi energy from an explicit report path.

    The parser intentionally accepts the historical ``Fermi energy: value``
    record only. Missing, malformed, or non-finite values are reported as a
    structured parser exception instead of being silently ignored.
    """

    path = Path(report_path)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise FermiEnergyParseError(
            f"could not read Fermi-energy report {path}: {exc}"
        ) from exc

    matches = [match for line in lines if (match := _FERMI_PATTERN.match(line))]
    if not matches:
        raise FermiEnergyParseError(
            f"Fermi energy not found in report file: {path}"
        )
    if len(matches) > 1:
        raise FermiEnergyParseError(
            f"multiple Fermi-energy records found in report file: {path}"
        )
    try:
        value = float(matches[0].group("value"))
    except ValueError as exc:
        raise FermiEnergyParseError(
            f"malformed Fermi energy in report file: {path}"
        ) from exc
    if not value == value or value in (float("inf"), float("-inf")):
        raise FermiEnergyParseError(
            f"Fermi energy is not finite in report file: {path}"
        )
    return value


# Descriptive alias for callers that use a noun-based parser name.
read_fermi_energy = parse_fermi_energy


__all__ = ["FermiEnergyParseError", "parse_fermi_energy", "read_fermi_energy"]
