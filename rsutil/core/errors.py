"""Exceptions raised by the reusable RS-utils core."""

from __future__ import annotations


class RsutilError(Exception):
    """Base class for expected RS-utils operation errors."""


class ValidationError(RsutilError, ValueError):
    """An input or typed contract failed validation."""


class SetupError(RsutilError):
    """A setup operation could not be completed."""


class ParserError(RsutilError):
    """An input or result file could not be parsed unambiguously."""


class ExchangeError(RsutilError):
    """An exchange mapping or exchange operation failed."""
