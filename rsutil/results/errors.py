"""Errors raised by the versioned RS-LMTO result parsers."""

from __future__ import annotations

from rsutil.core.errors import ParserError


class ResultParserError(ParserError):
    """Base class for an unambiguous result parsing failure."""

    code = "result_parse_error"
    error_code = code


class MissingResultFileError(ResultParserError):
    """A required result file is not present."""

    code = "missing_result_file"
    error_code = code


class UnsupportedResultFormatError(ResultParserError):
    """The available files do not identify a supported result format."""

    code = "unsupported_result_format"
    error_code = code


class MalformedResultError(ResultParserError):
    """A recognized result record is malformed or ambiguous."""

    code = "malformed_result"
    error_code = code


class SourceUnitMismatchError(MalformedResultError):
    """A known source format contradicts its writer-defined unit contract."""

    code = "SOURCE_UNIT_MISMATCH"
    error_code = code


class InteractionMappingError(MalformedResultError):
    """An interaction row cannot be mapped to exactly one canonical pair."""

    code = "interaction_mapping_error"
    error_code = code


__all__ = [
    "MalformedResultError",
    "InteractionMappingError",
    "MissingResultFileError",
    "ResultParserError",
    "SourceUnitMismatchError",
    "UnsupportedResultFormatError",
]
