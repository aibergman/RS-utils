"""Deterministic content and specification hashing helpers."""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
from enum import Enum
from pathlib import Path
from typing import Any, Mapping


_TEMPORAL_KEYS = frozenset(
    {
        "timestamp",
        "timestamps",
        "created_at",
        "updated_at",
        "started_at",
        "finished_at",
        "completed_at",
    }
)


def _json_compatible(value: Any, *, drop_temporal: bool = False) -> Any:
    """Convert common typed values to JSON-compatible values.

    Scientific identity helpers reject non-finite floats instead of allowing
    non-standard JSON values such as ``NaN`` into a hash.
    """

    if isinstance(value, Enum):
        return _json_compatible(value.value, drop_temporal=drop_temporal)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        value = dataclasses.asdict(value)
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            key_string = str(key)
            if drop_temporal and key_string.lower() in _TEMPORAL_KEYS:
                continue
            result[key_string] = _json_compatible(item, drop_temporal=drop_temporal)
        return result
    if isinstance(value, (list, tuple)):
        return [_json_compatible(item, drop_temporal=drop_temporal) for item in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("scientific specifications cannot contain non-finite floats")
        return value
    if isinstance(value, Path):
        return value.as_posix()
    raise TypeError(f"value of type {type(value).__name__} is not JSON-compatible")


def canonical_json(value: Any, *, scientific: bool = False) -> str:
    """Serialize a value deterministically for storage or hashing.

    Mapping insertion order does not affect the result.  Set-like values are
    intentionally unsupported because their order is not a scientific
    contract; callers should provide an ordered sequence instead.
    """

    normalized = _json_compatible(value, drop_temporal=scientific)
    return json.dumps(
        normalized,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def hash_json(value: Any, *, algorithm: str = "sha256") -> str:
    """Return a hash of canonical JSON, including all supplied fields."""

    digest = hashlib.new(algorithm)
    digest.update(canonical_json(value).encode("utf-8"))
    return digest.hexdigest()


def scientific_hash(value: Any, *, algorithm: str = "sha256") -> str:
    """Return a deterministic identity hash for scientific specifications.

    Provenance timestamps are excluded recursively.  They belong in run
    records, not in scientific identity.
    """

    digest = hashlib.new(algorithm)
    digest.update(canonical_json(value, scientific=True).encode("utf-8"))
    return digest.hexdigest()


def hash_file(path: str | Path, *, algorithm: str = "sha256", chunk_size: int = 1024 * 1024) -> str:
    """Hash file content without incorporating path metadata or timestamps."""

    file_path = Path(path)
    digest = hashlib.new(algorithm)
    with file_path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


file_hash = hash_file
specification_hash = scientific_hash
