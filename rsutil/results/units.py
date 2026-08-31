"""Small, explicit unit vocabulary used by parsed result quantities.

The parser intentionally does not convert values.  A unit is attached to the
value exactly as established by the source contract, or set to ``unknown``
when the source does not establish one.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Sequence


class Unit(str, Enum):
    RY = "Ry"
    ELECTRON = "electron"
    MU_B = "mu_B"
    UNKNOWN = "unknown"


RYDBERG = Unit.RY
RY = RYDBERG
ELECTRONS = Unit.ELECTRON
MU_B = Unit.MU_B
UNKNOWN_UNIT = Unit.UNKNOWN


@dataclass(frozen=True)
class Quantity:
    """A scalar or finite JSON-compatible array with an explicit unit."""

    value: Any
    unit: Unit | str

    def __post_init__(self) -> None:
        unit = self.unit if isinstance(self.unit, Unit) else Unit(str(self.unit))
        object.__setattr__(self, "unit", unit)

    def to_dict(self) -> dict[str, Any]:
        return {"value": _json_value(self.value), "unit": self.unit.value}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Quantity":
        return cls(value=data["value"], unit=data["unit"])


def scalar(value: float, unit: Unit | str) -> Quantity:
    """Construct a scalar quantity without performing unit conversion."""

    return Quantity(float(value), unit)


def vector(values: Sequence[float], unit: Unit | str) -> Quantity:
    """Construct a three-component quantity with an explicit unit."""

    return Quantity(tuple(float(value) for value in values), unit)


def _json_value(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Quantity):
        return value.to_dict()
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return value


__all__ = [
    "ELECTRONS",
    "MU_B",
    "Quantity",
    "RY",
    "RYDBERG",
    "UNKNOWN_UNIT",
    "Unit",
    "scalar",
    "vector",
]
