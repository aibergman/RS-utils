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
    MILLI_RY = "mRy"
    ELECTRONVOLT = "eV"
    ELECTRON = "electron"
    MU_B = "mu_B"
    SIEMENS_PER_METER = "S/m"
    SIEMENS_PER_CENTIMETER = "S/cm"
    PER_SECOND = "1/s"
    DIMENSIONLESS = "dimensionless"
    UNKNOWN = "unknown"


RYDBERG = Unit.RY
RY = RYDBERG
MILLI_RYDBERG = Unit.MILLI_RY
MILLI_RY = MILLI_RYDBERG
ELECTRONVOLT = Unit.ELECTRONVOLT
ELECTRONS = Unit.ELECTRON
MU_B = Unit.MU_B
SIEMENS_PER_METER = Unit.SIEMENS_PER_METER
SIEMENS_PER_CENTIMETER = Unit.SIEMENS_PER_CENTIMETER
PER_SECOND = Unit.PER_SECOND
DIMENSIONLESS = Unit.DIMENSIONLESS
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

    def convert(self, unit: Unit | str) -> "Quantity":
        """Return an explicitly converted quantity.

        Parsed values retain their source unit.  Conversion is opt-in and is
        currently limited to the Ry/mRy pair used by interaction outputs.
        """

        target = unit if isinstance(unit, Unit) else Unit(str(unit))
        if target == self.unit:
            return self
        try:
            factor = _CONVERSION_FACTORS[(self.unit, target)]
        except KeyError as exc:
            raise ValueError(f"no conversion is defined from {self.unit.value} to {target.value}") from exc
        return Quantity(_scale_value(self.value, factor), target)

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


_CONVERSION_FACTORS: dict[tuple[Unit, Unit], float] = {
    (Unit.RY, Unit.MILLI_RY): 1000.0,
    (Unit.MILLI_RY, Unit.RY): 1.0e-3,
}


def _scale_value(value: Any, factor: float) -> Any:
    if isinstance(value, tuple):
        return tuple(_scale_value(item, factor) for item in value)
    if isinstance(value, list):
        return [_scale_value(item, factor) for item in value]
    return float(value) * factor


__all__ = [
    "ELECTRONS",
    "ELECTRONVOLT",
    "DIMENSIONLESS",
    "MU_B",
    "MILLI_RY",
    "MILLI_RYDBERG",
    "PER_SECOND",
    "Quantity",
    "RY",
    "RYDBERG",
    "SIEMENS_PER_CENTIMETER",
    "SIEMENS_PER_METER",
    "UNKNOWN_UNIT",
    "Unit",
    "scalar",
    "vector",
]
