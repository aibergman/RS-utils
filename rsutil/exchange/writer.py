"""Compatibility serializers for exchange mapping artifacts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence


def write_Rij_mappings(
    mappings: Sequence[Mapping[str, Any]] | Mapping[Any, Sequence[Mapping[str, Any]]],
    filename: str | Path,
    fmt: str = "json",
) -> None:
    """Write legacy Cartesian ``R_ij`` records to an explicit path."""

    if isinstance(mappings, Mapping):
        records: list[dict[str, Any]] = []
        for key in sorted(mappings):
            for entry in mappings[key]:
                records.append(
                    {
                        "i": int(key),
                        "j": int(entry.get("j")),
                        "R_ij": [float(value) for value in entry.get("R_ij", ())],
                        "R_ij_kind": entry.get("R_ij_kind", "cartesian_displacement"),
                        "pair_key": entry.get("pair_key"),
                        "distance": float(entry.get("distance", 0.0)),
                        "species_i": str(entry.get("species_i", "")),
                        "species_j": str(entry.get("species_j", entry.get("species", ""))),
                    }
                )
    else:
        records = [dict(record) for record in mappings]

    path = Path(filename).expanduser()
    format_name = fmt.lower()
    if format_name == "json":
        with path.open("w", encoding="utf-8") as handle:
            json.dump(records, handle, indent=2)
        return
    if format_name == "txt":
        with path.open("w", encoding="utf-8") as handle:
            handle.write("# R_ij is Cartesian displacement; canonical identity is a PairKey\n")
            handle.write("# i j R_x R_y R_z distance species_i species_j\n")
            line_format = "{:4d} {:4d} {:12.8f} {:12.8f} {:12.8f} {:12.8f} {:10s} {:10s}\n"
            for record in records:
                rx, ry, rz = record.get("R_ij", (0.0, 0.0, 0.0))
                handle.write(
                    line_format.format(
                        int(record.get("i")),
                        int(record.get("j")),
                        float(rx),
                        float(ry),
                        float(rz),
                        float(record.get("distance", 0.0)),
                        str(record.get("species_i", ""))[:10],
                        str(record.get("species_j", ""))[:10],
                    )
                )
        return
    raise ValueError(f"Unsupported format: {fmt}")


__all__ = ["write_Rij_mappings"]
