"""Typed contracts for the lightweight campaign layer."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from ..core.contracts import SetupSpec
from ..core.hashing import scientific_hash


CAMPAIGN_SCHEMA_VERSION = "rsutil-campaign-v1"
CAMPAIGN_INDEX_SCHEMA_VERSION = "rsutil-campaign-index-v1"


def _mapping(value: Mapping[str, Any] | None) -> dict[str, Any]:
    return dict(value or {})


@dataclass
class CampaignEntry:
    """One independently generated setup run.

    ``source`` is either a directory containing ``lattice.nml`` and
    ``input.nml``, a named structure key, or a mapping with explicit
    ``lattice`` and ``input`` files.  ``setup_spec`` is the normalized
    scientific setup contract; ``setup`` is accepted as a JSON spelling for
    convenient construction and is normalized into the same field.  Entry
    ``metadata`` is intentionally part of scientific identity; campaign-wide
    provenance belongs on ``Campaign.metadata``.
    """

    name: str
    source: str | Mapping[str, Any] | None = None
    setup_spec: SetupSpec | Mapping[str, Any] | None = None
    calculation: Mapping[str, Any] = field(default_factory=dict)
    exchange: Mapping[str, Any] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)
    destination: str | None = None
    structure: str | None = None
    seed: int | None = None
    setup: SetupSpec | Mapping[str, Any] | None = None
    entry_id: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("campaign entry name must be non-empty")
        if self.source is None and self.structure is None:
            raise ValueError(f"campaign entry {self.name!r} needs source or structure")
        if self.source is not None and self.structure is not None:
            raise ValueError(f"campaign entry {self.name!r} cannot set both source and structure")
        supplied_setup = self.setup_spec if self.setup_spec is not None else self.setup
        self.setup_spec = supplied_setup if isinstance(supplied_setup, SetupSpec) else SetupSpec.from_dict(supplied_setup or {})
        if self.seed is not None:
            if self.setup_spec.seed is not None and self.setup_spec.seed != self.seed:
                raise ValueError(f"campaign entry {self.name!r} has conflicting seeds")
            self.setup_spec = SetupSpec.from_dict({**self.setup_spec.to_dict(), "seed": self.seed})
        self.setup = self.setup_spec
        self.calculation = _mapping(self.calculation)
        self.exchange = None if self.exchange is None else _mapping(self.exchange)
        self.metadata = _mapping(self.metadata)
        if self.destination is not None:
            if not isinstance(self.destination, str) or not self.destination.strip():
                raise ValueError("entry destination must be a non-empty string")

    def to_dict(self, *, include_id: bool = True) -> dict[str, Any]:
        result: dict[str, Any] = {
            "name": self.name,
            "setup": self.setup_spec.to_dict(),
            "calculation": dict(self.calculation),
            "metadata": dict(self.metadata),
        }
        if self.source is not None:
            result["source"] = self.source if isinstance(self.source, str) else dict(self.source)
        if self.structure is not None:
            result["structure"] = self.structure
        if self.exchange is not None:
            result["exchange"] = dict(self.exchange)
        if self.destination is not None:
            result["destination"] = self.destination
        if self.seed is not None:
            result["seed"] = self.seed
        if include_id and self.entry_id is not None:
            result["entry_id"] = self.entry_id
        return result

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CampaignEntry":
        known = dict(data)
        return cls(
            name=str(known["name"]),
            source=known.get("source"),
            structure=known.get("structure"),
            setup_spec=known.get("setup_spec", known.get("setup")),
            calculation=known.get("calculation", known.get("run", {})),
            exchange=known.get("exchange"),
            metadata=known.get("metadata"),
            destination=known.get("destination"),
            seed=known.get("seed"),
            entry_id=known.get("entry_id"),
        )

    def identity_payload(
        self,
        *,
        source_hashes: Mapping[str, Any],
        setup_template: Mapping[str, Any] | None = None,
        calculation_template: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Return fields participating in the stable scientific entry ID.

        Entry metadata remains identity-bearing by design.  Non-scientific
        scaffold provenance, such as a historical CIF snapshot, must be
        stored in campaign-level metadata instead.
        """

        setup = dict(setup_template or {})
        setup.update(self.setup_spec.to_dict())
        calculation = dict(calculation_template or {})
        calculation.update(self.calculation)
        return {
            "schema_version": CAMPAIGN_SCHEMA_VERSION,
            "source_hashes": dict(source_hashes),
            "setup": setup,
            "calculation": calculation,
            "exchange": None if self.exchange is None else dict(self.exchange),
            "metadata": dict(self.metadata),
        }

    def identity_hash(self, *, source_hashes: Mapping[str, Any], setup_template: Mapping[str, Any] | None = None, calculation_template: Mapping[str, Any] | None = None) -> str:
        return scientific_hash(self.identity_payload(
            source_hashes=source_hashes,
            setup_template=setup_template,
            calculation_template=calculation_template,
        ))


@dataclass
class Campaign:
    """Versioned campaign specification.

    ``metadata`` is campaign-level metadata and is included in the campaign
    hash.  Scaffolded CIF snapshots live here under ``structure_sources``;
    they do not participate in any entry's scientific identity.
    """

    name: str
    entries: tuple[CampaignEntry, ...]
    structures: Mapping[str, Any] = field(default_factory=dict)
    setup_template: Mapping[str, Any] = field(default_factory=dict)
    calculation_template: Mapping[str, Any] = field(default_factory=dict)
    destination_naming: str = "runs/{name}--{entry_id}"
    destination: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)
    schema_version: str = CAMPAIGN_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != CAMPAIGN_SCHEMA_VERSION:
            raise ValueError(f"unsupported campaign schema: {self.schema_version!r}")
        if not self.name.strip():
            raise ValueError("campaign name must be non-empty")
        self.entries = tuple(item if isinstance(item, CampaignEntry) else CampaignEntry.from_dict(item) for item in self.entries)
        if not self.entries:
            raise ValueError("campaign must contain at least one entry")
        names = [entry.name for entry in self.entries]
        if len(names) != len(set(names)):
            raise ValueError("campaign entry names must be unique")
        self.structures = _mapping(self.structures)
        self.setup_template = _mapping(self.setup_template)
        self.calculation_template = _mapping(self.calculation_template)
        if not isinstance(self.destination_naming, str) or not self.destination_naming.strip():
            raise ValueError("destination_naming must be a non-empty string")
        self.metadata = _mapping(self.metadata)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Campaign":
        setup_template = data.get("setup_template", data.get("setup", {})) or {}
        entries = []
        for item in data.get("entries", ()):
            raw = dict(item)
            entry_setup = raw.get("setup_spec", raw.get("setup", {})) or {}
            raw["setup_spec"] = {**dict(setup_template), **dict(entry_setup)}
            entries.append(CampaignEntry.from_dict(raw))
        return cls(
            name=str(data.get("name", "campaign")),
            entries=tuple(entries),
            structures=data.get("structures", {}),
            setup_template=setup_template,
            calculation_template=data.get("calculation_template", {}),
            destination_naming=str(data.get("destination_naming", data.get("naming_policy", "runs/{name}--{entry_id}"))),
            destination=data.get("destination"),
            metadata=data.get("metadata"),
            schema_version=str(data.get("schema_version", CAMPAIGN_SCHEMA_VERSION)),
        )

    def to_dict(self) -> dict[str, Any]:
        result = {
            "schema_version": self.schema_version,
            "name": self.name,
            "structures": dict(self.structures),
            "setup_template": dict(self.setup_template),
            "calculation_template": dict(self.calculation_template),
            "destination_naming": self.destination_naming,
            "entries": [entry.to_dict(include_id=False) for entry in self.entries],
            "metadata": dict(self.metadata),
        }
        if self.destination is not None:
            result["destination"] = self.destination
        return result

    def campaign_hash(self, entry_ids: list[str]) -> str:
        return scientific_hash({
            "schema_version": self.schema_version,
            "entries": sorted(entry_ids),
            "metadata": dict(self.metadata),
        })


@dataclass(frozen=True)
class CampaignIndexEntry:
    entry_id: str
    name: str
    run_path: str
    manifest_path: str
    manifest_hash: str
    status: str = "generated"

    def to_dict(self) -> dict[str, Any]:
        return {
            "entry_id": self.entry_id,
            "name": self.name,
            "run_path": self.run_path,
            "manifest_path": self.manifest_path,
            "manifest_hash": self.manifest_hash,
            "status": self.status,
        }


@dataclass(frozen=True)
class CampaignIndex:
    """Relocatable top-level index for generated campaign runs."""

    campaign_id: str
    campaign_name: str
    campaign_spec_hash: str
    entries: tuple[CampaignIndexEntry, ...]
    schema_version: str = CAMPAIGN_INDEX_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "campaign_id": self.campaign_id,
            "campaign_name": self.campaign_name,
            "campaign_spec_hash": self.campaign_spec_hash,
            "entries": [entry.to_dict() for entry in self.entries],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CampaignIndex":
        if data.get("schema_version") != CAMPAIGN_INDEX_SCHEMA_VERSION:
            raise ValueError(f"unsupported campaign index schema: {data.get('schema_version')!r}")
        return cls(
            campaign_id=str(data["campaign_id"]),
            campaign_name=str(data["campaign_name"]),
            campaign_spec_hash=str(data["campaign_spec_hash"]),
            entries=tuple(CampaignIndexEntry(**dict(item)) for item in data.get("entries", ())),
        )


__all__ = [
    "CAMPAIGN_INDEX_SCHEMA_VERSION",
    "CAMPAIGN_SCHEMA_VERSION",
    "Campaign",
    "CampaignEntry",
    "CampaignIndex",
    "CampaignIndexEntry",
]
