# Transactional exchange jobs

Exchange generation keeps the legacy `Atom_<site>` directory layout, but the
directory is now a generated serialization of a canonical exchange job rather
than its scientific identity.

Each generated directory contains:

- `lattice.nml` and `input.nml`;
- the explicitly required `<label>_out.nml` files from the selected source run;
- `rsutil_exchange_manifest.json`.

The manifest uses schema `rsutil-exchange-v2` and records the deterministic job
specification hash, source-run hash, cutoff/reduction/SOC policy, Fermi-energy
source, one-based legacy `ijpair` values, and the canonical `(i, j, R)` pair
list. `R` is an integer reference-lattice translation, not a Cartesian vector.
The authoritative pair identity is `spec.pair_keys`, which is covered by
`spec_hash`. The top-level `canonical_pair_list` and `pair_keys` fields are
compatibility views only; readers reject a manifest if either differs from
`spec.pair_keys`. Unambiguous v1 manifests remain readable with the same
consistency checks.

The Python API accepts these collision policies:

| Policy | Behavior |
|---|---|
| `error` | Refuse any existing destination. |
| `replace_generated` | Replace only an existing rsutil-owned directory. The old directory is held as a rollback backup until the full commit succeeds. |
| `reuse_if_identical` | Reuse an owned directory only when its specification hash, compatibility views, expected file set, and current generated-file hashes match exactly. |

Unmarked legacy or user directories are never deleted by default. The explicit
`force_legacy_collision=True` API option, or the legacy CLI option
`--force-legacy-collision`, is required to replace one. The legacy `-J` path
uses `replace_generated`, so clean `Atom_*` creation remains familiar and
previously generated rsutil jobs can be regenerated safely.

All job directories are populated and semantically checked in a same-filesystem
staging area before any destination is changed. A staging or validation error
therefore leaves existing job directories untouched; a commit error rolls back
already-replaced directories where the filesystem permits.
