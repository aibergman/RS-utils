# Transactional exchange jobs

Exchange generation keeps the legacy `Atom_<site>` directory layout, but the
directory is now a generated serialization of a canonical exchange job rather
than its scientific identity.

Each generated directory contains:

- `lattice.nml` and `input.nml`;
- the explicitly required `<label>_out.nml` files from the selected source run;
- `rsutil_exchange_manifest.json`.

The manifest uses schema `rsutil-exchange-v3` and records the deterministic job
specification hash, source-run hash, cutoff/reduction/SOC policy, Fermi-energy
source, one-based legacy `ijpair` values, and the canonical `(i, j, R)` pair
list. Each v3 pair record also carries endpoint `rs_type_i`/`rs_type_j` and
Cartesian displacement evidence. `R` is an integer reference-lattice
translation, not a Cartesian vector.
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
| `reuse_if_identical` | Stage the current job afresh and reuse an owned directory only when its specification hash, expected immutable file set, and fresh candidate hashes match exactly. |

Unmarked legacy or user directories are never deleted by default. The explicit
`force_legacy_collision=True` API option, or the legacy CLI option
`--force-legacy-collision`, is required to replace one. The legacy `-J` path
uses `replace_generated`, so clean `Atom_*` creation remains familiar and
previously generated rsutil jobs can be regenerated safely.

The manifest distinguishes primitive-cell site identity from RS electronic
type identity. The RS input's `izp(site) -> atoms.label(type)` mapping remains
authoritative, including when multiple sites share one type or chemically
identical sites have distinct types.

All job directories are populated and semantically checked in a same-filesystem
staging area before any destination is changed. That fresh staged candidate is
the authority for `reuse_if_identical`; the old directory manifest is only an
ownership/specification marker, so changing both an existing scientific file
and its historical manifest hash cannot make reuse succeed. Runtime outputs
(`*.out`, `*.log`, and `rsutil_execution.json`) are allowed to accumulate
without changing generated-job identity. A staging or validation error
therefore leaves existing job directories untouched; a commit error rolls back
already-replaced directories where the filesystem permits.
