# Exchange result parsing

RS-utils represents magnetic interactions as an `ExchangeSet` of
`InteractionRecord` objects keyed by the canonical
`PairKey(i, j, R)`.  `i` and `j` are zero-based reference-cell sites and `R`
is an integer translation in the reference lattice basis.  The output table's
`itype` and `jtype` columns are one-based and are converted only at this
boundary.

## Supported table evidence

The checked-in RS-LMTO fixture `examples/export/MnSi` contains the following
stable, headerless tables.  The column names below are also the names used in
the checked-in format note `docs/notes/RSLMTO_format.md`.

| File family | Columns, in order | Parsed field | Unit |
| --- | --- | --- | --- |
| `jijso.out`, `jijfo.out` | `itype jtype r_ij[3] j_ij |r_ij|` | isotropic exchange `J` | Ry |
| `dijso.out`, `dijfo.out`, `dmi*.out` | `itype jtype r_ij[3] d_ij[3] |r_ij|` | DMI vector `D` | Ry |
| explicit tensor file/header | `itype jtype r_ij[3] J_ij[3x3] [|r_ij|]` | full exchange tensor | Ry |

The first five fields after the two indices are Cartesian displacement
coordinates in the same scale as the supplied `RSStructure`; `|r_ij|` is a
redundant distance check.  The result values are retained in the source
energy unit established by the RS-LMTO fixture convention (Ry); no conversion
is performed.  Tensor-shaped rows are accepted only when the filename or
header explicitly identifies a tensor, so a widened scalar table cannot be
silently reinterpreted.

The suffix/channel (`so` or `fo`) is retained as metadata.  SOC is retained
only when a header explicitly states it; it is not inferred from a filename.

## Pair mapping and validation

Rows are mapped using the one-based indices plus Cartesian displacement against
the canonical pair list in `rsutil_exchange_manifest.json` and the supplied
reference structure.  A missing or mismatched/ambiguous mapping is fatal.  If
no manifest is supplied, a reference structure is still required so the
integer translation can be reconstructed from the displacement.  Row order
and `Atom_*` directory names are never scientific identity.

Parser failures are raised as result-parser errors.  Physics validation is
stored separately in `ExchangeSet.validation`: when both directions are
present, isotropic exchange is checked as symmetric, DMI as antisymmetric,
and a full tensor against the transpose relation.  These are warnings by
default and can be made errors with `strict_physics=True`.  No inversion or
crystal-symmetry constraint is applied without explicit structure/symmetry
proof.

## Export

`ExchangeSet.to_json()` and `write_exchange_json()` emit
`rsutil-exchange-result-v1`.  Every populated magnetic quantity contains its
value and unit, and the export retains source hashes, source rows, format
metadata, the canonical pair key, and optional derived geometry.

## Evidence boundary

This RS-utils checkout does not contain the requested `fable_v4/source/exchange.f90`
or a fable_v4 tensor fixture.  The supported scalar and DMI schemas above are
therefore documented from the checked-in real MnSi outputs and format note;
the tensor grammar is intentionally explicit and fail-closed pending an
upstream fable_v4 tensor fixture.  No fable_v4 sign convention beyond the
fixture's Ry-valued columns is inferred here.
