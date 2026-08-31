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
| `jijso.out`, `jijfo.out`, `jij.out` | `itype jtype r_ij[3] j_ij |r_ij|` | isotropic exchange `J` | mRy |
| `dijso.out`, `dijfo.out`, `dij.out` | `itype jtype r_ij[3] d_ij[3] |r_ij|` | DMI vector `D` | mRy |
| `aijso.out`, `aijfo.out`, `aij.out` | `itype jtype r_ij[3] A_ij[3x3] |r_ij|` | anisotropic component `A` | mRy |
| generic explicitly identified tensor file | `itype jtype r_ij[3] A_ij[3x3] [|r_ij|]` | anisotropic tensor only when the file contract says so | unknown unless an explicit unit header is present |

The first five fields after the two indices are Cartesian displacement
coordinates in the same scale as the supplied `RSStructure`; `|r_ij|` is a
redundant distance check. The current `fable_v4` writer multiplies scalar
exchange, DMI, and `aij` tensor values by `1.0d3/(4.0d0*pi)` before writing
the corresponding files, establishing the source output unit as mRy. See
[`exchange.f90` in `fable_v4`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/exchange.f90#L941-L1005)
for the scaling and [`jijso`/`dijso`/`aijso` writers](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/exchange.f90#L1038-L1071)
for the file mapping. Parsed values retain the numeric source value and
source unit; no conversion is performed during parsing. For a known
`fable_v4` filename, a contradictory `# unit:` header is rejected with
`SOURCE_UNIT_MISMATCH`; the source writer contract remains authoritative.

The three source observables are kept separate:

```text
J_iso          scalar isotropic exchange from jij*.out
D_vector       DMI vector from dij*.out
A_anisotropic  raw anisotropic tensor from aij*.out
J_full         derived complete tensor, only when J_iso, D_vector, and
               A_anisotropic all have compatible known energy units
```

For `fable_v4`, the derived tensor follows the source convention:

\[
\mathsf{J}^{full}_{ij} = J_{ij}\mathsf{I} +
\begin{pmatrix}
0 & D_z & -D_y\\
-D_z & 0 & D_x\\
D_y & -D_x & 0
\end{pmatrix} + \mathsf{A}_{ij}.
\]

The Python operation is `InteractionRecord.full_exchange_tensor()` (or its
raising counterpart `require_full_exchange_tensor()`). A raw `aij*.out` row
is not a full exchange tensor and missing components are not zero-filled.

The nine `aij` values are serialized by Fortran in column-major array order.
Thus a row containing `1 2 3 4 5 6 7 8 9` reconstructs as
`((1,4,7),(2,5,8),(3,6,9))`.

Tensor-shaped rows are accepted only when the filename or header explicitly
identifies a tensor, so a widened scalar table cannot be silently
reinterpreted. A generic tensor filename without a source-unit declaration
is retained with `unknown` units rather than being assigned mRy by shape or
filename similarity.

The suffix/channel (`so` or `fo`) is retained as metadata.  SOC is retained
only when a header explicitly states it; it is not inferred from a filename.

## Geometry units

`ExchangePair` carries `length_unit` alongside its Cartesian displacement and
distance. Exchange JSON and legacy `R_ij` mapping artifacts preserve this
field; legacy data with no recoverable physical scale is labeled `unknown` or
`lattice`, never assumed to be Angstrom.

## Pair mapping and validation

Rows are mapped using the one-based indices plus Cartesian displacement against
the authoritative `spec.pair_keys` list in `rsutil_exchange_manifest.json` and
the supplied reference structure.  The manifest specification hash and any
retained top-level compatibility pair views are checked before mapping. A
missing or mismatched/ambiguous mapping is fatal. If no manifest is supplied,
a reference structure is still required so the integer translation can be
reconstructed from the displacement. Row order and `Atom_*` directory names
are never scientific identity.

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

The raw tensor is serialized under `anisotropic_exchange_tensor`; the old
`exchange_tensor` input key remains readable as a deprecated compatibility
alias.

## Evidence boundary

The checked-in real MnSi fixture is headerless, so its source unit is supplied
by the audited `fable_v4` writer contract above. No raw fixture values were
changed. The result schema remains `rsutil-exchange-result-v1` because units
were already serialized explicitly per quantity; old serialized records that
say `Ry` remain `Ry` when loaded and are not silently migrated. Reparse the
raw source file to obtain the corrected mRy interpretation.
