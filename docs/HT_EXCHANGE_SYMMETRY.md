# HT exchange symmetry reduction

RS-utils uses the canonical periodic pair identity

```text
PairKey(i, j, R)
```

where `i` and `j` are reference-cell site indices and `R` is an integer
translation of the reference lattice. Cartesian displacement and distance are
derived values; shell labels are presentation metadata.

For a crystallographic operation in fractional coordinates,

```text
x' = W x + t
```

the site mapping is determined exactly (within the mapping tolerance) by

```text
W r[k] + t = r[p(k)] + n(k)
```

where `p(k)` is a reference-cell site and `n(k)` is an integer lattice shift.
The induced action on a pair is therefore

```text
(i, j, R) -> (p(i), p(j), n(j) + W R - n(i))
```

Only pair keys in the requested candidate set are reduced. A representative
is the deterministic lexicographically smallest member of each exact orbit.
This permits a subset of recursion centers to remain safe: an operation that
maps outside the requested candidate set cannot remove the only available
calculation.

## Supported magnetic domains

The reduction API now accepts an explicit `MagneticSymmetryContext` for
domains beyond the legacy default. It records the reference class, SOC flag,
and optional global axis. Only `collinear_fm` is currently proven for
reduction. Other classes are accepted as explicit input but resolve to
identity-only reduction, so no pair is eliminated on unverified magnetic
equivalence.

| Magnetic reference | SOC | Reduction status |
|---|---:|---|
| collinear FM | off | supported; exact chemical periodic operations |
| collinear FM + explicit global axis | on | supported/restricted; axial-axis-preserving operations only |
| AFM | off/on | conservative identity-only; site-resolved magnetic symmetry is not implemented |
| ferrimagnetic | off/on | conservative identity-only; site-resolved magnetic symmetry is not implemented |
| generic noncollinear | off/on | conservative identity-only |
| magnetic space group / time reversal | off/on | not implemented; time reversal is never assumed |
| unknown | off/on | conservative identity-only |

The RS-LMTO `fable_v4` control contract distinguishes collinear `nsp=1/2`
from noncollinear `nsp=3/4`, and identifies SOC in `nsp=2/4`; it does not
provide a machine-readable magnetic-space-group contract for this reducer.
See [`control.f90`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/control.f90#L51-L66)
and the exchange routine's explicit pair-list loop in
[`exchange.f90`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/exchange.f90#L148-L166).

## Legacy SOC convention

The Boolean `--soc` option does not claim generic noncollinear magnetic-space-
group support. It denotes the historical collinear-FM reference state with a
Cartesian magnetization axis of `+z`. The axis is treated as an axial vector:

```text
axis' = det(Q) Q axis
```

where `Q` is the Cartesian rotation corresponding to the fractional operation
`W`. An operation is used for reduction only if it preserves `+z`, preserves
the chemical structure/site species, and can be mapped to exact site and
translation identities. Time reversal is not assumed, so an operation that
reverses the magnetization is excluded. An uncertain mapping or magnetic
validity excludes the operation and leaves the pair unreduced.

Use `--full-pairs` (or `full_pairs=True` in the exchange service) to disable
reduction for validation against the complete candidate pair set. The legacy
`--shell-width` argument remains accepted but no longer defines scientific
equivalence.
