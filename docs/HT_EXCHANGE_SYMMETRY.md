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

The reduction API accepts an explicit `MagneticSymmetryContext` and the
exchange service derives its spin/SOC fields from the reference input's
`control.nsp`. The default service context is `unknown`, not FM. Only an
explicit `collinear_fm` assertion is currently proven for reduction. Other
classes resolve to identity-only reduction, so no pair is eliminated on
unverified magnetic equivalence.

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
Therefore `nsp` encodes spin representation/SOC, not FM versus AFM order.
Missing or unsupported `nsp` is explicit unknown. When order is unknown,
`unique=True` retains all candidate pairs and emits a warning. An explicit
`collinear_fm` assertion is required before reduction; `nsp=2` SOC reduction
also requires an explicit global magnetization axis.
See [`control.f90`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/control.f90#L51-L66)
and the exchange routine's explicit pair-list loop in
[`exchange.f90`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/exchange.f90#L148-L166).

## Legacy SOC convention

The Boolean `--soc` option does not claim generic noncollinear magnetic-space-
group support. It requests SOC and must agree with source `control.nsp`; it
does not assert FM order. For the restricted SOC reduction domain, use an
explicit FM reference and axis (for example,
`--magnetic-reference collinear_fm --magnetic-axis 0 0 1`). The axis is
treated as an axial vector:

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
