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

## Legacy SOC convention

The Boolean `--soc` option does not claim generic noncollinear magnetic-space-
group support. It denotes the historical collinear reference state with a
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
