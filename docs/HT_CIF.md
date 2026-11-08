# CIF to RS-LMTO source files

`rsutil setup cif` converts an ordered CIF into the `lattice.nml` and
`input.nml` pair that `rsutil setup bulk` consumes. It replaces the external
`cif2cell` step. Parsing and space-group expansion use pymatgen, which is
already a core dependency.

```bash
rsutil setup cif Ni2MnGa.cif --ct 4.3 -o runs/Ni2MnGa
```

```python
from rsutil.cif import convert_cif

files, conversion = convert_cif("Ni2MnGa.cif", "runs/Ni2MnGa", ct=4.3)
```

## Conventions

- `alat` is in Å and defaults to the length of the first lattice vector.
  `a(:, i)` and `crd(:, j)` are written in units of `alat`.
- Symmetry is not used. Every atom of the cell is its own electronic type:
  `izp = no = iu = ib = irec = atom index`, and
  `ntot = nbas = nrec = nbulk_bulk = ntype = number of atoms`.
- Atoms are ordered by descending atomic number, then by CIF site label, and
  labelled `<Element><running count>` (`Ga1`, ..., `Ni1`, ..., `Mn1`).
- Fractional coordinates are wrapped into `[0, 1)`. pymatgen rounds values
  such as `0.3333` to `1/3`; this matches `cif2cell` and is reported as a
  warning.
- `--supercell N1 N2 N3` replicates the cell as an ordered supercell;
  `alat` keeps the unreplicated cell length. `--primitive` reduces to the
  primitive cell.

## What is *not* derived from the CIF

`ct` has no safe structure-derived default and must be supplied; `r2`
defaults to `ct**2`. `rc`, `wav`, and the `&calculation`, `&self`,
`&energy`, `&control`, `&mix`, and `&hamiltonian` groups are settings to
review. Pass `--template input.nml` to reuse the non-structural groups from an
existing input instead of the built-in defaults.

## Limitations

- Partially occupied (disordered) CIFs are rejected with the offending sites
  listed. SQS or ordered-supercell handling of disorder is not implemented.
- Existing `lattice.nml`/`input.nml` are never replaced unless `--overwrite`
  is given.

Validated against `cif2cell 2.1.0` output for the pre-martensitic Ni2MnGa
structure (ICSD 190387, 24 atoms) in `tests/data/cif`.
