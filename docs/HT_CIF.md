# CIF to RS-LMTO source files

`rsutil setup cif` converts an ordered CIF into the `lattice.nml` and
`input.nml` pair that `rsutil setup bulk` consumes. It replaces the external
`cif2cell` step. Parsing and space-group expansion use pymatgen, which is
already a core dependency.

```bash
# convert only
rsutil setup cif Ni2MnGa.cif -o runs/Ni2MnGa

# convert into the current directory, then run the normal validated setup
rsutil setup bulk --cif Ni2MnGa.cif
```

```python
from rsutil.cif import convert_cif

files, conversion = convert_cif("Ni2MnGa.cif", "runs/Ni2MnGa")
print(files.ct, files.neighbors_min)
```

`setup bulk --cif` writes `lattice.nml` and `input.nml` into the working
directory and refuses to replace existing ones unless `--overwrite-source` is
given. The CIF is converted as one unit cell; the setup's own `--supercell`
and `--sqs` then apply as usual, so the cell is never replicated twice. The
legacy `-S` flag is unchanged.

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
- **Cell choice.** By default a centred CIF cell (F, I, A, B, C, R) is reduced
  to its **standard primitive cell**, as `cif2cell` does: for example the
  8-atom cubic zinc-blende FeN cell becomes 2 atoms with
  `a = (½,½,0), (½,0,½), (0,½,½)`. `--cell conventional` keeps the cell stored
  in the CIF. Primitive lattices (P) are never changed. `alat` is the
  conventional cell length in both cases. The F reduction is checked against
  `cif2cell`; the I, A, B, C and R reductions use the usual crystallographic
  primitive vectors and have not been compared with `cif2cell` output.
  If the standard reduction cannot be applied, pymatgen's primitive cell is
  used and a warning says that the vector orientation may differ.
- `--supercell N1 N2 N3` replicates the (primitive or conventional) cell as an
  ordered supercell; `alat` keeps the unreplicated length.

## Neighbour cutoff `ct`

`ct` is the neighbour cutoff in Å. Unless `--ct` is given it is chosen so that
**every** atom has at least `--ct-neighbors` atoms within it, the atom itself
included (default 15, which is the bcc value; fcc gives 19):

1. take the smallest radius at which every atom has that many atoms;
2. extend it to the end of that neighbour shell, treating distances within
   0.1 Å as one shell so a distorted, split shell is never cut through;
3. place `ct` 0.3 Å beyond it, or half the gap to the next shell if smaller.

For a bcc cell this encloses 8 + 6 neighbours plus the atom (15); for fcc
12 + 6 plus the atom (19); for B2 FeCo with `alat = 3.0` it gives `ct = 3.3`.
The chosen `ct` and the atoms per site are printed. Use `--ct-neighbors 19`
(or an explicit `--ct`) for a larger cluster; complex or distorted structures
are often run with more than 15. `r2` defaults to `ct**2`.

## What is *not* derived from the CIF

`rc`, `wav`, and the `&calculation`, `&self`, `&energy`, `&control`, `&mix`,
and `&hamiltonian` groups are settings to review. Pass `--template input.nml` to reuse the non-structural groups from an
existing input instead of the built-in defaults.

## Limitations

- Partially occupied (disordered) CIFs are rejected with the offending sites
  listed. SQS or ordered-supercell handling of disorder is not implemented.
- Existing `lattice.nml`/`input.nml` are never replaced unless `--overwrite`
  is given.

Validated against `cif2cell 2.1.0` output for the pre-martensitic Ni2MnGa
structure (ICSD 190387, 24 atoms) in `tests/data/cif`.
