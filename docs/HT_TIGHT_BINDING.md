# HT-14: versioned RS-LMTO tight-binding data

RS-utils supports the versioned real-space export written by the
`fable_v4` `hamiltonian_paoflow_io` module. The source reference is the
[`fable_v4` writer](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/hamiltonian_paoflow_io.f90).

## Supported fable_v4 files

The exporter writes a pair of files, normally `rs_tb.meta` and `rs_tb.tb`.

`rs_tb.meta` starts with `# rs-tb-meta-v1` and contains:

```text
energy_unit eV
length_unit Angstrom
n_atoms <integer>
norb_scalar <integer>
block_size <integer>
index_base 1
lattice_vectors_cart
<three Cartesian vectors>
positions_cart
<one-based site index and three Cartesian coordinates per site>
```

The current writer may also provide `fermi_energy <value> [unit]`,
`spin_representation`, `soc_representation`, or `control_nsp`. They are
optional because the current upstream writer does not emit a Fermi level or
the `control%nsp` value in the TB metadata file.

`rs_tb.tb` starts with `# rs-tb-hoppings-v1`. Data records have twelve fields:

```text
R1 R2 R3 ia jbasis k_neigh i_local j_local ipao jpao real_eV imag_eV
```

The site and orbital fields are one-based in the text format. RS-utils stores
them canonically as zero-based indices. `R=(R1,R2,R3)` is an integer lattice
translation locating the periodic image of `jbasis`; it is not a Cartesian
displacement. The hopping convention is therefore:

```text
H(k)[site_i, site_j] += H(R) * exp(+2π i k_fractional · R)
```

The complex value is preserved as native `complex128` in ML export. The
versioned writer has no overlap-matrix output, so `TightBindingModel.overlap`
is absent for this format.

`norb_scalar` and `block_size` establish a two-channel spinor layout when
`block_size == 2*norb_scalar`, but do not by themselves establish SOC.
When an `input.nml` is supplied, RS-utils applies the RS-LMTO control rule:

- `control%nsp = 1`: magnetic, collinear, no SOC; compatible with GBT;
- `control%nsp = 2` or `4`: SOC enabled;
- other or absent values remain explicit/unknown rather than inferred.

Use `read_fable_v4_tb(..., input_path=...)` when this distinction is needed.

## ML dataset

`export_ml_dataset(model, "dataset")` writes `dataset.npz` and
`dataset.json`. The NPZ contains deterministic, sorted arrays:

| Array | Meaning |
|---|---|
| `R` | integer target-cell translations, shape `(N, 3)` |
| `site_i`, `site_j` | zero-based reference-cell row/column sites |
| `orbital_i`, `orbital_j` | zero-based local row/column orbitals |
| `value` | native complex hopping values, dtype `complex128` |

The JSON carries schema, units, structure, basis, spin/SOC metadata, source
hashes and other provenance. No HDF5/Zarr dependency is required. A future
large-scale storage backend can use the same explicit identity arrays.

Legacy `rs2txt.dat` remains readable by the existing bands reader. Its
conversion method, `HamiltonianReader.to_tight_binding_model(lattice)`, is
fail-closed unless its Cartesian displacement maps to exactly one periodic
site/translation identity; legacy band/DOS analysis does not require this
conversion.
