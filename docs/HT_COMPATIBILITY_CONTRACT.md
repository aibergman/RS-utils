# RS-utils HT Compatibility Contract

This contract defines the external behavior protected during the
high-throughput refactor. It is subordinate to the governing foundation
blueprint and is based on the checked-out baseline documented in
[`HT_BASELINE_INVENTORY.md`](HT_BASELINE_INVENTORY.md).

## Preserved interfaces

The following legacy flag interfaces remain supported:

```text
rsutil -S ...
rsutil --setup ...
rsutil -J ...
rsutil --exchange ...
```

`-S`/`--setup` processes `lattice.nml` and `input.nml` from the working
directory and produces setup outputs. `-J`/`--exchange` builds exchange input
directories from the conventional RS-LMTO files. New APIs and aliases may be
added behind these flags, but the flags are not removed or made dependent on
a new subcommand design.

## Ordinary names and layouts

When the operation is scientifically valid, implementations retain the
ordinary names and recognizable layout expected by existing workflows:

| Workflow | Existing names/layout to preserve |
|---|---|
| Setup input | `lattice.nml`, `input.nml` |
| Setup generated namelists | `lattice_new.nml`, `input_new.nml` |
| Setup sphere/potential artifacts | `radii.out`, `spheres.log`, `<label>.nml`, `<label>_out.nml` |
| Exchange inputs | `Atom_<index>/lattice.nml`, `Atom_<index>/input.nml`, copied `<label>_out.nml` files |
| Exchange source data | `clust`, `report.out`, and root `lattice.nml`/`input.nml` |
| Optional exchange mapping | `R_ij` mapping artifact when the option is available and requested |

`Atom_<index>` is a compatibility disk/job layout, not the scientific pair
identity. Future services may use an isolated staging area and commit these
familiar artifacts on success.

## Compatibility meaning

Compatibility is semantic and workflow-oriented. Harmless changes to
whitespace, namelist ordering, diagnostics formatting, and additional
manifest/report files are allowed. Existing scientifically meaningful
parameters, site ordering, output meaning, and ordinary downstream workflow
conventions must remain stable unless a documented safety or correctness
change applies.

Byte-for-byte identity is not required. Generated files must remain readable
by the legacy RS-LMTO workflows for the supported modes.

## Permitted behavior changes

A behavior may change only when the old behavior is demonstrably one of the
following:

- scientifically incorrect or based on an invalid pair/symmetry identity;
- ambiguous or nondeterministically unsafe;
- vulnerable to stale-file or stale-directory reuse;
- capable of leaving partially modified source inputs;
- incapable of propagating a worker, mapping, or parsing failure;
- dependent on undocumented implementation accidents, private third-party
  attributes, or unsafe in-process process termination.

Each such change requires an explicit note in the master blueprint and a
focused regression or correctness test in the task that introduces it. The
protected Fortran numerical algorithms and their scientifically meaningful
defaults are not changed merely for code style or HT performance.

## Specific HT invariants

- Setup source inputs remain unchanged on failed operations.
- Future HT setup execution isolates fixed Fortran filenames and fatal-stop
  behavior in a worker workspace.
- Deterministic setup and campaign operations record their seed and relevant
  provenance.
- Exchange scientific identity is the periodic `(i, j, R)` key, where `R` is
  an integer lattice translation; `Atom_<index>` remains serialization only.
- Symmetry reduction is conservative for SOC/magnetic states and does not
  mutate third-party private symmetry state.
- Explicitly requested output failures propagate as operation failures.
- Raw calculation evidence is retained when results are parsed.
