# RS-utils HT Baseline Inventory

Generated from the checked-out `HT` branch on 2026-08-30 before functional
HT refactoring. This document records observations, commands, and source
locations rather than intended future behavior.

## Repository and runtime identity

| Item | Observed value |
|---|---|
| Package version | `0.10.0` (`setup.py:51`, `rsutil/__init__.py:55`) |
| Baseline git SHA | `e4a536d9660e3172ee5e6b4168eacf7b58331b4b` |
| Latest commit | `Updated docs` |
| External RS-LMTO reference | branch `fable_v4` (not modified) |
| Python used for source probe | `/Users/andersb/envs/p311/bin/python`, Python 3.11.16 |
| Source shared library | loaded: `/Users/andersb/Jobb/RS-utils/build/lib.macosx-14.0-arm64-cpython-311/rsutil/rsutil.so` |

Version and loading probe:

```text
version= 0.10.0
lib_filename= rsutil.so
lib_path= /Users/andersb/Jobb/RS-utils/build/lib.macosx-14.0-arm64-cpython-311/rsutil/rsutil.so
shared_library_loaded= True
```

The package initializer first looks for `rsutil*.so` in the package, then in
build/wheel locations, and loads the first discovered candidate with
`ctypes.CDLL`. It configures the `process_lattice` argument types and leaves
the package importable with `_rsutil = None` if loading fails
(`rsutil/__init__.py:57-96`).

## CLI baseline

The installed `rsutil` launcher is `/opt/homebrew/bin/rsutil` and points at a
different checkout, `/Users/andersb/Jobb/RS_tools`. Therefore the direct
probe failed before CLI parsing:

```text
ModuleNotFoundError: No module named 'packaging.version'
```

The equivalent source entry point succeeded with:

```bash
MPLCONFIGDIR=/tmp/rsutil-mpl PYTHONPATH=. python -m rsutil.cli --help
```

Relevant archived help text from that command:

```text
usage: cli.py [-h] [--setup] [--lattice] [--clust] [--exchange] [--data DATA]
              [--cutoff CUTOFF] [--unique] [--shell-width SHELL_WIDTH] [--soc]
              [--export] [--jij JIJ] [--dmi DMI] [--out OUT] [--lmto47]
              [--nspin NSPIN] [--lmx LMX]
              [--rec-method {block,chebyshev,lanczos}] [--rec-steps REC_STEPS]
              [--lattice-dist LATTICE_DIST] [--basis-dist BASIS_DIST]
              [--alat-scale ALAT_SCALE] [--seed SEED] [--2d]
              [--z-iso-factor Z_ISO_FACTOR] [--sqs SQS]
              [--sqs-shells SQS_SHELLS] [--sqs-steps SQS_STEPS]
              [--sqs-seed SQS_SEED] [--cluster {periodic,free}]
              [--empty_spheres] [--min_empty_radius MIN_EMPTY_RADIUS]
              [--max_empty_radius MAX_EMPTY_RADIUS] [--verbose] [--summary]
              [--validate] [--backup] [--in-place] [--supercell NX NY NZ]
              {bands,dos} ...

options:
  -h, --help            show this help message and exit
  --setup, -S           Process lattice
  --lattice, -L         Run visualizer for lattice input
  --clust, -C           Run visualizer for clust file
  --exchange, -J        Create exchange input files
  --data DATA           Directory for data files to be visualized
  --cutoff CUTOFF       Cutoff for exchange interactions, in lattice units
  --unique              Find only unique neighbors per shell (one
                        representative per unique species-vector combination)
  --shell-width SHELL_WIDTH
                        Width of shells for grouping neighbors when using
                        --unique (default: 0.1)
  --soc                 Account for spin-orbit coupling in symmetry
                        considerations
  --export, -E          Export RS-LMTO inputs to UppASD format
  --jij JIJ             Path to jijso.out file to include in export
  --dmi DMI             Path to dijso.out / dmi file to include in export (DM
                        vectors)
  --out OUT             Output directory for exported UppASD files (relative
                        to --data)
  --lmto47              Export RS-LMTO lattice/input to LMTO47 CTRL format
  --nspin NSPIN         NSPIN value for LMTO47 export (default: 2)
  --lmx LMX             Maximum L channel override (default: 2; spd)
  --rec-method {block,chebyshev,lanczos}
                        Recursion method for &control recur in setup output
                        (default: block)
  --rec-steps REC_STEPS
                        Recursion steps for &control lld in setup output
                        (default: 20, or 200 for chebyshev)
  --lattice-dist LATTICE_DIST
                        Std.dev of Gaussian lattice-vector distortion in units
                        of alat (default: 0.0)
  --basis-dist BASIS_DIST
                        Std.dev of Gaussian basis-vector distortion in units
                        of alat (default: 0.0)
  --alat-scale ALAT_SCALE
                        Scale factor applied to alat before setup/export
                        transforms (default: 1.0)
  --seed SEED           Random seed for reproducible lattice/basis distortions
                        (default: random)
  --2d                  Enable 2D setup mode and derive z-axis length before
                        WS/ASA processing
  --z-iso-factor Z_ISO_FACTOR
                        Isotropic-density surface margin factor for 2D z-axis
                        derivation (requires --2d)
  --sqs SQS             Enable SQS assignment on all sites, e.g. --sqs
                        "Fe=0.5,Ni=0.5"
  --sqs-shells SQS_SHELLS
                        Number of pair distance shells in SQS objective
                        (default: 3)
  --sqs-steps SQS_STEPS
                        Number of MC swap steps for SQS optimization (default:
                        5000)
  --sqs-seed SQS_SEED   Random seed for SQS optimizer (default: random)
  --cluster {periodic,free}
                        Cluster type: periodic (default) or free
  --empty_spheres       Enable empty sphere calculations
  --min_empty_radius MIN_EMPTY_RADIUS
                        Minimum empty sphere radius (default: 1.2)
  --max_empty_radius MAX_EMPTY_RADIUS
                        Maximum empty sphere radius (default: 2.5)
  --verbose, -v         Show detailed progress during processing
  --summary             Print summary of changes after processing
  --validate            Validate output files after generation
  --backup              Create backup files before any operations
  --in-place            Update input files in place (creates .bak backups)
  --supercell NX NY NZ  Create supercell by repeating unit cell NX×NY×NZ times
```

### `-S` / `--setup` options and behavior

The parser exposes setup-relevant controls for:

```text
--cluster {periodic,free}
--empty_spheres --min_empty_radius --max_empty_radius
--lmx --rec-method {block,chebyshev,lanczos} --rec-steps
--lattice-dist --basis-dist --alat-scale --seed
--2d --z-iso-factor
--supercell NX NY NZ
--sqs --sqs-shells --sqs-steps --sqs-seed
--validate --summary --backup --in-place --verbose
```

The CLI requires `lattice.nml` and `input.nml` in the current working
directory. It validates parameter ranges, optionally creates backups, and
calls the Fortran wrapper. Normal setup output is expected as
`lattice_new.nml` and `input_new.nml`; `--in-place` moves those over the
inputs after creating `.bak` copies.

### `-J` / `--exchange` options and current probe

The intended exchange controls are:

```text
--cutoff CUTOFF       default 1.0, in lattice units
--unique              one representative per shell/species/vector grouping
--shell-width WIDTH   default 0.1
--soc                 alter symmetry handling for SOC
--data DATA           accepted globally, but not passed into exchange_map()
```

The implementation at `rsutil/cli.py:2343-2350` also reads
`args.rij_out` and `args.rij_format`, but no corresponding parser options are
defined. A source probe therefore fails immediately:

```text
PYTHONPATH=. python -m rsutil.cli -J
AttributeError: 'Namespace' object has no attribute 'rij_out'
```

This is a baseline finding for a later focused task, not an HT-00 code fix.

## Setup call path: CLI -> ctypes -> Fortran

```text
rsutil.cli.main()
  -> process_lattice(...) imported from rsutil.__init__
  -> rsutil._rsutil.process_lattice(...), a ctypes.CDLL symbol
  -> source/rsutil_c_interface.f90:subroutine process_lattice
  -> source/rsutil_lib.f90:process_lattice_impl
  -> read_lattice_namelist("lattice.nml", ...)
  -> read_input_namelists("input.nml", ...)
  -> transforms / sphere calculations / potential parameters
  -> write_lattice_namelist("lattice_new.nml", ...)
  -> update_namelists("input.nml", "input_new.nml", ...)
```

The C interface sets cluster type, empty-sphere controls, distortion and
recursion configuration in module state before calling the implementation.
The implementation applies `alat`/lattice/basis/2D transforms, calls the
protected sphere routines, and writes the fixed-name outputs.

## Fortran `STOP` / `ERROR STOP` inventory

The following is the executable-source inventory. Commented-out occurrences
were excluded from the list. These paths are why future HT calls must be
isolated from the parent Python interpreter.

| File | Lines / paths |
|---|---|
| `source/io/atom_io.f90` | `77`, `148`: file open/write failure stops |
| `source/io/input_io.f90` | `60`, `273`, `283`, `379`, `455`: namelist open/read/write failures |
| `source/io/lattice_io.f90` | `76`, `87`: lattice file failures |
| `source/rsutil_lib.f90` | `111`, `112`, `113`, `114`, `119`, `131`, `151`, `164`, `166`, `174`, `178`: invalid parameters, singular geometry, and invalid 2D geometry (`ERROR STOP`) |
| `source/spheres/empty.f90` | `300`, `320`: empty-sphere failure paths |
| `source/spheres/radii.f90` | `389`: unknown chemical element |
| `source/spheres/utils.f90` | `69`, `151`, `818`, `1898`, `1917`: numerical/initial-data failure paths |
| `source/types/element_info.f90` | `76`, `93`: element file/read failures |
| `source/potpars/potpar_wrapper.f90` | `64`: unsupported `lmax` (`ERROR STOP`) |

`source/spheres/sphere.f90` has no executable `STOP`/`ERROR STOP` occurrence
in the current source, but it participates in the same fixed-file setup call
path.

## Fixed filenames, working-directory assumptions, and temporary behavior

### Fortran fixed names

- `source/rsutil_lib.f90` reads `lattice.nml` and `input.nml` from the
  process working directory.
- It writes `lattice_new.nml` and `input_new.nml`.
- `source/spheres/sphere.f90` and `source/spheres/radii.f90` write
  `spheres.log` and `radii.out`.
- `source/potpars/potpar_wrapper.f90:27` reads `<label>.nml` and writes
  `<label>_out.nml` through the potential-parameter path.

### Setup staging and supercells

Without `--supercell`, setup operates directly in the current directory and
the kernel sees the real `lattice.nml`/`input.nml` names. With a non-in-place
supercell, `rsutil/cli.py:2181-2250` creates
`.rsutil_supercell_lattice.nml` and `.rsutil_supercell_input.nml`, moves the
real inputs to `.rsutil_original_lattice.nml` and
`.rsutil_original_input.nml`, copies staged files into the fixed names,
invokes Fortran, and restores the originals in `finally` before removing the
staged files. With `--in-place`, `create_supercell()` edits the real inputs,
creates `lattice_original.nml.bak` / `input_original.nml.bak`, and may create
per-label `.nml` files. This is current legacy behavior and a known HT risk;
HT-03/04/05 must replace it in the core path with isolated workspace staging.

## Exchange directory-generation behavior

`rsutil/exchange/exchange_map.py` currently:

1. reads `input.nml` and `lattice.nml` from the current directory;
2. builds a pymatgen `Structure` and `SpacegroupAnalyzer`;
3. reads `clust` into a pymatgen `Molecule` by parsing two sites per ten-token
   cluster record;
4. finds all neighbors or shell/species/vector representatives within the
   cutoff;
5. uses `neighbor.index + 1` as the target index;
6. reads `report.out` for the Fermi energy;
7. for each recursion site, creates `Atom_<index>` only if absent, copies
   `lattice.nml` and every root `*_out.nml`, and force-writes
   `Atom_<index>/input.nml` with exchange post-processing data.

Existing `Atom_*` directories are reused rather than transactionally cleaned,
so stale files can remain. The CLI's global `--data` value is not passed to
the exchange function; exchange therefore remains current-working-directory
based. Requested R-mapping write errors are caught and printed by
`exchange_map()` rather than propagated.

## Current `R_ij` representation and index meaning

`get_Rij_mappings()` emits records with `i`, `j`, `R_ij`, `distance`,
`species_i`, and `species_j`. Current meanings are:

- `i`: the enumerated center position in the parsed cluster, normally
  one-based when serialized;
- `j`: `neighbor.index + 1`, an index in the parsed cluster `Molecule`, not a
  reference-cell basis index;
- `R_ij`: `neighbor.coords - center.coords`, a Cartesian displacement vector;
- `distance`: Cartesian neighbor distance;
- species fields: labels associated with those cluster sites.

Thus current `R_ij` is not an integer periodic lattice translation and cannot
serve as the canonical periodic identity. The required future identity is
`(i, j, R)` with `R=(R1,R2,R3)` integer and reverse operation
`(i,j,R) -> (j,i,-R)`.

## Current SOC and unique-symmetry behavior

`create_structure_and_analyzer(..., soc=True)` obtains pymatgen symmetry
operations and heuristically filters paired z-preserving/z-mirroring
operations. It assigns the filtered list to the private
`SpacegroupAnalyzer._symmetry_operations` attribute. This is not a safe HT
correctness boundary.

When `--unique` is selected, neighbors are grouped by rounded shell distance
(`round(distance / shell_width)`) and species. Candidate Cartesian vectors are
converted to fractional coordinates, rotated by the analyzer operations, and
compared modulo one. This is a prototype vector grouping, not exact periodic
pair mapping. It also does not establish whether an operation preserves a
supplied magnetic reference state.

## Tests and CI status

### Commands run

```text
python -m pytest --collect-only -q
no tests collected in 0.02s
```

The command exits with pytest's no-tests status (`5`). No `tests/` directory,
pytest configuration, or test modules are tracked in this branch.

### CI observations

- `.github/workflows/build.yaml` builds the Fortran/CMake project on Linux,
  macOS, and Windows.
- `.github/workflows/build_cibuildwheel.yaml` builds wheels and uploads them;
  it does not run a test suite.
- `.github/workflows/build_pip.yaml` builds and installs a wheel, then runs
  `python -c "import rsutil; rsutil.process_lattice()"` from
  `examples/FeCo`; it does not invoke pytest.
- `.github/workflows/deploy_docs.yaml` builds documentation.
- `docs/notes/TEST_RESULTS.md` contains a historical successful band/DOS
  report for v0.10.0, but its commands and generated files are not an active
  repository test suite.

## HT-00 evidence and scope check

Only the following repository files are introduced by HT-00:

- `docs/HT_FOUNDATION_BLUEPRINT.md`
- `docs/HT_COMPATIBILITY_CONTRACT.md`
- `docs/HT_BASELINE_INVENTORY.md`

No Python, Fortran, build, test, or existing documentation code was changed.
