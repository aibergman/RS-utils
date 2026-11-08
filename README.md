# RS-utils

RS-utils is a Python + Fortran toolkit for RS-LMTO workflows. It now includes a pure-Python band structure and DOS analysis module alongside the existing Fortran utilities for lattice setup and exchange/export workflows.

> **HT pilot release:** RS-utils `0.11.0rc1` is ready for controlled
> high-throughput RS-LMTO use within the documented support matrix. Start with
> the [HT quick start](docs/HT_QUICKSTART.md) and review the [HT readiness
> report](docs/HT_READINESS_REPORT.md).

## Components

- **Band/DOS analysis (Python):**
    - High-symmetry band structures with projections (orbital/site/spin)
    - DOS with tetrahedron or Gaussian broadening
    - DOS moments and Fermi level determination
    - Requires only Python + NumPy/SciPy/Matplotlib/Spglib (no Fortran runtime)

- **Legacy RS-LMTO utilities (Fortran):**
    - Namelist I/O and lattice processing
    - Exchange interaction mapping and export to UppASD
    - Export to LMTO47 CTRL format
    - Supercell generation with optional all-site SQS assignment
    - CMake-based build of the Fortran library and tools

## Features

- **Modular Design:** Separated modules for types, I/O, and sphere calculations.
- **Namelist Handling:** Robust routines to read and update namelist blocks without disturbing unmodified sections.
- **Atomic Data Lookups:** Efficient lookup tables for atomic numbers and electron counts.
- **CMake Build System:** A modern CMake-based build framework to compile the Fortran sources.

## Build and Installation

### Release wheel (recommended for users)

Install the `0.11.0rc1` release wheel in a dedicated virtual environment;
see the [HT quick start](docs/HT_QUICKSTART.md) for the cluster deployment
model.

### Developer checkout

```bash
pip install -e '.[test]'
```

The editable install builds and exposes the native Fortran library in the
rsutil/ package automatically. No manual shared-library copying or
build_ext --inplace step is required. For the complete clean-checkout
procedure, see [HT_TESTING.md](docs/HT_TESTING.md).

Core runtime dependencies are `numpy`, `f90nml`, and `pymatgen`. Install
`.[bands]` for band/DOS analysis and plotting, `.[visualization]` for the
PyVista structure viewer, or `.[full]` for all optional runtime features.

### Fortran build (optional; legacy utilities)

```bash
cmake -S . -B build
cmake --build build
```

This builds the Fortran library and executable targets (for example `rsutil.x` in the build tree).

### Tests

The active HT regression suite uses the test extra:

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```

Release wheels contain the ABI-qualified native library under rsutil/ and
load it from that installed package directory. For a clean-checkout
procedure, including the required Linux Fortran toolchain, see
[HT_TESTING.md](docs/HT_TESTING.md).

## CLI Usage

The installed CLI entrypoint is:

```bash
rsutil --help
```

### Single-run setup workflows (`setup`)

```bash
# Basic bulk setup (reads lattice.nml and input.nml in the current directory)
rsutil setup bulk

# Empty-sphere setup with the optional Questaal FINDES backend
rsutil setup bulk --empty_spheres --empty-sphere-backend questaal

# Supercell setup
rsutil setup bulk --supercell 1 1 2

# In-place replacement (otherwise setup is non-destructive)
rsutil setup bulk --supercell 1 1 2 --in-place

# 2D mode (default c-axis rule)
rsutil setup bulk --2d

# 2D mode with isotropic-density margins
rsutil setup bulk --2d --z-iso-factor 0.5

# SQS on all supercell sites (requires --supercell)
rsutil setup bulk --supercell 2 2 2 --sqs "Fe=0.5,Ni=0.3,Co=0.2" \
       --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234

# Recursion control override in &control
rsutil setup bulk --rec-method chebyshev           # default rec-steps becomes 200
rsutil setup bulk --rec-method lanczos --rec-steps 64

# Prepare an RS-LMTO-ASA surface from a converged bulk calculation
rsutil setup surface --bulk tmp/fccCu --output tmp/fccCu001 \
       --hkl 0 0 1 --layers 6 --empty-spheres 1

# Expand an existing surface using its type restarts
rsutil setup surface --from-surface tmp/fccCu001 --output tmp/fccCu001_more \
       --layers 7 --empty-spheres 2

# Set up an impurity in a converged bulk or surface host
rsutil setup impurity --host tmp/fccCu --xyz impurity.xyz --output tmp/Cu_imp

# Reproducible distortions and scaling
rsutil setup bulk --lattice-dist 0.01 --basis-dist 0.02 --alat-scale 1.01 --seed 1234
```

``rsutil -S`` and ``rsutil --setup`` remain supported as shortcuts for
``rsutil setup bulk``. The campaign workflow is separate and creates multiple
setup directories from a campaign specification.

The empty-sphere backend can also be selected through the Python setup API
with ``SetupSpec(empty_spheres=True, empty_sphere_backend="questaal")``.
The default backend remains ``xband``.

### Export workflows (`--export`, `--lmto47`)

```bash
# Export to UppASD format (requires built exchange files)
rsutil --export --jij jijso.out --dmi dmi --out uppasd_export

# Export to LMTO47 CTRL format
rsutil --lmto47 --nspin 2 --lmx 3 --out lmto47_out

# LMTO47 export from staged supercell (non-destructive by default)
rsutil --lmto47 --supercell 2 2 2 --out lmto47_out
```

The UppASD `momfile` moments (muB) are read from `report.out` in `--data`
(`Spin moment of atom N`), or from the file given with `--report`. If the
report or an atom's entry is missing, that moment falls back to 1.0 muB and a
warning is printed.

Existing UppASD export functionality remains available, but automated UppASD
integration is not currently part of the validated `0.11.0rc1` HT campaign
workflow. RS-utils does not claim an automatic `SCF -> exchange -> UppASD`
pipeline.

### Visualization workflows (`viz`)

```bash
# Plot DOS files from an RS-LMTO calculation (filled overlay by default)
rsutil viz dos ./calc --fermi 0.7 --output dos.png

# Select species, change layout, or disable filling
rsutil viz dos ./calc --species Ni-S Ni-S-1 --layout stack --no-fill

# Interactive structure viewers
rsutil viz lattice ./calc
rsutil viz clust ./calc

# Legacy structure-viewer flags remain supported
rsutil --lattice --data .
rsutil --clust --data .
```

``rsutil dos`` calculates DOS from a Hamiltonian; ``rsutil viz dos`` plots
DOS output files already present in a calculation directory. DOS plotting
requires the ``bands`` extra (``pip install 'rsutil[bands]'``); interactive
structure viewing requires the ``visualization`` extra.
The plotter reads ``&energy fermi`` from ``input.nml`` in the calculation
folder by default; ``--fermi`` overrides it.
In interactive overlay plots, click a legend entry to hide or show that
species' curves.

### Exchange mapping workflow (`--exchange`)

```bash
# Exchange map with readable R_ij output
rsutil --exchange --cutoff 3.0 --unique --rij-out rij_map.txt --rij-format txt
```

### Recommended high-throughput workflow

The established `-S` and `-J` flags remain supported for individual legacy
runs. For reproducible batch preparation, validate and generate a versioned
campaign first, then render a site-specific Slurm array for review:

```bash
rsutil setup bulk --validate --summary
rsutil -J --cutoff 3.0 --unique                # legacy exchange jobs

rsutil campaign validate examples/campaigns/ht18_pilot.json
rsutil campaign create examples/campaigns/ht18_pilot.json --destination ht18-runs
rsutil campaign render-slurm ht18-runs/campaign_index.json \
    --executable /path/to/rs-lmto --launcher 'srun --exclusive' \
    --time 02:00:00 --resource cpus_per_task=8
rsutil campaign status ht18-runs
rsutil campaign collect ht18-runs --csv ht18-runs/results.csv
```

Exchange is a separate, explicit second campaign stage. Inspect convergence
first; by default only successfully executed, parse-valid, explicitly
converged parents are eligible, and full pairs are the conservative starting
policy:

```bash
rsutil campaign exchange prepare ht18-runs/campaign_index.json \
    --destination ht18-runs/exchange-stage --cutoff 3.0 --full-pairs
rsutil campaign exchange render-slurm ht18-runs/exchange-stage/<stage-id> \
    --executable /path/to/rs-lmto-exchange \
    --launcher 'srun --exclusive' --time 02:00:00 \
    --resource cpus_per_task=8
cd ht18-runs/exchange-stage/<stage-id> && sbatch exchange_submit_array.sbatch && cd ../../..
rsutil campaign exchange status ht18-runs/exchange-stage/<stage-id>
rsutil campaign exchange collect ht18-runs/exchange-stage/<stage-id>
```

Each stage has its own persistent manifest and artifacts. Use `--unique` only
with an explicitly supported magnetic context. RS-utils does not yet
implement a generic automatic SCF→exchange DAG or retry engine.

Campaign generation is setup-only and does not submit jobs. Inspect the
rendered script and supply the executable, launcher, modules, and resources
for the target site. `status` and `collect` use filesystem evidence only;
they do not query Slurm or retry jobs. See [HT_CAMPAIGNS.md](docs/HT_CAMPAIGNS.md)
and [HT_READINESS_REPORT.md](docs/HT_READINESS_REPORT.md) for the support
matrix and closeout evidence.

Entry-level Slurm resources are optional operational metadata. When overrides
are present, rendering groups entries into homogeneous resource classes and
writes `resource_plan.json` plus `submit_all.sh`; `--resource-aware` requests
that layout explicitly. Resource allocation is operational policy, not
scientific calculation identity, and RS-utils does not yet claim an optimal
automatic MPI sizing model.

### Bands and DOS workflows (pure Python)

```bash
# Band structure
rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

# Projected bands
rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
       --projections --projection orbital

# DOS with Fermi level and moments
rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 \
       --electrons 16 --moments
```

LMAX notes
----------

- ``--lmx`` controls two paths:
  - Setup path (``--setup``): overrides RS-LMTO generation ``lmax`` (default ``2`` = spd).
  - LMTO47 path (``--lmto47``): writes ``LMX`` in LMTO47 ``CLASS`` entries.
- Current backend support is ``0 <= lmax <= 3`` (spdf); larger values are rejected.
- ``--lattice-dist`` and ``--basis-dist`` are Gaussian std.dev values in units of ``alat``.
- ``--alat-scale`` scales ``alat`` before applying distortions.
- ``--seed`` makes distortions reproducible across runs.
- ``--rec-method`` overrides ``&control recur`` in ``input_new.nml``.
- ``--rec-steps`` overrides ``&control lld`` in ``input_new.nml``.
- Recursion defaults: ``rec-method=block`` and ``rec-steps=20``.
- If ``rec-method=chebyshev`` and ``--rec-steps`` is omitted, default ``rec-steps`` is ``200``.

2D mode notes
-------------

- ``--2d`` enables slab-oriented setup/export geometry handling.
- If ``--2d`` is given without ``--z-iso-factor``:
  - ``Lz = (zmax - zmin) + 1.0*alat``.
- If ``--2d --z-iso-factor X`` is given:
  - ``Lz = (zmax - zmin) + X*sqrt(A/n_bottom) + X*sqrt(A/n_top)``.
- In both cases, the third lattice vector is aligned with the normal to ``a1,a2``:
  - ``a3 || cross(a1,a2)``, sign chosen to keep positive global ``z`` projection.
- In setup output (``input_new.nml``), 2D mode enforces ``b3=F`` and ``n3=1``.

Caveats:

- ``--z-iso-factor`` requires ``--2d``.
- Surface-layer counting for isotropic mode uses a small z-tolerance; strongly rumpled surfaces may need manual verification.
- This is a practical ASA setup heuristic, not a full electrostatic vacuum convergence protocol; verify sensitivity for production runs.

Surface preparation
-------------------

``rsutil setup surface`` creates a new surface input directory from a
converged bulk calculation. It stages the host ``*_out.nml`` potential files,
host radial ground-state data, and the bulk ``lattice.nml`` when
``crystal_sym='file'``. The generated ``input.nml`` uses
``pre_processing='buildsurf'``, ``calctype='S'``, the selected ``surftype``,
``fix_fermi=.true.``, and the bulk Fermi level. RS-LMTO creates its cluster and
surface maps when its ``buildsurf`` preprocessor runs.

``--layers`` sets ``nlay`` and the added surface/interface type count. In
RS-LMTO, ``nlay`` counts the near-surface planes assigned self-consistent
representative types. The output enforces ``ntype_surface = ntype_bulk +
nlay``. ``--empty-spheres`` chooses how many of those added types use the
empty-sphere label; it does not insert extra geometric sites. The
``buildsurf`` preprocessor constructs the cluster and maps its selected
planes to the labels. RS-utils clones the converged bulk host atom potential
into the interface-layer starting files. It generates empty-sphere radial
parameters with the existing Fortran setup kernel, then assigns a small
starting occupancy while preserving the kernel's native namelist array
indexing. The initial s, p, and d logarithmic-derivative parameters use the
lowest-branch free-electron values ``1.5``, ``2.25``, and ``3.1476`` for both
spins. No empty-sphere search in the bulk geometry is required. By default, each
empty-sphere type starts with a 0.2-electron occupancy;
when interface types are present, that charge is removed from their guesses
with an erf-derived profile concentrated near the outermost interface layer.
The profile is normalized so the total transfer equals the empty-sphere seed
charge. ``--charge-sigma`` controls its width in layer spacings and defaults
to 0.5. This is an initial guess, not a predicted surface charge. Use ``--empty-sphere-charge``
to change the occupancy or set it to zero. The workflow writes a
``surface_setup.json`` manifest alongside the namelist.
Defaults are ``n1=n2=n3=32``, ``rc=100``, ``nstep=15``, and ``beta=0.01``;
each can be overridden on the command line.

Crystal symmetry is inherited from the bulk input or can be selected with
``--crystal-sym fcc|bcc|hcp|file``. ``file`` requires the bulk
``lattice.nml``. ``--hkl`` accepts any nonzero integer surface normal and
reduces common factors before writing ``surftype``. For hcp, a three-value
``h k l`` input is expanded to Miller-Bravais ``h k -(h+k) l``; four values
may also be given directly, with the required ``i=-(h+k)`` relation.

``--orientation preserve`` passes the requested normal through as
``surftype``. ``--orientation global-z`` rotates the lattice and basis
together, writes the rotated cell to ``lattice.nml``, and uses
``crystal_sym='file'`` with ``surftype='0 0 1'``. This rigid rotation retains
the bulk stacking. For built-in fcc, bcc, or hcp symmetry, RS-utils creates
the same primitive cell and basis used by RS-LMTO when no source
``lattice.nml`` is present; hcp uses the bulk ``celldm`` value when set and
otherwise defaults to ``c/a=1.633``.

For ``--orientation global-z``, ``--align-bond X Y Z`` together with
``--align-bond-to U V W`` fixes the remaining in-plane rotation. The first is
a Cartesian bond direction in the original bulk frame; the second is its
desired direction in the rotated global XY plane. The bond must lie in the
selected surface plane. For example, fcc(111) can align the original
``(0,0.5,-0.5)`` bond with global ``(1,-1,0)`` using
``--align-bond 0 0.5 -0.5 --align-bond-to 1 -1 0``. The rotation matrix and
alignment are recorded in ``surface_setup.json``.

An existing RS-utils surface can be expanded without restarting from bulk
atomic guesses. For example, five interface types plus two empty-sphere types
use ``--layers 7 --empty-spheres 2`` because ``nlay`` includes ES types.
Matching per-type restarts are retained; additional interface types use the
deepest existing interface restart, and new ES types use an existing ES
restart. The surface geometry, fixed Fermi level, and other settings are
inherited from its manifest. The manifest's original bulk directory must
remain accessible, and the output must be a new directory:

```bash
rsutil setup surface --from-surface surface_3plus1 \
       --output surface_5plus2 --layers 7 --empty-spheres 2
```

Impurity preparation
--------------------

``rsutil setup impurity`` supports converged bulk and surface hosts. The XYZ
coordinates are in Å, in the host calculation's Cartesian frame, and must
match sites in its ``clust`` file. The host's converged restart potentials and
fixed Fermi level are carried into the new setup. A positive ``--scf-radius``
also selects nearby host sites as independent self-consistent types; empty
spheres are skipped. Each selected site adds a type, so larger radii increase
the number of types and atom files. The available growth region is limited by
the host's existing ``clust``.

```bash
rsutil setup impurity --host bulk_or_surface --xyz impurity.xyz \
       --output impurity_run --scf-radius 4.0
```

Setup safety workflow
---------------------

Recommended workflow for destructive or batch updates:

```bash
# 1) Generate outputs and validate
rsutil setup bulk --validate --summary

# 2) Optional automatic replacement with backups
rsutil setup bulk --in-place --backup
```

SQS workflow (all-site substitution)
------------------------------------

SQS can be applied to all sites in a generated supercell:

```bash
rsutil setup bulk --supercell 2 2 2 \
       --sqs "Fe=0.5,Ni=0.3,Co=0.2" \
       --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234
```

Notes:

- ``--sqs`` currently requires ``--supercell``.
- Fractions are normalized internally if they do not sum exactly to 1.
- Integer site counts are chosen to match the requested fractions as closely as possible for the given supercell size.
- Distortions/scaling can be combined and are applied after SQS assignment.

## Notes

Caution: The f-electron treatment is not fully optimal yet so please double check occupations and core/valence electrons in the input is properly
described.

## License

This project is licensed under the GNU Affero General Public License, version 3 or later.
See [LICENSE](LICENSE) for the license identifier and official license text link.
