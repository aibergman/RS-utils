# RS-utils

RS-utils is a Python + Fortran toolkit for RS-LMTO workflows. It now includes a pure-Python band structure and DOS analysis module alongside the existing Fortran utilities for lattice setup and exchange/export workflows.

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

### Python package (recommended)

```bash
pip install -e .
```

Key runtime dependencies: `numpy`, `scipy`, `matplotlib`, `spglib` (with `pymatgen/seekpath` for k-path workflows).

### Fortran build (optional; legacy utilities)

```bash
cmake -S . -B build
cmake --build build
```

This builds the Fortran library and executable targets (for example `rsutil.x` in the build tree).

## CLI Usage

The installed CLI entrypoint is:

```bash
rsutil --help
```

### Setup workflows (`--setup`)

```bash
# Basic setup
rsutil --setup

# Supercell setup
rsutil --setup --supercell 1 1 2

# In-place replacement (otherwise setup is non-destructive)
rsutil --setup --supercell 1 1 2 --in-place

# 2D mode (default c-axis rule)
rsutil --setup --2d

# 2D mode with isotropic-density margins
rsutil --setup --2d --z-iso-factor 0.5

# SQS on all supercell sites (requires --supercell)
rsutil --setup --supercell 2 2 2 --sqs "Fe=0.5,Ni=0.3,Co=0.2" \
       --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234

# Recursion control override in &control
rsutil --setup --rec-method chebyshev           # default rec-steps becomes 200
rsutil --setup --rec-method lanczos --rec-steps 64

# Reproducible distortions and scaling
rsutil --setup --lattice-dist 0.01 --basis-dist 0.02 --alat-scale 1.01 --seed 1234
```

### Export workflows (`--export`, `--lmto47`)

```bash
# Export to UppASD format (requires built exchange files)
rsutil --export --jij jijso.out --dmi dmi --out uppasd_export

# Export to LMTO47 CTRL format
rsutil --lmto47 --nspin 2 --lmx 3 --out lmto47_out

# LMTO47 export from staged supercell (non-destructive by default)
rsutil --lmto47 --supercell 2 2 2 --out lmto47_out
```

### Visualization workflows (`--lattice`, `--clust`)

```bash
# Lattice visualization
rsutil --lattice --data .

# Cluster visualization
rsutil --clust --data .
```

### Exchange mapping workflow (`--exchange`)

```bash
# Exchange map with readable R_ij output
rsutil --exchange --cutoff 3.0 --unique --rij-out rij_map.txt --rij-format txt
```

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

Setup safety workflow
---------------------

Recommended workflow for destructive or batch updates:

```bash
# 1) Generate outputs and validate
rsutil --setup --validate --summary

# 2) Optional automatic replacement with backups
rsutil --setup --in-place --backup
```

SQS workflow (all-site substitution)
------------------------------------

SQS can be applied to all sites in a generated supercell:

```bash
rsutil --setup --supercell 2 2 2 \
       --sqs "Fe=0.5,Ni=0.3,Co=0.2" \
       --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234
```

Notes:

- ``--sqs`` currently requires ``--supercell``.
- Fractions are normalized internally if they do not sum exactly to 1.
- Integer site counts are chosen to match the requested fractions as closely as possible for the given supercell size.
- Distortions/scaling can be combined and are applied after SQS assignment.

Namelist compatibility notes
----------------------------

The setup pipeline preserves selected extended ``&lattice`` fields when present in input:

- ``strux_backend``
- ``strux_want_sdot``
- ``strux_solve_scale``

These keys are only written to ``input_new.nml`` when they were present in the original input.

Caution: The f-electron treatment is not fully optimal yet so please double check occupations and core/valence electrons in the input is properly 
described.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
