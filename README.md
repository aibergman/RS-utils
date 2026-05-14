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
    - CMake-based build of the Fortran library and tools

## Features

- **Modular Design:** Separated modules for types, I/O, and sphere calculations.
- **Namelist Handling:** Robust routines to read and update namelist blocks without disturbing unmodified sections.
- **Atomic Data Lookups:** Efficient lookup tables for atomic numbers and electron counts.
- **CMake Build System:** A modern CMake-based build framework to compile the Fortran sources.

## Installation

### Python (recommended)

```bash
pip install -e .
```

Key runtime dependencies: numpy, scipy, matplotlib, spglib (pymatgen/seekpath used for k-paths).

### Fortran build (optional, for legacy utilities)

```bash
cmake -S . -B build
cmake --build build
```

## CLI Usage

Band structure and DOS (pure Python):

```bash
# Band structure (high-symmetry path)
rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

# Projected bands
rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --projections --projection orbital

# DOS with Fermi level and moments
rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 --electrons 16 --moments
```

Legacy workflow (Fortran-backed):

    ```bash
    cmake --build .
    ```

This will create a binary `process_lattice.x` in the `./build/bin` directory of your repository.

## Usage

After building, an executable (for example, `rsutil.x`) will be available in your build directory. To run the program you need
RS-LMTO-ASA compliant `input.nml` and `lattice.nml` files present. Based on the structural information in `lattice.nml`
and the chemical information in `input.nml`, the binary will calculate optimal Wigner–Seitz radii and reasonable starting
guesses for occupations, moments, and potential parameters. Results are written to atomic `*_out.nml` files and updated
`input.nml` / `lattice.nml` copies placed in per-atom output directories when post-processing is requested.

Alternatively, the Python package exposes a `rsutil` console script and higher-level tools (visualizer, exchange mapping,
export utilities). Install the package into your Python environment to use the CLI:

```bash
pip install .
# or for editable development install:
pip install -e .
```

Python CLI examples
-------------------

After installing the Python package, the `rsutil` console script is available. Example workflows:

```bash
# Show help and available options
rsutil --help

# Process lattice inputs (same as the Fortran driver)
rsutil --setup

# Create a 1x1x2 supercell (unique labels and per-species .nml files)
rsutil --setup --supercell 1 1 2

# Explicitly replace input files (otherwise setup is non-destructive)
rsutil --setup --supercell 1 1 2 --in-place

# Run exchange mapping and write a human-readable R_ij table
rsutil --exchange --cutoff 3.0 --unique --rij-out rij_map.txt --rij-format txt

# Export to UppASD format (requires built exchange files)
rsutil --export --jij jijso.out --dmi dmi --out uppasd_export

# Export to LMTO47 CTRL format
rsutil --lmto47 --nspin 2 --lmx 3 --out lmto47_out

# Apply reproducible random distortions and scale alat
rsutil --setup --lattice-dist 0.01 --basis-dist 0.02 --alat-scale 1.01 --seed 1234
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

Setup safety workflow
---------------------

Recommended workflow for destructive or batch updates:

```bash
# 1) Generate outputs and validate
rsutil --setup --validate --summary

# 2) Optional automatic replacement with backups
rsutil --setup --in-place --backup
```

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
