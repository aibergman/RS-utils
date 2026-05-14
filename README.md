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

# Run exchange mapping and write a human-readable R_ij table
rsutil --exchange --cutoff 3.0 --unique --rij-out rij_map.txt --rij-format txt

# Export to UppASD format (requires built exchange files)
rsutil --export --jij jijso.out --dmi dmi --out uppasd_export
```

Caution: The f-electron treatment is not fully optimal yet so please double check occupations and core/valence electrons in the input is properly 
described.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


