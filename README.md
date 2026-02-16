# RS-utils

RS-utils is a collection of modular Fortran utilities designed to support the RS‐LMTO‐ASA electronic structure code. This repository provides routines for reading and writing namelist files, processing lattice and atomic data, performing sphere-based calculations, and other related tasks used in LMTO calculations.

## Description

The RS‐LMTO‐ASA method requires a variety of input and processing routines to handle configuration files, atomic data, and lattice information. This project gathers those utilities into a well-structured, modern Fortran codebase that includes:

- **Namelist I/O:** Subroutines for reading and updating input files (e.g., lattice and atoms namelists).
- **Atomic Utilities:** Functions to determine atomic numbers, electron counts, and other properties.
- **Lattice Processing:** Routines to read, update, and write lattice configurations.
- **Sphere Calculations:** Modules for initializing spheres, parsing structures, performing calculations, and writing output as needed for RS‐LMTO‐ASA.

## Features

- **Modular Design:** Separated modules for types, I/O, and sphere calculations.
- **Namelist Handling:** Robust routines to read and update namelist blocks without disturbing unmodified sections.
- **Atomic Data Lookups:** Efficient lookup tables for atomic numbers and electron counts.
- **CMake Build System:** A modern CMake-based build framework to compile the Fortran sources.

## Installation

### Prerequisites

- A Fortran compiler that supports Fortran 2008 (e.g., gfortran, ifort)
- CMake (version 3.10 or higher)

### Building the Project

1. **Clone the repository:**

    ```bash
    git clone https://github.com/anbergman/RS-utils.git
    ```

2. **Navigate to the project directory:**

    ```bash
    cd RS-utils
    ```

3. **Create a build directory and configure the project with CMake:**

    ```bash
    cmake -S . -B build
    ```

4. **Build the project:**

    ```bash
    cmake --build .
    ```

This will create a binary `rsutil.x` in the `./build` directory of your repository.

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

# Run exchange mapping and write a human-readable R_ij table
rsutil --exchange --cutoff 3.0 --unique --rij-out rij_map.txt --rij-format txt

# Export to UppASD format (requires built exchange files)
rsutil --export --jij jijso.out --dmi dmi --out uppasd_export
```

Caution: The f-electron treatment is not fully optimal yet so please double check occupations and core/valence electrons in the input is properly 
described.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


