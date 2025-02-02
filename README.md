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
    cd LMTO-utils
    ```

3. **Create a build directory and configure the project with CMake:**

    ```bash
    cmake -S . -B build
    ```

4. **Build the project:**

    ```bash
    cmake --build .
    ```

## Usage

After building, an executable (for example, `expand_lattice`) will be available in your build directory. To run the program:

```bash
./expand_lattice
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


