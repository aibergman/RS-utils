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
# Lattice setup
rsutil --setup --empty_spheres

# Exchange interactions
rsutil --exchange --cutoff 2.5 --unique --shell-width 0.1

# Export to UppASD
rsutil --export --data . --jij jijso.out --out uppasd_out
```

## Programmatic Examples (Python)

```python
from rsutil.bands import LatticeReader, HamiltonianReader, BandStructureCalculator, DOSCalculator, FullBZCalculator

lattice = LatticeReader('lattice.nml')
hamiltonian = HamiltonianReader('rs2txt.dat')

# Bands
calc = BandStructureCalculator(lattice, hamiltonian)
kpoints, klabels, kdist = calc.setup_kpath(npoints=100)
bands = calc.calculate_bands(kpoints)
calc.plot_bands(kdist, bands, klabels, 'bands.png')

# DOS
bz = FullBZCalculator(lattice, hamiltonian)
dos_calc = DOSCalculator(bz)
dos = dos_calc.calculate_dos(mesh_size=(6, 6, 6))
dos_calc.plot_dos(dos, fermi_energy=0.0, output_file='dos.png')
```

## Usage

After building, an executable (for example, `process_lattics`) will be available in your build directory. To run the program 
you need to have RS-LMTO-ASA compliant `input.nml` and `lattice.nml` files present. Based on the structural information in 
`lattice.nml`  and the chemical information in `input.nml`. the binary `process_lattice.x`  will calculate optimal Wigner-Seitz radii
and reasonable starting guesses for occupations, moments, and potential parameters. This will be stored in atomic `.nml` files as well
as in the created `input_new.nml` and `lattice_new.nml`.

Caution: The f-electron treatment is not fully optimal yet so please double check occupations and core/valence electrons in the input is properly 
described.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


