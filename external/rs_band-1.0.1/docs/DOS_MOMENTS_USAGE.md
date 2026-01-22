## DOS Moments Analysis

This document describes how to use the DOS moments analysis functionality in the rs_band package.

### Overview

The DOS moments analysis calculates:
- **Zeroth moment (occupation)**: Number of electrons in each orbital/spin channel
- **First moment (band center)**: Average energy of occupied states
- **Second moment (band width)**: RMS width of the occupied states
- **Fermi level**: Automatically determined from total electron count

### Quick Start

#### Using the CLI
```bash
# Basic usage with default parameters
python -c "from rs_band.cli import dos_moments_analysis; dos_moments_analysis()"

# With custom parameters
python -c "from rs_band.cli import dos_moments_analysis; dos_moments_analysis()" \
    --mesh 8 8 8 \
    --electrons 8.0 \
    --temperature 300 \
    --output-prefix my_analysis

# View help
python -c "from rs_band.cli import dos_moments_analysis; dos_moments_analysis()" --help
```

#### Using Python directly
```python
from rs_band.band_structure_calculator import LatticeReader, HamiltonianReader
from rs_band.full_bz_dos import FullBZCalculator, DOSCalculator
from rs_band.dos_moments import DOSMomentsCalculator

# Read input files
lattice = LatticeReader('lattice.nml')
hamiltonian = HamiltonianReader('rs2txt.dat')

# Create calculators
bz_calc = FullBZCalculator(lattice, hamiltonian)
dos_calc = DOSCalculator(bz_calc)

# Calculate DOS
dos_data = dos_calc.calculate_dos(mesh_size=(6, 6, 6))

# Prepare projected DOS
projected_dos = {}
for key, dos in dos_data['orbital_dos'].items():
    projected_dos[f'orbital_{key}'] = dos
for key, dos in dos_data['spin_dos'].items():
    projected_dos[f'spin_{key}'] = dos

# Calculate moments
moments_calc = DOSMomentsCalculator(
    dos_data['energies'],
    dos_data['total_dos'],
    projected_dos,
    temperature=300.0
)

# Find Fermi level and calculate moments
fermi_level = moments_calc.find_fermi_level(target_electrons=8.0)
moments = moments_calc.calculate_moments()

# Print results
moments_calc.print_summary()

# Save results
moments_calc.save_moments("my_results")
```

### Parameters

#### CLI Parameters
- `--mesh X Y Z`: k-point mesh size (default: 8 8 8)
- `--electrons N`: Number of electrons (auto-detected if not specified)
- `--temperature T`: Temperature in Kelvin (default: 300)
- `--output-prefix PREFIX`: Output file prefix (default: dos_moments)
- `--no-symmetry`: Disable symmetry reduction
- `--gaussian`: Use Gaussian broadening instead of tetrahedron method

#### Python Parameters
- `energies`: Energy grid (numpy array)
- `total_dos`: Total DOS values (numpy array)
- `projected_dos`: Dictionary of projected DOS components
- `temperature`: Temperature in Kelvin (default: 300.0)
- `target_electrons`: Number of electrons for Fermi level determination

### Output

#### Console Output
```
============================================================
DOS MOMENTS SUMMARY
============================================================
Fermi Level: -0.154 eV
Temperature: 300.0 K
============================================================

total:
----------------------------------------
  Occupation:     8.000000 electrons
  Band Center:    -0.310422 eV
  Band Width:     0.136077 eV
  Avg E (rel EF): -0.156213 eV

orbital_s:
----------------------------------------
  Occupation:     1.463558 electrons
  Band Center:    -0.485898 eV
  Band Width:     0.110488 eV
  Avg E (rel EF): -0.331689 eV
...
```

#### File Output
Results are saved to an `.npz` file containing:
- `fermi_level`: Fermi level in eV
- `temperature`: Temperature in K
- `moments`: Dictionary with moments for each component
- `energies`: Energy grid
- `total_dos`: Total DOS
- `projected_dos`: All projected DOS components

### Example: BCC Iron
```python
# For BCC Fe with 8 valence electrons
python -c "from rs_band.cli import dos_moments_analysis; dos_moments_analysis()" \
    --mesh 6 6 6 \
    --electrons 8.0 \
    --temperature 300 \
    --output-prefix bccFe_analysis
```

Expected results for BCC Fe:
- Total occupation: 8.0 electrons
- s orbital: ~1.5 electrons
- p orbital: ~0.9 electrons  
- d orbital: ~5.6 electrons
- Fermi level: ~-0.15 eV

### Convergence Testing

The package includes a comprehensive test script for convergence analysis:
```bash
python test_dos_moments.py
```

This script:
- Tests multiple k-point mesh sizes (4³, 6³)
- Tests different temperatures (0K, 300K)
- Generates summary plots showing convergence
- Provides timing analysis
- Saves detailed results for each configuration

### Advanced Usage

#### Temperature Effects
- Use `temperature=0` for T=0K (step function)
- Use `temperature=300` for room temperature Fermi-Dirac distribution
- Higher temperatures smear the Fermi-Dirac distribution

#### Convergence
- Start with small meshes (4³) for testing
- Use larger meshes (8³ or higher) for production calculations
- Monitor convergence of Fermi level and orbital occupations

#### Error Handling
The code includes robust error handling for:
- No solution in the specified energy range
- Division by zero in moments calculation
- Missing or invalid DOS data
- Numerical integration failures

### Integration with Existing Workflow

The DOS moments analysis integrates seamlessly with the existing rs_band workflow:
1. Use existing lattice and Hamiltonian files
2. Leverages the same DOS calculation infrastructure
3. Compatible with all symmetry and k-point mesh options
4. Can be combined with band structure calculations
