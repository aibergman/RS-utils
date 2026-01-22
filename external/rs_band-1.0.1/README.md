# rs_band - Real-space to Band structure analysis

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![PyPI](https://img.shields.io/badge/PyPI-rs--band-orange.svg)](https://pypi.org/project/rs-band/)

A comprehensive Python package for analyzing TB-LMTO (Tight-Binding Linear Muffin-Tin Orbital) electronic structure data, including band structure calculations, density of states, orbital projections, and non-collinear spin analysis.

## 🚀 Features

- **Band Structure Calculations**: High-symmetry k-path generation and band plotting
- **Orbital & Site Projections**: s, p, d orbital contributions and site-specific analysis
- **Non-collinear Spin Analysis**: Full spin texture and magnetic structure analysis
- **Full Brillouin Zone DOS**: Tetrahedron method with k-space symmetry reduction
- **DOS Moments Analysis**: Fermi level determination and orbital occupations (0th, 1st, 2nd moments)
- **Professional Visualization**: Publication-quality plots and comprehensive analysis

## 📦 Installation

### From PyPI (recommended)
```bash
pip install rs-band
```

### From source
```bash
git clone https://github.com/username/rs_band.git
cd rs_band
pip install -e .
```

### Dependencies

- numpy >= 1.20.0
- matplotlib >= 3.5.0
- scipy >= 1.7.0
- spglib >= 2.0.0 (for symmetry detection)

#### Optional
- pymatgen >= 2022.0.0 (enhanced symmetry features)

## 🎯 Quick Start

### Command Line Interface

```bash
# Complete analysis (recommended)
rs-band

# DOS moments analysis with orbital occupations
rs-band --moments --mesh 6 6 6 --electrons 8.0 --temperature 300

# High-resolution DOS calculation
rs-band-dos --mesh 8 8 8

# Orbital projected bands
rs-band-projections --projection orbital
```

### DOS Moments Analysis

Calculate Fermi level and orbital occupations automatically:

```bash
# Basic DOS moments analysis
rs-band --moments

# With custom parameters for BCC Fe
rs-band --moments --mesh 6 6 6 --electrons 8.0 --temperature 300 --output-prefix bccfe_analysis
```

This provides:

- **Fermi level determination** from total electron count
- **Orbital occupations** (s, p, d electrons)
- **Band centers and widths** for each orbital
- **Spin-resolved analysis** for magnetic systems
- **Temperature-dependent** Fermi-Dirac distribution

### Python API

```python
import rs_band

# Load your TB-LMTO data
lattice = rs_band.LatticeReader('lattice.nml')
hamiltonian = rs_band.HamiltonianReader('rs2txt.dat')

# Calculate band structure
calculator = rs_band.BandStructureCalculator(lattice, hamiltonian)
kpoints, klabels, kdistances = calculator.setup_kpath(100)
band_energies = calculator.calculate_bands(kpoints)

# Plot results
calculator.plot_bands(kdistances, band_energies, klabels, 'bands.png')

# DOS moments analysis
from rs_band.dos_moments import DOSMomentsCalculator
from rs_band.full_bz_dos import DOSCalculator

dos_calc = DOSCalculator(calculator.bz_calculator)
dos_data = dos_calc.calculate_dos_with_moments(
    mesh_size=(6, 6, 6),
    target_electrons=8.0,
    temperature=300.0
)
```

## 📚 Documentation

- [DOS Moments Usage Guide](examples/bccFe/DOS_MOMENTS_USAGE.md) - Complete guide for DOS moments analysis
- [Examples](examples/) - BCC Fe and other test cases
- [API Reference](rs_band/) - Full Python API documentation

## 📞 Support

- � Issues: [GitHub Issues](https://github.com/username/rs_band/issues)
- � Discussions: [GitHub Discussions](https://github.com/username/rs_band/discussions)
