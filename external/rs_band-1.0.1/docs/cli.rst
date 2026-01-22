CLI Reference
=============

The ``rs-band`` command provides a unified interface to all analysis tools in the rslmto-band package.

Overview
--------

The CLI uses a subcommand structure: ``rs-band <subcommand> [options]``. This replaces multiple separate entry points with a single, consistent interface.

Available subcommands:

- ``dos``: Calculate Density of States (DOS)
- ``moments``: Analyze DOS moments (occupations, centers, widths)
- ``complete``: Run full combined analysis (bands, DOS, projections)
- ``projections``: Calculate projected band structures
- ``structure``: Analyze crystal structure and generate k-paths
- ``quick``: Run quick analysis (if available)

Common Options
--------------

These options are available for most subcommands:

- ``--lattice LATTICE``: Path to the lattice input file (default: ``lattice.nml``). This file contains crystal structure information including lattice vectors, atomic positions, and atomic numbers.
- ``--hamiltonian HAMILTONIAN``: Path to the Hamiltonian input file (default: ``rs2txt.dat``). This file contains the electronic Hamiltonian matrix elements from TB-LMTO calculations.

DOS Subcommand
--------------

Calculate the Density of States using either the tetrahedron method (accurate, default) or Gaussian broadening (fast).

Usage
~~~~~

.. code-block:: bash

   rs-band dos [options]

Options
~~~~~~~

- ``--mesh NX NY NZ``: k-point mesh size for Brillouin zone sampling (default: 6 6 6). Larger meshes provide better convergence but increase computation time. Example: ``--mesh 8 8 8`` for higher accuracy.
- ``--energy-range EMIN EMAX``: Energy range for DOS calculation in eV (no default). Specify the minimum and maximum energies to sample. Example: ``--energy-range -5 5`` for energies from -5 to 5 eV.
- ``--num-energies N``: Number of energy points in the DOS (default: 1000). More points give smoother DOS curves but increase memory usage.
- ``--fermi-energy EF``: Fermi energy in eV for plotting reference lines (default: 0.0). This draws a vertical line at the Fermi level on DOS plots.
- ``--electrons N``: Number of electrons for automatic Fermi level determination (optional). When specified, the Fermi level is calculated from the DOS to match the target electron count.
- ``--no-symmetry``: Disable symmetry reduction in k-point sampling. Forces full uniform mesh, which is required for tetrahedron method but slower for large meshes.
- ``--no-tetrahedron``: Use Gaussian broadening instead of tetrahedron method. Faster for large meshes but less accurate near band edges.
- ``--gaussian-smearing WIDTH``: Gaussian smearing width in eV when using Gaussian broadening (default: 0.1). Controls the width of the Gaussian functions used for DOS broadening.
- ``--output-prefix PREFIX``: Prefix for output files (default: ``dos``). Files will be named ``PREFIX.png`` (plot) and ``PREFIX_data.npz`` (data).

Examples
~~~~~~~~

.. code-block:: bash

   # Basic DOS calculation with tetrahedron method
   rs-band dos --mesh 6 6 6 --energy-range -3 3

   # DOS with automatic Fermi level determination
   rs-band dos --electrons 8.0 --mesh 6 6 6 --energy-range -5 5

   # High-resolution DOS with custom Gaussian smearing
   rs-band dos --mesh 10 10 10 --energy-range -5 5 --num-energies 2000 --no-tetrahedron --gaussian-smearing 0.05 --output-prefix high_res_dos

   # DOS without symmetry reduction
   rs-band dos --mesh 4 4 4 --energy-range -2 2 --no-symmetry

Moments Subcommand
------------------

Analyze DOS moments including electron occupations, band centers, and bandwidths. Useful for understanding electronic structure properties.

Usage
~~~~~

.. code-block:: bash

   rs-band moments [options]

Options
~~~~~~~

- ``--mesh NX NY NZ``: k-point mesh size (default: 8 8 8). Higher values improve accuracy of moment calculations.
- ``--electrons N``: Target number of electrons for Fermi level determination (optional). If not specified, the code attempts auto-detection from the DOS.
- ``--temperature T``: Temperature in Kelvin for Fermi-Dirac distribution (default: 300.0). Affects occupation calculations at finite temperature.
- ``--no-symmetry``: Disable symmetry reduction (forces uniform mesh).
- ``--gaussian``: Use Gaussian broadening instead of tetrahedron method.
- ``--num-energies N``: Number of energy points (default: 1000).
- ``--energy-range EMIN EMAX``: Energy range for DOS calculation (optional).
- ``--output-prefix PREFIX``: Prefix for output files (default: ``dos_moments``).

Examples
~~~~~~~~

.. code-block:: bash

   # Moments analysis with auto-detected electrons
   rs-band moments --mesh 8 8 8 --temperature 0

   # Specify electron count and energy range
   rs-band moments --electrons 10.5 --energy-range -10 5 --mesh 6 6 6

Complete Subcommand
-------------------

Run a comprehensive analysis including band structure, DOS, and projections in a single command.

Usage
~~~~~

.. code-block:: bash

   rs-band complete [options]

Options
~~~~~~~

- ``--nk-line N``: Number of k-points along each segment of the band path (default: 100). More points give smoother band plots.
- ``--dos-mesh NX NY NZ``: k-point mesh for DOS calculation (default: 6 6 6).
- ``--energy-range EMIN EMAX``: Energy range for DOS (optional).
- ``--num-energies N``: Number of energy points for DOS (default: 500).
- ``--fermi-energy EF``: Fermi energy for plots (default: 0.0).
- ``--skip-dos``: Skip DOS calculation to speed up analysis.
- ``--skip-projections``: Skip projection calculations.
- ``--no-symmetry``: Disable symmetry reduction for DOS.

Examples
~~~~~~~~

.. code-block:: bash

   # Full analysis with default settings
   rs-band complete

   # Skip DOS for faster band structure only
   rs-band complete --skip-dos --nk-line 200

   # High-resolution DOS with custom energy range
   rs-band complete --dos-mesh 8 8 8 --energy-range -5 5 --num-energies 1000

Projections Subcommand
----------------------

Calculate projected band structures showing orbital, site, or spin character.

Usage
~~~~~

.. code-block:: bash

   rs-band projections [options]

Options
~~~~~~~

- ``--projection TYPE``: Type of projection to calculate (choices: ``orbital``, ``site``, ``spin``, ``all``; default: ``all``). ``all`` computes all projection types.
- ``--output-prefix PREFIX``: Prefix for output files (default: ``projected``).

Examples
~~~~~~~~

.. code-block:: bash

   # All projections
   rs-band projections

   # Orbital projections only
   rs-band projections --projection orbital

   # Site-resolved projections
   rs-band projections --projection site

Structure Subcommand
--------------------

Analyze crystal structure and generate high-symmetry k-paths for band structure calculations.

Usage
~~~~~

.. code-block:: bash

   rs-band structure [options]

Note: This subcommand has limited options and uses default file paths.

Quick Subcommand
----------------

Run a quick analysis (if available in the installation).

Usage
~~~~~

.. code-block:: bash

   rs-band quick [options]

Note: Availability depends on the specific installation.

Output Files
------------

Depending on the subcommand, the following files are generated:

- **DOS plots**: ``PREFIX.png`` - Matplotlib plot of the density of states
- **DOS data**: ``PREFIX_data.npz`` - NumPy archive containing:

  - ``energies``: Energy grid array
  - ``total_dos``: Total DOS values
  - ``orbital_dos``: Orbital-projected DOS (if calculated)
  - ``site_dos``: Site-projected DOS (if calculated)
  - ``spin_dos``: Spin-projected DOS (if calculated)

- **Band structure plots**: ``PREFIX_bands.png`` (for complete analysis)
- **Projection data**: Additional arrays in the NPZ file for projections

Tips and Best Practices
-----------------------

1. **Mesh Size**: Start with smaller meshes (4x4x4) for testing, then increase to 8x8x8 or higher for production calculations.
2. **Energy Range**: Choose based on your system's band structure. For metals, include the Fermi level; for insulators, focus on the band gap region.
3. **Tetrahedron vs Gaussian**: Tetrahedron method is more accurate but slower. Use Gaussian for quick previews or large meshes.
4. **Symmetry**: Enable symmetry reduction for faster calculations unless using tetrahedron method (which requires uniform mesh).
5. **Memory**: Large meshes and high energy resolution increase memory usage. Monitor system resources for very large calculations.

Troubleshooting
---------------

- **No tetrahedron method**: If tetrahedron calculation fails, try ``--no-tetrahedron`` for Gaussian broadening.
- **Memory errors**: Reduce mesh size or number of energy points.
- **Slow calculations**: Use ``--no-symmetry`` only when necessary, and consider Gaussian broadening for large meshes.
- **File not found**: Ensure ``lattice.nml`` and ``rs2txt.dat`` are in the current directory or specify full paths.

For more advanced usage or custom analysis, see the API documentation or source code.
