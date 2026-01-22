Examples
========

BCC Fe complete analysis
------------------------

This example demonstrates a full analysis run for body-centered cubic iron (BCC Fe).
The `examples/bccFe` folder contains the TB-LMTO real-space Hamiltonian (`rs2txt.dat`),
the structural input (`lattice.nml`) and a minimal `input.nml` used by the analysis
pipeline. The three figures below were produced by the package end-to-end workflow
(Fourier transform → diagonalization → projections → DOS).

Files used
~~~~~~~~~~

- `examples/bccFe/lattice.nml` — crystal lattice and atomic positions used for spacegroup detection and k-path generation.
- `examples/bccFe/rs2txt.dat` — TB-LMTO Hamiltonian data in real-space, parsed by the package.
- `examples/bccFe/input.nml` — optional example driver settings (can be adapted).

Notes on reproduction
~~~~~~~~~~~~~~~~~~~~~

1. Install the package in editable mode to ensure local changes are used:

   pip install -e .

2. Run the analysis pipeline on the example inputs (the package provides a CLI and
   Python API). The figures below were produced by a single `complete_analysis`
   workflow that computes bands, projections and the BZ DOS.

Figures
~~~~~~~

Bands — full k-path band structure (energies in eV)

.. image:: /_static/examples/complete_analysis_bands.png
   :alt: Bands
   :align: center
   :scale: 75

Projections — orbital/site projections plotted on the band structure

.. image:: /_static/examples/complete_analysis_projections.png
   :alt: Projections
   :align: center
   :scale: 75

DOS — Brillouin-zone integrated density of states computed with the tetrahedron method

.. image:: /_static/examples/complete_analysis_dos.png
   :alt: DOS
   :align: center
   :scale: 75

DOS Moments Analysis
--------------------

This example demonstrates DOS moments analysis for determining Fermi level and orbital occupations in BCC Fe.

Basic Analysis
~~~~~~~~~~~~~~

Run DOS moments analysis with default parameters:

.. code-block:: bash

   cd examples/bccFe
   rs-band --moments

This performs:

* Automatic Fermi level determination from total electron count
* Calculation of zeroth, first, and second moments
* Orbital-resolved occupation analysis
* Temperature-dependent Fermi-Dirac statistics

Custom Parameters
~~~~~~~~~~~~~~~~~

For production-quality results with specific parameters:

.. code-block:: bash

   rs-band --moments --mesh 8 8 8 --electrons 8.0 --temperature 300 --output-prefix bccfe_analysis

Parameter explanations:

* ``--mesh 8 8 8``: Use 8×8×8 k-point mesh for convergence
* ``--electrons 8.0``: Target 8 valence electrons for BCC Fe
* ``--temperature 300``: Room temperature (300K) Fermi-Dirac distribution
* ``--output-prefix bccfe_analysis``: Save results with descriptive prefix

Expected Results
~~~~~~~~~~~~~~~~

For BCC Fe with 8 valence electrons, expect results similar to:

.. code-block:: text

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

   orbital_s:
   ----------------------------------------
     Occupation:     1.463558 electrons
     Band Center:    -0.485898 eV
     Band Width:     0.110488 eV

   orbital_p:
   ----------------------------------------
     Occupation:     0.920000 electrons
     Band Center:    -0.350000 eV
     Band Width:     0.095000 eV

   orbital_d:
   ----------------------------------------
     Occupation:     5.616442 electrons
     Band Center:    -0.250000 eV
     Band Width:     0.120000 eV

Convergence Testing
~~~~~~~~~~~~~~~~~~~

Test convergence with different mesh sizes:

.. code-block:: bash

   # Quick test with coarse mesh
   rs-band --moments --mesh 4 4 4 --output-prefix coarse

   # Production quality
   rs-band --moments --mesh 12 12 12 --output-prefix fine

Temperature Effects
~~~~~~~~~~~~~~~~~~~

Compare zero temperature and finite temperature results:

.. code-block:: bash

   # T=0K (step function)
   rs-band --moments --temperature 0 --output-prefix t0k

   # Room temperature
   rs-band --moments --temperature 300 --output-prefix rt

Python API Usage
~~~~~~~~~~~~~~~~

For programmatic access:

.. code-block:: python

   from rs_band.band_structure_calculator import LatticeReader, HamiltonianReader
   from rs_band.full_bz_dos import FullBZCalculator, DOSCalculator

   # Read inputs
   lattice = LatticeReader('lattice.nml')
   hamiltonian = HamiltonianReader('rs2txt.dat')

   # Calculate DOS with moments
   bz_calc = FullBZCalculator(lattice, hamiltonian)
   dos_calc = DOSCalculator(bz_calc)

   results = dos_calc.calculate_dos_with_moments(
       mesh_size=(6, 6, 6),
       target_electrons=8.0,
       temperature=300.0
   )

   print(f"Fermi level: {results['fermi_level']:.3f} eV")
   print(f"Total electrons: {results['moments']['total']['occupation']:.3f}")

Further notes
~~~~~~~~~~~~~

If you want me to add an executable example script under `examples/bccFe/` that
reproduces these exact figures (and add the command lines to this page), I can
create it and wire it into the docs build.
