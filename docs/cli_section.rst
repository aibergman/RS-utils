.. _cli_section:

Command-Line Interface (CLI)
=============================

The **RS-utils CLI** provides two main Python-only subcommands for band/DOS analysis plus legacy flags for lattice, exchange, export, and visualization.

Top-level usage:

.. code-block:: bash

    rsutil {bands,dos} [options]
    rsutil [--setup|--exchange|--export|--lmto47|--lattice|--clust] [options]

Common Tasks
------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Task
     - Command
   * - Setup from ``input.nml`` + ``lattice.nml``
     - ``rsutil --setup``
   * - Setup and validate generated files
     - ``rsutil --setup --validate --summary``
   * - Supercell setup (non-destructive default)
     - ``rsutil --setup --supercell 2 2 2``
   * - 2D slab setup with isotropic margin
     - ``rsutil --setup --2d --z-iso-factor 0.5``
   * - Override recursion controls in ``&control``
     - ``rsutil --setup --rec-method chebyshev --rec-steps 200``
   * - Build exchange map
     - ``rsutil --exchange --cutoff 2.5 --unique``
   * - Export UppASD input set
     - ``rsutil --export --jij jijso.out --out uppasd_out``
   * - Export LMTO47 CTRL
     - ``rsutil --lmto47 --nspin 2 --lmx 3 --out lmto47_out``
   * - Plot band structure
     - ``rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png``
   * - Compute DOS + moments
     - ``rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 --electrons 16 --moments``

Band and DOS
------------

**Bands**

.. code-block:: bash

    # Band structure (high-symmetry path)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

    # Projected bands
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
                 --projections --projection orbital

**DOS**

.. code-block:: bash

    # DOS with tetrahedron method
    rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6

    # DOS moments + Fermi level from electron count
    rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 \
               --electrons 16 --moments


Advanced Band Structure Features
---------------------------------

**K-Path Generation Modes**

The ``--kpath-mode`` option controls how high-symmetry k-paths are generated:

.. code-block:: bash

    # Auto mode (simple, fast; default)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200

    # Seekpath mode (explicit high-symmetry path)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode seekpath

    # Custom path with specific labels
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 300 \
                 --kpath-mode custom --custom-path "K-Gamma-K"

    # Interactive mode: prompts for seekpath vs custom
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode interactive

**Seekpath Configuration**

Fine-tune seekpath behavior:

.. code-block:: bash

    # Stricter symmetry tolerance
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode seekpath --seekpath-symprec 1e-6

    # Finer Brillouin zone sampling
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode seekpath --seekpath-distance 0.01

**Custom K-Paths with Prime Convention**

Define custom paths and optionally apply the prime convention (K' = -K):

.. code-block:: bash

    # Simple custom path (K → Γ → M)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 300 \
                 --kpath-mode custom --custom-path "K-Gamma-M"

    # With prime points (K' is K with in-plane flip)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 300 \
                 --kpath-mode custom --custom-path "K-Prime-Gamma-K" --prime-points

**Degenerate Band Stabilization**

For materials with degenerate bands, choose a stabilization strategy for smooth color projections:

.. code-block:: bash

    # Block average (recommended; default)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --projections --projection spin \
                 --stabilize-degenerate block_average

    # Pair cancellation (heuristic; useful for many degeneracies)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --projections --projection spin \
                 --stabilize-degenerate pair_cancel

    # No stabilization (raw values; for debugging)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --projections --projection spin \
                 --stabilize-degenerate none

**Brillouin Zone Visualization**

Plot the 3D Brillouin zone polyhedron with high-symmetry points labeled:

.. code-block:: bash

    # Plot BZ alongside bands
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --plot-bz

    # With prime points labeled
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode custom --custom-path "K-Gamma-M" \
                 --plot-bz --prime-points

**K-Path Export for Reproducibility**

Save the exact k-path used for future reference or use with other tools:

.. code-block:: bash

    # Export k-path to file
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --export-kpath kpath_band_structure.dat

Output file format (space-separated columns):

::

    x_linearcoord  kx  ky  kz
    0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
    1.234567e-02  -1.234567e-02   2.345678e-03   0.000000e+00
    ...

**Complete Workflow Example**

Combine seekpath, spin projections, degenerate stabilization, BZ plot, and k-path export:

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 250 \
                 --projections --projection spin \
                 --projections-mode expectation \
                 --kpath-mode seekpath \
                 --stabilize-degenerate block_average \
                 --plot-bz \
                 --export-kpath kpath_seekpath.dat \
                 --output-prefix band_analysis_spin

This command:
1. Generates a seekpath k-path with 250 points
2. Computes spin projections with ⟨S_z⟩ expectation values
3. Stabilizes degenerate manifolds using block averaging
4. Plots the Brillouin zone
5. Saves the k-path to ``kpath_seekpath.dat``
6. Outputs plots to ``band_analysis_spin_*.png``

**Available Options Summary**

=================================  ============================================  ================================================================
Option                             Values                                        Description
=================================  ============================================  ================================================================
``--kpath-mode``                   auto, seekpath, custom, interactive           K-path generation method (default: auto)
``--seekpath-symprec``             float [0.00001, 0.1]                          Symmetry tolerance for seekpath (default: 1e-5)
``--seekpath-distance``            float [0.001, 0.1]                            Reference distance for BZ sampling in Å^-1 (default: 0.025)
``--custom-path``                  string (e.g., "K-Gamma-K")                    Labels for custom k-path (required with --kpath-mode custom)
``--prime-points``                 flag                                          Apply prime convention to k-points
``--stabilize-degenerate``         none, block_average, pair_cancel               Degenerate band handling (default: block_average)
``--plot-bz``                      flag                                          Plot Brillouin zone polyhedron (requires pymatgen)
``--export-kpath``                 filename                                      Export k-path to file for reproducibility
=================================  ============================================  ================================================================

Lattice Processing
-------------------

The ``--setup`` flag processes lattice structures and calculates optimal Wigner-Seitz radii and starting parameters.

**Flags:**

- ``--setup, -S``: Process lattice and generate optimized atomic parameters.
- ``--cluster``: Cluster type: ``periodic`` (default) or ``free``.
- ``--empty_spheres``: Enable empty sphere calculations.
- ``--min_empty_radius``: Minimum empty sphere radius (default: 1.2).
- ``--max_empty_radius``: Maximum empty sphere radius (default: 2.5).
- ``--supercell NX NY NZ``: Repeat the input unit cell along the lattice vectors.
- ``--validate``: Validate output files after processing.
- ``--summary``: Print a summary of file changes after processing.
- ``--backup``: Create .bak backups of input files before processing.
- ``--in-place``: Replace input files with generated outputs (creates .bak backups).
- ``--lmx``: Override maximum l-channel used in setup path (default: ``2``; supported range ``0..3``).
- ``--lattice-dist``: Gaussian std.dev for lattice-vector distortion in units of ``alat``.
- ``--basis-dist``: Gaussian std.dev for basis-vector distortion in units of ``alat``.
- ``--alat-scale``: Scalar factor applied to ``alat`` before distortion.
- ``--seed``: Integer seed for reproducible random distortions.
- ``--2d``: Enable 2D slab mode for setup/export geometry handling.
- ``--z-iso-factor``: Surface-margin factor for density-aware 2D ``Lz`` (requires ``--2d``).
- ``--rec-method``: Override ``&control recur`` in ``input_new.nml`` (``block``, ``chebyshev``, ``lanczos``).
- ``--rec-steps``: Override ``&control lld`` in ``input_new.nml``.
- ``--sqs``: Enable SQS assignment on all sites, e.g. ``"Fe=0.5,Ni=0.5"`` (requires ``--supercell``).
- ``--sqs-shells``: Number of pair-distance shells in the SQS objective.
- ``--sqs-steps``: Number of MC swap steps for SQS optimization.
- ``--sqs-seed``: Integer seed for reproducible SQS optimization.

**Example:**

.. code-block:: bash

    # Basic lattice processing
    rsutil --setup

    # With empty spheres
    rsutil --setup --empty_spheres --min_empty_radius 1.0 --max_empty_radius 3.0

    # Free cluster calculation
    rsutil --setup --cluster free

    # Create a 1x1x2 supercell with broken symmetry labels
    rsutil --setup --supercell 1 1 2

    # Supercell plus empty spheres (use in-place updates and validation)
    rsutil --setup --supercell 1 1 2 --empty_spheres --in-place --validate

    # Distort lattice/basis vectors reproducibly and scale alat
    rsutil --setup --lattice-dist 0.01 --basis-dist 0.02 --alat-scale 1.01 --seed 1234

    # 2D mode with simple z-rule
    rsutil --setup --2d

    # 2D mode with isotropic-density margins
    rsutil --setup --2d --z-iso-factor 0.5

    # Override recursion controls in &control
    rsutil --setup --rec-method chebyshev
    rsutil --setup --rec-method lanczos --rec-steps 64

    # SQS on all supercell sites
    rsutil --setup --supercell 2 2 2 --sqs "Fe=0.5,Ni=0.3,Co=0.2" \
           --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234

Supercell Notes
~~~~~~~~~~~~~~~

- Supercells are generated with broken symmetry indexing (arrays run 1..N).
- Labels in input.nml are expanded into unique species (e.g., Fe1, N1, Fe2, N2).
- ntype and ct(i) are expanded to match the unique labels.
- Per-species .nml files are created by copying from the original label files.
- By default, ``--supercell`` is staged internally and original input files are restored after setup.
- Use ``--in-place`` if you explicitly want to replace ``lattice.nml``/``input.nml``.
- In ``--in-place`` mode, backups with ``.bak`` extension are created first.
- The same staged supercell behavior also applies when using ``--lmto47``.

2D Notes
~~~~~~~~

- ``--2d`` modifies setup/export geometry before WS/ASA steps.
- Default 2D z-length (when ``--z-iso-factor`` is omitted):
  - ``Lz = (zmax-zmin) + 1.0*alat``.
- Isotropic-density mode:
  - ``Lz = (zmax-zmin) + X*sqrt(A/n_bottom) + X*sqrt(A/n_top)`` with ``X = --z-iso-factor``.
- The third lattice vector is built from the surface normal, not forced to global z:
  - ``a3 || cross(a1,a2)``, sign chosen so ``a3·z >= 0``.
- In setup output, 2D mode enforces ``b3=F`` and ``n3=1`` in ``input_new.nml``.
- ``--z-iso-factor`` is only valid together with ``--2d``.

Recursion Control Notes
~~~~~~~~~~~~~~~~~~~~~~~

- ``--rec-method`` maps to ``recur`` in the ``&control`` namelist.
- ``--rec-steps`` maps to ``lld`` in the ``&control`` namelist.
- Defaults:
  - method: ``block``
  - steps: ``20``
- Special default:
  - if ``--rec-method chebyshev`` is selected and ``--rec-steps`` is omitted, steps default to ``200``.

Setup Safety Workflow
~~~~~~~~~~~~~~~~~~~~~

For robust setup runs, prefer:

.. code-block:: bash

    # dry-style run with checks
    rsutil --setup --validate --summary

    # if happy with output, replace files and keep backups
    rsutil --setup --in-place --backup

Extended Lattice Key Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The setup reader/writer supports selected extended ``&lattice`` keys used by RS-LMTO variants:

- ``strux_backend``
- ``strux_want_sdot``
- ``strux_solve_scale``

These keys are only written to output if present in the input.


Visualization
--------------

The ``--lattice`` and ``--clust`` flags launch an interactive 3D visualizer for structures.

**Flags:**

- ``--lattice, -L``: Run visualizer for lattice input.
- ``--clust, -C``: Run visualizer for clust file.
- ``--data``: Directory containing data files (default: current directory).

**Examples:**

.. code-block:: bash

    # Visualize lattice structure
    rsutil --lattice --data /path/to/data

    # Visualize clust file
    rsutil --clust --data /path/to/data

Visualizer Controls
~~~~~~~~~~~~~~~~~~~

- Mouse/trackpad: rotate, pan, and zoom the scene (PyVista defaults)
- Toggles: show/hide spins, bonds, and unit cell guides
- Colors: CPK element colors are used by default (see ``rsutil/viz/cpk_colors.json``)
- Snapshots: save screenshots from the window menu (platform dependent)

Data Conventions
~~~~~~~~~~~~~~~~

- Expected files in ``--data``: ``lattice.nml`` and ``input.nml`` for lattice view; ``clust`` for cluster view
- Magnetic moments and bonds are derived from RS-LMTO inputs when available
- Periodic vs free cluster visualization follows the underlying input setup


Exchange Interactions
----------------------

The ``--exchange`` flag creates exchange interaction input files (J_ij) for spin dynamics simulations.

**Flags:**

- ``--exchange, -J``: Create exchange input files.
- ``--cutoff``: Cutoff distance for exchange interactions in lattice units (default: 1.0).
- ``--unique``: Find only unique neighbors per shell (one representative per unique species-vector combination).
- ``--shell-width``: Width of shells for grouping neighbors when using ``--unique`` (default: 0.1).
- ``--soc``: Account for spin-orbit coupling in symmetry considerations.
- ``--data``: Directory for input files (default: current directory).
- ``--rij-out``: Path to write R_ij neighbor mappings (optional). Use with ``--exchange``.
- ``--rij-format``: Format for R_ij output; supported values: ``json`` (machine readable) and
    ``txt`` (human readable, fixed-width columns). Default: ``json``.

**Example:**

.. code-block:: bash

    # Basic exchange calculation
    rsutil --exchange --cutoff 2.5

    # With unique neighbor selection and SOC
    rsutil --exchange --cutoff 3.0 --unique --soc --shell-width 0.15

    # From a specific directory
    rsutil --exchange --data /path/to/data --cutoff 2.0

    # Write R_ij mapping to a JSON file
    rsutil --exchange --cutoff 3.0 --rij-out rij_map.json --rij-format json

    # Write human-readable fixed-width R_ij table
    rsutil --exchange --cutoff 3.0 --rij-out rij_map.txt --rij-format txt --unique


Export to UppASD
----------------

The ``--export`` flag exports RS-LMTO inputs to the UppASD (Uppsala Atomistic Spin Dynamics) simulation format. This includes structural information, exchange interactions, and optionally DM vectors and spin-orbit coupling.

**Flags:**

- ``--export, -E``: Export RS-LMTO inputs to UppASD format.
- ``--data``: Directory containing input files and exchange data (default: current directory).
- ``--out``: Output directory for exported files relative to ``--data`` (default: ``uppasd_out``).
- ``--jij``: Path to ``jijso.out`` file containing exchange interactions.
- ``--dmi``: Path to ``dijso.out`` or ``dmi`` file containing DM vector interactions.

**Example:**

.. code-block:: bash

    # Basic export
    rsutil --export --data /path/to/simulation

    # With exchange and DM data
    rsutil --export --data /path/to/simulation \\
           --jij /path/to/jijso.out \\
           --dmi /path/to/dmi \\
           --out my_uppasd_output


Export to LMTO47 CTRL
---------------------

The ``--lmto47`` flag exports RS-LMTO structure/input data into an LMTO47 ``CTRL`` file.

**Flags:**

- ``--lmto47``: Enable LMTO47 CTRL export.
- ``--data``: Directory containing ``lattice.nml`` and ``input.nml``.
- ``--out``: Output directory for generated ``CTRL`` file.
- ``--nspin``: ``NSPIN`` value in LMTO47 output (default: ``2``).
- ``--lmx``: ``LMX`` value written for each LMTO47 ``CLASS`` entry.
- ``--lattice-dist``: Gaussian std.dev for lattice-vector distortion in units of ``alat``.
- ``--basis-dist``: Gaussian std.dev for basis-vector distortion in units of ``alat``.
- ``--alat-scale``: Scalar factor applied to ``alat`` before distortion.
- ``--seed``: Integer seed for reproducible random distortions.
- ``--supercell NX NY NZ``: Build CTRL from a repeated supercell input.
- ``--sqs``: Enable SQS assignment on all sites (requires ``--supercell``).
- ``--sqs-shells`` / ``--sqs-steps`` / ``--sqs-seed``: SQS optimization controls.

**Example:**

.. code-block:: bash

    # Basic LMTO47 export
    rsutil --lmto47 --data /path/to/simulation --out lmto47_out

    # Explicit NSPIN/LMX
    rsutil --lmto47 --data /path/to/simulation --nspin 2 --lmx 3 --out lmto47_out

    # Distorted/scaled LMTO47 export with reproducible seed
    rsutil --lmto47 --data /path/to/simulation --lattice-dist 0.01 --basis-dist 0.02 \
           --alat-scale 1.01 --seed 1234 --out lmto47_out

    # Supercell-aware LMTO47 export (non-destructive unless --in-place is set)
    rsutil --lmto47 --data /path/to/simulation --supercell 2 2 2 --out lmto47_out

    # LMTO47 export from SQS supercell
    rsutil --lmto47 --data /path/to/simulation --supercell 2 2 2 \
           --sqs "Fe=0.5,Ni=0.5" --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234 \
           --out lmto47_out


LMAX Semantics
--------------

The ``--lmx`` option has two contexts:

- Setup context (``--setup``): overrides RS-LMTO generation ``lmax``.
- LMTO47 context (``--lmto47``): writes LMTO47 ``LMX`` in ``CLASS`` entries.

Current RS-LMTO backend support in this repository is ``0 <= lmax <= 3``.


Quick Reference
----------------

.. code-block:: bash

    # Show all available options
    rsutil --help

    # Typical workflow: setup → exchange → export
    rsutil --setup --empty_spheres
    rsutil --exchange --unique --cutoff 2.5
    rsutil --export --jij jijso.out
