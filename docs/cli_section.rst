.. _cli_section:

Command-Line Interface (CLI)
=============================

The **RS-utils CLI** provides a powerful command-line interface for interacting with the `rsutil` package. It allows users to process lattice structures, visualize data, create exchange input files, and export to other simulation formats with ease.

You can use the CLI to run various tasks by providing the corresponding flags. Below is an overview of the available options and usage examples.

**Usage:**

To use the CLI, run the script with the desired flags from your terminal:

.. code-block:: bash

    rsutil [--flag] [options]


Lattice Processing
-------------------

The ``--setup`` flag processes lattice structures and calculates optimal Wigner-Seitz radii and starting parameters.

**Flags:**

- ``--setup, -S``: Process lattice and generate optimized atomic parameters.
- ``--cluster``: Cluster type: ``periodic`` (default) or ``free``.
- ``--empty_spheres``: Enable empty sphere calculations.
- ``--min_empty_radius``: Minimum empty sphere radius (default: 1.2).
- ``--max_empty_radius``: Maximum empty sphere radius (default: 2.5).

**Example:**

.. code-block:: bash

    # Basic lattice processing
    rsutil --setup

    # With empty spheres
    rsutil --setup --empty_spheres --min_empty_radius 1.0 --max_empty_radius 3.0

    # Free cluster calculation
    rsutil --setup --cluster free


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


Quick Reference
----------------

.. code-block:: bash

    # Show all available options
    rsutil --help

    # Typical workflow: setup → exchange → export
    rsutil --setup --empty_spheres
    rsutil --exchange --unique --cutoff 2.5
    rsutil --export --jij jijso.out

