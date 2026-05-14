.. _rsutil_module:

The `rsutil` module
-------------------

The **`rsutil`** module is the core of the RS-utils package. It provides utilities for processing lattice structures in RS-LMTO-ASA calculations, including:

- **Lattice Processing**: The main ``process_lattice()`` function calculates optimal Wigner-Seitz radii and generates starting parameters for atomic calculations.
- **Disorder Workflow**: Supercell construction plus optional all-site SQS assignment for substitutional alloys.
- **Data Visualization**: Interactive 3D visualization of atomic structures and magnetic properties.
- **Exchange Interactions**: Generation of exchange interaction files (J_ij) for spin dynamics simulations.
- **Format Conversion**: Export to UppASD format and other simulation engines.

Core Functions
~~~~~~~~~~~~~~

.. automodule:: rsutil
   :members:

Notes:

- ``process_lattice(..., lmx=2)`` exposes setup-side ``lmax`` override from Python.
- ``process_lattice`` also supports ``lattice_dist``, ``basis_dist``, ``alat_scale``, and ``seed`` for reproducible geometric perturbations.
- ``process_lattice`` supports 2D slab mode via ``mode_2d`` and optional ``z_iso_factor``.
- Current backend support in this repository is ``0 <= lmax <= 3``.
- Extended lattice keys ``strux_backend``, ``strux_want_sdot``, and ``strux_solve_scale`` are preserved when present.
- CLI-level SQS controls (``--sqs``, ``--sqs-shells``, ``--sqs-steps``, ``--sqs-seed``) can be combined with ``--supercell`` for both ``--setup`` and ``--lmto47`` workflows.

Submodules
~~~~~~~~~~

The package is organized into several submodules:

- **``rsutil.bands``**: Band structure, DOS, projections, DOS moments
- **``rsutil.viz``**: Visualization tools for structures and simulations
- **``rsutil.exchange``**: Exchange interaction calculation and mapping
- **``rsutil.export``**: Export functionality for other simulation formats

- **``rsutil.cli``**: Command-line entry points and parsers for `rsutil` console commands


Command-line Interface
~~~~~~~~~~~~~~~~~~~~~~

The package exposes console-oriented entrypoints and helpers in `rsutil.cli` for running `bands`, `dos` and visualization workflows from the shell.

.. automodule:: rsutil.cli
   :members:

Exporters
~~~~~~~~~

Export utilities (e.g., UppASD exporters) are available under `rsutil.export`.

.. automodule:: rsutil.export.uppasd
   :members:
   :noindex:
