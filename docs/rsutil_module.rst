.. _rsutil_module:

The ``rsutil`` module
=====================

The top-level ``rsutil`` module is the compatibility-facing entry point for
RS-utils.  It retains the historical Fortran setup wrapper and exposes the
high-throughput services and analysis packages:

- **Lattice processing**: ``process_lattice()`` remains available for legacy
  in-process setup calls.
- **Safe setup**: :mod:`rsutil.setup` provides isolated execution, validation,
  manifests, supercells, SQS, distortions, and 2D handling.
- **HT campaigns**: :mod:`rsutil.campaign` validates and materializes
  deterministic campaigns, renders Slurm helpers, and collects evidence.
- **Results**: :mod:`rsutil.results` parses SCF, exchange, transport/damping,
  and tight-binding data into versioned contracts.
- **Analysis and interoperability**: bands/DOS, visualization, exchange
  mapping, ASE conversion, and UppASD/LMTO47 export.
- **Surface preparation**: :func:`rsutil.prepare_surface` builds an RS-LMTO
  ``buildsurf`` input bundle from a converged bulk calculation.

Core Functions
~~~~~~~~~~~~~~

.. autofunction:: rsutil.process_lattice

Notes:

- ``process_lattice(..., lmx=2)`` exposes setup-side ``lmax`` override from Python.
- ``process_lattice`` also supports ``lattice_dist``, ``basis_dist``, ``alat_scale``, and ``seed`` for reproducible geometric perturbations.
- ``process_lattice`` supports 2D slab mode via ``mode_2d`` and optional ``z_iso_factor``.
- ``process_lattice`` supports recursion overrides via ``rec_method`` and ``rec_steps`` (written to ``&control`` as ``recur`` and ``lld``).
- Current backend support in this repository is ``0 <= lmax <= 3``.
- Extended lattice keys ``strux_backend``, ``strux_want_sdot``, and ``strux_solve_scale`` are preserved when present.
- CLI-level SQS controls (``--sqs``, ``--sqs-shells``, ``--sqs-steps``, ``--sqs-seed``) can be combined with ``--supercell`` for both ``--setup`` and ``--lmto47`` workflows.

Surface Preparation API
~~~~~~~~~~~~~~~~~~~~~~~~

The Python API mirrors ``rsutil surface prepare``. It writes to a new output
directory and returns a manifest describing the staged files and settings:

.. autoclass:: rsutil.SurfaceSpec
   :members:

.. autofunction:: rsutil.prepare_surface

Submodules
~~~~~~~~~~

The package is organized into the following submodules:

- **``rsutil.core``**: Shared data contracts, hashing, and safe paths
- **``rsutil.setup``**: Setup runner, service, and validation APIs
- **``rsutil.campaign``**: Deterministic HT campaign APIs
- **``rsutil.results``**: Versioned result models and parsers
- **``rsutil.ase``**: Optional ASE interoperability
- **``rsutil.bands``**: Band structure, DOS, projections, DOS moments
- **``rsutil.viz``**: Visualization tools for structures and simulations
- **``rsutil.exchange``**: Exchange interaction calculation and mapping
- **``rsutil.export``**: Export functionality for other simulation formats

- **``rsutil.cli``**: Command-line entry points and parsers for ``rsutil``
  console commands


Command-line Interface
~~~~~~~~~~~~~~~~~~~~~~

The package exposes console-oriented entrypoints and helpers in `rsutil.cli` for running `bands`, `dos` and visualization workflows from the shell.

.. automodule:: rsutil.cli
   :members:

Exporters
~~~~~~~~~

Export utilities (e.g., UppASD exporters) are available under
``rsutil.export``.

.. automodule:: rsutil.export.uppasd
   :members:
   :noindex:
