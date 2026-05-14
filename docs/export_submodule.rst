.. _export_submodule:

The `export` submodule
----------------------

The **`export`** submodule provides functionality to convert RS-LMTO-ASA simulation data into formats compatible with other simulation packages. It currently supports export to:

- UppASD (Uppsala Atomistic Spin Dynamics)
- LMTO47 CTRL format

**Features:**

- Export RS-LMTO lattice and atomic structures to UppASD
- Include exchange interactions (J_ij) in the exported files
- Optional DM vector (Dzyaloshinsky-Moriya) interactions
- Spin-orbit coupling support
- Automatic handling of different cluster types (periodic and free)
- LMTO47 CTRL generation with configurable ``NSPIN`` and ``LMX``

**Submodules:**

.. _uppasd_export:

UppASD Export
~~~~~~~~~~~~~

The ``uppasd`` module handles conversion to the UppASD format, which is widely used for spin dynamics simulations.

.. automodule:: rsutil.export.uppasd
   :members:

**Usage:**

.. code-block:: python

    from rsutil.export.uppasd import export_from_files

    # Export with lattice data and exchange interactions
    export_from_files(
        lattice_nml='path/to/lattice.nml',
        input_nml='path/to/input.nml',
        clust='path/to/clust',
        out_dir='output',
        jij_file='path/to/jijso.out',
        dij_file='path/to/dmi'
    )

**From the Command Line:**

.. code-block:: bash

    # Basic export
    rsutil --export --data /path/to/simulation

    # With exchange and DM interactions
    rsutil --export --data /path/to/simulation \
           --jij /path/to/jijso.out \
           --dmi /path/to/dmi \
           --out uppasd_output

.. _lmto47_export:

LMTO47 CTRL Export
~~~~~~~~~~~~~~~~~~

The ``lmto47`` module exports RS-LMTO lattice/input data to an LMTO47 ``CTRL`` file.

.. automodule:: rsutil.export.lmto47
   :members:

**From the Command Line:**

.. code-block:: bash

    # Basic CTRL export
    rsutil --lmto47 --data /path/to/simulation --out lmto47_out

    # Explicit spin and l-channel settings
    rsutil --lmto47 --data /path/to/simulation --nspin 2 --lmx 3 --out lmto47_out

    # 2D mode with default z-rule
    rsutil --lmto47 --data /path/to/simulation --2d --out lmto47_out

    # 2D mode with isotropic-density margins
    rsutil --lmto47 --data /path/to/simulation --2d --z-iso-factor 0.5 --out lmto47_out

Notes:

- ``--nspin`` controls ``NSPIN`` in the LMTO47 ``OPTIONS`` block.
- ``--lmx`` writes ``LMX`` for all ``CLASS`` entries in CTRL export.
- ``--lattice-dist`` and ``--basis-dist`` apply Gaussian distortions in units of ``alat``.
- ``--alat-scale`` scales ``alat`` before distortion.
- ``--seed`` makes distortions reproducible.
- ``--2d`` enables slab-oriented z-length derivation.
- ``--z-iso-factor`` switches 2D z-length from fixed ``+1.0*alat`` to density-aware margins.
