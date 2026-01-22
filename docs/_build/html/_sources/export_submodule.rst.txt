.. _export_submodule:

The `export` submodule
----------------------

The **`export`** submodule provides functionality to convert RS-LMTO-ASA simulation data into formats compatible with other simulation packages. Currently, it supports export to UppASD (Uppsala Atomistic Spin Dynamics) format.

**Features:**

- Export RS-LMTO lattice and atomic structures to UppASD
- Include exchange interactions (J_ij) in the exported files
- Optional DM vector (Dzyaloshinsky-Moriya) interactions
- Spin-orbit coupling support
- Automatic handling of different cluster types (periodic and free)

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
