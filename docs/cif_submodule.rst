.. _cif_submodule:

CIF conversion
==============

``rsutil.cif`` converts ordered CIF files into RS-LMTO ``lattice.nml`` and
``input.nml``. See :doc:`HT_CIF` for the conventions and :doc:`cli_section` for
the command-line interface.

.. code-block:: python

    from rsutil.cif import convert_cif

    files, conversion = convert_cif("FeN.cif", "run")   # primitive cell, automatic ct
    print(files.ct, files.neighbors_min, conversion.structure.metadata["centering"])

.. automodule:: rsutil.cif
   :members: read_cif, write_rs_input, convert_cif, suggest_ct, neighbor_counts,
             format_lattice_nml, format_input_nml, CifConversion, RSInputFiles,
             CifConversionError
   :member-order: bysource
