.. _ase_module:

ASE interoperability
====================

The optional :mod:`rsutil.ase` adapter converts between ASE ``Atoms`` and the
canonical :class:`rsutil.core.RSStructure`.  Install the extra with
``pip install rsutil[ase]``.  Importing :mod:`rsutil.core` itself does not
require ASE.

The boundary preserves atom order, Cartesian positions, cell vectors, PBC
flags, chemical symbols, and JSON-compatible ``Atoms.info``/per-atom arrays.
RS labels, RS types, and numeric RS type indices use the explicit metadata keys
``rsutil_rs_labels``, ``rsutil_rs_types``, and ``rsutil_rs_type_indices``.
ASE lengths are represented in Å; conversion from ``angstrom``, ``bohr``,
``nm``, and ``pm`` is supported when converting an ``RSStructure`` to ASE.

.. code-block:: python

    from ase import Atoms
    from rsutil.ase import from_ase, to_ase

    atoms = Atoms("Fe2", cell=[2.8, 2.8, 2.8], pbc=True)
    structure = from_ase(atoms)
    restored = to_ase(structure)

ASE-to-setup source generation is an explicitly experimental convenience path.
It writes placeholder LMTO namelists and requires ``experimental=True`` after
reviewing the generated policy; it is not an HT-ready automatic setup policy.
The supported HT ASE contract is the structure conversion pair above.

.. automodule:: rsutil.ase
   :members:
   :show-inheritance:
