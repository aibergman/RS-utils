.. _core_submodule:

The ``core`` submodule
======================

:mod:`rsutil.core` contains the dependency-light contracts shared by setup,
exchange, result parsing, campaign generation, and ASE interoperability.  The
models use stable zero-based site and pair identities and provide explicit
serialization for machine-readable artifacts.

Core contracts
--------------

``RSStructure`` is the canonical structure representation.  It stores three
lattice vectors, site positions, chemical identity, optional RS labels/type
indices, coordinate convention, periodic-boundary flags, and metadata.
``SetupSpec`` is the normalized parameter contract for the legacy ``-S`` setup
workflow, including supercells, SQS, 2D geometry, distortions, and recursion
controls.  ``PairKey(i, j, R)`` identifies a periodic interaction from site
``i`` to site ``j`` in translation ``R``.

``OperationResult``, ``SetupResult``, ``ValidationResult``, ``SetupReport``,
and ``RunManifest`` carry diagnostics, output inventory, hashes, status, and
provenance without requiring a scheduler or database.  ``ExchangePair`` is
the canonical exchange-pair record used by job generation and result parsing.

Hashing and filesystem primitives
---------------------------------

The hashing helpers provide canonical JSON, ordinary content hashes, and
scientific hashes.  The path helpers create owned temporary/staging workspaces
and perform atomic replacement; they are the basis of the non-destructive HT
setup and campaign contracts.

Public API
----------

.. automodule:: rsutil.core.contracts
   :members:
   :show-inheritance:

.. automodule:: rsutil.core.hashing
   :members:

.. automodule:: rsutil.core.paths
   :members:
   :show-inheritance:

.. automodule:: rsutil.core.errors
   :members:
   :show-inheritance:
