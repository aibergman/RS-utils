.. _setup_submodule:

Setup services and validation
=============================

RS-utils keeps the numerical Fortran setup kernel behind a small, path-explicit
Python service.  :class:`rsutil.setup_service.SetupService` stages inputs,
invokes :class:`rsutil.setup_runner.FortranSetupRunner` in an isolated worker,
validates generated namelists and outputs, and commits results atomically.
The caller's working directory is not changed.

The compatibility function :func:`rsutil.process_lattice` remains available
for historical in-process use.  It is unsafe for HT orchestration because the
legacy Fortran call graph can terminate the Python interpreter; use
``SetupService``/``run_setup`` for batch and campaign work.

Setup features
--------------

``SetupSpec`` and the CLI support:

* periodic or free clusters and optional empty spheres, using the default
  ``xband`` backend or the revisioned ``questaal`` FINDES backend
  (``empty_sphere_backend="questaal"`` / ``--empty-sphere-backend questaal``);
* supercell generation and seeded all-site SQS assignment;
* reproducible lattice/basis Gaussian distortions and ``alat`` scaling;
* 2D slab geometry with fixed or density-aware z margins;
* ``block``, ``chebyshev``, and ``lanczos`` recursion controls;
* ``lmx`` values from 0 through 3 in the current backend;
* validation reports, output summaries, backups, and explicit in-place commit.

Typical service use
-------------------

.. code-block:: python

    from rsutil.core import SetupSpec
    from rsutil.setup import run_setup

    result = run_setup(
        "source-run",
        "prepared-run",
        SetupSpec(supercell=(2, 2, 2), seed=1234),
    )
    if not result.success:
        for diagnostic in result.diagnostics:
            print(diagnostic.code, diagnostic.message)

The setup manifest records input/output hashes, normalized parameters,
structure identity, software/environment metadata, timestamps, and stochastic
provenance.  It is an evidence record and should be retained with a run.

The Questaal adapter uses identity symmetry, matching the symmetry information
currently passed by RS-utils' xband path.  The vendored source revision,
adaptation boundary, and AGPL provenance are recorded in
``vendor/questaal_findes/README.md``.

Validation API
--------------

The validation functions separate source namelist validation from generated
output validation.  Reports can be written and read as JSON-compatible
objects, making them suitable for campaign status and collection.

.. automodule:: rsutil.setup_service
   :members:
   :show-inheritance:

.. automodule:: rsutil.setup_runner
   :members:
   :show-inheritance:

.. automodule:: rsutil.setup_validation
   :members:
   :show-inheritance:
