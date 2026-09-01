.. _results_submodule:

The ``results`` submodule
=========================

:mod:`rsutil.results` provides versioned, unit-aware result contracts for
RS-LMTO runs.  Parsers retain source provenance and diagnostics instead of
silently guessing when a file is malformed, ambiguous, unsupported, or uses
an unestablished unit.

SCF and exchange results
------------------------

Use :func:`parse_run` or :func:`parse_report` for SCF evidence.  An
:class:`SCFResult` contains convergence state, energies, Fermi energy, site
moments, validation diagnostics, and source-file hashes.  Exchange tables are
parsed with :func:`parse_exchange_results` and represented as an
:class:`ExchangeSet` of :class:`InteractionRecord` objects keyed by canonical
``PairKey`` values.  Scalar ``J``, DMI ``D``, anisotropic ``Aij``, and the
derived full tensor remain separate fields; source units are preserved.

Transport and damping
---------------------

The transport boundary is explicit.  Stable supported conductivity and
damping records use :class:`ConductivityResult` and :class:`DampingResult`;
ambiguous or unsupported source formats return an
:class:`UnsupportedTransportResult` with the reason and provenance rather than
inventing a parser contract.  See :doc:`HT_TRANSPORT_DAMPING_FORMATS` for the
format audit and support boundary.

Tight binding and ML export
---------------------------

:class:`TightBindingModel` is a typed sparse/list representation of a real-space
Hamiltonian.  It supports the versioned ``fable_v4`` metadata/hopping files,
legacy reader conversion, ``H(k)`` reconstruction, eigenvalues, and a
versioned NPZ/JSON ML export.  Spin/SOC metadata and ``control%nsp`` semantics
are retained; unsupported representations fail explicitly.

Units
-----

:class:`Quantity` attaches an explicit :class:`Unit` to scalar and structured
values.  Parsing does not convert source values.  Conversion is opt-in and is
currently defined for the Ry/mRy pair used by interaction outputs.

Public API
----------

.. automodule:: rsutil.results.models
   :members:
   :show-inheritance:

.. automodule:: rsutil.results.parser
   :members:

.. automodule:: rsutil.results.exchange
   :members:

.. automodule:: rsutil.results.transport
   :members:
   :show-inheritance:

.. automodule:: rsutil.results.tight_binding
   :members:
   :show-inheritance:

.. automodule:: rsutil.results.units
   :members:
   :show-inheritance:

.. automodule:: rsutil.results.errors
   :members:
   :show-inheritance:
