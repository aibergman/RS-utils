.. _rsutil_module:

The `rsutil` module
-------------------

The **`rsutil`** module is the core of the RS-utils package. It provides utilities for processing lattice structures in RS-LMTO-ASA calculations, including:

- **Lattice Processing**: The main ``process_lattice()`` function calculates optimal Wigner-Seitz radii and generates starting parameters for atomic calculations.
- **Data Visualization**: Interactive 3D visualization of atomic structures and magnetic properties.
- **Exchange Interactions**: Generation of exchange interaction files (J_ij) for spin dynamics simulations.
- **Format Conversion**: Export to UppASD format and other simulation engines.

Core Functions
~~~~~~~~~~~~~~

.. automodule:: rsutil
   :members:

Submodules
~~~~~~~~~~

The package is organized into several submodules:

- **``rsutil.bands``**: Band structure, DOS, projections, DOS moments
- **``rsutil.viz``**: Visualization tools for structures and simulations
- **``rsutil.exchange``**: Exchange interaction calculation and mapping
- **``rsutil.export``**: Export functionality for other simulation formats

