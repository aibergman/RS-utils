.. _exchange_map_module:

The `exchange_map` module
-------------------------

The **`exchange_map`** module provides the core functionality for mapping and calculating magnetic exchange interactions between atoms in a crystal structure.

**Overview:**

The exchange mapping process:

1. Reads the lattice structure from ``lattice.nml``
2. Identifies all atomic neighbors within the specified cutoff distance
3. Calculates exchange interaction strengths (J_ij values)
4. Optionally selects exact periodic-symmetry representatives
5. Applies conservative axial-vector filtering when ``--soc`` is used
6. Generates output files suitable for spin dynamics simulations

**Key Parameters:**

- **cutoff**: Maximum distance (in lattice units) for considering exchange interactions
- **unique**: If ``True``, identifies one representative per exact periodic pair orbit
- **shell_width**: Retained for compatibility; it does not define pair equivalence
- **soc**: If ``True``, applies the documented legacy collinear +z magnetic-reference policy

**Functions:**

.. automodule:: rsutil.exchange.exchange_map
   :members:

**Example Usage:**

.. code-block:: python

    from rsutil.exchange.exchange_map import exchange_map

    # Standard calculation: all neighbors up to 2.5 lattice units
    exchange_map(cutoff=2.5)

    # Advanced: find exact periodic-symmetry representatives
    exchange_map(
        cutoff=3.0,
        unique=True,
        shell_width=0.1  # Accepted for legacy compatibility only
    )

    # With spin-orbit coupling considerations
    exchange_map(
        cutoff=2.5,
        soc=True  # Use the conservative legacy +z reference policy
    )

**Output:**

The module generates exchange interaction files in formats compatible with:

- RS-LMTO-ASA simulations
- UppASD (Uppsala Atomistic Spin Dynamics)
- Other spin dynamics packages
