.. _exchange_map_module:

The `exchange_map` module
-------------------------

The **`exchange_map`** module provides the core functionality for mapping and calculating magnetic exchange interactions between atoms in a crystal structure.

**Overview:**

The exchange mapping process:

1. Reads the lattice structure from ``lattice.nml``
2. Identifies all atomic neighbors within the specified cutoff distance
3. Calculates exchange interaction strengths (J_ij values)
4. Optionally groups neighbors into shells and selects unique representatives
5. Accounts for symmetry operations when ``--soc`` flag is used
6. Generates output files suitable for spin dynamics simulations

**Key Parameters:**

- **cutoff**: Maximum distance (in lattice units) for considering exchange interactions
- **unique**: If ``True``, identifies only unique neighbors per shell (reduces output size)
- **shell_width**: Tolerance for grouping neighbors into shells by distance
- **soc**: If ``True``, includes spin-orbit coupling effects in symmetry considerations

**Functions:**

.. automodule:: rsutil.exchange.exchange_map
   :members:

**Example Usage:**

.. code-block:: python

    from rsutil.exchange.exchange_map import exchange_map

    # Standard calculation: all neighbors up to 2.5 lattice units
    exchange_map(cutoff=2.5)

    # Advanced: find unique representatives of each neighbor shell
    exchange_map(
        cutoff=3.0,
        unique=True,
        shell_width=0.1  # Group neighbors within 0.1 lattice units
    )

    # With spin-orbit coupling considerations
    exchange_map(
        cutoff=2.5,
        soc=True  # Account for SOC in symmetry
    )

**Output:**

The module generates exchange interaction files in formats compatible with:

- RS-LMTO-ASA simulations
- UppASD (Uppsala Atomistic Spin Dynamics)
- Other spin dynamics packages

