.. _exchange_submodule:

The `exchange` submodule
------------------------

The **`exchange`** submodule provides utilities for calculating and managing exchange interactions (magnetic couplings) between atoms in a crystal structure. This is essential for setting up spin dynamics simulations.

For most users, this page is the primary exchange documentation.  
If you need implementation-level API details, see :ref:`exchange_map_module`.

**Key Features:**

- **Exchange Mapping**: Identifies and calculates magnetic exchange interactions (J_ij) between atomic sites
- **Neighbor Analysis**: Automatically finds interacting neighbors within a specified cutoff distance
- **Symmetry Handling**: Supports symmetry operations including spin-orbit coupling (SOC) effects
- **Unique Neighbor Selection**: Option to identify only unique magnetic neighbors per shell to reduce computational complexity
- **Shell Grouping**: Groups neighbors by distance with configurable shell width for organized output

**Main Functions:**

.. automodule:: rsutil.exchange
   :members:

**Usage:**

Basic exchange interaction calculation:

.. code-block:: python

    from rsutil.exchange.exchange_map import exchange_map

    # Calculate exchange interactions with 2.5 Å cutoff
    exchange_map(cutoff=2.5)

    # With unique neighbor selection and SOC
    exchange_map(
        cutoff=3.0,
        unique=True,
        shell_width=0.15,
        soc=True
    )

**From the Command Line:**

.. code-block:: bash

    # Basic exchange calculation
    rsutil --exchange --cutoff 2.5

    # Advanced: unique neighbors with SOC
    rsutil --exchange --cutoff 3.0 --unique --soc --shell-width 0.15

    # From a specific directory
    rsutil --exchange --data /path/to/data --cutoff 2.0
