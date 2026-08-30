.. _exchange_submodule:

The `exchange` submodule
------------------------

The **`exchange`** submodule provides utilities for calculating and managing exchange interactions (magnetic couplings) between atoms in a crystal structure. This is essential for setting up spin dynamics simulations.

For most users, this page is the primary exchange documentation.  
If you need implementation-level API details, see :ref:`exchange_map_module`.

**Key Features:**

- **Exchange Mapping**: Identifies and calculates magnetic exchange interactions (J_ij) between atomic sites
- **Neighbor Analysis**: Automatically finds interacting neighbors within a specified cutoff distance
- **Symmetry Handling**: Uses exact periodic pair-key operations and conservative SOC filtering
- **Unique Neighbor Selection**: Option to identify one representative per exact periodic orbit
- **Shell Grouping**: Shell width is retained as a compatibility/presentation option, not an equivalence rule

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

    # Advanced: exact unique neighbors with the legacy SOC reference
    rsutil --exchange --cutoff 3.0 --unique --soc --shell-width 0.15

    # Full pair inventory for validation
    rsutil --exchange --cutoff 3.0 --unique --full-pairs

    # From a specific directory
    rsutil --exchange --data /path/to/data --cutoff 2.0
