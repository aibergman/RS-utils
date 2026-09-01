.. _viz_submodule:

The `viz` submodule
-------------------

The **`viz`** submodule contains all visualization-related functionality for RS-utils. It provides interactive 3D visualization of atomic structures, magnetic properties, and crystal lattices using PyVista for advanced rendering.

Key Features
~~~~~~~~~~~~

- Interactive 3D visualization of crystal structures
- Magnetic moment (spin) vectors and glyphs
- Bond visualization with adjustable styling
- Periodic lattice rendering and cluster views
- RS-LMTO input support: ``lattice.nml``, ``input.nml``, and ``clust`` files
- CPK color schemes for element visualization

Core Classes
~~~~~~~~~~~~

.. autoclass:: rsutil.viz.glyphs.Atom
   :members:
   :undoc-members:

.. autoclass:: rsutil.viz.glyphs.Bond
   :members:
   :undoc-members:

.. autoclass:: rsutil.viz.glyphs.Spin
   :members:
   :undoc-members:

Input and Structure
~~~~~~~~~~~~~~~~~~~

.. autoclass:: rsutil.viz.inputreader.InputReader
   :members:
   :undoc-members:

.. autoclass:: rsutil.viz.structure.Structure
   :members:
   :undoc-members:

Main Visualizer
~~~~~~~~~~~~~~~

.. autoclass:: rsutil.viz.visualizer.Visualizer
   :members:
   :undoc-members:

Usage Examples
~~~~~~~~~~~~~~

From the command line:

.. code-block:: bash

    # Visualize lattice structure
    rsutil --lattice --data /path/to/data

    # Visualize cluster definition
    rsutil --clust --data /path/to/data

From Python code:

.. code-block:: python

    from rsutil.viz.rs_visualizer import rs_visualizer

    # Launch interactive visualizer for lattice
    rs_visualizer(data_dir='/path/to/data', clust=False)

    # Launch visualizer for clust file
    rs_visualizer(data_dir='/path/to/data', clust=True)
