.. _vis_submodule:

The `viz` submodule
-------------------
The **`viz`** submodule contains all the visualization-related functionality. It helps to visualize atomic structures, spins, bonds, and other scientific data in a clear and intuitive way. The `viz` module uses libraries like PyVista to provide interactive 3D visualizations.

.. automodule:: rsutil.viz
   :members:

Classes in `viz`
----------------
This section documents the individual classes within the `viz` submodule, which are used to represent atoms, bonds, spins, and other visual components.

- **`Atom`**: Represents an atom in a structure, storing its position and properties for visualization.
- **`Bond`**: Represents a bond between two atoms, including its visual attributes like radius.
- **`Spin`**: Represents a magnetic moment or spin, visualized as an arrow or vector.
- **`InputReader`**: A utility for reading input files and parsing structure data.
- **`Structure`**: Represents a structure consisting of atoms and bonds, and provides methods to manipulate and visualize the structure.
- **`Visualizer`**: The main class for creating interactive visualizations of atomic structures and their properties.

.. autoclass:: rsutil.viz.glyphs.Atom
   :members:

.. autoclass:: rsutil.viz.glyphs.Bond
   :members:

.. autoclass:: rsutil.viz.glyphs.Spin
   :members:

.. autoclass:: rsutil.viz.inputreader.InputReader
   :members:

.. autoclass:: rsutil.viz.structure.Structure
   :members:

.. autoclass:: rsutil.viz.visualizer.Visualizer
   :members:
