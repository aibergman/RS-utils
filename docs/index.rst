
.. RS-utils documentation master file, created by
   sphinx-quickstart on Tue Feb 11 14:55:38 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RS-utils's documentation!
====================================

RS-utils is a Python package designed to assist with various utilities and tools for scientific computing, particularly in the context of simulations, structure visualization, and data processing. This documentation will guide you through the features and modules provided by the package.

Contents
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. _cli_section:

Command-Line Interface (CLI)
=============================

The **RS-utils CLI** provides a powerful command-line interface for interacting with the `rsutil` package. It allows users to process lattice structures, visualize data, and create exchange input files with ease.

You can use the CLI to run various tasks by providing the corresponding flags. Below is an overview of the available options and usage examples.

**Usage:**
To use the CLI, run the script with the desired flags from your terminal:

  ```bash
  rsutil <--flag1> <--flag2> ...
  ```



**Available Flags:**
- `--setup, -S`: Process lattice.

- `--lattice, -L`: Run visualizer for lattice input.

- `--clust, -C`: Run visualizer for clust file.

- `--exchange, -J`: Create exchange input files.

- `--data`: Directory for data files to be visualized (default: current directory).

- `--cutoff`: Cutoff for exchange interactions, in lattice units.


**Example:**
- To process a lattice:

  ```bash
  rsutil --setup
  ```

- To visualize a lattice input:
  
  ```bash
  rsutil --lattice --data /path/to/data
  ```

- To create exchange input files:
  
  ```bash
  rsutil --exchange --cutoff 1.5
  ```

For more detailed descriptions of each option, you can refer to the `cli.py` module directly.


************************

.. _modules:

Modules
=============

The RS-utils package consists of several modules that provide different functionalities for scientific computing and data visualization. This section provides an overview of the main modules and their features.

The `rsutil` module
-------------------
The **`rsutil`** module is the core of the package. It provides various utility functions and base classes that can be used to process and analyze atomic structures, perform data processing, and handle other aspects of scientific computing.

.. automodule:: rsutil
   :members:


************************

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

************************

The `exchange` submodule
-------------------------
The **`exchange`** submodule provides utilities for handling data exchange between different parts of a simulation or processing pipeline. This includes data mapping, input/output handling, and other exchange-related functions.

.. automodule:: rsutil.exchange
   :members:

The `exchange_map` module
-------------------------
The **`exchange_map`** module provides mapping functions that facilitate the transfer of data between different representations or formats. This is especially useful for handling large datasets or converting between different simulation outputs.

.. automodule:: rsutil.exchange.exchange_map
   :members:


************************

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
