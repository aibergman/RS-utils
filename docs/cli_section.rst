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
- To preprocess a lattice:

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

