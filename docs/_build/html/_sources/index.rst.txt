
.. RS-utils documentation master file, created by
   sphinx-quickstart on Tue Feb 11 14:55:38 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RS-utils's documentation!
====================================

RS-utils is a toolkit for setting up RS-LMTO-ASA calculations. It processes lattice and atomic inputs to compute optimal Wigner–Seitz radii and reasonable starting parameters, builds exchange interaction maps for spin dynamics, provides interactive 3D visualization of structures and magnetic moments, and exports results to UppASD for further simulation. It now also includes a pure-Python band structure and DOS workflow with projections and DOS moments.

Getting Started
===============

Typical workflow from RS-LMTO inputs to UppASD:

.. code-block:: bash

    # 1) Setup lattice and initial atomic parameters
    rsutil --setup --empty_spheres

    # 2) Build exchange interactions (J_ij)
    rsutil --exchange --cutoff 2.5 --unique --shell-width 0.1

    # 3) Export to UppASD
    rsutil --export --data . --jij jijso.out --out uppasd_out

Band/DOS workflow:

.. code-block:: bash

   # Band structure (high-symmetry path)
   rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

   # Projected bands (orbital/site/spin)
   rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --projections --projection orbital

   # DOS with Fermi level and moments
   rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 --electrons 16 --moments

Visualization:

.. code-block:: bash

    # Lattice view
    rsutil --lattice --data .

    # Cluster view
    rsutil --clust --data .

.. toctree::
   :maxdepth: 1
   :caption: Command-Line Interface:

   cli_section

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   rsutil_module
   viz_submodule
   exchange_submodule
   exchange_map_module
   export_submodule
   bands_submodule

************************

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
