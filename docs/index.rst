
.. RS-utils documentation master file, created by
   sphinx-quickstart on Tue Feb 11 14:55:38 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RS-utils's documentation!
====================================

RS-utils is a comprehensive toolkit for RS-LMTO-ASA calculations. It processes lattice and atomic inputs to compute optimal Wigner–Seitz radii and reasonable starting parameters, builds exchange interaction maps for spin dynamics, provides interactive 3D visualization of structures and magnetic moments, and exports results to UppASD for further simulation.

It now includes an advanced pure-Python band structure and DOS workflow with orbital/site/spin projections, gauge fixing for smooth eigenvectors, degenerate band stabilization, Brillouin zone visualization, and high-order quantities like Berry curvature and effective mass.

Getting Started
===============

**RS-LMTO setup workflow:**

.. code-block:: bash

    # 1) Setup lattice and initial atomic parameters
    rsutil --setup --empty_spheres

    # 2) Build exchange interactions (J_ij)
    rsutil --exchange --cutoff 2.5 --unique --shell-width 0.1

    # 3) Export to UppASD
    rsutil --export --data . --jij jijso.out --out uppasd_out

    # Optional: all-site SQS on supercell before setup
    rsutil --setup --supercell 2 2 2 --sqs "Fe=0.5,Ni=0.5" \
        --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234

    # Optional: 2D slab mode
    rsutil --setup --2d --z-iso-factor 0.5

    # Optional: recursion control overrides
    rsutil --setup --rec-method chebyshev

**Band structure & DOS workflow:**

.. code-block:: bash

   # Band structure with automatic k-path (fast)
   rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

   # Spin-projected bands with Sz expectation values
   rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
       --projections --projection spin --projections-mode expectation --output bands_spin.png

   # Orbital-projected bands (s, p, d character)
   rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
       --projections --projection orbital --output bands_orbital.png

   # Full Brillouin zone DOS with Fermi level determination
   rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat \
       --mesh 6 6 6 --electrons 16 --moments --output dos.png

**Visualization:**

.. code-block:: bash

    # Lattice view  
    rsutil --lattice --data .

    # Cluster view
    rsutil --clust --data .

.. toctree::
   :maxdepth: 1
   :caption: Command-Line Interface:

   cli_section
   band_structure_guide

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   rsutil_module
   bands_submodule
   viz_submodule
   exchange_submodule
   export_submodule

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
