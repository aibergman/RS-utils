Usage
=====

Quick start
-----------

Install editable:

.. code-block:: bash

   pip install -e .

Run the CLI:

.. code-block:: bash

   rs-band --help

DOS Moments Analysis
--------------------

The package includes specialized functionality for DOS moments analysis, including Fermi level determination and orbital occupation calculations.

Basic DOS moments:

.. code-block:: bash

   rs-band --moments

With custom parameters:

.. code-block:: bash

   rs-band --moments --mesh 8 8 8 --electrons 8.0 --temperature 300 --output-prefix analysis

Complete Analysis Pipeline
--------------------------

For comprehensive analysis including bands, DOS, and projections:

.. code-block:: bash

   rs-band

This runs the full analysis pipeline with default parameters.

Example
-------

Use the example data in `examples/bccFe/` to run a quick test.

The example includes:

* ``lattice.nml`` — Crystal structure definition
* ``rs2txt.dat`` — TB-LMTO Hamiltonian data
* ``input.nml`` — Optional configuration file

Run analysis on the example:

.. code-block:: bash

   cd examples/bccFe
   rs-band --moments --mesh 6 6 6 --electrons 8.0

Expected output for BCC Fe:

* Fermi level: ~-0.15 eV
* s electrons: ~1.5
* p electrons: ~0.9
* d electrons: ~5.6
* Total: 8.0 electrons
