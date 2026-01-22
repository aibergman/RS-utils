DOS Moments Analysis
====================

The DOS moments analysis functionality provides comprehensive analysis of electronic structure moments, including Fermi level determination and orbital occupation calculations.

Overview
--------

The DOS moments calculator computes:

* **Zeroth moment (occupation)**: Number of electrons in each orbital/spin channel
* **First moment (band center)**: Average energy of occupied states
* **Second moment (band width)**: RMS width of the occupied states
* **Fermi level**: Automatically determined from total electron count using Fermi-Dirac statistics

Key Features
------------

* **Automatic Fermi Level Determination**: Uses scipy.optimize.brentq for robust root finding
* **Temperature-Dependent Analysis**: Supports both T=0K and finite temperature calculations
* **Orbital and Spin Resolution**: Separate analysis for s, p, d orbitals and spin channels
* **Numerical Stability**: Robust handling of degenerate tetrahedra and numerical edge cases
* **Comprehensive Output**: Detailed results with convergence information and error bounds

Command Line Interface
----------------------

Basic Usage
~~~~~~~~~~~

Run DOS moments analysis with default parameters:

.. code-block:: bash

   rs-band --moments

Custom Parameters
~~~~~~~~~~~~~~~~~

Specify mesh size, electron count, temperature, and output options:

.. code-block:: bash

   rs-band --moments --mesh 8 8 8 --electrons 8.0 --temperature 300 --output-prefix analysis

Parameter Options
~~~~~~~~~~~~~~~~~

.. list-table:: CLI Parameters
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Description
   * - ``--mesh X Y Z``
     - 8 8 8
     - k-point mesh size for Brillouin zone integration
   * - ``--electrons N``
     - auto
     - Total number of electrons (auto-detected if not specified)
   * - ``--temperature T``
     - 300
     - Temperature in Kelvin for Fermi-Dirac distribution
   * - ``--output-prefix PREFIX``
     - dos_moments
     - Prefix for output files
   * - ``--no-symmetry``
     - False
     - Disable symmetry reduction in k-point generation
   * - ``--gaussian``
     - False
     - Use Gaussian broadening instead of tetrahedron method

Help System
~~~~~~~~~~~

Get detailed help for the moments command:

.. code-block:: bash

   rs-band --moments --help

Python API
----------

Direct Python Usage
~~~~~~~~~~~~~~~~~~~

For advanced users, the DOS moments functionality can be accessed directly through Python:

.. code-block:: python

   from rs_band.band_structure_calculator import LatticeReader, HamiltonianReader
   from rs_band.full_bz_dos import FullBZCalculator, DOSCalculator
   from rs_band.dos_moments import DOSMomentsCalculator

   # Read input files
   lattice = LatticeReader('lattice.nml')
   hamiltonian = HamiltonianReader('rs2txt.dat')

   # Create calculators
   bz_calc = FullBZCalculator(lattice, hamiltonian)
   dos_calc = DOSCalculator(bz_calc)

   # Calculate DOS with moments
   dos_data = dos_calc.calculate_dos_with_moments(
       mesh_size=(6, 6, 6),
       target_electrons=8.0,
       temperature=300.0
   )

   # Access results
   fermi_level = dos_data['fermi_level']
   moments = dos_data['moments']

DOSMomentsCalculator Class
~~~~~~~~~~~~~~~~~~~~~~~~~~

The core analysis is performed by the :class:`rs_band.dos_moments.DOSMomentsCalculator` class:

.. autoclass:: rs_band.dos_moments.DOSMomentsCalculator
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

Key Methods
~~~~~~~~~~~

.. automethod:: rs_band.dos_moments.DOSMomentsCalculator.find_fermi_level
   :no-index:

.. automethod:: rs_band.dos_moments.DOSMomentsCalculator.calculate_moments
   :no-index:

.. automethod:: rs_band.dos_moments.DOSMomentsCalculator.print_summary
   :no-index:

.. automethod:: rs_band.dos_moments.DOSMomentsCalculator.save_moments
   :no-index:

Analysis Functions
~~~~~~~~~~~~~~~~~~

The DOS moments analysis is performed using the :class:`rs_band.dos_moments.DOSMomentsCalculator` class described above.

Mathematical Background
-----------------------

Fermi Level Determination
~~~~~~~~~~~~~~~~~~~~~~~~~

The Fermi level is determined by solving the equation:

.. math::

   N = \int_{-\infty}^{\mu} D(\epsilon) f(\epsilon, T) d\epsilon

where:
* :math:`N` is the target number of electrons
* :math:`D(\epsilon)` is the density of states
* :math:`f(\epsilon, T)` is the Fermi-Dirac distribution
* :math:`\mu` is the Fermi level

For finite temperature, the Fermi-Dirac distribution is:

.. math::

   f(\epsilon, T) = \frac{1}{1 + \exp\left(\frac{\epsilon - \mu}{k_B T}\right)}

DOS Moments Definition
~~~~~~~~~~~~~~~~~~~~~~

The moments are calculated as weighted integrals:

**Zeroth moment (occupation):**

.. math::

   M_0 = \int_{-\infty}^{\mu} D(\epsilon) f(\epsilon, T) d\epsilon

**First moment (band center):**

.. math::

   M_1 = \frac{\int_{-\infty}^{\mu} \epsilon D(\epsilon) f(\epsilon, T) d\epsilon}{M_0}

**Second moment (band width):**

.. math::

   M_2 = \sqrt{\frac{\int_{-\infty}^{\mu} (\epsilon - M_1)^2 D(\epsilon) f(\epsilon, T) d\epsilon}{M_0}}

Numerical Implementation
~~~~~~~~~~~~~~~~~~~~~~~~

* **Integration**: Uses scipy.integrate.simpson for accurate numerical integration
* **Root Finding**: Employs scipy.optimize.brentq for robust Fermi level determination
* **Tetrahedron Method**: Handles Brillouin zone integration with proper degeneracy handling
* **Overflow Protection**: Implements numerical stability for extreme temperature ranges

Output and Results
------------------

Console Output
~~~~~~~~~~~~~~

The analysis produces detailed console output:

.. code-block:: text

   ============================================================
   DOS MOMENTS SUMMARY
   ============================================================
   Fermi Level: -0.154 eV
   Temperature: 300.0 K
   ============================================================

   total:
   ----------------------------------------
     Occupation:     8.000000 electrons
     Band Center:    -0.310422 eV
     Band Width:     0.136077 eV
     Avg E (rel EF): -0.156213 eV

   orbital_s:
   ----------------------------------------
     Occupation:     1.463558 electrons
     Band Center:    -0.485898 eV
     Band Width:     0.110488 eV
     Avg E (rel EF): -0.331689 eV

File Output
~~~~~~~~~~~

Results are saved to compressed NumPy files (.npz) containing:

* ``fermi_level``: Fermi level in eV
* ``temperature``: Temperature in K
* ``moments``: Dictionary with moments for each component
* ``energies``: Energy grid used for calculations
* ``total_dos``: Total density of states
* ``projected_dos``: All projected DOS components
* ``convergence_info``: Details about numerical convergence

Example Analysis
----------------

BCC Iron Case Study
~~~~~~~~~~~~~~~~~~~

For body-centered cubic iron with 8 valence electrons:

.. code-block:: bash

   # Run analysis
   rs-band --moments --mesh 6 6 6 --electrons 8.0 --temperature 300 --output-prefix bccfe

Expected results for BCC Fe:

* **Total occupation**: 8.0 electrons
* **s orbital**: ~1.5 electrons (band center: ~-0.49 eV)
* **p orbital**: ~0.9 electrons (band center: ~-0.35 eV)
* **d orbital**: ~5.6 electrons (band center: ~-0.25 eV)
* **Fermi level**: ~-0.15 eV

Temperature Effects
~~~~~~~~~~~~~~~~~~~

Compare T=0K and finite temperature results:

.. code-block:: bash

   # T=0K calculation (step function)
   rs-band --moments --temperature 0 --output-prefix t0k

   # Room temperature
   rs-band --moments --temperature 300 --output-prefix rt

Convergence Analysis
~~~~~~~~~~~~~~~~~~~~

Test convergence with different mesh sizes:

.. code-block:: bash

   # Coarse mesh for testing
   rs-band --moments --mesh 4 4 4 --output-prefix coarse

   # Production quality
   rs-band --moments --mesh 12 12 12 --output-prefix fine

Advanced Usage
--------------

Custom DOS Data
~~~~~~~~~~~~~~~

Use pre-computed DOS data:

.. code-block:: python

   import numpy as np
   from rs_band.dos_moments import DOSMomentsCalculator

   # Load your DOS data
   energies = np.load('energies.npy')
   total_dos = np.load('total_dos.npy')
   projected_dos = np.load('projected_dos.npy').item()

   # Create calculator
   calc = DOSMomentsCalculator(energies, total_dos, projected_dos, temperature=100)

   # Find Fermi level
   ef = calc.find_fermi_level(target_electrons=10.5)

   # Calculate moments
   moments = calc.calculate_moments(ef)

Integration with Existing Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DOS moments analysis integrates seamlessly with the existing rs_band workflow:

1. **Compatible Files**: Uses the same lattice.nml and rs2txt.dat files
2. **Shared Infrastructure**: Leverages existing DOS calculation methods
3. **Consistent Parameters**: Same mesh sizes and symmetry options
4. **Combined Analysis**: Can be run alongside band structure calculations

Error Handling
~~~~~~~~~~~~~~

The implementation includes robust error handling for:

* **No Solution**: When target electron count cannot be achieved in the energy range
* **Numerical Instability**: Division by zero in degenerate tetrahedra
* **Invalid Data**: Missing or malformed DOS/projected DOS data
* **Integration Failures**: Numerical issues in moment calculations

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Fermi level not found:**

* Check energy range covers the relevant states
* Verify electron count is reasonable for the system
* Try adjusting temperature (lower T can help convergence)

**Unrealistic occupations:**

* Check mesh convergence (increase k-point density)
* Verify input files are correct
* Compare with literature values for similar systems

**Numerical warnings:**

* These are usually harmless but indicate edge cases
* Check for degenerate energy levels in DOS
* Consider using Gaussian broadening for smoother DOS

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

For large calculations:

* Start with coarse meshes (4³) for parameter testing
* Use symmetry reduction when possible (--no-symmetry to disable)
* Consider temperature=0 for faster convergence in simple cases
* Monitor memory usage for very fine meshes (>16³)

References
----------

The DOS moments implementation is based on:

* Linear tetrahedron method for Brillouin zone integration
* Fermi-Dirac statistics for occupation calculations
* Numerical root finding for Fermi level determination
* Standard definitions of DOS moments in electronic structure theory
