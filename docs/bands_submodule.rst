.. _bands_submodule:

The `bands` submodule
---------------------

The **`bands`** submodule provides pure-Python band structure and DOS analysis for RS-LMTO real-space Hamiltonians.

**Required input files:**

- ``lattice.nml`` – Crystal structure (lattice vectors, atomic positions)
- ``rs2txt.dat`` – Real-space Hamiltonian matrix elements (pre-computed by RS-LMTO)

Key capabilities
~~~~~~~~~~~~~~~~
- High-symmetry band structures with k-paths from spglib/seekpath
- Projected bands (orbital/site/spin) with publication-quality plots
- Full Brillouin zone DOS via tetrahedron or Gaussian broadening
- DOS moments and Fermi level determination from electron count

Core classes
~~~~~~~~~~~~
.. autoclass:: rsutil.bands.LatticeReader
   :members:

.. autoclass:: rsutil.bands.HamiltonianReader
   :members:

.. autoclass:: rsutil.bands.BandStructureCalculator
   :members:

.. autoclass:: rsutil.bands.ProjectedBandCalculator
   :members:

.. autoclass:: rsutil.bands.FullBZCalculator
   :members:

.. autoclass:: rsutil.bands.DOSCalculator
   :members:

.. autoclass:: rsutil.bands.DOSMomentsCalculator
   :members:

CLI examples
~~~~~~~~~~~~

.. code-block:: bash

    # Band structure (high-symmetry path)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

    # Projected bands (orbital)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --projections --projection orbital

    # DOS with moments and Fermi level
    rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 --electrons 16 --moments

Programmatic examples
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from rsutil.bands import (
        LatticeReader, HamiltonianReader,
        BandStructureCalculator, ProjectedBandCalculator,
        FullBZCalculator, DOSCalculator, DOSMomentsCalculator,
    )

    lattice = LatticeReader('lattice.nml')
    ham = HamiltonianReader('rs2txt.dat')

    # Bands
    bcalc = BandStructureCalculator(lattice, ham)
    kpts, klabels, kdist = bcalc.setup_kpath(npoints=100)
    bands = bcalc.calculate_bands(kpts)
    bcalc.plot_bands(kdist, bands, klabels, 'bands.png')

    # DOS + moments
    bz = FullBZCalculator(lattice, ham)
    dcalc = DOSCalculator(bz)
    dos = dcalc.calculate_dos(mesh_size=(6, 6, 6))
    dcalc.plot_dos(dos, fermi_energy=0.0, output_file='dos.png')

    moments_calc = DOSMomentsCalculator(dos['energies'], dos['total_dos'])
    ef = moments_calc.find_fermi_level(target_electrons=16)
    moments = moments_calc.calculate_moments(ef)
