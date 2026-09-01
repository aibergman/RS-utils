.. _bands_submodule:

The `bands` submodule: Band Structure and DOS Analysis
=======================================================

The **`bands`** submodule provides comprehensive pure-Python band structure and DOS analysis for RS-LMTO real-space Hamiltonians with advanced features for publication-quality visualization and analysis.

**Required input files:**

- ``lattice.nml`` – Crystal structure (lattice vectors, atomic positions, Fermi energy)
- ``rs2txt.dat`` – Real-space Hamiltonian matrix elements (complex H(R) blocks)

**Optional dependencies:**

- ``seekpath`` – For advanced k-path generation (recommended for automatic paths)
- ``pymatgen`` – For Brillouin zone visualization
- ``spglib`` – For automatic space group detection

Key capabilities
~~~~~~~~~~~~~~~~~

**Band Structure:**
- High-symmetry band structures with automatic k-paths from spglib, pymatgen, or seekpath
- Interactive k-path selection with custom path builder
- Prime convention support for symmetry-related k-points
- Non-collinear and collinear spin structures
- Eigenvector storage for advanced analysis

**Projections & Analysis:**
- Orbital-projected bands (s, p, d character)
- Site-projected bands (atom-specific contributions)
- Spin projections with 4 visualization modes (expectation, weight_up/down, polarization)
- Gauge fixing for smooth eigenvector phase continuity
- Degenerate band stabilization (3 modes: none, block_average, pair_cancel)

**DOS & Integration:**
- Full Brillouin zone DOS with tetrahedron method
- Site and orbital projected DOS
- Fermi level determination from electron count
- DOS moments calculation
- Symmetry reduction for efficient BZ sampling

**Visualization:**
- Publication-quality band structure plots with projections
- 3D Brillouin zone polyhedron with high-symmetry point labels
- Colored scatter plots with projection-weighted marker sizing
- Automated k-path export for reproducibility
- TeX-formatted labels for high-symmetry points

**Magnetic Analysis:**
- Non-collinear spin texture analysis
- Effective mass calculation
- Berry curvature and topological invariants
- Altermagnetic property detection (spin splitting without net magnetization)
- Magnetic moment calculations

Core I/O Classes
~~~~~~~~~~~~~~~~

.. autoclass:: rsutil.bands.LatticeReader
   :members:
   :undoc-members:

.. autoclass:: rsutil.bands.HamiltonianReader
   :members:
   :undoc-members:

Core Band Structure Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rsutil.bands.BandStructureCalculator
   :members: setup_kpath, setup_kpath_seekpath, calculate_bands, fourier_transform_hamiltonian
   :undoc-members:

.. autoclass:: rsutil.bands.ProjectedBandCalculator
   :members: calculate_projected_bands
   :undoc-members:

Core DOS Classes
~~~~~~~~~~~~~~~~

.. autoclass:: rsutil.bands.FullBZCalculator
   :members:
   :undoc-members:

.. autoclass:: rsutil.bands.DOSCalculator
   :members: calculate_dos, plot_dos
   :undoc-members:

.. autoclass:: rsutil.bands.DOSMomentsCalculator
   :members: find_fermi_level, calculate_moments
   :undoc-members:

Advanced Analysis Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rsutil.bands.SpinMomentAnalyzer
    :members:
    :undoc-members:

.. autoclass:: rsutil.bands.SpinAnalyzer
    :members:
    :undoc-members:
    :noindex:

.. autoclass:: rsutil.bands.BandAnalyzer
   :members:
   :undoc-members:

.. autoclass:: rsutil.bands.MagneticMomentCalculator
   :members:
   :undoc-members:

.. autoclass:: rsutil.bands.EffectiveMassCalculator
   :members:
   :undoc-members:

.. autoclass:: rsutil.bands.BerryCurvatureCalculator
   :members:
   :undoc-members:

.. autoclass:: rsutil.bands.AltermagneticAnalyzer
   :members:
   :undoc-members:

Spin Plotting Helpers
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: rsutil.bands.analysis.spin.plot_spin_texture

.. autofunction:: rsutil.bands.analysis.spin.plot_orbital_magnetism

.. autofunction:: rsutil.bands.analysis.spin.plot_spin_texture_2d

.. autofunction:: rsutil.bands.analysis.spin.plot_fermi_surface_spin_texture

K-Path Utilities (Advanced)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These utilities provide advanced k-path generation, management, and manipulation.

.. code-block:: python

    from rsutil.bands import (
        seekpath_kpath, build_custom_kpath,
        apply_prime_convention, allocate_points_per_segment,
        normalize_label, tex_format_label,
        get_reciprocal_lattice, fractional_to_cartesian
    )

**seekpath integration:**

.. autofunction:: rsutil.bands.seekpath_kpath
   :noindex:

**Custom k-path builder:**

.. autofunction:: rsutil.bands.build_custom_kpath
   :noindex:

**Prime convention:**

.. autofunction:: rsutil.bands.apply_prime_convention
   :noindex:

**Coordinate utilities:**

.. autofunction:: rsutil.bands.fractional_to_cartesian
   :noindex:

.. autofunction:: rsutil.bands.get_reciprocal_lattice
   :noindex:

Projection Utilities
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: rsutil.bands.stabilize_degenerate_projections
   :noindex:

.. autofunction:: rsutil.bands.compute_quantity_from_sexp
   :noindex:

.. autofunction:: rsutil.bands.get_projection_vrange
   :noindex:

.. autofunction:: rsutil.bands.gauge_fix_evecs
   :noindex:

Visualization Utilities
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: rsutil.bands.plot_brillouin_zone
   :noindex:

.. autofunction:: rsutil.bands.plot_bands_with_projections
   :noindex:

.. autofunction:: rsutil.bands.export_kpath
   :noindex:

CLI examples
~~~~~~~~~~~~

**Basic band structure (automatic k-path):**

.. code-block:: bash

    # Default path with 100 k-points
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

**Projected bands (orbital contributions):**

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
        --projections --projection orbital --output bands_orbital.png

**Spin-projected bands (expectation values):**

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
        --projections --projection spin --projections-mode expectation \
        --output bands_spin.png

**DOS with Fermi level and moments:**

.. code-block:: bash

    rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat \
        --mesh 6 6 6 --electrons 16 --moments --output dos.png

Bands CLI option reference (synchronized with ``rsutil bands --help``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input and format options:

- ``--lattice`` (default: ``lattice.nml``)
- ``--hamiltonian`` (default: ``rs2txt.dat``)
- ``--lattice-format``: ``rs-lmto``, ``paoflow``, ``kite``
- ``--hamiltonian-format``: ``rs-lmto``, ``paoflow``, ``kite``

Core band/projection options:

- ``--nk`` number of k-points
- ``--output`` output plot filename
- ``--projections`` enable projected bands
- ``--projection``: ``orbital``, ``site``, ``spin``, ``all``
- ``--projections-mode``: ``expectation``, ``weight_up``, ``weight_down``, ``polarization``
- ``--output-prefix``
- ``--energy-range emin emax``

Altermagnetic options:

- ``--altermagnetic``
- ``--mag-threshold``
- ``--plot-splitting``
- ``--plot-character``

Berry/topology options:

- ``--berry``
- ``--fermi-energy``
- ``--temperature``
- ``--chern-bands``
- ``--chern-tol``
- ``--ahc``
- ``--orbital-mag``
- ``--num-bands``
- ``--plot-berry``
- ``--plot-chern``
- ``--plot-bands``

K-path and Brillouin-zone options:

- ``--kpath-mode``: ``auto``, ``seekpath``, ``custom``, ``interactive``
- ``--seekpath-symprec``
- ``--seekpath-distance``
- ``--custom-path``
- ``--prime-points``
- ``--stabilize-degenerate``: ``none``, ``block_average``, ``pair_cancel``
- ``--plot-bz``
- ``--export-kpath``

Spin texture options:

- ``--spin-texture``
- ``--spin-band``
- ``--spin-plane``: ``xy``, ``xz``, ``yz``
- ``--spin-component``: ``x``, ``y``, ``z``, ``magnitude``
- ``--rashba-analysis``
- ``--plot-spin-fermi``

Fat-band options:

- ``--fat-bands``
- ``--fat-type``: ``linewidth``, ``color``, ``spin``
- ``--fat-orbitals``
- ``--fat-orbital-groups`` (JSON mapping)
- ``--fat-scale``

Additional focused CLI examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Berry + Chern workflow**

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 150 \
        --berry --chern-bands 0 1 2 --plot-chern --plot-bands

**Spin texture workflow**

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 180 \
        --spin-texture --spin-band 3 --spin-plane xy --spin-component z \
        --rashba-analysis --plot-spin-fermi

**Fat-band workflow**

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 160 \
        --fat-bands --fat-type spin --output-prefix fatbands

Programmatic Examples
~~~~~~~~~~~~~~~~~~~~~

**Example 1: Basic band structure with automatic k-path**

.. code-block:: python

    from rsutil.bands import LatticeReader, HamiltonianReader, BandStructureCalculator

    # Load data
    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    # Create calculator
    calc = BandStructureCalculator(lattice, hamiltonian)

    # Automatic k-path from structure analysis
    kpoints, klabels, kdistances = calc.setup_kpath(npoints=100)

    # Calculate band structure
    band_energies = calc.calculate_bands(kpoints)

    # Plot
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for band in band_energies.T:
        ax.plot(kdistances, band - lattice.fermi_energy, 'k-', lw=0.5)
    ax.set_ylabel('Energy (eV)')
    ax.axhline(0, color='r', linestyle='--', alpha=0.5)
    plt.savefig('bands.png', dpi=300, bbox_inches='tight')

**Example 2: Interactive k-path selection with seekpath**

.. code-block:: python

    from rsutil.bands import BandStructureCalculator, LatticeReader, HamiltonianReader

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')
    calc = BandStructureCalculator(lattice, hamiltonian)

    # Interactive mode: user chooses standard or custom path
    kpoints, klabels, kdistances = calc.setup_kpath_seekpath(
        num_points=200, 
        interactive=True
    )
    # Prompts: 
    #   "Use standard seekpath k-path? [Y/n]: "
    #   If 'n': "Enter k-path labels separated by '-', e.g K'-Γ-K : "

    band_energies = calc.calculate_bands(kpoints)


**Example 3: Spin-projected band structure with degenerate stabilization**

.. code-block:: python

    from rsutil.bands import (
        ProjectedBandCalculator, LatticeReader, HamiltonianReader,
        stabilize_degenerate_projections, plot_bands_with_projections
    )
    import numpy as np

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    # Create projected calculator
    calc = ProjectedBandCalculator(
        lattice, hamiltonian,
        spin_projection_mode='expectation'  # or 'weight_up', 'weight_down', 'polarization'
    )

    # K-path
    kpoints, klabels, kdistances = calc.setup_kpath_seekpath(num_points=150)

    # Calculate bands with info
    band_energies = calc.calculate_bands(kpoints)

    # For each k-point, stabilize degenerate band projections
    spin_proj = {}
    for k_idx in range(len(kpoints)):
        # Get raw numerator and denominator
        numerator = calc.get_spin_numerator(k_idx, component='z')
        denominator = calc.get_projection_weight(k_idx)
        
        # Stabilize to avoid oscillations in degenerate blocks
        spin_proj[k_idx] = stabilize_degenerate_projections(
            energies=band_energies[k_idx],
            numerators=numerator,
            denominators=denominator,
            mode='block_average',  # Subspace averaging
            energy_tol=1e-6
        )

    # Plot with colored projections
    fig, ax = plot_bands_with_projections(
        kpath_x=kdistances,
        band_energies=band_energies,
        projections=np.array([spin_proj[i] for i in range(len(kpoints))]),
        xticks=xticks,
        xticklabels=xticklabels,
        fermi_energy=lattice.fermi_energy,
        projection_mode='expectation',
        marker_sizes=projection_weights,
        ylim=(-3, 2),
        title='Spin-Projected Band Structure (Sz expectation)'
    )
    plt.savefig('bands_spin.png', dpi=300, bbox_inches='tight')
    plt.show()

**Example 4: 3D Brillouin zone visualization**

.. code-block:: python

    from rsutil.bands import (
        seekpath_kpath, plot_brillouin_zone,
        LatticeReader, HamiltonianReader
    )

    lattice = LatticeReader('lattice.nml')

    # Get BZ information from seekpath
    kpath_data = seekpath_kpath(
        lattice_vectors=lattice.lattice_vectors,
        atomic_positions=lattice.atomic_positions,
        atomic_numbers=lattice.atomic_numbers
    )

    # Plot 3D BZ with high-symmetry points
    fig, ax = plot_brillouin_zone(
        reciprocal_lattice=kpath_data['reciprocal_lattice'],
        point_coords=kpath_data['point_coords'],
        label_hs_points=True,
        label_primes=True,  # Show K', M', etc.
        figsize=(10, 10)
    )
    plt.savefig('brillouin_zone.png', dpi=300, bbox_inches='tight')
    plt.show()

**Example 5: Custom k-path with prime convention**

.. code-block:: python

    from rsutil.bands import (
        BandStructureCalculator, seekpath_kpath,
        build_custom_kpath, LatticeReader, HamiltonianReader
    )

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    # Get seekpath data for high-symmetry points
    kpath_data = seekpath_kpath(
        lattice_vectors=lattice.lattice_vectors,
        atomic_positions=lattice.atomic_positions,
        atomic_numbers=lattice.atomic_numbers
    )

    # Build custom path: K' → Γ → K → M
    # (K' is K with in-plane flip via prime convention)
    kpoints, xvals, xticks, xticklabels, metadata = build_custom_kpath(
        point_coords=kpath_data['point_coords'],
        reciprocal_lattice=kpath_data['reciprocal_lattice'],
        labels_list=["K'", "Gamma", "K", "M"],
        total_kpoints=400
    )

    # Which points are primed?
    print("Prime flags:", metadata['node_isprime'])  # [True, False, False, False]
    print("Segment lengths:", metadata['seg_lengths'])  # Physical lengths in |k|

    # Calculate bands
    calc = BandStructureCalculator(lattice, hamiltonian)
    band_energies = calc.calculate_bands(kpoints)

**Example 6: DOS and Fermi level determination**

.. code-block:: python

    from rsutil.bands import (
        FullBZCalculator, DOSCalculator, DOSMomentsCalculator,
        LatticeReader, HamiltonianReader
    )

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    # Full BZ calculator
    bz_calc = FullBZCalculator(lattice, hamiltonian)

    # DOS calculator
    dos_calc = DOSCalculator(bz_calc)

    # Calculate DOS on full BZ mesh
    dos_result = dos_calc.calculate_dos(
        mesh_size=(8, 8, 8),  # k-mesh
        energy_grid=np.linspace(-5, 5, 500),
        method='tetrahedron'  # Accurate integration method
    )

    # Plot
    dos_calc.plot_dos(dos_result, output_file='dos.png')

    # Determine Fermi level from electron count
    moments_calc = DOSMomentsCalculator(
        dos_result['energies'],
        dos_result['total_dos']
    )

    fermi_level = moments_calc.find_fermi_level(
        target_electrons=16,  # Valence electrons
        temperature=0.0
    )
    print(f"Fermi level: {fermi_level:.4f} eV")

    # Calculate moments at Fermi level
    moments = moments_calc.calculate_moments(fermi_level)
    print(f"Band energy: {moments['band_energy']:.4f} eV")
    print(f"DOS at EF: {moments['dos_at_ef']:.4f} states/eV")

**Example 7: Orbital-projected bands**

.. code-block:: python

    from rsutil.bands import ProjectedBandCalculator, LatticeReader, HamiltonianReader

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    calc = ProjectedBandCalculator(lattice, hamiltonian)
    kpoints, klabels, kdistances = calc.setup_kpath(npoints=100)

    # Calculate projected bands
    band_energies, projections = calc.calculate_projected_bands(kpoints)

    # projections dict contains:
    # - 'orbital': {s, p, d} contributions
    # - 'site': per-atom contributions
    # - 'spin_expectation': Sz values
    # - 'orbital_site': combined

    # Plot orbital character
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot bands
    for band in band_energies.T:
        ax.plot(kdistances, band - lattice.fermi_energy, 'k-', lw=0.5, alpha=0.3)

    # Color by d-orbital character
    d_proj = projections['orbital']['d']
    sc = ax.scatter(
        np.tile(kdistances, band_energies.shape[1]),
        (band_energies - lattice.fermi_energy).ravel(),
        c=d_proj.ravel(),
        s=20, cmap='Reds', alpha=0.7
    )
    plt.colorbar(sc, label='d-orbital character')
    ax.set_ylabel('Energy (eV)')
    ax.set_ylim(-3, 2)
    plt.savefig('bands_d_character.png', dpi=300, bbox_inches='tight')

**Example 8: Effective mass calculation**

.. code-block:: python

    from rsutil.bands import EffectiveMassCalculator, BandStructureCalculator
    from rsutil.bands import LatticeReader, HamiltonianReader

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    calc = BandStructureCalculator(lattice, hamiltonian)
    kpoints, klabels, kdistances = calc.setup_kpath(npoints=50)
    band_energies, eigenvectors = calc.calculate_bands_with_eigenvectors(kpoints)

    # Effective mass around band extrema
    em_calc = EffectiveMassCalculator(
        band_energies, eigenvectors, kpoints,
        lattice.lattice_vectors
    )

    # Calculate around conduction band minimum (band 0 at k=0)
    em = em_calc.effective_mass(band_index=0, k_index=0)
    print(f"Effective mass tensor:\n{em}")

**Example 9: Berry curvature and topological properties**

.. code-block:: python

    from rsutil.bands import BerryCurvatureCalculator, BandStructureCalculator
    from rsutil.bands import LatticeReader, HamiltonianReader

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    calc = BandStructureCalculator(lattice, hamiltonian)
    
    # Calculate on dense k-mesh
    kpts_mesh = calc.get_dense_kpoints(mesh_size=(20, 20, 20))
    energies, eigenvectors = calc.calculate_bands_with_eigenvectors(kpts_mesh)

    # Berry curvature
    bc_calc = BerryCurvatureCalculator(
        energies, eigenvectors, kpts_mesh,
        reciprocal_lattice=calc.get_reciprocal_lattice(),
        fermi_energy=lattice.fermi_energy
    )

    # Berry curvature Ω_z(k)
    bc_z = bc_calc.calculate_berry_curvature(component='z')

    # Chern number (integral over BZ)
    chern = bc_calc.calculate_chern_number(component='z')
    print(f"Chern number (Ω_z): {chern:.4f}")

**Example 10: K-path export for reproducibility**

.. code-block:: python

    from rsutil.bands import export_kpath, BandStructureCalculator
    from rsutil.bands import LatticeReader, HamiltonianReader

    lattice = LatticeReader('lattice.nml')
    hamiltonian = HamiltonianReader('rs2txt.dat')

    calc = BandStructureCalculator(lattice, hamiltonian)
    kpoints, klabels, kdistances = calc.setup_kpath(npoints=100)

    # Save k-path to file for reproducibility
    export_kpath(
        filename='band_kpath.dat',
        kpath_x=kdistances,
        kpoints_abs=kpoints,
        kpoint_labels=klabels,
        header='x_linearcoord  kx  ky  kz  label'
    )

    # Other researchers can now use this exact k-path

Complete public API
~~~~~~~~~~~~~~~~~~~

The sections above focus on the most common workflows.  The following module
references expose the remaining public reader, numerical, path, projection,
spin-texture, and fat-band helpers for users building their own workflows.

.. automodule:: rsutil.bands.readers
   :members:
   :noindex:

.. automodule:: rsutil.bands.core.band_calculations
   :members:
   :noindex:

.. automodule:: rsutil.bands.core.hamiltonian
   :members:
   :noindex:

.. automodule:: rsutil.bands.core.kpath_generator
   :members:
   :noindex:

.. automodule:: rsutil.bands.dos.integration
   :members:
   :noindex:

.. automodule:: rsutil.bands.analysis.structure_analysis
   :members:
   :noindex:

.. automodule:: rsutil.bands.analysis.spin.spin_texture
   :members:
   :noindex:

.. automodule:: rsutil.bands.visualization.fat_bands
   :members:
   :noindex:

.. automodule:: rsutil.bands.paths.coordinates
   :members:
   :noindex:

.. automodule:: rsutil.bands.paths.custom_paths
   :members:
   :noindex:

.. automodule:: rsutil.bands.paths.labels
   :members:
   :noindex:

.. automodule:: rsutil.bands.paths.point_allocation
   :members:
   :noindex:

.. automodule:: rsutil.bands.projections.projection_utils
   :members:
   :noindex:

Configuration & Spin Projection Modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Spin projection modes** (for colored band plots):

- ``"expectation"``: Direct spin expectation value ⟨S_z⟩ in [-0.5, +0.5]
- ``"weight_up"``: Spin-up probability = 0.5 + ⟨S_z⟩ in [0, 1]
- ``"weight_down"``: Spin-down probability = 0.5 - ⟨S_z⟩ in [0, 1]  
- ``"polarization"``: Signed polarization = 2⟨S_z⟩ in [-1, +1]

**Degenerate band stabilization modes**:

- ``"none"``: No stabilization (inspect raw values)
- ``"block_average"``: Replace degenerate manifolds with subspace average ⟨O⟩ = (Σ num)/(Σ den)
- ``"pair_cancel"``: Heuristic to find and cancel near-opposite pairs

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **K-path generation**: seekpath (explicit k-path) is faster than pymatgen (path optimization)
- **DOS calculation**: Tetrahedron method is more accurate but slower; Gaussian broadening is approximate but faster
- **Symmetry reduction**: Use spglib for large systems to reduce k-mesh size
- **Eigenvector storage**: Only enable for projection analysis to save memory

Installation Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Minimal installation** (basic band structure):

.. code-block:: bash

    pip install rsutil

**With k-path enhancements**:

.. code-block:: bash

    pip install rsutil seekpath

**With full visualization**:

.. code-block:: bash

    pip install rsutil seekpath pymatgen

**For developers**:

.. code-block:: bash

    git clone <repo>
    cd RS_tools
    pip install -e . [seekpath,pymatgen]

See Also
~~~~~~~~~

See the command examples in this page for band/DOS workflows. High-throughput
deployment and support boundaries are covered in :doc:`HT_QUICKSTART` and
:doc:`HT_READINESS_REPORT`.
