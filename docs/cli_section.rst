.. _cli_section:

Command-Line Interface (CLI)
=============================

The **RS-utils CLI** provides single-run ``setup`` workflows, band/DOS
analysis, visualization, and campaign workflows, plus legacy flags for bulk
setup, exchange, and export.

Top-level usage:

.. code-block:: bash

    rsutil {setup,bands,dos,campaign,viz} [options]
    rsutil [--setup|--exchange|--export|--lmto47|--lattice|--clust] [options]

Common Tasks
------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Task
     - Command
   * - Setup from ``input.nml`` + ``lattice.nml``
     - ``rsutil setup bulk``
   * - Setup and validate generated files
     - ``rsutil setup bulk --validate --summary``
   * - Supercell setup (non-destructive default)
     - ``rsutil setup bulk --supercell 2 2 2``
   * - 2D slab setup with isotropic margin
     - ``rsutil setup bulk --2d --z-iso-factor 0.5``
   * - Override recursion controls in ``&control``
     - ``rsutil setup bulk --rec-method chebyshev --rec-steps 200``
   * - Build exchange map
     - ``rsutil --exchange --cutoff 2.5 --unique``
   * - Export UppASD input set
     - ``rsutil --export --jij jijso.out --out uppasd_out``
   * - Export LMTO47 CTRL
     - ``rsutil --lmto47 --nspin 2 --lmx 3 --out lmto47_out``
   * - Plot band structure
     - ``rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png``
   * - Prepare a surface setup from converged bulk files
     - ``rsutil setup surface --bulk bulk_run --output surface_run``
   * - Prepare an impurity in a converged bulk or surface host
     - ``rsutil setup impurity --host host_run --xyz impurity.xyz --output impurity_run``
   * - Compute DOS + moments
     - ``rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 --electrons 16 --moments``

Surface preparation
--------------------

Prepare an RS-LMTO-ASA ``buildsurf`` calculation from an existing converged
bulk directory:

.. code-block:: bash

    rsutil setup surface --bulk tmp/fccCu --output tmp/fccCu001 \
        --hkl 0 0 1 --layers 6 --empty-spheres 1

The command creates a new output directory, copies the bulk host potential and
radial ground-state files, carries over ``lattice.nml`` for ``crystal_sym='file'``,
and writes ``input.nml`` plus ``surface_setup.json``. It sets
``pre_processing='buildsurf'``, ``calctype='S'``, and the bulk Fermi level with
``fix_fermi=.true.``. RS-LMTO generates the surface cluster when it runs.
RS-utils clones the converged host atom file for both the self-consistent
surface types. It generates empty-sphere radial parameters with the existing
Fortran setup kernel, then assigns a small starting occupancy while preserving
the kernel's native namelist array indexing. It does not need to discover
empty spheres in the bulk geometry. The ES seed uses the lowest-branch
free-electron PL values ``1.5``, ``2.25``, and ``3.1476`` for s, p, and d,
respectively, in both spin channels.

Options:

* ``--crystal-sym fcc|bcc|hcp|file`` selects or overrides the bulk symmetry.
  ``file`` requires ``lattice.nml`` in the bulk directory.
* ``--hkl H K L`` selects a nonzero integer surface normal. For hcp, this is
  expanded to Miller-Bravais ``H K -(H+K) L``; four values can also be supplied
  directly as ``H K I L`` with ``I=-(H+K)``.
* ``--layers`` sets ``nlay``. RS-LMTO uses these near-surface planes as
  self-consistent representatives, and the output sets
  ``ntype_surface = ntype_bulk + nlay``.
* ``--empty-spheres`` assigns that many of the added representative types the
  empty-sphere label. It controls type labels, not extra geometric sites;
  ``buildsurf`` constructs the geometry and maps its planes to those labels.
* ``--empty-sphere-charge`` sets the initial occupancy per empty-sphere
  type. The default is 0.2 electrons; the same total amount is removed from
  interface-layer starting potentials when interface types are present.
* ``--charge-sigma`` sets the width, in layer spacings, of the erf-derived
  transfer profile. The default is 0.5; larger values spread the compensation
  farther into the interface. The weights are normalized over the available
  interface layers, so the net initial transfer is charge balanced.
* ``--orientation preserve`` keeps the input orientation and writes the
  requested normal to ``surftype``. ``--orientation global-z`` rigidly rotates
  the lattice and basis so the normal is global Z, writes ``lattice.nml``, and
  sets ``crystal_sym='file'`` and ``surftype='0 0 1'``.
* With ``--orientation global-z``, ``--align-bond X Y Z`` and
  ``--align-bond-to U V W`` remove the remaining in-plane rotation ambiguity.
  The first vector is a Cartesian bond direction in the original bulk frame
  (the units do not matter); the second is its desired Cartesian direction in
  the rotated frame and must lie in global XY. The source bond must lie in the
  selected surface plane. For example, for fcc(111),
  ``--align-bond 0 0.5 -0.5 --align-bond-to 1 -1 0`` aligns that bond with the
  global ``(1,-1,0)`` direction. The resulting Cartesian rotation matrix and
  alignment are recorded in ``surface_setup.json``. If omitted, the existing
  minimal normal-to-Z rotation is used.
* ``--n1``, ``--n2``, ``--n3``, ``--rc``, ``--self-steps``, and
  ``--mixing-beta`` set surface calculation controls. Defaults are 32, 32, 32,
  100, 15, and 0.01 respectively.

To expand an existing RS-utils surface setup, use it as the restart source:

.. code-block:: bash

    rsutil setup surface --from-surface tmp/fccCu001_3plus1 \
        --output tmp/fccCu001_5plus2 --layers 7 --empty-spheres 2

This example expands from three interface types plus one empty-sphere type to
five interface types plus two empty-sphere types. ``--layers`` is total ``nlay``
and includes empty-sphere types, so the expanded setup has ``nlay=7`` and
``ntype = ntype_bulk + 7``. Existing type restarts are preserved by label;
additional interface types are seeded from the deepest existing interface
restart, and additional empty-sphere types from the existing ES restart. The
surface normal, rotation, lattice, fixed Fermi level, charge-seed settings,
and other surface controls are inherited from ``surface_setup.json``. Pass
``--self-steps`` or ``--mixing-beta`` only when you want to override those
inherited SCF controls. The source manifest's original bulk directory must
remain available for rebuilding the input bundle. The source surface is left
unchanged and RS-LMTO regenerates the larger surface cluster when the new
calculation runs.

The source and destination must be separate directories; the destination must
not already exist. Built-in fcc/bcc/hcp inputs can use ``global-z`` without a
bulk ``lattice.nml``. Hcp uses the bulk ``celldm`` value when available and
otherwise uses ``c/a=1.633``.

Impurity preparation
--------------------

Prepare an embedded impurity calculation from a converged bulk or surface
calculation. The XYZ file uses standard XYZ structure: atom count, comment,
then one element symbol and Cartesian position in Å per row. Each listed atom
must coincide with a site in the host calculation's ``clust`` file. Coordinates
use the host calculation's Cartesian frame and origin.

.. code-block:: bash

    rsutil setup impurity --host tmp/fccCu --xyz impurity.xyz \
        --output tmp/Cu_imp --scf-radius 4.0

Without ``--scf-radius``, only the XYZ impurity sites are self-consistent.
With a positive radius, host sites within that distance of any listed impurity
are added to the ``inclu`` list and each receives its own host-like atom type.
Empty spheres are excluded from the grown host region. The command sets
``calctype='I'``, chooses ``newclubulk`` or ``newclusurf`` from the host, fixes
the host Fermi level, stages converged host potentials, and generates seed atom
files for the added types. The output includes ``impurity_sites.xyz`` with the
resolved site positions and ``impurity_setup.json`` with a setup manifest.

The host ``clust`` file defines the available sites and limits the largest
usable growth radius. Every selected site becomes a distinct self-consistent
type, so larger radii increase the number of types and atom files.

Band and DOS
------------

**Bands**

.. code-block:: bash

    # Band structure (high-symmetry path)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 --output bands.png

    # Projected bands
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 100 \
                 --projections --projection orbital

**DOS**

.. code-block:: bash

    # DOS with tetrahedron method
    rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6

    # DOS moments + Fermi level from electron count
    rsutil dos --lattice lattice.nml --hamiltonian rs2txt.dat --mesh 6 6 6 \
               --electrons 16 --moments


Advanced Band Structure Features
---------------------------------

**K-Path Generation Modes**

The ``--kpath-mode`` option controls how high-symmetry k-paths are generated:

.. code-block:: bash

    # Auto mode (simple, fast; default)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200

    # Seekpath mode (explicit high-symmetry path)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode seekpath

    # Custom path with specific labels
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 300 \
                 --kpath-mode custom --custom-path "K-Gamma-K"

    # Interactive mode: prompts for seekpath vs custom
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode interactive

**Seekpath Configuration**

Fine-tune seekpath behavior:

.. code-block:: bash

    # Stricter symmetry tolerance
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode seekpath --seekpath-symprec 1e-6

    # Finer Brillouin zone sampling
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode seekpath --seekpath-distance 0.01

**Custom K-Paths with Prime Convention**

Define custom paths and optionally apply the prime convention (K' = -K):

.. code-block:: bash

    # Simple custom path (K → Γ → M)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 300 \
                 --kpath-mode custom --custom-path "K-Gamma-M"

    # With prime points (K' is K with in-plane flip)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 300 \
                 --kpath-mode custom --custom-path "K-Prime-Gamma-K" --prime-points

**Degenerate Band Stabilization**

For materials with degenerate bands, choose a stabilization strategy for smooth color projections:

.. code-block:: bash

    # Block average (recommended; default)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --projections --projection spin \
                 --stabilize-degenerate block_average

    # Pair cancellation (heuristic; useful for many degeneracies)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --projections --projection spin \
                 --stabilize-degenerate pair_cancel

    # No stabilization (raw values; for debugging)
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --projections --projection spin \
                 --stabilize-degenerate none

**Brillouin Zone Visualization**

Plot the 3D Brillouin zone polyhedron with high-symmetry points labeled:

.. code-block:: bash

    # Plot BZ alongside bands
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --plot-bz

    # With prime points labeled
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --kpath-mode custom --custom-path "K-Gamma-M" \
                 --plot-bz --prime-points

**K-Path Export for Reproducibility**

Save the exact k-path used for future reference or use with other tools:

.. code-block:: bash

    # Export k-path to file
    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 200 \
                 --export-kpath kpath_band_structure.dat

Output file format (space-separated columns):

::

    x_linearcoord  kx  ky  kz
    0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
    1.234567e-02  -1.234567e-02   2.345678e-03   0.000000e+00
    ...

**Complete Workflow Example**

Combine seekpath, spin projections, degenerate stabilization, BZ plot, and k-path export:

.. code-block:: bash

    rsutil bands --lattice lattice.nml --hamiltonian rs2txt.dat --nk 250 \
                 --projections --projection spin \
                 --projections-mode expectation \
                 --kpath-mode seekpath \
                 --stabilize-degenerate block_average \
                 --plot-bz \
                 --export-kpath kpath_seekpath.dat \
                 --output-prefix band_analysis_spin

This command:
1. Generates a seekpath k-path with 250 points
2. Computes spin projections with ⟨S_z⟩ expectation values
3. Stabilizes degenerate manifolds using block averaging
4. Plots the Brillouin zone
5. Saves the k-path to ``kpath_seekpath.dat``
6. Outputs plots to ``band_analysis_spin_*.png``

**Available Options Summary**

=================================  ============================================  ================================================================
Option                             Values                                        Description
=================================  ============================================  ================================================================
``--kpath-mode``                   auto, seekpath, custom, interactive           K-path generation method (default: auto)
``--seekpath-symprec``             float [0.00001, 0.1]                          Symmetry tolerance for seekpath (default: 1e-5)
``--seekpath-distance``            float [0.001, 0.1]                            Reference distance for BZ sampling in Å^-1 (default: 0.025)
``--custom-path``                  string (e.g., "K-Gamma-K")                    Labels for custom k-path (required with --kpath-mode custom)
``--prime-points``                 flag                                          Apply prime convention to k-points
``--stabilize-degenerate``         none, block_average, pair_cancel               Degenerate band handling (default: block_average)
``--plot-bz``                      flag                                          Plot Brillouin zone polyhedron (requires pymatgen)
``--export-kpath``                 filename                                      Export k-path to file for reproducibility
=================================  ============================================  ================================================================

Lattice Processing
-------------------

The ``setup bulk`` command processes lattice structures and calculates
Wigner-Seitz radii and starting parameters. It reads ``lattice.nml`` and
``input.nml`` from the current directory. ``rsutil -S`` and ``rsutil --setup``
remain as compatibility shortcuts for this command.

**Flags:**

- ``rsutil setup bulk``: Process lattice and generate optimized atomic parameters.
- ``--cluster``: Cluster type: ``periodic`` (default) or ``free``.
- ``--empty_spheres``: Enable empty sphere calculations.
- ``--min_empty_radius``: Minimum empty sphere radius (default: 1.2).
- ``--max_empty_radius``: Maximum empty sphere radius (default: 2.5).
- ``--supercell NX NY NZ``: Repeat the input unit cell along the lattice vectors.
- ``--validate``: Validate output files after processing.
- ``--summary``: Print a summary of file changes after processing.
- ``--backup``: Create .bak backups of input files before processing.
- ``--in-place``: Replace input files with generated outputs (creates .bak backups).
- ``--lmx``: Override maximum l-channel used in setup path (default: ``2``; supported range ``0..3``).
- ``--preset KS``: Enable ``use_kspace = .true.`` in the ``&self`` namelist.
- ``--lattice-dist``: Gaussian std.dev for lattice-vector distortion in units of ``alat``.
- ``--basis-dist``: Gaussian std.dev for basis-vector distortion in units of ``alat``.
- ``--alat-scale``: Scalar factor applied to ``alat`` before distortion.
- ``--seed``: Integer seed for reproducible random distortions.
- ``--2d``: Enable 2D slab mode for setup/export geometry handling.
- ``--z-iso-factor``: Surface-margin factor for density-aware 2D ``Lz`` (requires ``--2d``).
- ``--rec-method``: Override ``&control recur`` in ``input_new.nml`` (``block``, ``chebyshev``, ``lanczos``).
- ``--rec-steps``: Override ``&control lld`` in ``input_new.nml``.
- ``--sqs``: Enable SQS assignment on all sites, e.g. ``"Fe=0.5,Ni=0.5"`` (requires ``--supercell``).
- ``--sqs-shells``: Number of pair-distance shells in the SQS objective.
- ``--sqs-steps``: Number of MC swap steps for SQS optimization.
- ``--sqs-seed``: Integer seed for reproducible SQS optimization.

Generated ``input_new.nml`` uses ``strux_backend = 'strux_lib'`` in its
``&lattice`` namelist unless the input specifies a backend explicitly.

**Example:**

.. code-block:: bash

    # Basic lattice processing
    rsutil setup bulk

    # Enable the KS setup preset
    rsutil setup bulk --preset KS

    # With empty spheres
    rsutil setup bulk --empty_spheres --min_empty_radius 1.0 --max_empty_radius 3.0

    # Free cluster calculation
    rsutil setup bulk --cluster free

    # Create a 1x1x2 supercell with broken symmetry labels
    rsutil setup bulk --supercell 1 1 2

    # Supercell plus empty spheres (use in-place updates and validation)
    rsutil setup bulk --supercell 1 1 2 --empty_spheres --in-place --validate

    # Distort lattice/basis vectors reproducibly and scale alat
    rsutil setup bulk --lattice-dist 0.01 --basis-dist 0.02 --alat-scale 1.01 --seed 1234

    # 2D mode with simple z-rule
    rsutil setup bulk --2d

    # 2D mode with isotropic-density margins
    rsutil setup bulk --2d --z-iso-factor 0.5

    # Override recursion controls in &control
    rsutil setup bulk --rec-method chebyshev
    rsutil setup bulk --rec-method lanczos --rec-steps 64

    # SQS on all supercell sites
    rsutil setup bulk --supercell 2 2 2 --sqs "Fe=0.5,Ni=0.3,Co=0.2" \
           --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234

Supercell Notes
~~~~~~~~~~~~~~~

- Supercells are generated with broken symmetry indexing (arrays run 1..N).
- Labels in input.nml are expanded into unique species (e.g., Fe1, N1, Fe2, N2).
- ntype and ct(i) are expanded to match the unique labels.
- Per-species .nml files are created by copying from the original label files.
- By default, ``--supercell`` is staged internally and original input files are restored after setup.
- Use ``--in-place`` if you explicitly want to replace ``lattice.nml``/``input.nml``.
- In ``--in-place`` mode, backups with ``.bak`` extension are created first.
- The same staged supercell behavior also applies when using ``--lmto47``.

2D Notes
~~~~~~~~

- ``--2d`` modifies setup/export geometry before WS/ASA steps.
- Default 2D z-length (when ``--z-iso-factor`` is omitted):
  - ``Lz = (zmax-zmin) + 1.0*alat``.
- Isotropic-density mode:
  - ``Lz = (zmax-zmin) + X*sqrt(A/n_bottom) + X*sqrt(A/n_top)`` with ``X = --z-iso-factor``.
- The third lattice vector is built from the surface normal, not forced to global z:
  - ``a3 || cross(a1,a2)``, sign chosen so ``a3·z >= 0``.
- In setup output, 2D mode enforces ``b3=F`` and ``n3=1`` in ``input_new.nml``.
- ``--z-iso-factor`` is only valid together with ``--2d``.

Recursion Control Notes
~~~~~~~~~~~~~~~~~~~~~~~

- ``--rec-method`` maps to ``recur`` in the ``&control`` namelist.
- ``--rec-steps`` maps to ``lld`` in the ``&control`` namelist.
- Defaults:
  - method: ``block``
  - steps: ``20``
- Special default:
  - if ``--rec-method chebyshev`` is selected and ``--rec-steps`` is omitted, steps default to ``200``.

Setup Safety Workflow
~~~~~~~~~~~~~~~~~~~~~

For robust setup runs, prefer:

.. code-block:: bash

    # dry-style run with checks
    rsutil setup bulk --validate --summary

    # if happy with output, replace files and keep backups
    rsutil setup bulk --in-place --backup

Extended Lattice Key Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The setup reader/writer supports selected extended ``&lattice`` keys used by RS-LMTO variants:

- ``strux_backend``
- ``strux_want_sdot``
- ``strux_solve_scale``

These keys are only written to output if present in the input.


Visualization
--------------

The ``viz`` command group provides DOS plotting and structure visualization.

**Commands:**

- ``rsutil viz dos [FOLDER]``: Plot existing RS-LMTO DOS output files.
- ``rsutil viz lattice [FOLDER]``: Open the interactive lattice viewer.
- ``rsutil viz clust [FOLDER]``: Open the interactive cluster viewer.

DOS plotting defaults to a filled overlay plot. Use ``--layout stack`` for one
panel per DOS file, or ``--no-fill`` to turn off curve filling. The energy
values are in Ry. By default, ``--fermi`` is read from ``&energy fermi`` in
``FOLDER/input.nml``; supplying ``--fermi VALUE`` overrides it. In an
interactive overlay plot, click a legend entry to hide or show that DOS
file's curves. This interaction is unavailable in saved images.

**DOS plotting examples:**

.. code-block:: bash

    rsutil viz dos .
    rsutil viz dos ./calc --fermi 0.7 --species Ni-S Ni-S-1 --output dos.png
    rsutil viz dos ./calc --layout stack --no-fill --xlim -0.5 0.2

**Structure visualization:**

The older flags remain supported for compatibility. Both structure viewers
require the ``visualization`` extra; DOS plotting requires the ``bands`` extra.

**Flags:**

- ``--lattice, -L``: Run visualizer for lattice input.
- ``--clust, -C``: Run visualizer for clust file.
- ``--data``: Directory containing data files (default: current directory).

**Examples:**

.. code-block:: bash

    # Visualize lattice structure
    rsutil --lattice --data /path/to/data
    rsutil viz lattice /path/to/data

    # Visualize clust file
    rsutil --clust --data /path/to/data
    rsutil viz clust /path/to/data

Visualizer Controls
~~~~~~~~~~~~~~~~~~~

- Mouse/trackpad: rotate, pan, and zoom the scene (PyVista defaults)
- Toggles: show/hide spins, bonds, and unit cell guides
- Colors: CPK element colors are used by default (see ``rsutil/viz/cpk_colors.json``)
- Snapshots: save screenshots from the window menu (platform dependent)

Data Conventions
~~~~~~~~~~~~~~~~

- Expected files in ``--data``: ``lattice.nml`` and ``input.nml`` for lattice view; ``clust`` for cluster view
- Magnetic moments and bonds are derived from RS-LMTO inputs when available
- Periodic vs free cluster visualization follows the underlying input setup


Exchange Interactions
----------------------

The ``--exchange`` flag creates exchange interaction input files (J_ij) for spin dynamics simulations.

**Flags:**

- ``--exchange, -J``: Create exchange input files.
- ``--cutoff``: Cutoff distance for exchange interactions in lattice units (default: 1.0).
- ``--unique``: Find one representative per exact periodic pair orbit.
- ``--shell-width``: Retained for compatibility; it does not define pair equivalence (default: 0.1).
- ``--soc``: Use conservative symmetry for the legacy collinear +z magnetic reference.
- ``--full-pairs``: Disable symmetry reduction and generate the complete pair set.
- ``--max-pairs-per-job N`` (alias ``--max-pairs``): Split the exchange pairs of each atom into
  jobs of at most ``N`` pairs, written as ``Atom_1a``, ``Atom_1b``, ... instead of one ``Atom_1``.
  Use this to keep a single RS-LMTO-ASA run within memory limits. ``0`` (default) keeps all pairs
  of an atom in one ``Atom_<n>`` directory. Pairs are sliced in neighbor order, so the slices
  together contain exactly the pairs of the unsliced job. Remove or ignore any earlier
  ``Atom_<n>`` directory when switching layouts.
- ``--data``: Directory for input files (default: current directory).
- ``--rij-out``: Path to write R_ij neighbor mappings (optional). Use with ``--exchange``.
- ``--rij-format``: Format for R_ij output; supported values: ``json`` (machine readable) and
    ``txt`` (human readable, fixed-width columns). Default: ``json``.

**Example:**

.. code-block:: bash

    # Basic exchange calculation
    rsutil --exchange --cutoff 2.5

    # With unique neighbor selection and SOC
    rsutil --exchange --cutoff 3.0 --unique --soc --shell-width 0.15

    # From a specific directory
    rsutil --exchange --data /path/to/data --cutoff 2.0

    # Write R_ij mapping to a JSON file
    rsutil --exchange --cutoff 3.0 --rij-out rij_map.json --rij-format json

    # Write human-readable fixed-width R_ij table
    rsutil --exchange --cutoff 3.0 --rij-out rij_map.txt --rij-format txt --unique

    # Disable reduction for validation
    rsutil --exchange --cutoff 3.0 --unique --full-pairs


Export to UppASD
----------------

The ``--export`` flag exports RS-LMTO inputs to the UppASD (Uppsala Atomistic Spin Dynamics) simulation format. This includes structural information, exchange interactions, and optionally DM vectors and spin-orbit coupling.

**Flags:**

- ``--export, -E``: Export RS-LMTO inputs to UppASD format.
- ``--data``: Directory containing input files and exchange data (default: current directory).
- ``--out``: Output directory for exported files relative to ``--data`` (default: ``uppasd_out``).
- ``--jij``: Path to ``jijso.out`` file containing exchange interactions.
- ``--report``: Path to ``report.out`` providing the atomic spin moments (muB) for the ``momfile`` (default: ``report.out`` in ``--data``). Moments fall back to 1.0 with a warning if the file or an atom entry is missing.
- ``--dmi``: Path to ``dijso.out`` or ``dmi`` file containing DM vector interactions.

**Example:**

.. code-block:: bash

    # Basic export
    rsutil --export --data /path/to/simulation

    # With exchange and DM data
    rsutil --export --data /path/to/simulation \\
           --jij /path/to/jijso.out \\
           --dmi /path/to/dmi \\
           --out my_uppasd_output


Export to LMTO47 CTRL
---------------------

The ``--lmto47`` flag exports RS-LMTO structure/input data into an LMTO47 ``CTRL`` file.

**Flags:**

- ``--lmto47``: Enable LMTO47 CTRL export.
- ``--data``: Directory containing ``lattice.nml`` and ``input.nml``.
- ``--out``: Output directory for generated ``CTRL`` file.
- ``--nspin``: ``NSPIN`` value in LMTO47 output (default: ``2``).
- ``--lmx``: ``LMX`` value written for each LMTO47 ``CLASS`` entry.
- ``--lattice-dist``: Gaussian std.dev for lattice-vector distortion in units of ``alat``.
- ``--basis-dist``: Gaussian std.dev for basis-vector distortion in units of ``alat``.
- ``--alat-scale``: Scalar factor applied to ``alat`` before distortion.
- ``--seed``: Integer seed for reproducible random distortions.
- ``--supercell NX NY NZ``: Build CTRL from a repeated supercell input.
- ``--sqs``: Enable SQS assignment on all sites (requires ``--supercell``).
- ``--sqs-shells`` / ``--sqs-steps`` / ``--sqs-seed``: SQS optimization controls.

**Example:**

.. code-block:: bash

    # Basic LMTO47 export
    rsutil --lmto47 --data /path/to/simulation --out lmto47_out

    # Explicit NSPIN/LMX
    rsutil --lmto47 --data /path/to/simulation --nspin 2 --lmx 3 --out lmto47_out

    # Distorted/scaled LMTO47 export with reproducible seed
    rsutil --lmto47 --data /path/to/simulation --lattice-dist 0.01 --basis-dist 0.02 \
           --alat-scale 1.01 --seed 1234 --out lmto47_out

    # Supercell-aware LMTO47 export (non-destructive unless --in-place is set)
    rsutil --lmto47 --data /path/to/simulation --supercell 2 2 2 --out lmto47_out

    # LMTO47 export from SQS supercell
    rsutil --lmto47 --data /path/to/simulation --supercell 2 2 2 \
           --sqs "Fe=0.5,Ni=0.5" --sqs-shells 3 --sqs-steps 5000 --sqs-seed 1234 \
           --out lmto47_out


LMAX Semantics
--------------

The ``--lmx`` option has two contexts:

- Setup context (``--setup``): overrides RS-LMTO generation ``lmax``.
- LMTO47 context (``--lmto47``): writes LMTO47 ``LMX`` in ``CLASS`` entries.

Current RS-LMTO backend support in this repository is ``0 <= lmax <= 3``.


High-throughput campaigns
-------------------------

The ``campaign`` command group prepares deterministic, setup-owned run trees
for review and external execution.  It is independent of the legacy ``-S``
and ``-J`` flags and does not submit jobs or query Slurm.

.. code-block:: bash

    # Validate without touching the destination
    rsutil campaign validate examples/campaigns/ht18_pilot.json

    # Materialize isolated setup runs and manifests
    rsutil campaign create examples/campaigns/ht18_pilot.json \
        --destination ht18-runs

    # Render, but do not submit, a site-specific array helper
    rsutil campaign render-slurm ht18-runs/campaign_index.json \
        --executable /opt/rs-lmto/bin/rs-lmto \
        --launcher 'srun --exclusive' \
        --partition compute --account project --time 02:00:00 \
        --resource cpus_per_task=8 --env-preamble 'module load gcc'

    # Inspect local evidence and collect structured results
    rsutil campaign status ht18-runs
    rsutil campaign collect ht18-runs \
        --csv ht18-runs/results.csv --jsonl ht18-runs/results.jsonl

The generated ``campaign_index.json``, ``run_manifest.json``,
``slurm_run_list.tsv``, ``submit_array.sbatch``, ``rsutil_execution.json``,
and ``campaign_results.json`` artifacts are versioned and relocatable.  Status
keeps execution, parser/output validity, and scientific convergence separate;
``COMPLETED_VALID`` requires all three success conditions.  Collection keeps
full parsed models in JSON sidecars and records per-entry errors without
stopping the campaign.


Quick Reference
----------------

.. code-block:: bash

    # Show all available options
    rsutil --help

    # Typical workflow: setup → exchange → export
    rsutil setup bulk --empty_spheres
    rsutil --exchange --unique --cutoff 2.5
    rsutil --export --jij jijso.out
