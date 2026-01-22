#!/usr/bin/env python3
"""Single-entry CLI for rs_band providing subcommands:

  rs-band dos ...
  rs-band moments ...
  rs-band complete ...
  rs-band projections ...
  rs-band dump-hamiltonian ...

The top-level help prints aggregated help for all subcommands.
"""

import argparse
import sys
import numpy as np


def _make_common_parser():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument('--lattice', default='lattice.nml', help='Lattice file (default: lattice.nml)')
    p.add_argument('--hamiltonian', default='rs2txt.dat', help='Hamiltonian file (default: rs2txt.dat)')
    return p


def run_dos_subcommand(args):
    # If electrons are specified, use programmatic API for Fermi level determination
    if args.electrons is not None:
        from .band_structure_calculator import LatticeReader, HamiltonianReader
        from .full_bz_dos import FullBZCalculator, DOSCalculator

        lattice = LatticeReader(args.lattice)
        hamiltonian = HamiltonianReader(args.hamiltonian)
        bz_calc = FullBZCalculator(lattice, hamiltonian)
        dos_calc = DOSCalculator(bz_calc)

        energy_range_val = tuple(args.energy_range) if args.energy_range else None  # type: ignore
        dos_data, moments_calc = dos_calc.calculate_dos_with_moments(
            mesh_size=tuple(args.mesh),
            target_electrons=args.electrons,
            temperature=300.0,  # Default temperature
            energy_range=energy_range_val,  # type: ignore
            num_energies=args.num_energies,
            use_symmetry=not args.no_symmetry,
            use_tetrahedron=not args.no_tetrahedron,
            gaussian_smearing=args.gaussian_smearing,
            output_prefix=(args.output_prefix or 'dos')
        )

        # Print Fermi level information
        fermi_level = moments_calc.fermi_level
        if fermi_level is not None:
            print(f"Fermi level determined: {fermi_level:.6f} eV")
            print(f"Total electrons: {args.electrons}")

            # Plot with determined Fermi level
            dos_calc.plot_dos(dos_data, fermi_energy=fermi_level,  # type: ignore
                             output_file=f"{args.output_prefix or 'dos'}.png")
        else:
            print("Warning: Could not determine Fermi level")
            dos_calc.plot_dos(dos_data, fermi_energy=0.0,
                             output_file=f"{args.output_prefix or 'dos'}.png")

    else:
        # Reuse the existing full_bz_dos CLI entry by building sys.argv
        from . import full_bz_dos

        sys_argv = ['rs-band',]
        if args.lattice:
            sys_argv += ['--lattice', args.lattice]
        if args.hamiltonian:
            sys_argv += ['--hamiltonian', args.hamiltonian]
        if args.mesh:
            sys_argv += ['--mesh', str(args.mesh[0]), str(args.mesh[1]), str(args.mesh[2])]
        if args.energy_range:
            sys_argv += ['--energy-range', str(args.energy_range[0]), str(args.energy_range[1])]
        if args.num_energies is not None:
            sys_argv += ['--num-energies', str(args.num_energies)]
        if args.fermi_energy is not None:
            sys_argv += ['--fermi-energy', str(args.fermi_energy)]
        if args.no_symmetry:
            sys_argv += ['--no-symmetry']
        if args.no_tetrahedron:
            sys_argv += ['--no-tetrahedron']
        if args.gaussian_smearing != 0.1:  # Only pass if different from default
            sys_argv += ['--gaussian-smearing', str(args.gaussian_smearing)]
        if args.output_prefix:
            sys_argv += ['--output-prefix', args.output_prefix]

        sys.argv = sys_argv
        full_bz_dos.main()


def run_moments_subcommand(args):
    # Perform moments analysis using programmatic API
    from .band_structure_calculator import LatticeReader, HamiltonianReader
    from .full_bz_dos import FullBZCalculator, DOSCalculator

    lattice = LatticeReader(args.lattice)
    hamiltonian = HamiltonianReader(args.hamiltonian)
    bz_calc = FullBZCalculator(lattice, hamiltonian)
    dos_calc = DOSCalculator(bz_calc)

    energy_range_val = tuple(args.energy_range) if args.energy_range else None
    dos_calc.calculate_dos_with_moments(
        mesh_size=tuple(args.mesh),
        target_electrons=args.electrons,
        temperature=args.temperature,
        energy_range=energy_range_val,
        num_energies=args.num_energies,
        use_symmetry=not args.no_symmetry,
        use_tetrahedron=not args.gaussian,
        output_prefix=(args.output_prefix or 'dos_moments')
    )


def run_complete_subcommand(args):
    # Forward to existing complete_analysis main function
    from . import complete_analysis

    sys_argv = ['rs-band-complete']
    if args.lattice:
        sys_argv += ['--lattice', args.lattice]
    if args.hamiltonian:
        sys_argv += ['--hamiltonian', args.hamiltonian]
    if args.dos_mesh:
        sys_argv += ['--dos-mesh', str(args.dos_mesh[0]), str(args.dos_mesh[1]), str(args.dos_mesh[2])]
    if args.energy_range:
        sys_argv += ['--energy-range', str(args.energy_range[0]), str(args.energy_range[1])]
    if args.num_energies is not None:
        sys_argv += ['--num-energies', str(args.num_energies)]
    if args.fermi_energy is not None:
        sys_argv += ['--fermi-energy', str(args.fermi_energy)]
    if args.skip_dos:
        sys_argv += ['--skip-dos']
    if args.skip_projections:
        sys_argv += ['--skip-projections']
    if args.no_symmetry:
        sys_argv += ['--no-symmetry']
    if args.output_prefix:
        sys_argv += ['--output-prefix', args.output_prefix]

    sys.argv = sys_argv
    complete_analysis.main()


def run_projections_subcommand(args):
    from . import projected_bands

    sys_argv = ['rs-band-projections']
    if args.lattice:
        sys_argv += ['--lattice', args.lattice]
    if args.hamiltonian:
        sys_argv += ['--hamiltonian', args.hamiltonian]
    if args.projection:
        sys_argv += ['--projection', args.projection]
    if args.output_prefix:
        sys_argv += ['--output-prefix', args.output_prefix]

    sys.argv = sys_argv
    projected_bands.main()


def run_dump_hamiltonian_subcommand(args):
    """Dump H(k) matrices at Gamma and H points for BCC structure"""
    from .band_structure_calculator import LatticeReader, HamiltonianReader, BandStructureCalculator
    
    lattice = LatticeReader(args.lattice)
    hamiltonian = HamiltonianReader(args.hamiltonian)
    calculator = BandStructureCalculator(lattice, hamiltonian)
    
    # First dump the standard Gamma and H points
    calculator.dump_hamiltonian_at_points(args.output_prefix)
    
    # If custom points are specified, dump those too
    if args.custom_points:
        custom_points = []
        for point_str in args.custom_points:
            parts = point_str.split(':')
            if len(parts) == 2:
                name = parts[0].strip()
                coords = [float(x.strip()) for x in parts[1].split(',')]
                if len(coords) == 3:
                    custom_points.append((name, np.array(coords)))
        
        if custom_points:
            # print(f"\nDumping {len(custom_points)} custom k-point(s)...")
            calculator.dump_hamiltonian_at_custom_points(custom_points, args.output_prefix)


def main():
    parent = _make_common_parser()

    parser = argparse.ArgumentParser(
        prog='rs-band',
        description='Unified CLI for rs_band — subcommands: dos, moments, complete, projections, dump-hamiltonian, structure',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Top-level flags to support single-command mode (e.g. `rsband --flags`)
    parser.add_argument('--dos', action='store_true', help='Run DOS calculation (shortcut for subcommand dos)')
    parser.add_argument('--moments', action='store_true', help='Run DOS moments analysis (shortcut)')
    parser.add_argument('--complete', action='store_true', help='Run complete analysis (shortcut)')
    parser.add_argument('--projections', action='store_true', help='Run projected bands (shortcut)')
    parser.add_argument('--dump-hamiltonian', action='store_true', help='Dump H(k) matrices (shortcut)')

    subparsers = parser.add_subparsers(dest='subcommand', metavar='SUBCOMMAND')

    # DOS subcommand
    p_dos = subparsers.add_parser('dos', parents=[parent], help='Calculate DOS (tetrahedron or Gaussian)')
    p_dos.add_argument('--mesh', nargs=3, type=int, default=[6, 6, 6], help='k-point mesh size (default: 6 6 6)')
    p_dos.add_argument('--energy-range', nargs=2, type=float, help='Energy range (emin emax)')
    p_dos.add_argument('--num-energies', type=int, default=1000, help='Number of energy points (default: 1000)')
    p_dos.add_argument('--fermi-energy', type=float, default=0.0, help='Fermi energy for plots (default: 0.0)')
    p_dos.add_argument('--electrons', type=float, help='Number of electrons (auto-determine Fermi level)')
    p_dos.add_argument('--no-symmetry', action='store_true', help='Disable symmetry reduction')
    p_dos.add_argument('--no-tetrahedron', action='store_true', help='Use Gaussian broadening instead of tetrahedron')
    p_dos.add_argument('--gaussian-smearing', type=float, default=0.1, help='Gaussian smearing width in eV (default: 0.1)')
    p_dos.add_argument('--output-prefix', default='dos', help='Output file prefix (default: dos)')
    p_dos.set_defaults(func=run_dos_subcommand)

    # Moments subcommand
    p_mom = subparsers.add_parser('moments', parents=[parent], help='DOS moments analysis (occupations, centers, widths)')
    p_mom.add_argument('--mesh', nargs=3, type=int, default=[8, 8, 8], help='k-point mesh size (default: 8 8 8)')
    p_mom.add_argument('--electrons', type=float, help='Number of electrons (auto-detect if not specified)')
    p_mom.add_argument('--temperature', type=float, default=300.0, help='Temperature in K (default: 300.0)')
    p_mom.add_argument('--no-symmetry', action='store_true', help='Disable symmetry reduction')
    p_mom.add_argument('--gaussian', action='store_true', help='Use Gaussian broadening instead of tetrahedron')
    p_mom.add_argument('--num-energies', type=int, default=1000, help='Number of energy points (default: 1000)')
    p_mom.add_argument('--energy-range', nargs=2, type=float, help='Energy range (emin emax)')
    p_mom.add_argument('--output-prefix', default='dos_moments', help='Output file prefix (default: dos_moments)')
    p_mom.set_defaults(func=run_moments_subcommand)

    # Complete analysis subcommand
    p_cmp = subparsers.add_parser('complete', parents=[parent], help='Run complete combined analysis (bands, DOS, projections)')
    p_cmp.add_argument('--nk-line', type=int, default=100, help='Number of k-points for band path (default: 100)')
    p_cmp.add_argument('--dos-mesh', nargs=3, type=int, default=[6, 6, 6], help='k-point mesh for DOS (default: 6 6 6)')
    p_cmp.add_argument('--energy-range', nargs=2, type=float, help='Energy range for DOS (emin emax)')
    p_cmp.add_argument('--num-energies', type=int, default=500, help='Number of energy points for DOS (default: 500)')
    p_cmp.add_argument('--fermi-energy', type=float, default=0.0, help='Fermi energy (default: 0.0)')
    p_cmp.add_argument('--skip-dos', action='store_true', help='Skip DOS calculation')
    p_cmp.add_argument('--skip-projections', action='store_true', help='Skip projection calculations')
    p_cmp.add_argument('--no-symmetry', action='store_true', help='Disable symmetry reduction for DOS')
    p_cmp.set_defaults(func=run_complete_subcommand)

    # Projections subcommand
    p_proj = subparsers.add_parser('projections', parents=[parent], help='Calculate projected bands')
    p_proj.add_argument('--projection', choices=['orbital', 'site', 'spin', 'all'], default='all', help='Type of projection (default: all)')
    p_proj.add_argument('--output-prefix', default='projected', help='Output file prefix (default: projected)')
    p_proj.set_defaults(func=run_projections_subcommand)

    # Dump Hamiltonian subcommand
    p_dump = subparsers.add_parser('dump-hamiltonian', parents=[parent], help='Dump H(k) matrices at Gamma and H points (BCC)')
    p_dump.add_argument('--output-prefix', default='H_k', help='Output file prefix (default: H_k)')
    p_dump.add_argument('--custom-points', nargs='*', 
                       help='Custom k-points in format "Name:x,y,z" (e.g., "K_half:0.375,0.125,0.125")')
    p_dump.set_defaults(func=run_dump_hamiltonian_subcommand)

    # structure and quick can be forwarded to existing modules if present
    subparsers.add_parser('structure', help='Analyze crystal structure (k-paths)')
    subparsers.add_parser('quick', help='Run quick analysis')

    # If no args provided, show aggregated help
    if len(sys.argv) == 1:
        parser.print_help()
        print('\nAvailable subcommands and their options:\n')
        for name, sub in [('dos', p_dos), ('moments', p_mom), ('complete', p_cmp), ('projections', p_proj), ('dump-hamiltonian', p_dump)]:
            print(f"--- {name} ---")
            sub.print_help()
            print()
        return

    args = parser.parse_args()

    # If a top-level shortcut flag was used, dispatch to the corresponding handler
    # by building a minimal args namespace and calling the existing functions.
    from types import SimpleNamespace

    if getattr(args, 'dos', False):
        ns = SimpleNamespace(
            lattice=args.lattice,
            hamiltonian=args.hamiltonian,
            mesh=[6, 6, 6],
            energy_range=None,
            num_energies=1000,
            fermi_energy=0.0,
            electrons=None,
            no_symmetry=False,
            no_tetrahedron=False,
            gaussian_smearing=0.1,
            output_prefix='dos'
        )
        return run_dos_subcommand(ns)

    if getattr(args, 'moments', False):
        ns = SimpleNamespace(
            lattice=args.lattice,
            hamiltonian=args.hamiltonian,
            mesh=[8, 8, 8],
            electrons=None,
            temperature=300.0,
            no_symmetry=False,
            gaussian=False,
            num_energies=1000,
            energy_range=None,
            output_prefix='dos_moments'
        )
        return run_moments_subcommand(ns)

    if getattr(args, 'complete', False):
        ns = SimpleNamespace(
            lattice=args.lattice,
            hamiltonian=args.hamiltonian,
            nk_line=100,
            dos_mesh=[6, 6, 6],
            energy_range=None,
            num_energies=500,
            fermi_energy=0.0,
            skip_dos=False,
            skip_projections=False,
            no_symmetry=False,
            output_prefix='complete_analysis'
        )
        return run_complete_subcommand(ns)

    if getattr(args, 'projections', False):
        ns = SimpleNamespace(
            lattice=args.lattice,
            hamiltonian=args.hamiltonian,
            projection='all',
            output_prefix='projected'
        )
        return run_projections_subcommand(ns)

    if getattr(args, 'dump_hamiltonian', False):
        ns = SimpleNamespace(
            lattice=args.lattice,
            hamiltonian=args.hamiltonian,
            output_prefix='H_k',
            custom_points=None
        )
        return run_dump_hamiltonian_subcommand(ns)

    # Dispatch to subcommand
    if hasattr(args, 'func'):
        return args.func(args)

    # Handle simple forwards for structure/quick
    if args.subcommand == 'structure':
        from .structure_analysis import analyze_structure
        # Use defaults if not provided via earlier parent
        analyze_structure('lattice.nml', True, 'structure_analysis.png')
        return
    if args.subcommand == 'quick':
        try:
            from .quick_analysis import main as quick_main
            quick_main()
        except (ImportError, RuntimeError):
            print('Quick analysis not available in this installation')
        return


if __name__ == '__main__':
    main()
