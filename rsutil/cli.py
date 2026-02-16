"""
This module provides a command-line interface (CLI) for the rsutil package.
It allows users to process lattice structures, visualize data, and create
exchange input files through various command-line arguments.

Functions
---------
    Sets up the command-line argument parser and handles the execution of
    different functionalities based on the provided arguments.

Usage
-----
To use this CLI, run the script with the desired flags:
    python cli.py --setup
    python cli.py --lattice --data <data_directory>
    python cli.py --clust --data <data_directory>
    python cli.py --exchange

Flags
-----
--setup, -S
    Process lattice.
--lattice, -L
    Run visualizer for lattice input.
--clust, -C
    Run visualizer for clust file.
--exchange, -J
    Create exchange input files.
--data
    Directory for data files to be visualized (default: current directory).
"""

import argparse
import sys
from pathlib import Path
from . import process_lattice
from .viz.rs_visualizer import rs_visualizer
from .exchange.exchange_map import exchange_map
from .export.uppasd import export_from_files
import os


def validate_setup_inputs():
    """
    Validate that required input files exist for --setup workflow.
    
    Returns
    -------
    bool
        True if all required files exist, False otherwise
    """
    required_files = ['lattice.nml', 'input.nml']
    missing_files = []
    
    for filename in required_files:
        if not Path(filename).exists():
            missing_files.append(filename)
    
    if missing_files:
        print(f"ERROR: Required input files not found: {', '.join(missing_files)}", file=sys.stderr)
        print("\nThe --setup workflow requires the following files in the current directory:", file=sys.stderr)
        print("  - lattice.nml   : Lattice structure data (lattice_t namelist)", file=sys.stderr)
        print("  - input.nml     : Atomic positions and parameters (&lattice, &atoms namelists)", file=sys.stderr)
        print("\nPlease ensure these files exist before running --setup.", file=sys.stderr)
        return False
    
    return True


def validate_setup_parameters(min_radius, max_radius):
    """
    Validate setup parameter ranges.
    
    Parameters
    ----------
    min_radius : float
        Minimum empty sphere radius
    max_radius : float
        Maximum empty sphere radius
    
    Returns
    -------
    bool
        True if parameters are valid, False otherwise
    """
    if min_radius <= 0:
        print(f"ERROR: --min_empty_radius must be positive, got {min_radius}", file=sys.stderr)
        return False
    
    if max_radius <= 0:
        print(f"ERROR: --max_empty_radius must be positive, got {max_radius}", file=sys.stderr)
        return False
    
    if min_radius >= max_radius:
        print(f"ERROR: --min_empty_radius ({min_radius}) must be less than --max_empty_radius ({max_radius})", file=sys.stderr)
        return False
    
    if max_radius > 10.0:
        print(f"WARNING: --max_empty_radius ({max_radius}) is very large. Typical values are 1.5-3.0.", file=sys.stderr)
    
    if min_radius < 0.5:
        print(f"WARNING: --min_empty_radius ({min_radius}) is very small. Typical values are 1.0-1.5.", file=sys.stderr)
    
    return True


def validate_output_files():
    """
    Validate that output files were created and are readable.
    
    Returns
    -------
    dict
        Dictionary with validation results for each file
    """
    validation_results = {}
    expected_files = ['lattice_new.nml', 'input_new.nml']
    
    for filename in expected_files:
        filepath = Path(filename)
        result = {
            'exists': filepath.exists(),
            'readable': False,
            'size': 0,
            'valid_format': False
        }
        
        if result['exists']:
            result['size'] = filepath.stat().st_size
            try:
                with open(filepath, 'r') as f:
                    content = f.read()
                    result['readable'] = True
                    # Basic format check - should contain namelist markers
                    result['valid_format'] = '&' in content and '/' in content
            except Exception:
                pass
        
        validation_results[filename] = result
    
    return validation_results


def print_validation_report(validation_results):
    """Print validation report for output files."""
    print("\nValidation Report:")
    print("-" * 70)
    
    all_valid = True
    for filename, result in validation_results.items():
        status = "✓ PASS" if result['exists'] and result['readable'] and result['valid_format'] else "✗ FAIL"
        print(f"{status} {filename}")
        
        if not result['exists']:
            print(f"     ERROR: File not created")
            all_valid = False
        elif not result['readable']:
            print(f"     ERROR: File not readable")
            all_valid = False
        elif not result['valid_format']:
            print(f"     WARNING: File may not be valid namelist format")
        else:
            print(f"     Size: {result['size']} bytes")
    
    print("-" * 70)
    if all_valid:
        print("All validation checks passed ✓")
    else:
        print("Some validation checks failed ✗")
    
    return all_valid


def create_backup_files():
    """Create backup copies of input files."""
    input_files = ['lattice.nml', 'input.nml']
    backed_up = []
    
    for filename in input_files:
        if Path(filename).exists():
            backup_path = Path(f"{filename}.bak")
            try:
                import shutil
                shutil.copy2(filename, backup_path)
                backed_up.append(str(backup_path))
            except Exception as e:
                print(f"WARNING: Could not create backup of {filename}: {e}", file=sys.stderr)
    
    if backed_up:
        print(f"Created backups: {', '.join(backed_up)}")
    
    return backed_up


def apply_in_place_updates():
    """Apply in-place updates by renaming _new files to originals."""
    updates = [
        ('lattice_new.nml', 'lattice.nml'),
        ('input_new.nml', 'input.nml')
    ]
    
    updated = []
    for new_file, orig_file in updates:
        if Path(new_file).exists():
            # Create backup first
            if Path(orig_file).exists():
                backup_path = Path(f"{orig_file}.bak")
                try:
                    import shutil
                    shutil.copy2(orig_file, backup_path)
                    print(f"Backup created: {backup_path}")
                except Exception as e:
                    print(f"WARNING: Could not create backup of {orig_file}: {e}", file=sys.stderr)
            
            # Replace original with new
            try:
                import shutil
                shutil.move(new_file, orig_file)
                updated.append(orig_file)
                print(f"Updated: {orig_file}")
            except Exception as e:
                print(f"ERROR: Could not update {orig_file}: {e}", file=sys.stderr)
    
    return updated


def create_supercell(nx, ny, nz):
    """
    Create a supercell by repeating the unit cell nx×ny×nz times.
    
    Reads lattice.nml, creates supercell, and writes updated file.
    Coordinates in lattice.nml are Cartesian (normalized, lattice parameter = 1).
    
    Parameters
    ----------
    nx, ny, nz : int
        Repetitions along each lattice vector
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    import re
    import shutil
    import numpy as np
    
    print(f"\nCreating {nx}×{ny}×{nz} supercell...")
    
    lattice_file = Path('lattice.nml')
    if not lattice_file.exists():
        print("ERROR: lattice.nml not found", file=sys.stderr)
        return False
    
    try:
        # === Process lattice.nml ===
        with open(lattice_file, 'r') as f:
            content = f.read()
        
        # Extract lattice vectors (a(:,1), a(:,2), a(:,3))
        lattice_vectors = []
        for i in [1, 2, 3]:
            pattern = rf'a\s*\(\s*:\s*,\s*{i}\s*\)\s*=\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)'
            match = re.search(pattern, content)
            if match:
                vec = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                lattice_vectors.append(vec)
        
        if len(lattice_vectors) != 3:
            print(f"ERROR: Could not parse all 3 lattice vectors (found {len(lattice_vectors)})", file=sys.stderr)
            return False
        
        a_vec, b_vec, c_vec = np.array(lattice_vectors)
        
        # Extract atomic positions (Cartesian coordinates, normalized)
        positions = []
        crd_pattern = r'crd\s*\(\s*:\s*,\s*(\d+)\s*\)\s*=\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)'
        for match in re.finditer(crd_pattern, content):
            atom_idx = int(match.group(1))
            x, y, z = float(match.group(2)), float(match.group(3)), float(match.group(4))
            positions.append(np.array([x, y, z]))
        
        if not positions:
            print("ERROR: Could not parse atomic positions from lattice.nml", file=sys.stderr)
            return False
        
        n_atoms_orig = len(positions)
        print(f"  Original unit cell: {n_atoms_orig} atoms")
        
        # Create supercell positions (Cartesian)
        # For each repetition (ix, iy, iz), shift by ix*a + iy*b + iz*c
        supercell_positions = []
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for pos in positions:
                        new_pos = pos + ix * a_vec + iy * b_vec + iz * c_vec
                        supercell_positions.append(new_pos)
        
        n_atoms_super = len(supercell_positions)
        print(f"  Supercell: {n_atoms_super} atoms ({nx}×{ny}×{nz} = {nx*ny*nz} unit cells)")
        
        # Update nbas and related counts
        content = re.sub(r'(nbulk_bulk\s*=\s*)\d+', f'nbulk_bulk = {n_atoms_super}', content, flags=re.IGNORECASE)
        content = re.sub(r'(ntot\s*=\s*)\d+', f'ntot = {n_atoms_super}', content, flags=re.IGNORECASE)
        content = re.sub(r'(nbas\s*=\s*)\d+', f'nbas = {n_atoms_super}', content, flags=re.IGNORECASE)
        content = re.sub(r'(nrec\s*=\s*)\d+', f'nrec = {n_atoms_super}', content, flags=re.IGNORECASE)
        
        # Scale lattice vectors
        content = re.sub(
            r'(a\s*\(\s*:\s*,\s*1\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
            lambda m: f"{m.group(1)}{float(m.group(2))*nx:.10f}, {float(m.group(3))*nx:.10f}, {float(m.group(4))*nx:.10f}",
            content,
            flags=re.IGNORECASE
        )
        
        content = re.sub(
            r'(a\s*\(\s*:\s*,\s*2\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
            lambda m: f"{m.group(1)}{float(m.group(2))*ny:.10f}, {float(m.group(3))*ny:.10f}, {float(m.group(4))*ny:.10f}",
            content,
            flags=re.IGNORECASE
        )
        
        content = re.sub(
            r'(a\s*\(\s*:\s*,\s*3\s*\)\s*=\s*)([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
            lambda m: f"{m.group(1)}{float(m.group(2))*nz:.10f}, {float(m.group(3))*nz:.10f}, {float(m.group(4))*nz:.10f}",
            content,
            flags=re.IGNORECASE
        )
        
        # Replace atomic positions
        crd_start = re.search(r'crd\s*\(\s*:\s*,\s*1\s*\)', content)
        if not crd_start:
            print("ERROR: Could not find coordinate section", file=sys.stderr)
            return False
        
        # Find end of crd section
        crd_section_start = crd_start.start()
        lines = content[crd_section_start:].split('\n')
        crd_lines_count = sum(1 for line in lines if 'crd(' in line)
        
        temp_pos = crd_section_start
        for i in range(crd_lines_count):
            temp_pos = content.find('\n', temp_pos) + 1
        crd_section_end = temp_pos
        
        # Build new coordinate section
        new_crd_section = ""
        for i, pos in enumerate(supercell_positions, start=1):
            new_crd_section += f"crd(:, {i}) = {pos[0]:.10f}, {pos[1]:.10f}, {pos[2]:.10f}\n"
        
        new_content = content[:crd_section_start] + new_crd_section + content[crd_section_end:]
        
        # Update arrays that depend on atom count (broken symmetry: 1..N)
        for array_name in ['izp', 'no', 'iu', 'ib', 'irec']:
            array_start = re.search(rf'{array_name}\s*\(\s*1\s*\)', new_content)
            if array_start:
                array_section_start = array_start.start()
                lines = new_content[array_section_start:].split('\n')
                array_lines = sum(1 for line in lines if f'{array_name}(' in line)
                
                temp_pos = array_section_start
                for i in range(array_lines):
                    temp_pos = new_content.find('\n', temp_pos) + 1
                array_section_end = temp_pos
                
                new_array_section = ""
                for i in range(1, n_atoms_super + 1):
                    new_array_section += f"{array_name}({i}) = {i}\n"
                
                new_content = new_content[:array_section_start] + new_array_section + new_content[array_section_end:]
        
        # Backup and write lattice.nml
        backup_file = Path('lattice_original.nml.bak')
        shutil.copy2('lattice.nml', backup_file)
        
        with open('lattice.nml', 'w') as f:
            f.write(new_content)
        
        print(f"  ✓ Updated lattice.nml with supercell (backup: {backup_file})")
        
        # === Update input.nml to replicate labels ===
        input_file = Path('input.nml')
        if input_file.exists():
            with open(input_file, 'r') as f:
                input_content = f.read()
            
            # Extract labels
            label_pattern = r"label\s*\(\s*(\d+)\s*\)\s*=\s*['\"]([^'\"]+)['\"]"
            orig_labels = {}
            for match in re.finditer(label_pattern, input_content):
                idx = int(match.group(1))
                label = match.group(2)
                if idx <= n_atoms_orig:
                    orig_labels[idx] = label
            
            # Extract ct values (match original label ordering)
            ct_pattern = r'ct\s*\(\s*(\d+)\s*\)\s*=\s*([\-\d.eE+]+)'
            orig_ct = {}
            for match in re.finditer(ct_pattern, input_content):
                idx = int(match.group(1))
                val = match.group(2)
                orig_ct[idx] = val
            
            if orig_labels:
                # Replicate labels for supercell (unique labels, broken symmetry)
                orig_label_keys = sorted(orig_labels.keys())
                label_meta = []
                for key in orig_label_keys:
                    label = orig_labels[key]
                    match = re.match(r'^(.*?)(\d+)$', label)
                    if match:
                        base = match.group(1)
                        start_num = int(match.group(2))
                    else:
                        base = label
                        start_num = 0
                    label_meta.append((base, start_num, label))
                
                orig_ct_values = []
                if orig_ct:
                    for key in orig_label_keys:
                        if key in orig_ct:
                            orig_ct_values.append(orig_ct[key])
                
                label_counts = {}
                super_labels = []
                super_ct = []
                source_labels = []
                for _ in range(nx * ny * nz):
                    for idx, key in enumerate(orig_label_keys):
                        base, start_num, src_label = label_meta[idx]
                        if base not in label_counts:
                            label_counts[base] = start_num - 1 if start_num > 0 else 0
                        label_counts[base] += 1
                        super_labels.append(f"{base}{label_counts[base]}")
                        source_labels.append(src_label)
                        if orig_ct_values:
                            super_ct.append(orig_ct_values[idx])
                new_ntype = len(super_labels)
                
                # Find label section and replace
                label_start = re.search(r"label\s*\(\s*1\s*\)", input_content)
                if label_start:
                    # Find end of label section
                    label_section_start = label_start.start()
                    lines = input_content[label_section_start:].split('\n')
                    label_lines = sum(1 for line in lines if 'label(' in line)
                    
                    temp_pos = label_section_start
                    for i in range(label_lines):
                        temp_pos = input_content.find('\n', temp_pos) + 1
                    label_section_end = temp_pos
                    
                    # Build new label section
                    new_label_section = ""
                    for i, label in enumerate(super_labels, start=1):
                        new_label_section += f"label({i}) = '{label}'\n"
                    
                    input_content = input_content[:label_section_start] + new_label_section + input_content[label_section_end:]
                    
                    # Update ntype in &lattice
                    input_content = re.sub(
                        r'(ntype\s*=\s*)\d+',
                        f'ntype       = {new_ntype}',
                        input_content,
                        flags=re.IGNORECASE
                    )
                    
                    # Expand ct(i) definitions to match ntype
                    if orig_ct and super_ct:
                        ct_start = re.search(r'ct\s*\(\s*1\s*\)', input_content)
                        if ct_start:
                            ct_section_start = ct_start.start()
                            lines = input_content[ct_section_start:].split('\n')
                            ct_lines = sum(1 for line in lines if 'ct(' in line)
                            
                            temp_pos = ct_section_start
                            for i in range(ct_lines):
                                temp_pos = input_content.find('\n', temp_pos) + 1
                            ct_section_end = temp_pos
                            
                            new_ct_section = ""
                            for i, val in enumerate(super_ct, start=1):
                                new_ct_section += f"ct({i}) = {val}\n"
                            
                            input_content = input_content[:ct_section_start] + new_ct_section + input_content[ct_section_end:]
                    
                    # Backup and write
                    input_backup = Path('input_original.nml.bak')
                    shutil.copy2('input.nml', input_backup)
                    
                    with open('input.nml', 'w') as f:
                        f.write(input_content)
                    
                    print(f"  ✓ Updated input.nml with replicated labels (backup: {input_backup})")
                    
                    # Create per-species .nml files for unique labels
                    for new_label, src_label in zip(super_labels, source_labels):
                        src_file = Path(f"{src_label}.nml")
                        dest_file = Path(f"{new_label}.nml")
                        if src_file.exists() and not dest_file.exists():
                            shutil.copy2(src_file, dest_file)
        
        print(f"  ✓ Supercell creation complete")
        
        return True
        
    except Exception as e:
        print(f"ERROR creating supercell: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False


def print_summary_diff():
    """Print a summary of changes between old and new files."""
    print("\nSummary of Changes:")
    print("=" * 70)
    
    # Compare lattice files
    old_lattice = Path('lattice.nml')
    new_lattice = Path('lattice_new.nml')
    
    if old_lattice.exists() and new_lattice.exists():
        old_size = old_lattice.stat().st_size
        new_size = new_lattice.stat().st_size
        print(f"\nlattice.nml:")
        print(f"  Old size: {old_size} bytes")
        print(f"  New size: {new_size} bytes")
        print(f"  Change:   {new_size - old_size:+d} bytes")
    
    # Compare input files  
    old_input = Path('input.nml')
    new_input = Path('input_new.nml')
    
    if old_input.exists() and new_input.exists():
        old_size = old_input.stat().st_size
        new_size = new_input.stat().st_size
        print(f"\ninput.nml:")
        print(f"  Old size: {old_size} bytes")
        print(f"  New size: {new_size} bytes")
        print(f"  Change:   {new_size - old_size:+d} bytes")
        
        # Try to extract some key info
        try:
            with open(new_input, 'r') as f:
                content = f.read()
                if 'nbas' in content:
                    # Try to find nbas value
                    for line in content.split('\n'):
                        if 'nbas' in line.lower() and '=' in line:
                            print(f"  {line.strip()}")
                            break
        except Exception:
            pass
    
    print("=" * 70)


def main():
    """
    Command-line interface for rsutil.

    This function sets up the command-line argument parser and handles the
    execution of different functionalities based on the provided arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(
        description="RS Tools Command Line Interface"
    )
    parser.add_argument(
        '--setup', '-S', action='store_true', default=False,
        help='Process lattice'
    )
    parser.add_argument(
        '--lattice', '-L', action='store_true', default=False,
        help='Run visualizer for lattice input'
    )
    parser.add_argument(
        '--clust', '-C', action='store_true', default=False,
        help='Run visualizer for clust file'
    )
    parser.add_argument(
        '--exchange', '-J', action='store_true', default=False,
        help='Create exchange input files'
    )
    parser.add_argument(
        '--data', type=str, default='.',
        help='Directory for data files to be visualized'
    )
    parser.add_argument(
        '--cutoff', type=float, default='1.0',
        help='Cutoff for exchange interactions, in lattice units'
    )
    parser.add_argument(
        '--unique', action='store_true', default=False,
        help='Find only unique neighbors per shell (one representative per unique species-vector combination)'
    )
    parser.add_argument(
        '--shell-width', type=float, default=0.1,
        help='Width of shells for grouping neighbors when using --unique (default: 0.1)'
    )
    parser.add_argument(
        '--soc', action='store_true', default=False,
        help='Account for spin-orbit coupling in symmetry considerations'
    )
    parser.add_argument(
        '--rij-out', type=str, default=None,
        help='Write R_ij mappings to this file (JSON)'
    )
    parser.add_argument(
        '--rij-format', type=str, default='json', choices=['json', 'txt'],
        help='Format for R_ij output (json or txt; default: json)'
    )
    parser.add_argument(
        '--export', '-E', action='store_true', default=False,
        help='Export RS-LMTO inputs to UppASD format'
    )
    parser.add_argument(
        '--jij', type=str, default=None,
        help='Path to jijso.out file to include in export'
    )
    parser.add_argument(
        '--dmi', type=str, default=None,
        help='Path to dijso.out / dmi file to include in export (DM vectors)'
    )
    parser.add_argument(
        '--out', type=str, default='uppasd_out',
        help='Output directory for exported UppASD files (relative to --data)'
    )
    parser.add_argument(
        '--cluster', type=str, default='periodic', choices=['periodic', 'free'],
        help='Cluster type: periodic (default) or free'
    )
    parser.add_argument(
        '--empty_spheres', action='store_true', default=False,
        help='Enable empty sphere calculations'
    )
    parser.add_argument(
        '--min_empty_radius', type=float, default=1.2,
        help='Minimum empty sphere radius (default: 1.2)'
    )
    parser.add_argument(
        '--max_empty_radius', type=float, default=2.5,
        help='Maximum empty sphere radius (default: 2.5)'
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true', default=False,
        help='Show detailed progress during processing'
    )
    parser.add_argument(
        '--summary', action='store_true', default=False,
        help='Print summary of changes after processing'
    )
    parser.add_argument(
        '--validate', action='store_true', default=False,
        help='Validate output files after generation'
    )
    parser.add_argument(
        '--backup', action='store_true', default=False,
        help='Create backup files before any operations'
    )
    parser.add_argument(
        '--in-place', action='store_true', default=False,
        help='Update input files in place (creates .bak backups)'
    )
    parser.add_argument(
        '--supercell', nargs=3, type=int, metavar=('NX', 'NY', 'NZ'),
        help='Create supercell by repeating unit cell NX×NY×NZ times'
    )

    args = parser.parse_args()

    # Process lattice if --setup flag is provided
    if args.setup:
        print("=" * 70)
        print("RS-utils Lattice Setup")
        print("=" * 70)
        
        # Validate input files exist
        if not validate_setup_inputs():
            sys.exit(1)
        
        # Validate parameter ranges
        if not validate_setup_parameters(args.min_empty_radius, args.max_empty_radius):
            sys.exit(1)
        
        # Display configuration
        print(f"\nConfiguration:")
        print(f"  Cluster type:      {args.cluster}")
        print(f"  Empty spheres:     {'enabled' if args.empty_spheres else 'disabled'}")
        if args.empty_spheres:
            print(f"  Min radius:        {args.min_empty_radius}")
            print(f"  Max radius:        {args.max_empty_radius}")
        if args.supercell:
            nx, ny, nz = args.supercell
            print(f"  Supercell:         {nx}×{ny}×{nz}")
        if args.verbose:
            print(f"  Verbose mode:      enabled")
        if args.validate:
            print(f"  Validation:        enabled")
        if args.summary:
            print(f"  Summary:           enabled")
        if args.backup:
            print(f"  Backup:            enabled")
        if args.in_place:
            print(f"  In-place update:   enabled")
        print()
        
        # Create backups if requested
        if args.backup or args.in_place:
            print("Creating backup files...")
            create_backup_files()
            print()
        
        # Create supercell if requested (do this BEFORE processing)
        if args.supercell:
            nx, ny, nz = args.supercell
            if nx < 1 or ny < 1 or nz < 1:
                print("ERROR: Supercell dimensions must be positive integers", file=sys.stderr)
                sys.exit(1)
            
            if not create_supercell(nx, ny, nz):
                print("ERROR: Failed to create supercell", file=sys.stderr)
                sys.exit(1)
            print()
        
        # Process lattice with error handling
        try:
            if args.verbose:
                print("Starting lattice processing...")
            
            process_lattice(
                cluster_type=args.cluster, 
                empty_spheres=args.empty_spheres,
                min_empty_radius=args.min_empty_radius, 
                max_empty_radius=args.max_empty_radius
            )
            
            if args.verbose:
                print("Lattice processing completed.")
            
            # Check if output files were created
            output_files = ['lattice_new.nml', 'input_new.nml']
            created_files = [f for f in output_files if Path(f).exists()]
            
            if created_files:
                print("\n" + "=" * 70)
                print("SUCCESS: Lattice processing completed")
                print("=" * 70)
                print("\nGenerated files:")
                for filename in created_files:
                    file_size = Path(filename).stat().st_size
                    print(f"  ✓ {filename} ({file_size} bytes)")
                
                # Run validation if requested
                if args.validate:
                    validation_results = validate_output_files()
                    print_validation_report(validation_results)
                
                # Print summary if requested
                if args.summary:
                    print_summary_diff()
                
                # Apply in-place updates if requested
                if args.in_place:
                    print("\nApplying in-place updates...")
                    updated = apply_in_place_updates()
                    if updated:
                        print(f"\nUpdated {len(updated)} file(s) in place.")
                        print("Original files backed up with .bak extension.")
                        print("\nNext steps:")
                        print("  Continue with exchange calculations:")
                        print("     rsutil --exchange --cutoff 2.5")
                    else:
                        print("No files were updated.", file=sys.stderr)
                else:
                    # Standard workflow message
                    print("\nNext steps:")
                    print("  1. Review the generated files")
                    if not args.validate:
                        print("     (run with --validate to check file integrity)")
                    if not args.summary:
                        print("     (run with --summary to see detailed changes)")
                    print("  2. Rename them to replace original files:")
                    print("     mv lattice_new.nml lattice.nml")
                    print("     mv input_new.nml input.nml")
                    print("     (or use --in-place flag to do this automatically)")
                    print("  3. Continue with exchange calculations:")
                    print("     rsutil --exchange --cutoff 2.5")
            else:
                print("\nWARNING: Processing completed but expected output files were not found.", file=sys.stderr)
                print("Check for error messages above.", file=sys.stderr)
                
        except FileNotFoundError as e:
            print(f"\nERROR: File not found during processing: {e}", file=sys.stderr)
            print("This may indicate missing data files or incorrect file paths.", file=sys.stderr)
            sys.exit(1)
        except ValueError as e:
            print(f"\nERROR: Invalid parameter or data: {e}", file=sys.stderr)
            sys.exit(1)
        except RuntimeError as e:
            print(f"\nERROR: Processing failed: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"\nERROR: Unexpected error during lattice processing:", file=sys.stderr)
            print(f"  {type(e).__name__}: {e}", file=sys.stderr)
            print("\nPlease check your input files and try again.", file=sys.stderr)
            print("If the problem persists, please report this error.", file=sys.stderr)
            sys.exit(1)
    # Run visualizer for lattice or clust file if --lattice or --clust flag is provided
    elif args.lattice or args.clust:
        rs_visualizer(args.data, clust=args.clust)
    # Create exchange input files if --exchange flag is provided
    elif args.exchange:
        exchange_map(
            cutoff=args.cutoff,
            unique=args.unique,
            shell_width=args.shell_width,
            soc=args.soc,
            rij_out=args.rij_out,
            rij_fmt=args.rij_format,
        )
    # Export to UppASD if requested
    elif args.export:
        # default lattice path in data dir
        lattice_path = os.path.join(args.data, 'lattice.nml')
        out_dir = os.path.join(args.data, args.out) if args.data else args.out
        export_from_files(lattice_nml=lattice_path, input_nml=None, clust=None,
                          out_dir=out_dir, jij_file=args.jij, dij_file=args.dmi)
    # Print help message if no valid flag is provided
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
