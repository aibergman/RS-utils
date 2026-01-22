#!/usr/bin/env python3
"""
Comprehensive DOS moments calculation test for bccFe example

This script demonstrates the complete DOS moments analysis workflow with:
- Multiple k-point mesh sizes for convergence testing
- Multiple temperatures (0K and 300K)
- Detailed results comparison and plotting
- Performance timing analysis
"""

import sys
import os
import time
import numpy as np

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("⚠️  matplotlib not available - plotting will be skipped")

# Add the parent directory to Python path for local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# Local imports after path modification
from rs_band.band_structure_calculator import LatticeReader, HamiltonianReader  # noqa: E402
from rs_band.full_bz_dos import FullBZCalculator, DOSCalculator  # noqa: E402
from rs_band.dos_moments import DOSMomentsCalculator  # noqa: E402


def test_single_calculation(mesh_size, temperature, use_tetrahedron=True, label=""):
    """
    Perform single DOS moments calculation
    
    Returns:
        dict: Results containing timing, Fermi level, and moments
    """
    print(f"\n{'='*60}")
    print(f"Testing {label}")
    print(f"Mesh: {mesh_size}, Temperature: {temperature}K, Method: {'Tetrahedron' if use_tetrahedron else 'Gaussian'}")
    print(f"{'='*60}")
    
    # Read input files (only once at start)
    if not hasattr(test_single_calculation, '_lattice'):
        print("📖 Reading input files...")
        test_single_calculation._lattice = LatticeReader('lattice.nml')
        test_single_calculation._hamiltonian = HamiltonianReader('rs2txt.dat')
        print("✅ Input files loaded")
    
    # Create calculators
    bz_calc = FullBZCalculator(test_single_calculation._lattice, test_single_calculation._hamiltonian)
    dos_calc = DOSCalculator(bz_calc)
    
    # Time the DOS calculation
    start_time = time.time()
    print("⚡ Calculating DOS...")
    
    dos_data = dos_calc.calculate_dos(
        mesh_size=mesh_size,
        use_symmetry=True,
        use_tetrahedron=use_tetrahedron
    )
    
    dos_time = time.time() - start_time
    print(f"⏱️  DOS calculation completed in {dos_time:.2f} seconds")
    
    if dos_data['total_dos'] is None:
        raise RuntimeError("DOS calculation failed")
    
    # Prepare projected DOS
    projected_dos = {}
    for key, dos in dos_data['orbital_dos'].items():
        projected_dos[f'orbital_{key}'] = dos
    for key, dos in dos_data['spin_dos'].items():
        projected_dos[f'spin_{key}'] = dos
    
    # Time the moments calculation
    start_time = time.time()
    print("🧮 Calculating DOS moments...")
    
    # Create moments calculator
    moments_calc = DOSMomentsCalculator(
        dos_data['energies'],
        dos_data['total_dos'],
        projected_dos,
        temperature=temperature
    )
    
    # Find Fermi level
    fermi_level = moments_calc.find_fermi_level(target_electrons=8.0)
    
    # Calculate moments
    moments = moments_calc.calculate_moments()
    
    moments_time = time.time() - start_time
    print(f"⏱️  Moments calculation completed in {moments_time:.2f} seconds")
    
    # Print summary
    moments_calc.print_summary()
    
    # Save results
    output_prefix = f"bccFe_{mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}_T{temperature:.0f}K"
    moments_calc.save_moments(output_prefix)
    print(f"💾 Results saved with prefix: {output_prefix}")
    
    return {
        'mesh_size': mesh_size,
        'temperature': temperature,
        'method': 'tetrahedron' if use_tetrahedron else 'gaussian',
        'dos_time': dos_time,
        'moments_time': moments_time,
        'total_time': dos_time + moments_time,
        'fermi_level': fermi_level,
        'moments': moments,
        'energies': dos_data['energies'],
        'total_dos': dos_data['total_dos'],
        'projected_dos': projected_dos,
        'n_kpoints': dos_data.get('n_kpoints', 'unknown'),
        'n_irreducible': dos_data.get('n_irreducible', 'unknown')
    }


def create_summary_plot(results_list):
    """Create summary plots comparing different calculations"""
    if not HAS_MATPLOTLIB or not results_list:
        return
    
    print("\n📊 Creating summary plots...")
    
    # Create figure with subplots
    plt.figure(figsize=(15, 10))
    
    # Plot 1: DOS comparison
    ax1 = plt.subplot(2, 3, 1)
    for result in results_list:
        label = f"{result['mesh_size'][0]}x{result['mesh_size'][1]}x{result['mesh_size'][2]}, T={result['temperature']:.0f}K"
        plt.plot(result['energies'], result['total_dos'], label=label, alpha=0.7)
        plt.axvline(result['fermi_level'], linestyle='--', alpha=0.5)
    plt.xlabel('Energy (eV)')
    plt.ylabel('DOS (states/eV)')
    plt.title('Total DOS Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 2: Fermi level convergence
    ax2 = plt.subplot(2, 3, 2)
    mesh_sizes = [f"{r['mesh_size'][0]}x{r['mesh_size'][1]}x{r['mesh_size'][2]}" for r in results_list]
    fermi_levels = [r['fermi_level'] for r in results_list]
    colors = ['blue' if r['temperature'] == 0 else 'red' for r in results_list]
    plt.scatter(range(len(mesh_sizes)), fermi_levels, c=colors, s=60)
    plt.xlabel('Calculation')
    plt.ylabel('Fermi Level (eV)')
    plt.title('Fermi Level Convergence')
    plt.xticks(range(len(mesh_sizes)), mesh_sizes, rotation=45)
    plt.grid(True, alpha=0.3)
    
    # Plot 3: Orbital occupations
    ax3 = plt.subplot(2, 3, 3)
    orbital_types = ['s', 'p', 'd']
    width = 0.25
    x = np.arange(len(orbital_types))
    
    for i, result in enumerate(results_list):
        occupations = []
        for orb in orbital_types:
            key = f'orbital_{orb}'
            if key in result['moments']:
                occupations.append(result['moments'][key]['occupation'])
            else:
                occupations.append(0.0)
        
        label = f"{result['mesh_size'][0]}x{result['mesh_size'][1]}x{result['mesh_size'][2]}"
        plt.bar(x + i*width, occupations, width, label=label, alpha=0.7)
    
    plt.xlabel('Orbital Type')
    plt.ylabel('Occupation (electrons)')
    plt.title('Orbital Occupations')
    plt.xticks(x + width, orbital_types)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 4: Band centers (first moments)
    ax4 = plt.subplot(2, 3, 4)
    for i, result in enumerate(results_list):
        centers = []
        for orb in orbital_types:
            key = f'orbital_{orb}'
            if key in result['moments'] and result['moments'][key]['occupation'] > 0.01:
                centers.append(result['moments'][key]['band_center'])
            else:
                centers.append(0.0)
        
        label = f"{result['mesh_size'][0]}x{result['mesh_size'][1]}x{result['mesh_size'][2]}"
        plt.bar(x + i*width, centers, width, label=label, alpha=0.7)
    
    plt.xlabel('Orbital Type')
    plt.ylabel('Band Center (eV)')
    plt.title('Orbital Band Centers')
    plt.xticks(x + width, orbital_types)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 5: Timing comparison
    ax5 = plt.subplot(2, 3, 5)
    labels = [f"{r['mesh_size'][0]}³" for r in results_list]
    dos_times = [r['dos_time'] for r in results_list]
    moments_times = [r['moments_time'] for r in results_list]
    
    x = np.arange(len(labels))
    plt.bar(x, dos_times, label='DOS calculation', alpha=0.7)
    plt.bar(x, moments_times, bottom=dos_times, label='Moments calculation', alpha=0.7)
    
    plt.xlabel('Mesh Size')
    plt.ylabel('Time (seconds)')
    plt.title('Computational Time')
    plt.xticks(x, labels)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 6: k-point efficiency
    ax6 = plt.subplot(2, 3, 6)
    n_total = [r['mesh_size'][0]**3 for r in results_list]
    n_irreducible = [r['n_irreducible'] if r['n_irreducible'] != 'unknown' else 0 for r in results_list]
    efficiency = [n_irr/n_tot*100 if n_tot > 0 else 0 for n_tot, n_irr in zip(n_total, n_irreducible)]
    
    plt.bar(range(len(labels)), efficiency, alpha=0.7)
    plt.xlabel('Mesh Size')
    plt.ylabel('Symmetry Reduction (%)')
    plt.title('k-point Symmetry Efficiency')
    plt.xticks(range(len(labels)), labels)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('bccFe_dos_moments_summary.png', dpi=150, bbox_inches='tight')
    plt.savefig('bccFe_dos_moments_summary.pdf', bbox_inches='tight')
    print("📊 Summary plots saved as bccFe_dos_moments_summary.png/pdf")
    
    # Show plot if in interactive mode
    try:
        plt.show()
    except:
        pass

def main():
    """Run comprehensive DOS moments analysis"""
    print("🚀 Starting comprehensive BCC Fe DOS moments analysis...")
    print("This will test multiple mesh sizes and temperatures for convergence")
    
    # Test configurations
    test_configs = [
        {'mesh_size': (4, 4, 4), 'temperature': 300.0, 'label': 'Quick Test (4³, 300K)'},
        {'mesh_size': (4, 4, 4), 'temperature': 0.0, 'label': 'Quick Test (4³, 0K)'},
        {'mesh_size': (6, 6, 6), 'temperature': 300.0, 'label': 'Medium Test (6³, 300K)'},
        {'mesh_size': (6, 6, 6), 'temperature': 0.0, 'label': 'Medium Test (6³, 0K)'},
    ]
    
    results = []
    total_start_time = time.time()
    
    try:
        for config in test_configs:
            try:
                result = test_single_calculation(**config)
                results.append(result)
            except Exception as e:
                print(f"❌ Failed {config['label']}: {e}")
                continue
        
        total_time = time.time() - total_start_time
        
        # Print summary
        print(f"\n{'='*80}")
        print("📋 COMPREHENSIVE ANALYSIS SUMMARY")
        print(f"{'='*80}")
        print(f"⏱️  Total analysis time: {total_time:.1f} seconds")
        print(f"🧪 Successful calculations: {len(results)}/{len(test_configs)}")
        
        if results:
            print("📊 CONVERGENCE ANALYSIS:")
            print(f"{'Mesh':<8} {'T(K)':<6} {'EF(eV)':<10} {'Total e⁻':<10} {'s e⁻':<8} {'p e⁻':<8} {'d e⁻':<8}")
            print("-" * 70)
            
            for result in results:
                mesh_str = f"{result['mesh_size'][0]}³"
                temp_str = f"{result['temperature']:.0f}"
                ef_str = f"{result['fermi_level']:.3f}"
                
                total_e = result['moments']['total']['occupation']
                s_e = result['moments'].get('orbital_s', {}).get('occupation', 0.0)
                p_e = result['moments'].get('orbital_p', {}).get('occupation', 0.0)
                d_e = result['moments'].get('orbital_d', {}).get('occupation', 0.0)
                
                print(f"{mesh_str:<8} {temp_str:<6} {ef_str:<10} {total_e:<10.2f} {s_e:<8.2f} {p_e:<8.2f} {d_e:<8.2f}")
            
            # Create plots
            create_summary_plot(results)
            
            print(f"\n✅ Analysis completed successfully!")
            print(f"📁 Individual results saved with mesh-specific prefixes")
            if HAS_MATPLOTLIB:
                print(f"📊 Summary plots saved as bccFe_dos_moments_summary.png/pdf")
        else:
            print("❌ No successful calculations completed")
            return False
            
    except KeyboardInterrupt:
        print("\n⚠️  Analysis interrupted by user")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return len(results) > 0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
