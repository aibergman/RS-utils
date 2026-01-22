#!/usr/bin/env python3
"""
Quick TB-LMTO Analysis - Simple Interface

This script provides the simplest possible interface for TB-LMTO analysis.
Just run it in a directory with lattice.nml and rs2txt.dat files.
"""

import sys
import subprocess
from pathlib import Path

def check_files():
    """Check if required files exist"""
    required_files = ['lattice.nml', 'rs2txt.dat']
    missing_files = []
    
    for file in required_files:
        if not Path(file).exists():
            missing_files.append(file)
    
    if missing_files:
        print("❌ Missing required files:")
        for file in missing_files:
            print(f"   - {file}")
        print("\nPlease ensure you have:")
        print("   - lattice.nml (crystal structure)")
        print("   - rs2txt.dat (Hamiltonian data from Fortran code)")
        return False
    
    print("✅ All required files found!")
    return True

def run_analysis(analysis_type="high"):
    """Run the appropriate analysis"""
    
    if analysis_type == "quick":
        print("\n🚀 Running Quick analysis (4×4×4 k-mesh)...")
        # Run the package module as a subprocess so the script works when executed
        # directly (no package context required).
        cmd = [
            sys.executable, "-m", "rs_band.complete_analysis",
            "--dos-mesh", "4", "4", "4",
            "--num-energies", "200",
            "--output-prefix", "quick"
        ]
    
    elif analysis_type == "standard":
        print("\n🚀 Running standard analysis (6×6×6 k-mesh)...")
        cmd = [
            sys.executable, "-m", "rs_band.complete_analysis", 
            "--dos-mesh", "6", "6", "6",
            "--num-energies", "300",
            "--output-prefix", "standard"
        ]
    
    elif analysis_type == "high":
        print("\n🚀 Running high-resolution analysis (8×8×8 k-mesh)...")
        cmd = [
            sys.executable, "-m", "rs_band.complete_analysis",
            "--dos-mesh", "8", "8", "8", 
            "--num-energies", "500",
            "--nk-line", "200",
            "--output-prefix", "high_res"
        ]
    
    else:
        print(f"❌ Unknown analysis type: {analysis_type}")
        return False
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ Analysis completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print("❌ Analysis failed!")
        print(f"Error: {e}")
        print(f"Output: {e.stdout}")
        print(f"Error output: {e.stderr}")
        return False

def main():
    """Main function"""
    
    print("=" * 50)
    print("🎵 TB-LMTO Quick Analysis Tool")
    print("=" * 50)
    
    # Check command line arguments
    if len(sys.argv) > 1:
        analysis_type = sys.argv[1].lower()
        if analysis_type not in ["quick", "standard", "high"]:
            print("Usage: python quick_analysis.py [quick|standard|high]")
            print("\nAnalysis types:")
            print("  quick    - Fast analysis (4×4×4 mesh, ~30 seconds)")
            print("  standard - Standard analysis (6×6×6 mesh, ~1 minute)")  
            print("  high     - High-res analysis (8×8×8 mesh, ~3 minutes)")
            print("\nDefault: standard")
            return
    else:
        analysis_type = "standard"
    
    # Check files
    if not check_files():
        return
    
    # Run analysis
    success = run_analysis(analysis_type)
    
    if success:
        print(f"\n📊 Results saved with prefix '{analysis_type}_':")
        print(f"   - {analysis_type}_bands.png (band structure)")
        print(f"   - {analysis_type}_projections.png (orbital projections)")
        print(f"   - {analysis_type}_dos.png (density of states)")
        print(f"   - {analysis_type}_complete.png (combined analysis)")
        print(f"   - {analysis_type}_dos_data.npz (raw data)")
        
        print(f"\n🎉 TB-LMTO analysis complete!")
        print("📖 See ANALYSIS_SUITE_SUMMARY.md for detailed documentation")
    else:
        print("\n💡 Troubleshooting tips:")
        print("   - Check that lattice.nml and rs2txt.dat are in current directory")
        print("   - Verify Python environment has required packages")
        print("   - Try 'python -c \"import numpy, matplotlib, scipy, spglib\"'")

if __name__ == "__main__":
    main()
