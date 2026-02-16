"""
Structure analysis module for TB-LMTO calculations.

This module provides tools to analyze crystal structure, determine space groups,
and generate appropriate k-paths for band structure calculations.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
import warnings

# Try to import spglib
try:
    import spglib
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False
    warnings.warn("spglib not available. Structure analysis will be limited.")

from .readers import LatticeReader


class StructureAnalyzer:
    """Analyze crystal structure and generate appropriate k-paths"""
    
    def __init__(self, lattice_reader: LatticeReader):
        self.lattice = lattice_reader
        self.spacegroup_info = None
        self.conventional_cell = None
        self.special_points = None
        self.band_path = None
        
        if HAS_SPGLIB:
            self._analyze_structure()
    
    def _analyze_structure(self):
        """Analyze crystal structure using spglib"""
        if not HAS_SPGLIB:
            return
            
        # Check that lattice data is available
        if (self.lattice.lattice_vectors is None or 
            self.lattice.atomic_positions is None or 
            self.lattice.atomic_numbers is None):
            print("⚠️  Incomplete lattice data. Cannot analyze structure.")
            return
            
        # Convert to spglib format
        lattice_vectors = self.lattice.lattice_vectors.T  # spglib expects row vectors
        positions = self.lattice.atomic_positions
        numbers = self.lattice.atomic_numbers
        
        cell = (lattice_vectors, positions, numbers)
        
        # Get space group information
        spacegroup_string = spglib.get_spacegroup(cell, symprec=1e-5)
        if spacegroup_string:
            print(f"🔍 Raw spacegroup string: {spacegroup_string}")
            
            # Parse space group string like "Im-3m (229)"
            try:
                if '(' in spacegroup_string and ')' in spacegroup_string:
                    # Extract space group number from parentheses
                    sg_number = int(spacegroup_string.split('(')[-1].split(')')[0])
                    # Extract symbol (everything before the number in parentheses)
                    sg_symbol = spacegroup_string.split('(')[0].strip()
                    
                    # Create our own space group info dict since spglib.get_spacegroup_type uses hall numbers
                    self.spacegroup_info = {
                        'number': sg_number,
                        'international_short': sg_symbol,
                        'international_full': sg_symbol,  # We don't have the full symbol from the string
                        'crystal_system': self._determine_crystal_system(sg_number)
                    }
                    
                    print(f"✅ Space group identified: {sg_symbol} ({sg_number})")
                    
                else:
                    raise ValueError(f"Cannot parse spacegroup string: {spacegroup_string}")
                    
            except (ValueError, IndexError) as e:
                print(f"⚠️  Could not parse space group: {spacegroup_string}, error: {e}")
                self._fallback_spacegroup_info()
        else:
            print("⚠️  Could not determine space group")
            self._fallback_spacegroup_info()
        
        # Get conventional cell
        self.conventional_cell = spglib.standardize_cell(cell)
        
        # Print structure information
        self._print_structure_info()
        
        # Generate k-path
        self._generate_kpath()
    
    def _fallback_spacegroup_info(self):
        """Create fallback space group info when detection fails"""
        self.spacegroup_info = {
            'number': 1,
            'international_short': 'P1',
            'international_full': 'P 1',
            'crystal_system': 'triclinic'
        }
    
    def _print_structure_info(self):
        """Print detailed structure information"""
        if not self.spacegroup_info:
            return
            
        print(f"\n{'='*60}")
        print("🔬 Crystal Structure Analysis")
        print(f"{'='*60}")
        
        # Space group information
        try:
            sg_number = self.spacegroup_info.number if hasattr(self.spacegroup_info, 'number') else self.spacegroup_info['number']
            sg_symbol = self.spacegroup_info.international_short if hasattr(self.spacegroup_info, 'international_short') else self.spacegroup_info['international_short']
            sg_full = self.spacegroup_info.international_full if hasattr(self.spacegroup_info, 'international_full') else self.spacegroup_info['international_full']
            
            # Crystal system is not always available in newer spglib
            if hasattr(self.spacegroup_info, 'crystal_system'):
                crystal_system = self.spacegroup_info.crystal_system
            elif 'crystal_system' in self.spacegroup_info:
                crystal_system = self.spacegroup_info['crystal_system']
            else:
                # Determine crystal system from space group number
                crystal_system = self._determine_crystal_system(sg_number)
            
            print(f"Space Group: {sg_number} ({sg_symbol})")
            print(f"Full Symbol: {sg_full}")
            print(f"Crystal System: {crystal_system}")
            
        except Exception as e:
            print(f"⚠️  Error accessing space group info: {e}")
            return
        
        # Lattice parameters
        if self.conventional_cell:
            conv_lattice = self.conventional_cell[0]
            a, b, c = np.linalg.norm(conv_lattice, axis=1)
            
            # Calculate angles
            alpha = np.arccos(np.dot(conv_lattice[1], conv_lattice[2]) / (b * c)) * 180 / np.pi
            beta = np.arccos(np.dot(conv_lattice[0], conv_lattice[2]) / (a * c)) * 180 / np.pi
            gamma = np.arccos(np.dot(conv_lattice[0], conv_lattice[1]) / (a * b)) * 180 / np.pi
            
            print(f"\nLattice Parameters:")
            print(f"  a = {a:.6f} Å")
            print(f"  b = {b:.6f} Å") 
            print(f"  c = {c:.6f} Å")
            print(f"  α = {alpha:.2f}°")
            print(f"  β = {beta:.2f}°")
            print(f"  γ = {gamma:.2f}°")
        
        # Volume
        volume = np.abs(np.linalg.det(self.lattice.lattice_vectors))
        print(f"\nUnit Cell Volume: {volume:.6f} Å³")
        
        print(f"{'='*60}")
    
    def _determine_crystal_system(self, sg_number: int) -> str:
        """Determine crystal system from space group number"""
        if 1 <= sg_number <= 2:
            return 'triclinic'
        elif 3 <= sg_number <= 15:
            return 'monoclinic'
        elif 16 <= sg_number <= 74:
            return 'orthorhombic'
        elif 75 <= sg_number <= 142:
            return 'tetragonal'
        elif 143 <= sg_number <= 167:
            return 'trigonal'
        elif 168 <= sg_number <= 194:
            return 'hexagonal'
        elif 195 <= sg_number <= 230:
            return 'cubic'
        else:
            return 'unknown'
    
    def _generate_kpath(self):
        """Generate appropriate k-path based on crystal system"""
        # Prefer specialized k-path generators when available (seekpath -> pymatgen)
        if not HAS_SPGLIB or not self.spacegroup_info:
            self._fallback_kpath()
            return

        # Prepare cell in fractional coords expected by external libs
        lattice = self.lattice.lattice_vectors.T
        positions = (self.lattice.atomic_positions.tolist()
                     if isinstance(self.lattice.atomic_positions, np.ndarray)
                     else list(self.lattice.atomic_positions))
        numbers = (self.lattice.atomic_numbers.tolist()
                   if isinstance(self.lattice.atomic_numbers, np.ndarray)
                   else list(self.lattice.atomic_numbers))

        # Try seekpath first (produces standardized k-paths matching conventions)
        try:
            import seekpath
            print('\n🗺️  Generating k-path via seekpath')
            sp = seekpath.get_path((lattice, positions, numbers))
            # sp contains 'point_coords' (label->frac coords) and 'path' (list of tuples)
            self.special_points = {lab: list(coords) for lab, coords in sp.get('point_coords', {}).items()}
            # convert path (list of tuples) into list-of-lists for our code
            self.band_path = [list(seg) for seg in sp.get('path', [])]
            if self.band_path and self.special_points:
                print(f"✅ seekpath produced {len(self.band_path)} segments and {len(self.special_points)} special points")
                return
        except Exception:
            # seekpath not available or failed; fall through to pymatgen
            pass

        # Try pymatgen's HighSymmKpath if available
        try:
            from pymatgen.core import Lattice, Structure
            from pymatgen.symmetry.bandstructure import HighSymmKpath
            print('\n🗺️  Generating k-path via pymatgen HighSymmKpath')
            pmg_lat = Lattice(lattice)
            # build pymatgen Structure from fractional positions and atomic numbers
            species = []
            for z in numbers:
                # fall back to atomic number -> element symbol minimal mapping
                try:
                    from pymatgen.core.periodic_table import Element
                    species.append(Element(int(z)).symbol)
                except Exception:
                    species.append('X')
            structure = Structure(pmg_lat, species, positions)
            kpath = HighSymmKpath(structure)
            pts = kpath.kpath['kpoints']
            path = kpath.kpath['path']
            self.special_points = {lab: list(coord) for lab, coord in pts.items()}
            self.band_path = [list(seg) for seg in path]
            if self.band_path and self.special_points:
                print(f"✅ pymatgen produced {len(self.band_path)} segments and {len(self.special_points)} special points")
                return
        except Exception:
            pass

        # Fallback to older/manual generation if neither library is available
        print('\n⚠️  Falling back to manual k-path generation (install seekpath or pymatgen for standardized paths)')
        self._generate_manual_kpath()
    
    def _generate_spglib_kpath(self):
        """Generate k-path using spglib's built-in paths"""
        try:
            # For spglib >= 2.0, use seekpath-like functionality
            from spglib import get_ir_reciprocal_mesh
            
            # Get primitive cell
            lattice_vectors = self.lattice.lattice_vectors.T
            positions = self.lattice.atomic_positions
            numbers = self.lattice.atomic_numbers
            cell = (lattice_vectors, positions, numbers)
            
            primitive_cell = spglib.find_primitive(cell)
            if primitive_cell:
                self._generate_kpath_from_primitive(primitive_cell)
            else:
                self._generate_manual_kpath()
                
        except Exception as e:
            print(f"⚠️  Failed to use spglib k-path generation: {e}")
            self._generate_manual_kpath()
    
    def _generate_kpath_from_primitive(self, primitive_cell):
        """Generate k-path from primitive cell"""
        crystal_system = self.spacegroup_info['crystal_system']
        sg_number = self.spacegroup_info['number']
        
        # Get reciprocal lattice
        lattice = primitive_cell[0]
        recip_lattice = 2 * np.pi * np.linalg.inv(lattice).T
        
        if crystal_system == 'cubic':
            if sg_number == 229 or 'I' in str(self.spacegroup_info['international_short']):  # Body-centered cubic
                self._bcc_kpath(recip_lattice)
            elif 195 <= sg_number <= 230:  # Face-centered cubic  
                self._fcc_kpath(recip_lattice)
            else:  # Simple cubic
                self._sc_kpath(recip_lattice)
        elif crystal_system == 'hexagonal':
            self._hexagonal_kpath(recip_lattice)
        elif crystal_system == 'tetragonal':
            self._tetragonal_kpath(recip_lattice)
        else:
            self._generic_kpath(recip_lattice)
    
    def _bcc_kpath(self, recip_lattice):
        """Generate k-path for body-centered cubic lattice"""
        print("📍 Using BCC k-path: Γ-H-N-Γ-P-H | P-N")
        
        # BCC special points in reciprocal lattice coordinates
        self.special_points = {
            'Γ': [0.0, 0.0, 0.0],      # Gamma
            'H': [0.5, -0.5, 0.5],     # H  
            'P': [0.25, 0.25, 0.25],   # P
            'N': [0.0, 0.0, 0.5],      # N
        }
        
        # Define the path
        self.band_path = [
            ['Γ', 'H', 'N', 'Γ', 'P', 'H'],
            ['P', 'N']
        ]
        
        print("   Special points:")
        for label, coords in self.special_points.items():
            print(f"     {label}: {coords}")
    
    def _fcc_kpath(self, recip_lattice):
        """Generate k-path for face-centered cubic lattice"""
        print("📍 Using FCC k-path: Γ-X-W-K-Γ-L-U-W-L-K | U-X")
        
        self.special_points = {
            'Γ': [0.0, 0.0, 0.0],      # Gamma
            'X': [0.5, 0.0, 0.5],      # X
            'L': [0.5, 0.5, 0.5],      # L
            'W': [0.5, 0.25, 0.75],    # W
            'K': [0.375, 0.375, 0.75], # K
            'U': [0.625, 0.25, 0.625]  # U
        }
        
        self.band_path = [
            ['Γ', 'X', 'W', 'K', 'Γ', 'L', 'U', 'W', 'L', 'K'],
            ['U', 'X']
        ]
        
        print("   Special points:")
        for label, coords in self.special_points.items():
            print(f"     {label}: {coords}")
    
    def _sc_kpath(self, recip_lattice):
        """Generate k-path for simple cubic lattice"""
        print("📍 Using SC k-path: Γ-X-M-Γ-R-X | M-R")
        
        self.special_points = {
            'Γ': [0.0, 0.0, 0.0],      # Gamma
            'X': [0.5, 0.0, 0.0],      # X
            'M': [0.5, 0.5, 0.0],      # M
            'R': [0.5, 0.5, 0.5]       # R
        }
        
        self.band_path = [
            ['Γ', 'X', 'M', 'Γ', 'R', 'X'],
            ['M', 'R']
        ]
        
        print("   Special points:")
        for label, coords in self.special_points.items():
            print(f"     {label}: {coords}")
    
    def _hexagonal_kpath(self, recip_lattice):
        """Generate k-path for hexagonal lattice"""
        print("📍 Using Hexagonal k-path: Γ-M-K-Γ-A-L-H-A | L-M | K-H")
        
        self.special_points = {
            'Γ': [0.0, 0.0, 0.0],          # Gamma
            'M': [0.5, 0.0, 0.0],          # M
            'K': [1/3, 1/3, 0.0],          # K
            'A': [0.0, 0.0, 0.5],          # A
            'L': [0.5, 0.0, 0.5],          # L
            'H': [1/3, 1/3, 0.5]           # H
        }
        
        self.band_path = [
            ['Γ', 'M', 'K', 'Γ', 'A', 'L', 'H', 'A'],
            ['L', 'M'],
            ['K', 'H']
        ]
        
        print("   Special points:")
        for label, coords in self.special_points.items():
            print(f"     {label}: {coords}")
    
    def _tetragonal_kpath(self, recip_lattice):
        """Generate k-path for tetragonal lattice"""
        print("📍 Using Tetragonal k-path: Γ-X-M-Γ-Z-R-A-Z | X-R | M-A")
        
        self.special_points = {
            'Γ': [0.0, 0.0, 0.0],      # Gamma
            'X': [0.5, 0.0, 0.0],      # X
            'M': [0.5, 0.5, 0.0],      # M
            'Z': [0.0, 0.0, 0.5],      # Z
            'R': [0.5, 0.0, 0.5],      # R
            'A': [0.5, 0.5, 0.5]       # A
        }
        
        self.band_path = [
            ['Γ', 'X', 'M', 'Γ', 'Z', 'R', 'A', 'Z'],
            ['X', 'R'],
            ['M', 'A']
        ]
        
        print("   Special points:")
        for label, coords in self.special_points.items():
            print(f"     {label}: {coords}")
    
    def _generic_kpath(self, recip_lattice):
        """Generate generic k-path for other crystal systems"""
        print("📍 Using generic k-path: Γ-X-Y-Γ-Z")
        
        self.special_points = {
            'Γ': [0.0, 0.0, 0.0],      # Gamma
            'X': [0.5, 0.0, 0.0],      # X
            'Y': [0.0, 0.5, 0.0],      # Y
            'Z': [0.0, 0.0, 0.5]       # Z
        }
        
        self.band_path = [['Γ', 'X', 'Y', 'Γ', 'Z']]
        
        print("   Special points:")
        for label, coords in self.special_points.items():
            print(f"     {label}: {coords}")
    
    def _generate_manual_kpath(self):
        """Generate k-path manually based on space group"""
        if not self.spacegroup_info:
            self._fallback_kpath()
            return
            
        try:
            crystal_system = (self.spacegroup_info.crystal_system 
                            if hasattr(self.spacegroup_info, 'crystal_system') 
                            else self._determine_crystal_system(
                                self.spacegroup_info.number 
                                if hasattr(self.spacegroup_info, 'number') 
                                else self.spacegroup_info['number']))
            
            sg_number = (self.spacegroup_info.number 
                        if hasattr(self.spacegroup_info, 'number') 
                        else self.spacegroup_info['number'])
        except Exception as e:
            print(f"⚠️  Error accessing space group data: {e}")
            self._fallback_kpath()
            return
        
        # Identify lattice type from space group
        if crystal_system == 'cubic':
            if 195 <= sg_number <= 230:  # Could be FCC or BCC
                # Check for BCC specifically
                if sg_number == 229 or 'I' in str(self.spacegroup_info['international_short']):  # BCC
                    self._bcc_kpath(None)
                else:  # FCC
                    self._fcc_kpath(None)
            else:  # Simple cubic
                self._sc_kpath(None)
        elif crystal_system == 'hexagonal':
            self._hexagonal_kpath(None)
        elif crystal_system == 'tetragonal':
            self._tetragonal_kpath(None)
        else:
            self._generic_kpath(None)
    
    def _fallback_kpath(self):
        """Fallback k-path when spglib is not available"""
        print("⚠️  Using fallback generic k-path (install spglib for proper analysis)")
        self._generic_kpath(None)
    
    def generate_kpath_points(self, num_points: int = 100) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Generate k-points along the band path"""
        if not self.special_points or not self.band_path:
            raise ValueError("No k-path defined. Run structure analysis first.")
            
        if self.lattice.lattice_vectors is None:
            raise ValueError("No lattice vectors available.")
        
        all_kpoints = []
        all_labels = []
        label_positions = []
        first_point_added = False

        total_segments = sum(len(seg) - 1 for seg in self.band_path)
        points_per_segment = max(10, num_points // total_segments)

        for segment in self.band_path:
            for i in range(len(segment) - 1):
                start_label = segment[i]
                end_label = segment[i + 1]

                start_point = np.array(self.special_points[start_label])
                end_point = np.array(self.special_points[end_label])

                # Number of points for this segment (shared across segments)
                segment_points = points_per_segment

                # Generate points along segment
                for j in range(segment_points):
                    # Compute interpolation parameter
                    t = j / (segment_points - 1) if segment_points > 1 else 0
                    point = start_point + t * (end_point - start_point)

                    # If this is the very first point overall, append it and label
                    if not first_point_added and j == 0:
                        all_labels.append(start_label)
                        label_positions.append(len(all_kpoints))
                        all_kpoints.append(point)
                        first_point_added = True
                        continue

                    # For subsequent segments, skip the duplicate start point (j==0)
                    if first_point_added and j == 0:
                        # Do not append duplicate point
                        continue

                    # Normal interior or endpoint points
                    all_labels.append('')
                    all_kpoints.append(point)

                # Assign end label to the last k-point of this segment
                if len(all_kpoints) > 0:
                    all_labels[-1] = end_label
                    label_positions.append(len(all_kpoints) - 1)
        
        kpath = np.array(all_kpoints)
        
        # Keep k-path in fractional coordinates for the band calculator
        kpath_frac = kpath

        # Convert to Cartesian coordinates for distance calculations only
        recip_lattice = 2 * np.pi * np.linalg.inv(self.lattice.lattice_vectors)
        kpath_cart = np.dot(kpath_frac, recip_lattice)

        # Calculate cumulative distances along k-path (in Cartesian reciprocal units)
        kdistances = np.zeros(len(kpath_cart))
        for i in range(1, len(kpath_cart)):
            dk = np.linalg.norm(kpath_cart[i] - kpath_cart[i-1])
            kdistances[i] = kdistances[i-1] + dk

        # Ensure labels array matches kdistances length
        all_labels = all_labels[:len(kdistances)]  # Trim to match kdistances length
        while len(all_labels) < len(kdistances):  # Pad with empty strings if needed
            all_labels.append('')

        # Return fractional k-points (so the band calculator can convert them once),
        # the labels, and the Cartesian k-distances for plotting.
        return kpath_frac, all_labels, kdistances
    
    def plot_structure_info(self, save_path: str = "structure_analysis.png"):
        """Create a plot showing structure information"""
        if not HAS_SPGLIB or not self.spacegroup_info:
            print("⚠️  Cannot plot structure info without spglib")
            return
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Real space unit cell
        self._plot_unit_cell(ax1)
        
        # Plot 2: Reciprocal space and k-path
        self._plot_brillouin_zone(ax2)
        
        # Plot 3: Structure information table
        self._plot_structure_table(ax3)
        
        # Plot 4: k-path segments
        self._plot_kpath_segments(ax4)
        
        plt.suptitle(f"Structure Analysis: {self.spacegroup_info['international_short']} "
                    f"({self.spacegroup_info['crystal_system']})", fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📊 Structure analysis plot saved: {save_path}")
    
    def _plot_unit_cell(self, ax):
        """Plot the unit cell in real space"""
        # Plot lattice vectors
        origin = np.zeros(3)
        vectors = self.lattice.lattice_vectors
        
        # 2D projection (xy plane)
        ax.arrow(0, 0, vectors[0, 0], vectors[0, 1], 
                head_width=0.05, head_length=0.05, fc='red', ec='red', label='a₁')
        ax.arrow(0, 0, vectors[1, 0], vectors[1, 1], 
                head_width=0.05, head_length=0.05, fc='green', ec='green', label='a₂')
        ax.arrow(0, 0, vectors[2, 0], vectors[2, 1], 
                head_width=0.05, head_length=0.05, fc='blue', ec='blue', label='a₃')
        
        # Plot atomic positions
        for i in range(self.lattice.nbas):
            pos = self.lattice.atomic_positions[i]
            real_pos = np.dot(pos, self.lattice.lattice_vectors)
            ax.scatter(real_pos[0], real_pos[1], s=100, c='black', marker='o')
            ax.text(real_pos[0], real_pos[1] + 0.1, f'Atom {i+1}', 
                   ha='center', va='bottom', fontsize=8)
        
        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')
        ax.set_title('Unit Cell (xy projection)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
    
    def _plot_brillouin_zone(self, ax):
        """Plot the Brillouin zone and k-path"""
        if not self.special_points:
            ax.text(0.5, 0.5, 'No k-path data', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('Brillouin Zone')
            return
        
        # Plot special points
        for label, coords in self.special_points.items():
            ax.scatter(coords[0], coords[1], s=100, c='red', marker='o')
            ax.text(coords[0], coords[1] + 0.05, label, 
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # Plot k-path
        colors = ['blue', 'green', 'orange', 'purple']
        for i, segment in enumerate(self.band_path):
            color = colors[i % len(colors)]
            for j in range(len(segment) - 1):
                start = self.special_points[segment[j]]
                end = self.special_points[segment[j + 1]]
                ax.plot([start[0], end[0]], [start[1], end[1]], 
                       color=color, linewidth=2, alpha=0.7)
        
        ax.set_xlabel('kₓ (2π/a)')
        ax.set_ylabel('k_y (2π/a)')
        ax.set_title('Brillouin Zone & k-path')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
    
    def _plot_structure_table(self, ax):
        """Plot structure information as a table"""
        ax.axis('off')
        
        if not self.spacegroup_info:
            ax.text(0.5, 0.5, 'No structure data', ha='center', va='center')
            return
        
        # Create table data
        data = [
            ['Space Group', f"{self.spacegroup_info['number']} ({self.spacegroup_info['international_short']})"],
            ['Crystal System', self.spacegroup_info['crystal_system']],
            ['Full Symbol', self.spacegroup_info['international_full']],
            ['Number of atoms', str(self.lattice.nbas)],
        ]
        
        if self.conventional_cell:
            conv_lattice = self.conventional_cell[0]
            a, b, c = np.linalg.norm(conv_lattice, axis=1)
            data.extend([
                ['a (Å)', f"{a:.4f}"],
                ['b (Å)', f"{b:.4f}"],
                ['c (Å)', f"{c:.4f}"],
            ])
        
        volume = np.abs(np.linalg.det(self.lattice.plat))
        data.append(['Volume (Å³)', f"{volume:.4f}"])
        
        # Create table
        table = ax.table(cellText=data, cellLoc='left', loc='center',
                        colWidths=[0.4, 0.6])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        # Style the table
        for i in range(len(data)):
            table[(i, 0)].set_facecolor('#E8E8E8')
            table[(i, 0)].set_text_props(weight='bold')
        
        ax.set_title('Structure Information')
    
    def _plot_kpath_segments(self, ax):
        """Plot k-path segments information"""
        ax.axis('off')
        
        if not self.band_path:
            ax.text(0.5, 0.5, 'No k-path data', ha='center', va='center')
            return
        
        # Create k-path string
        path_strings = []
        for i, segment in enumerate(self.band_path):
            path_str = '-'.join(segment)
            path_strings.append(f"Segment {i+1}: {path_str}")
        
        full_path = ' | '.join(['-'.join(seg) for seg in self.band_path])
        
        text = f"K-path: {full_path}\n\n"
        text += "\n".join(path_strings)
        
        if self.special_points:
            text += "\n\nSpecial Points:\n"
            for label, coords in self.special_points.items():
                text += f"{label}: ({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f})\n"
        
        ax.text(0.05, 0.95, text, ha='left', va='top', fontsize=9, 
               family='monospace', transform=ax.transAxes)
        ax.set_title('K-path Details')


def analyze_structure(lattice_file: str = "lattice.nml", 
                     plot: bool = True, 
                     save_path: str = "structure_analysis.png") -> StructureAnalyzer:
    """
    Analyze crystal structure and generate appropriate k-path.
    
    Parameters:
    -----------
    lattice_file : str
        Path to lattice.nml file
    plot : bool
        Whether to generate structure analysis plots
    save_path : str
        Path to save the analysis plot
        
    Returns:
    --------
    StructureAnalyzer
        The analyzer object with structure information
    """
    print("🔬 Starting crystal structure analysis...")
    
    # Read lattice
    lattice_reader = LatticeReader(lattice_file)
    
    # Analyze structure
    analyzer = StructureAnalyzer(lattice_reader)
    
    if plot and HAS_SPGLIB:
        analyzer.plot_structure_info(save_path)
    
    return analyzer


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze crystal structure")
    parser.add_argument('--lattice', default='lattice.nml', 
                       help='Lattice file (default: lattice.nml)')
    parser.add_argument('--output', default='structure_analysis.png',
                       help='Output plot file (default: structure_analysis.png)')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip plotting')
    
    args = parser.parse_args()
    
    analyzer = analyze_structure(args.lattice, not args.no_plot, args.output)
