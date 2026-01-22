#!/usr/bin/env python3
"""
Full Brillouin Zone DOS Calculator with Tetrahedron Method

This module provides:
- Full BZ k-point sampling with symmetry reduction
- Tetrahedron method for accurate DOS calculation
- Total, site-projected, and spin-projected DOS
- Integration with existing band structure tools
"""

import numpy as np
import matplotlib.pyplot as plt
# ConvexHull not used here; left out to avoid unused import
from .projections import ProjectedBandCalculator
from .readers import LatticeReader, HamiltonianReader
from .dos_moments import DOSMomentsCalculator, analyze_dos_moments
from typing import Dict, Tuple, Optional
import time
import argparse

try:
    import spglib
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False


class TetrahedronDOS:
    """Tetrahedron method for DOS calculation"""
    
    def __init__(self, kpoints: np.ndarray, tetrahedra: np.ndarray, weights: np.ndarray):
        """
        Initialize tetrahedron DOS calculator
        
        Args:
            kpoints: Array of k-points (nk, 3)
            tetrahedra: Array of tetrahedron vertex indices (ntet, 4)
            weights: Weight of each k-point (nk,)
        """
        self.kpoints = kpoints
        self.tetrahedra = tetrahedra
        self.weights = weights
        self.num_kpoints = len(kpoints)
        self.num_tetrahedra = len(tetrahedra)
    
    def calculate_dos(self, energies: np.ndarray, energy_grid: np.ndarray,
                      projections: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Calculate DOS using tetrahedron method
        
        Args:
            energies: Band energies (nk, nbands)
            energy_grid: Energy points for DOS calculation
            projections: Optional projection weights (nk, nbands, nproj)
            
        Returns:
            dos: DOS on energy grid (if projections: (ne, nproj))
        """
        
        num_bands = energies.shape[1]
        num_energies = len(energy_grid)

        if projections is not None:
            num_proj = projections.shape[2]
            dos = np.zeros((num_energies, num_proj))
        else:
            dos = np.zeros(num_energies)

        # Process each tetrahedron
        for tet_idx in range(self.num_tetrahedra):
            tet_vertices = self.tetrahedra[tet_idx]

            # Determine tetrahedron weight from k-point weights
            tet_weight = np.mean(self.weights[tet_vertices])

            # Get energies at tetrahedron vertices for each band
            for band in range(num_bands):
                tet_energies = energies[tet_vertices, band]

                # Sort energies
                sort_idx = np.argsort(tet_energies)
                sorted_energies = tet_energies[sort_idx]

                # Get projections if provided
                if projections is not None:
                    tet_proj = projections[tet_vertices, band, :]
                    sorted_proj = tet_proj[sort_idx, :]

                # Vectorized contribution over the full energy grid
                contrib_array = self._tetrahedron_contribution_array(energy_grid, sorted_energies)
                if projections is not None:
                    avg_proj = np.mean(sorted_proj, axis=0)
                    # contrib_array shape (ne,), avg_proj (nproj,) -> outer product
                    dos += np.outer(contrib_array * tet_weight, avg_proj)
                else:
                    dos += contrib_array * tet_weight

        # No further normalization here; weights already incorporate k-point sampling
        return dos
    
    def _tetrahedron_contribution(self, energy: float, sorted_energies: np.ndarray) -> float:
        """Calculate tetrahedron contribution to DOS at given energy"""
        
        e1, e2, e3, e4 = sorted_energies
        
        # Check for degenerate cases (avoid division by zero)
        tolerance = 1e-12
        
        # If all energies are the same (completely degenerate)
        if abs(e4 - e1) < tolerance:
            # Use delta function approximation: contribution only if energy matches
            if abs(energy - e1) < tolerance:
                return 1.0  # This will be normalized by tetrahedron volume
            else:
                return 0.0
        
        # Check for partially degenerate cases
        if abs(e2 - e1) < tolerance:
            # e1 = e2, use simplified formula
            if energy < e1 or energy > e4:
                return 0.0
            elif energy < e3:
                if abs(e3 - e1) < tolerance or abs(e4 - e1) < tolerance:
                    return 0.0  # Too degenerate, skip
                return 3 * (energy - e1)**2 / ((e3 - e1) * (e4 - e1))
            else:
                if abs(e4 - e3) < tolerance:
                    return 0.0  # Too degenerate, skip
                return 3 * (e4 - energy)**2 / ((e4 - e1) * (e4 - e3))
        
        if abs(e3 - e2) < tolerance:
            # e2 = e3, use simplified formula
            if energy < e1 or energy > e4:
                return 0.0
            elif energy < e2:
                if abs(e2 - e1) < tolerance or abs(e4 - e1) < tolerance:
                    return 0.0  # Too degenerate, skip
                return 3 * (energy - e1)**2 / ((e2 - e1) * (e4 - e1))
            else:
                if abs(e4 - e2) < tolerance:
                    return 0.0  # Too degenerate, skip
                return 3 * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2))
        
        if abs(e4 - e3) < tolerance:
            # e3 = e4, use simplified formula
            if energy < e1 or energy > e3:
                return 0.0
            elif energy < e2:
                if abs(e2 - e1) < tolerance or abs(e3 - e1) < tolerance:
                    return 0.0  # Too degenerate, skip
                return 3 * (energy - e1)**2 / ((e2 - e1) * (e3 - e1))
            else:
                # energy in [e2, e3]
                if abs(e2 - e1) < tolerance or abs(e3 - e1) < tolerance:
                    return 0.0  # Too degenerate, skip
                A = (e2 - e1) * (e3 - e1)
                B = (energy - e1)**2 * (2 * e3 - energy - e1)
                C = (energy - e2)**2 * (energy - e1)
                return 3 * (B + C) / A
        
        # Standard non-degenerate case
        # Check which case we're in
        if energy < e1 or energy > e4:
            return 0.0
        elif energy < e2:
            # Case 1: energy in [e1, e2]
            denom = (e2 - e1) * (e3 - e1) * (e4 - e1)
            if abs(denom) < tolerance:
                return 0.0  # Avoid division by zero
            return 3 * (energy - e1)**2 / denom
        elif energy < e3:
            # Case 2: energy in [e2, e3]
            A = (e2 - e1) * (e3 - e1) * (e4 - e1)
            if abs(A) < tolerance:
                return 0.0  # Avoid division by zero
            B = (energy - e1)**2 * (e4 + e3 - energy - e1)
            C = (energy - e2)**2 * (energy - e1)
            D = (energy - e2)**2 * (e4 - energy)
            return 3 * (B + C + D) / A
        else:
            # Case 3: energy in [e3, e4]
            denom = (e4 - e1) * (e4 - e2) * (e4 - e3)
            if abs(denom) < tolerance:
                return 0.0  # Avoid division by zero
            return 3 * (e4 - energy)**2 / denom

    def _tetrahedron_contribution_array(self, energy_grid: np.ndarray, sorted_energies: np.ndarray) -> np.ndarray:
        """Vectorized tetrahedron contribution over an energy grid."""
        e1, e2, e3, e4 = sorted_energies
        tol = 1e-12
        contrib = np.zeros_like(energy_grid)

        # Degenerate full case
        if abs(e4 - e1) < tol:
            mask = np.abs(energy_grid - e1) < tol
            contrib[mask] = 1.0
            return contrib

        # Prepare masks
        mask1 = (energy_grid >= e1) & (energy_grid < e2)
        mask2 = (energy_grid >= e2) & (energy_grid < e3)
        mask3 = (energy_grid >= e3) & (energy_grid <= e4)

        # Case 1: [e1, e2]
        if np.any(mask1):
            denom = (e2 - e1) * (e3 - e1) * (e4 - e1)
            if abs(denom) > tol:
                e = energy_grid[mask1]
                contrib[mask1] = 3.0 * (e - e1)**2 / denom

        # Case 2: [e2, e3]
        if np.any(mask2):
            A = (e2 - e1) * (e3 - e1) * (e4 - e1)
            if abs(A) > tol:
                e = energy_grid[mask2]
                B = (e - e1)**2 * (e4 + e3 - e - e1)
                C = (e - e2)**2 * (e - e1)
                D = (e - e2)**2 * (e4 - e)
                contrib[mask2] = 3.0 * (B + C + D) / A

        # Case 3: [e3, e4]
        if np.any(mask3):
            denom = (e4 - e1) * (e4 - e2) * (e4 - e3)
            if abs(denom) > tol:
                e = energy_grid[mask3]
                contrib[mask3] = 3.0 * (e4 - e)**2 / denom

        # Zero outside [e1,e4] already covered by initialization
        return contrib


class FullBZCalculator:
    """Full Brillouin zone calculator with symmetry reduction"""
    
    def __init__(self, lattice_reader: LatticeReader, hamiltonian_reader: HamiltonianReader):
        self.lattice = lattice_reader
        self.hamiltonian = hamiltonian_reader
        self.calculator = ProjectedBandCalculator(lattice_reader, hamiltonian_reader)
        
    def generate_kpoint_mesh(self, mesh_size: Tuple[int, int, int], use_symmetry: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate k-point mesh with optional symmetry reduction
        
        Args:
            mesh_size: (nx, ny, nz) mesh dimensions
            use_symmetry: Whether to use symmetry reduction
            
        Returns:
            kpoints: Irreducible k-points
            weights: Weights of each k-point
            tetrahedra: Tetrahedron connectivity
        """
        
        if use_symmetry and HAS_SPGLIB:
            return self._generate_symmetric_mesh(mesh_size)
        else:
            return self._generate_uniform_mesh(mesh_size)
    
    def _generate_symmetric_mesh(self, mesh_size: Tuple[int, int, int]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate symmetry-reduced k-point mesh using spglib"""
        
        # Get crystal structure for spglib
        cell = (
            self.lattice.lattice_vectors.T,
            self.lattice.atomic_positions,
            self.lattice.atomic_numbers
        )
        
        # Generate irreducible k-points
        mapping, grid = spglib.get_ir_reciprocal_mesh(mesh_size, cell, is_shift=[0, 0, 0])
        
        # Get unique k-points and their weights
        unique_indices = np.unique(mapping)
        weights = np.bincount(mapping)[unique_indices]
        
        # Convert grid indices to fractional coordinates
        kpoints = []
        nx, ny, nz = mesh_size
        
        for idx in unique_indices:
            # Find first occurrence of this irreducible k-point index
            grid_indices = np.where(mapping == idx)[0]
            if len(grid_indices) > 0:
                grid_point = grid_indices[0]
                
                # Convert linear index to 3D grid coordinates
                iz = grid_point % nz
                iy = (grid_point // nz) % ny
                ix = grid_point // (ny * nz)
                
                # Convert to fractional coordinates
                kx = ix / nx
                ky = iy / ny
                kz = iz / nz
                
                kpoints.append([kx, ky, kz])
        
        kpoints = np.array(kpoints)
        weights = weights.astype(float) / np.sum(weights)
        
        # Generate tetrahedra (simplified - using Delaunay triangulation)
        tetrahedra = self._generate_tetrahedra(kpoints, mesh_size)
        
        print(f"Generated {len(kpoints)} irreducible k-points from {np.prod(mesh_size)} total points")
        print(f"Symmetry reduction factor: {np.prod(mesh_size) / len(kpoints):.1f}")
        
        return kpoints, weights, tetrahedra
    
    def _generate_uniform_mesh(self, mesh_size: Tuple[int, int, int]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate uniform k-point mesh without symmetry reduction"""
        nx, ny, nz = mesh_size

        # Generate uniform grid (fractional coordinates in [0,1))
        kx = np.linspace(0.0, 1.0 - 1.0 / nx, nx)
        ky = np.linspace(0.0, 1.0 - 1.0 / ny, ny)
        kz = np.linspace(0.0, 1.0 - 1.0 / nz, nz)

        # Build the list of k-points
        kpoints_list = []
        for x in kx:
            for y in ky:
                for z in kz:
                    kpoints_list.append([x, y, z])

        kpoints = np.array(kpoints_list)
        weights = np.ones(len(kpoints)) / float(len(kpoints))

        # Generate tetrahedra for the regular grid
        tetrahedra = self._generate_tetrahedra_uniform(mesh_size)

        print(f"Generated {len(kpoints)} k-points in uniform mesh")

        return kpoints, weights, tetrahedra
    
    def _generate_tetrahedra(self, kpoints: np.ndarray, mesh_size: Tuple[int, int, int]) -> np.ndarray:
        """Generate tetrahedra using Delaunay triangulation"""
        
        try:
            from scipy.spatial import Delaunay
            tri = Delaunay(kpoints)
            return tri.simplices
        except:
            # Fallback: simple tetrahedra for small meshes
            print("Warning: Using simplified tetrahedra")
            return self._generate_simple_tetrahedra(len(kpoints))
    
    def _generate_tetrahedra_uniform(self, mesh_size: Tuple[int, int, int]) -> np.ndarray:
        """Generate tetrahedra for uniform mesh"""
        
        nx, ny, nz = mesh_size
        tetrahedra = []
        
        # Create tetrahedra from cubic cells
        for i in range(nx - 1):
            for j in range(ny - 1):
                for k in range(nz - 1):
                    # Each cube is divided into 6 tetrahedra
                    base = i * ny * nz + j * nz + k
                    
                    # Define 8 corners of cube
                    corners = [
                        base,                           # (i,j,k)
                        base + nz,                      # (i,j+1,k)
                        base + ny * nz,                 # (i+1,j,k)
                        base + ny * nz + nz,            # (i+1,j+1,k)
                        base + 1,                       # (i,j,k+1)
                        base + nz + 1,                  # (i,j+1,k+1)
                        base + ny * nz + 1,             # (i+1,j,k+1)
                        base + ny * nz + nz + 1         # (i+1,j+1,k+1)
                    ]
                    
                    # Add 6 tetrahedra (standard decomposition)
                    cube_tetrahedra = [
                        [corners[0], corners[1], corners[2], corners[4]],
                        [corners[1], corners[2], corners[3], corners[7]],
                        [corners[1], corners[4], corners[5], corners[7]],
                        [corners[2], corners[4], corners[6], corners[7]],
                        [corners[1], corners[2], corners[4], corners[7]],
                        [corners[0], corners[2], corners[4], corners[6]]
                    ]
                    
                    tetrahedra.extend(cube_tetrahedra)
        
        return np.array(tetrahedra)
    
    def _generate_simple_tetrahedra(self, num_kpoints: int) -> np.ndarray:
        """Simple tetrahedra for fallback"""
        
        # Very simple approach: group every 4 points into tetrahedra
        tetrahedra = []
        for i in range(0, num_kpoints - 3, 4):
            tetrahedra.append([i, i+1, i+2, i+3])
        
        return np.array(tetrahedra)
    
    def calculate_full_bz_bands(self, kpoints: np.ndarray, 
                               include_projections: bool = True) -> Tuple[np.ndarray, Dict]:
        """Calculate band structure over full BZ"""
        
        print(f"Calculating bands for {len(kpoints)} k-points...")
        start_time = time.time()
        
        if include_projections:
            band_energies, projections = self.calculator.calculate_projected_bands(kpoints)
        else:
            band_energies = self.calculator.calculate_bands(kpoints)
            projections = {}
        
        elapsed_time = time.time() - start_time
        print(f"Band calculation completed in {elapsed_time:.2f} seconds")
        
        return band_energies, projections


class DOSCalculator:
    """Enhanced DOS calculator with tetrahedron method and projections"""
    
    def __init__(self, full_bz_calculator: FullBZCalculator):
        self.calculator = full_bz_calculator
    
    def calculate_dos(self, mesh_size: Tuple[int, int, int], 
                     energy_range: Optional[Tuple[float, float]] = None,
                     num_energies: int = 1000,
                     use_symmetry: bool = True,
                     use_tetrahedron: bool = True,
                     gaussian_smearing: float = 0.1) -> Dict:
        """
        Calculate total and projected DOS
        
        Returns:
            Dictionary with DOS data and metadata
        """
        
        print("Generating k-point mesh...")
        # For the tetrahedron method we must use the full uniform mesh so
        # the tetrahedra connectivity corresponds to the regular grid.
        # Symmetry reduction breaks the regular connectivity and leads to
        # invalid tetrahedra (and spurious peaks). Force use_symmetry=False
        # when using tetrahedron.
        mesh_symmetry = use_symmetry if not use_tetrahedron else False
        kpoints, weights, tetrahedra = self.calculator.generate_kpoint_mesh(
            mesh_size, use_symmetry=mesh_symmetry)
        
        print("Calculating band structure...")
        band_energies, projections = self.calculator.calculate_full_bz_bands(
            kpoints, include_projections=True)
        
        # Determine energy range
        if energy_range is None:
            emin = np.min(band_energies) - 1.0
            emax = np.max(band_energies) + 1.0
        else:
            emin, emax = energy_range
        
        energy_grid = np.linspace(emin, emax, num_energies)
        
        dos_data = {
            'energies': energy_grid,
            'total_dos': None,
            'orbital_dos': {},
            'site_dos': {},
            'spin_dos': {},
            'kpoints': kpoints,
            'weights': weights,
            'band_energies': band_energies
        }
        
        if use_tetrahedron and len(tetrahedra) > 0:
            print("Calculating DOS using tetrahedron method...")
            tet_dos = TetrahedronDOS(kpoints, tetrahedra, weights)
            
            # Total DOS
            dos_data['total_dos'] = tet_dos.calculate_dos(band_energies, energy_grid)
            
            # Projected DOS
            if projections:
                self._calculate_projected_dos_tetrahedron(
                    tet_dos, band_energies, energy_grid, projections, dos_data)
        
        else:
            print("Calculating DOS using Gaussian broadening...")
            self._calculate_dos_gaussian(band_energies, energy_grid, weights, 
                                       projections, dos_data, gaussian_smearing)
        
        return dos_data
    
    def calculate_dos_with_moments(self, mesh_size: Tuple[int, int, int],
                                  target_electrons: Optional[float] = None,
                                  temperature: float = 300.0,
                                  energy_range: Optional[Tuple[float, float]] = None,
                                  num_energies: int = 1000,
                                  use_symmetry: bool = True,
                                  use_tetrahedron: bool = True,
                                  gaussian_smearing: float = 0.1,
                                  output_prefix: str = "dos_analysis") -> Tuple[Dict, DOSMomentsCalculator]:
        """
        Calculate DOS and analyze moments (occupations, band centers, widths)
        
        Parameters:
        -----------
        mesh_size : tuple
            (nx, ny, nz) k-point mesh dimensions
        target_electrons : float, optional
            Number of electrons to determine Fermi level
        temperature : float
            Temperature in K for Fermi-Dirac distribution
        energy_range : tuple, optional
            (emin, emax) energy range for DOS
        num_energies : int
            Number of energy points
        use_symmetry : bool
            Use symmetry reduction
        use_tetrahedron : bool
            Use tetrahedron method
        output_prefix : str
            Prefix for output files
            
        Returns:
        --------
        dos_data : dict
            DOS data dictionary
        moments_calc : DOSMomentsCalculator
            Moments calculator with results
        """
        
        # Calculate DOS
        print("Calculating DOS...")
        dos_data = self.calculate_dos(
            mesh_size, energy_range, num_energies,
            use_symmetry, use_tetrahedron, gaussian_smearing)
        
        if dos_data['total_dos'] is None:
            raise ValueError("DOS calculation failed")
        
        # Prepare projected DOS for moments analysis
        projected_dos = {}
        
        # Add orbital projections
        for orbital_name, orbital_dos in dos_data['orbital_dos'].items():
            projected_dos[f'orbital_{orbital_name}'] = orbital_dos
        
        # Add site projections  
        for site_name, site_dos in dos_data['site_dos'].items():
            projected_dos[f'site_{site_name}'] = site_dos
        
        # Add spin projections
        for spin_name, spin_dos in dos_data['spin_dos'].items():
            projected_dos[f'spin_{spin_name}'] = spin_dos
        
        # Calculate moments
        print("Analyzing DOS moments...")
        moments_calc = analyze_dos_moments(
            dos_data['energies'], 
            dos_data['total_dos'],
            projected_dos,
            target_electrons,
            temperature,
            output_prefix
        )
        
        # Add moments data to dos_data
        dos_data['moments'] = moments_calc.moments
        dos_data['fermi_level'] = moments_calc.fermi_level
        dos_data['temperature'] = temperature
        
        return dos_data, moments_calc
    
    def _calculate_projected_dos_tetrahedron(self, tet_dos: TetrahedronDOS,
                                           band_energies: np.ndarray, energy_grid: np.ndarray,
                                           projections: Dict, dos_data: Dict):
        """Calculate projected DOS using tetrahedron method"""
        
        # Orbital projections
        for orbital in projections['orbital']:
            proj_weights = projections['orbital'][orbital][:, :, np.newaxis]
            dos_data['orbital_dos'][orbital] = tet_dos.calculate_dos(
                band_energies, energy_grid, proj_weights)[:, 0]
        
        # Site projections
        for site in projections['site']:
            proj_weights = projections['site'][site][:, :, np.newaxis]
            dos_data['site_dos'][site] = tet_dos.calculate_dos(
                band_energies, energy_grid, proj_weights)[:, 0]
        
        # Spin projections
        for spin in projections['spin']:
            proj_weights = projections['spin'][spin][:, :, np.newaxis]
            dos_data['spin_dos'][spin] = tet_dos.calculate_dos(
                band_energies, energy_grid, proj_weights)[:, 0]
    
    def _calculate_dos_gaussian(self, band_energies: np.ndarray, energy_grid: np.ndarray,
                              weights: np.ndarray, projections: Dict, dos_data: Dict,
                              broadening: float = 0.1):
        """Fallback Gaussian broadening method"""
        
        dos_data['total_dos'] = np.zeros(len(energy_grid))
        
        # Initialize projected DOS
        if projections:
            for orbital in projections['orbital']:
                dos_data['orbital_dos'][orbital] = np.zeros(len(energy_grid))
            for site in projections['site']:
                dos_data['site_dos'][site] = np.zeros(len(energy_grid))
            for spin in projections['spin']:
                dos_data['spin_dos'][spin] = np.zeros(len(energy_grid))
        
        # Calculate DOS
        for k_idx in range(len(band_energies)):
            weight = weights[k_idx]
            
            for band in range(band_energies.shape[1]):
                energy = band_energies[k_idx, band]
                
                # Gaussian contribution
                gaussian = np.exp(-((energy_grid - energy) / broadening)**2)
                gaussian /= (broadening * np.sqrt(np.pi))
                
                dos_data['total_dos'] += weight * gaussian
                
                # Projected contributions
                if projections:
                    for orbital in projections['orbital']:
                        proj_weight = projections['orbital'][orbital][k_idx, band]
                        dos_data['orbital_dos'][orbital] += weight * proj_weight * gaussian
                    
                    for site in projections['site']:
                        proj_weight = projections['site'][site][k_idx, band]
                        dos_data['site_dos'][site] += weight * proj_weight * gaussian
                    
                    for spin in projections['spin']:
                        proj_weight = projections['spin'][spin][k_idx, band]
                        dos_data['spin_dos'][spin] += weight * proj_weight * gaussian
    
    def plot_dos(self, dos_data: Dict, fermi_energy: float = 0.0, 
                output_file: str = None, show_projections: bool = True):
        """Plot comprehensive DOS analysis"""
        
        energies = dos_data['energies']
        
        if show_projections:
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            axes = [ax]
        
        # Plot total DOS
        ax = axes[0]
        ax.plot(energies, dos_data['total_dos'], 'k-', linewidth=2, label='Total DOS')
        ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7, label='E_F')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS (states/eV)')
        ax.set_title('Total Density of States')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if show_projections and len(axes) > 1:
            # Orbital DOS
            ax = axes[1]
            colors = ['red', 'blue', 'green']
            for i, (orbital, color) in enumerate(zip(['s', 'p', 'd'], colors)):
                if orbital in dos_data['orbital_dos']:
                    ax.plot(energies, dos_data['orbital_dos'][orbital], 
                           color=color, linewidth=2, label=f'{orbital}-orbital')
            ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7)
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('DOS (states/eV)')
            ax.set_title('Orbital Projected DOS')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Site DOS
            if len(axes) > 2:
                ax = axes[2]
                num_sites = len(dos_data['site_dos'])
                colors = plt.cm.Set1(np.linspace(0, 1, num_sites))
                for site, color in zip(dos_data['site_dos'], colors):
                    ax.plot(energies, dos_data['site_dos'][site], 
                           color=color, linewidth=2, label=f'Site {site+1}')
                ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7)
                ax.set_xlabel('Energy (eV)')
                ax.set_ylabel('DOS (states/eV)')
                ax.set_title('Site Projected DOS')
                ax.legend()
                ax.grid(True, alpha=0.3)
            
            # Spin DOS
            if len(axes) > 3:
                ax = axes[3]
                if 'up' in dos_data['spin_dos'] and 'down' in dos_data['spin_dos']:
                    ax.plot(energies, dos_data['spin_dos']['up'], 
                           'r-', linewidth=2, label='Spin up')
                    ax.plot(energies, -dos_data['spin_dos']['down'], 
                           'b-', linewidth=2, label='Spin down')
                    ax.axhline(0, color='k', linestyle='-', alpha=0.3)
                ax.axvline(fermi_energy, color='r', linestyle='--', alpha=0.7)
                ax.set_xlabel('Energy (eV)')
                ax.set_ylabel('DOS (states/eV)')
                ax.set_title('Spin Projected DOS')
                ax.legend()
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"DOS plot saved to {output_file}")
        
        plt.show()


def main():
    """Main function for full BZ DOS calculation"""
    
    parser = argparse.ArgumentParser(description='Calculate full BZ DOS with tetrahedron method')
    parser.add_argument('--lattice', default='lattice.nml', 
                       help='Lattice file (default: lattice.nml)')
    parser.add_argument('--hamiltonian', default='rs2txt.dat',
                       help='Hamiltonian file (default: rs2txt.dat)')
    parser.add_argument('--mesh', nargs=3, type=int, default=[8, 8, 8],
                       help='k-point mesh size (default: 8 8 8)')
    parser.add_argument('--energy-range', nargs=2, type=float,
                       help='Energy range for DOS (emin emax)')
    parser.add_argument('--num-energies', type=int, default=1000,
                       help='Number of energy points (default: 1000)')
    parser.add_argument('--no-symmetry', action='store_true',
                       help='Disable symmetry reduction')
    parser.add_argument('--no-tetrahedron', action='store_true',
                       help='Use Gaussian broadening instead of tetrahedron method')
    parser.add_argument('--gaussian-smearing', type=float, default=0.1,
                       help='Gaussian smearing width in eV (default: 0.1)')
    parser.add_argument('--output-prefix', default='full_bz_dos',
                       help='Output file prefix (default: full_bz_dos)')
    parser.add_argument('--fermi-energy', type=float, default=0.0,
                       help='Fermi energy (default: 0.0)')
    
    args = parser.parse_args()
    
    print("Setting up full BZ calculator...")
    lattice_reader = LatticeReader(args.lattice)
    hamiltonian_reader = HamiltonianReader(args.hamiltonian)
    full_bz_calc = FullBZCalculator(lattice_reader, hamiltonian_reader)
    
    print("Setting up DOS calculator...")
    dos_calc = DOSCalculator(full_bz_calc)
    
    energy_range = tuple(args.energy_range) if args.energy_range else None
    
    print("Calculating DOS...")
    dos_data = dos_calc.calculate_dos(
        mesh_size=tuple(args.mesh),
        energy_range=energy_range,  # type: ignore
        num_energies=args.num_energies,
        use_symmetry=not args.no_symmetry,
        use_tetrahedron=not args.no_tetrahedron,
        gaussian_smearing=args.gaussian_smearing
    )
    
    print("Plotting results...")
    dos_calc.plot_dos(dos_data, fermi_energy=args.fermi_energy,
                     output_file=f"{args.output_prefix}.png")
    
    # Save DOS data
    output_data = {
        'energies': dos_data['energies'],
        'total_dos': dos_data['total_dos'],
        'orbital_dos': dos_data['orbital_dos'],
        'site_dos': dos_data['site_dos'],
        'spin_dos': dos_data['spin_dos']
    }
    
    np.savez(f"{args.output_prefix}_data.npz", **output_data)
    print(f"DOS data saved to {args.output_prefix}_data.npz")
    
    print("Done!")


if __name__ == "__main__":
    main()
