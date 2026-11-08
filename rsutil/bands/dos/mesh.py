#!/usr/bin/env python3
"""
Full Brillouin zone k-point mesh generation with symmetry reduction.

Provides utilities for generating k-point meshes over the full Brillouin zone,
with optional symmetry reduction to the irreducible Brillouin zone, and
tetrahedron generation for integration methods.
"""

import numpy as np
from typing import Dict, Tuple, Optional

try:
    import spglib
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False

from ..readers import LatticeReader, HamiltonianReader
from ..projections import ProjectedBandCalculator
import time
import warnings


class FullBZCalculator:
    """Full Brillouin zone calculator with symmetry reduction
    
    Generates k-point meshes with or without symmetry reduction, and
    creates tetrahedra for integration purposes.
    
    Parameters
    ----------
    lattice_reader : LatticeReader
        Crystal structure information
    hamiltonian_reader : HamiltonianReader
        Hamiltonian data
    """
    
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
        except Exception as exc:  # ImportError or a Qhull/degenerate-input failure
            # Fallback: simple tetrahedra for small meshes
            warnings.warn(
                f"Delaunay tetrahedralisation failed ({type(exc).__name__}: {exc}); "
                "using simplified tetrahedra",
                RuntimeWarning,
                stacklevel=2,
            )
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
