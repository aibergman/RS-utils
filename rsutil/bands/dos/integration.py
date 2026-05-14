#!/usr/bin/env python3
"""
Tetrahedron method for accurate DOS integration.

Implements the tetrahedron method for computing density of states,
which is more accurate than Gaussian broadening, especially near
band edges and Van Hove singularities.

Theory:
  The tetrahedron method evaluates the density of states contribution
  from each band within each tetrahedron of the reciprocal space mesh.
  It uses analytical formulas for the integration within tetrahedra,
  providing exact integration for linear functions within each tetrahedron.

References:
  - Blochl et al., Phys. Rev. B 49, 16223 (1994)
  - Lehmann & Taut, Phys. Status Solidi B 54, 469 (1972)
"""

import numpy as np
from typing import Optional


class TetrahedronDOS:
    """Tetrahedron method for DOS calculation
    
    Parameters
    ----------
    kpoints : np.ndarray
        Array of k-points (nk, 3)
    tetrahedra : np.ndarray
        Array of tetrahedron vertex indices (ntet, 4)
    weights : np.ndarray
        Weight of each k-point (nk,)
    """
    
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
