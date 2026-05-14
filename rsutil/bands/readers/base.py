#!/usr/bin/env python3
"""
Abstract base classes for Hamiltonian and lattice readers.

All format-specific readers (RS-LMTO, PAOFLOW, KITE, etc.) inherit from these base classes.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, Optional


class BaseLatticeReader(ABC):
    """Abstract base class for lattice/structure readers"""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.lattice_vectors = None
        self.atomic_positions = None
        self.atomic_numbers = None
        self.nbas = None
        self.parse_lattice_file()
    
    @abstractmethod
    def parse_lattice_file(self):
        """Parse the lattice file (format-specific implementation)"""
        pass
    
    def get_structure_info(self) -> Dict:
        """Return structure information for spglib/pymatgen"""
        return {
            'lattice': self.lattice_vectors.T,  # spglib expects row vectors
            'positions': self.atomic_positions,
            'numbers': self.atomic_numbers
        }


class BaseHamiltonianReader(ABC):
    """Abstract base class for Hamiltonian readers"""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.hamiltonian_data = {}
        self.orbital_size = 18  # Default: 9 orbitals per spin channel (can be overridden)
        self.parse_hamiltonian_file()
    
    @abstractmethod
    def parse_hamiltonian_file(self):
        """Parse the Hamiltonian file (format-specific implementation)"""
        pass
    
    def fourier_transform_hamiltonian(self, k_point: np.ndarray) -> np.ndarray:
        """
        Fourier transform real-space Hamiltonian to k-space.
        
        This is a common operation for all formats that use real-space representation.
        Subclasses can override if needed for other representations.
        
        Parameters
        ----------
        k_point : array
            k-point in reciprocal space (units of 2π/a)
        
        Returns
        -------
        H_k : array
            k-space Hamiltonian matrix (complex)
        """
        if not self.hamiltonian_data:
            raise ValueError("No Hamiltonian data loaded")
        
        # Initialize with on-site (R=0) terms
        H_k = None
        
        for atom_type in self.hamiltonian_data:
            for neighbor_idx in self.hamiltonian_data[atom_type]:
                data = self.hamiltonian_data[atom_type][neighbor_idx]
                displacement = data['displacement']
                matrix = data['matrix']
                
                # Fourier phase: exp(i k·R)
                phase = np.exp(1j * np.dot(k_point, displacement) * 2 * np.pi)
                
                if H_k is None:
                    H_k = phase * matrix
                else:
                    H_k += phase * matrix
        
        if H_k is None:
            raise ValueError("Failed to construct k-space Hamiltonian")
        
        return H_k
