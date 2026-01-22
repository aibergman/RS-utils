#!/usr/bin/env python3
"""
Input file readers for TB-LMTO data

This module provides readers for:
- lattice.nml files (crystal structure)
- rs2txt.dat files (Hamiltonian matrix elements)
"""

import numpy as np
import re
from typing import Dict, List


class LatticeReader:
    """Read and parse lattice.nml file"""
    
    def __init__(self, filename: str = "lattice.nml"):
        self.filename = filename
        self.lattice_vectors = None
        self.atomic_positions = None
        self.atomic_numbers = None
        self.nbas = None
        self.parse_lattice_file()
    
    def parse_lattice_file(self):
        """Parse the Fortran namelist file"""
        with open(self.filename, 'r') as f:
            content = f.read()
        
        # Remove comments
        lines = []
        for line in content.split('\n'):
            # Remove comments starting with !
            if '!' in line:
                line = line[:line.index('!')]
            lines.append(line.strip())
        
        content = '\n'.join(lines)
        
        # Extract lattice vectors
        a_pattern = r'a\(\s*:\s*,\s*(\d+)\s*\)\s*=\s*([\d\.,\s-]+)'
        matches = re.findall(a_pattern, content)
        
        if len(matches) >= 3:
            self.lattice_vectors = np.zeros((3, 3))
            for match in matches:
                col = int(match[0]) - 1  # Convert to 0-based indexing
                values = [float(x.strip()) for x in match[1].split(',') if x.strip()]
                if len(values) == 3 and col < 3:
                    self.lattice_vectors[:, col] = values
        
        # Extract atomic positions
        crd_pattern = r'crd\(\s*:\s*,\s*(\d+)\s*\)\s*=\s*([\d\.,\s-]+)'
        matches = re.findall(crd_pattern, content)
        
        if matches:
            self.atomic_positions = []
            for match in matches:
                values = [float(x.strip()) for x in match[1].split(',') if x.strip()]
                if len(values) == 3:
                    self.atomic_positions.append(values)
            self.atomic_positions = np.array(self.atomic_positions)
        
        # Extract atomic numbers
        izp_pattern = r'izp\(\s*(\d+)\s*\)\s*=\s*(\d+)'
        matches = re.findall(izp_pattern, content)
        
        if matches:
            self.atomic_numbers = []
            for match in matches:
                self.atomic_numbers.append(int(match[1]))
            self.atomic_numbers = np.array(self.atomic_numbers)
        
        # Extract nbas
        nbas_pattern = r'nbas\s*=\s*(\d+)'
        match = re.search(nbas_pattern, content)
        if match:
            self.nbas = int(match.group(1))
    
    def get_structure_info(self):
        """Return structure information for spglib"""
        return {
            'lattice': self.lattice_vectors.T,  # spglib expects row vectors
            'positions': self.atomic_positions,
            'numbers': self.atomic_numbers
        }


class HamiltonianReader:
    """Read and parse rs2txt.dat file containing the TB-LMTO Hamiltonian"""
    
    def __init__(self, filename: str = "rs2txt.dat"):
        self.filename = filename
        self.hamiltonian_data = {}
        self.orbital_size = 18  # 9 orbitals per spin channel
        self.parse_hamiltonian_file()
    
    def parse_hamiltonian_file(self):
        """Parse the Hamiltonian data file"""
        data = []
        
        with open(self.filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 7:
                    ntype = int(parts[0])
                    k = int(parts[1])
                    rij = [float(parts[2]), float(parts[3]), float(parts[4])]
                    real_part = float(parts[5])
                    imag_part = float(parts[6])
                    
                    data.append({
                        'ntype': ntype,
                        'neighbor': k,
                        'displacement': np.array(rij),
                        'matrix_element': complex(real_part, imag_part)
                    })
        
        # Organize data by atom type and neighbor
        self.organize_hamiltonian_data(data)
    
    def organize_hamiltonian_data(self, data):
        """Organize the Hamiltonian data into matrices"""
        # Group by atom type and neighbor
        for entry in data:
            ntype = entry['ntype']
            neighbor = entry['neighbor']
            
            if ntype not in self.hamiltonian_data:
                self.hamiltonian_data[ntype] = {}
            
            if neighbor not in self.hamiltonian_data[ntype]:
                self.hamiltonian_data[ntype][neighbor] = {
                    'displacement': entry['displacement'],
                    'matrix': np.zeros((self.orbital_size, self.orbital_size), dtype=complex)
                }
        
        # Fill in matrix elements (this assumes data is ordered by i, j indices)
        matrix_idx = {}
        for entry in data:
            ntype = entry['ntype']
            neighbor = entry['neighbor']
            
            if (ntype, neighbor) not in matrix_idx:
                matrix_idx[(ntype, neighbor)] = [0, 0]  # [i, j] indices
            
            i, j = matrix_idx[(ntype, neighbor)]
            self.hamiltonian_data[ntype][neighbor]['matrix'][i, j] = entry['matrix_element']
            
            # Update indices
            j += 1
            if j >= self.orbital_size:
                j = 0
                i += 1
            matrix_idx[(ntype, neighbor)] = [i, j]
