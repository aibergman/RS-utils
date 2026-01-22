#!/usr/bin/env python3
"""
Reader factory and format detection for Hamiltonian and lattice files.

Supports multiple formats: RS-LMTO (native), PAOFLOW, KITE.
"""

from .base import BaseLatticeReader, BaseHamiltonianReader
from .rs_lmto import LatticeReader, HamiltonianReader
from typing import Union, Optional


def detect_lattice_format(filename: str) -> str:
    """
    Attempt to auto-detect lattice file format from filename/content.
    
    Parameters
    ----------
    filename : str
        Path to lattice file
    
    Returns
    -------
    format : str
        Detected format: 'rs-lmto', 'paoflow', 'kite', or raises ValueError
    """
    if filename.endswith('.nml') or 'lattice' in filename.lower():
        return 'rs-lmto'
    elif 'paoflow' in filename.lower():
        return 'paoflow'
    elif 'kite' in filename.lower():
        return 'kite'
    else:
        # Try to peek at file content for format signatures
        try:
            with open(filename, 'r') as f:
                first_line = f.readline()
            
            if '&LATTICE' in first_line or '&ATLATC' in first_line:
                return 'rs-lmto'
            elif 'PAOFLOW' in first_line:
                return 'paoflow'
            elif 'KITE' in first_line:
                return 'kite'
        except:
            pass
        
        # Default to RS-LMTO
        return 'rs-lmto'


def detect_hamiltonian_format(filename: str) -> str:
    """
    Attempt to auto-detect Hamiltonian file format from filename/content.
    
    Parameters
    ----------
    filename : str
        Path to Hamiltonian file
    
    Returns
    -------
    format : str
        Detected format: 'rs-lmto', 'paoflow', 'kite', or raises ValueError
    """
    if filename.endswith('.dat') or 'rs2txt' in filename.lower():
        return 'rs-lmto'
    elif 'paoflow' in filename.lower():
        return 'paoflow'
    elif 'kite' in filename.lower():
        return 'kite'
    else:
        # Try to peek at file content for format signatures
        try:
            with open(filename, 'r') as f:
                first_line = f.readline()
            
            # RS-LMTO format: lines start with integer ntype, then k, then floats
            parts = first_line.strip().split()
            if len(parts) >= 7:
                try:
                    int(parts[0])
                    int(parts[1])
                    float(parts[2])
                    return 'rs-lmto'
                except ValueError:
                    pass
            
            if 'PAOFLOW' in first_line:
                return 'paoflow'
            elif 'KITE' in first_line:
                return 'kite'
        except:
            pass
        
        # Default to RS-LMTO
        return 'rs-lmto'


def load_lattice(filename: str, format: Optional[str] = None) -> BaseLatticeReader:
    """
    Load lattice from file, auto-detecting or using specified format.
    
    Parameters
    ----------
    filename : str
        Path to lattice file
    format : str, optional
        File format: 'rs-lmto' (default), 'paoflow', 'kite'.
        If None, attempts auto-detection.
    
    Returns
    -------
    reader : BaseLatticeReader
        Appropriate reader instance for the format
    
    Examples
    --------
    >>> lattice = load_lattice('lattice.nml')
    >>> lattice = load_lattice('paoflow_lattice.txt', format='paoflow')
    """
    if format is None:
        format = detect_lattice_format(filename)
    
    format = format.lower().strip()
    
    if format == 'rs-lmto':
        return LatticeReader(filename)
    elif format == 'paoflow':
        # Will be implemented later
        raise NotImplementedError(f"PAOFLOW lattice reader not yet implemented")
    elif format == 'kite':
        # Will be implemented later
        raise NotImplementedError(f"KITE lattice reader not yet implemented")
    else:
        raise ValueError(f"Unknown lattice format: {format}. "
                        f"Supported: 'rs-lmto', 'paoflow', 'kite'")


def load_hamiltonian(filename: str, format: Optional[str] = None) -> BaseHamiltonianReader:
    """
    Load Hamiltonian from file, auto-detecting or using specified format.
    
    Parameters
    ----------
    filename : str
        Path to Hamiltonian file
    format : str, optional
        File format: 'rs-lmto' (default), 'paoflow', 'kite'.
        If None, attempts auto-detection.
    
    Returns
    -------
    reader : BaseHamiltonianReader
        Appropriate reader instance for the format
    
    Examples
    --------
    >>> ham = load_hamiltonian('rs2txt.dat')
    >>> ham = load_hamiltonian('paoflow_ham.txt', format='paoflow')
    """
    if format is None:
        format = detect_hamiltonian_format(filename)
    
    format = format.lower().strip()
    
    if format == 'rs-lmto':
        return HamiltonianReader(filename)
    elif format == 'paoflow':
        # Will be implemented later
        raise NotImplementedError(f"PAOFLOW Hamiltonian reader not yet implemented")
    elif format == 'kite':
        # Will be implemented later
        raise NotImplementedError(f"KITE Hamiltonian reader not yet implemented")
    else:
        raise ValueError(f"Unknown Hamiltonian format: {format}. "
                        f"Supported: 'rs-lmto', 'paoflow', 'kite'")


__all__ = [
    'BaseLatticeReader',
    'BaseHamiltonianReader',
    'LatticeReader',
    'HamiltonianReader',
    'load_lattice',
    'load_hamiltonian',
    'detect_lattice_format',
    'detect_hamiltonian_format',
]
