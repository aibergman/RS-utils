"""
K-point label formatting utilities.

Provides functions for parsing and formatting k-point labels, including 
support for Greek letters and prime convention notation.

Author: Ramon Cardias
"""

from typing import Tuple


def normalize_label(label: str) -> Tuple[str, bool]:
    """Parse k-point label and extract base name and prime flag.
    
    Parameters
    ----------
    label : str
        Label (e.g., "K", "Γ", "K'", "gamma'")
        
    Returns
    -------
    tuple
        (base_label_uppercase, is_prime)
        
    Author
    ------
    Ramon Cardias
    """
    s = label.strip().replace("Γ", "Gamma").replace("γ", "Gamma").replace("'", "'")
    
    if s.lower() in ("g", "gamma"):
        base = "GAMMA"
    elif s.endswith("'"):
        base = s[:-1].upper()
    else:
        base = s.upper()
    
    is_prime = s.endswith("'")
    return base, is_prime


def tex_format_label(label: str) -> str:
    """Format label for LaTeX rendering in plots.
    
    Parameters
    ----------
    label : str
        K-point label (e.g., "K", "Gamma", "K'")
        
    Returns
    -------
    str
        LaTeX-formatted label (e.g., r"$\\Gamma^{\\prime}$")
        
    Author
    ------
    Ramon Cardias
    """
    if label.upper() == "GAMMA":
        if label.endswith("'"):
            return r"$\Gamma^{\prime}$"
        return r"$\Gamma$"
    
    if label.endswith("'"):
        return r"$" + label[:-1] + r"^{\prime}$"
    
    return r"$" + label + r"$"
