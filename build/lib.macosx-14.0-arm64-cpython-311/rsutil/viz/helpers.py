"""
This module provides utility functions for visualizing elements, including
loading CPK colors from a JSON file, retrieving element colors, and obtaining
atomic radii.
Functions:
    load_cpk_colors(file_path: str) -> dict:
    get_element_color(element_symbol: str, cpk_colors: dict) -> str:
    get_atomic_radius(element_symbol: str) -> float:
"""

import os
import json
from pymatgen.core.periodic_table import Element

ANG_TO_BOHR = 1.88973


def load_cpk_colors(file_path: str = "cpk_colors.json") -> dict:
    """
    Loads CPK colors from an external JSON file.

    Args:
        file_path (str): Path to the JSON file containing element colors.

    Returns:
        dict: A dictionary mapping element symbols to hex color codes.
    """
    module_dir = os.path.dirname(os.path.abspath(__file__))
    with open(module_dir + "/" + file_path, "r", encoding="utf-8") as file:
        return json.load(file)


def get_element_color(element_symbol: str, cpk_colors: dict) -> str:
    """
    Returns the standard CPK color for the given element.

    Args:
        element_symbol (str): Atomic symbol (e.g., "H", "O", "C").
        cpk_colors (dict): Dictionary of element colors.

    Returns:
        str: Hex color code (e.g., "#FFFFFF"), or default gray if not found.
    """
    return cpk_colors.get(element_symbol, "#A0A0A0")  # Default to gray


def get_atomic_radius(element_symbol: str) -> float:
    """
    Returns the atomic radius of the given element symbol.

    Args:
        element_symbol (str): The atomic symbol (e.g., "Fe", "O").

    Returns:
        float: Atomic radius in Ångströms (or None if not available).
    """
    try:
        return Element(element_symbol).atomic_radius  # *ANG_TO_BOHR
    except ValueError:
        return 1.0  # Default to 1.0 Å if not found, assuming empty sphere
