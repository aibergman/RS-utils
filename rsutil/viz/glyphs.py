"""
Classes for visualizing atoms, bonds, and spins.

Classes:
    Atom: Atom with position, species, radius, and color.
    Bond: Bond as a cylinder with start and end positions.
    Spin: Magnetic moment as an arrow with direction and magnitude.
"""
import numpy as np
from .helpers import get_element_color, get_atomic_radius, load_cpk_colors


# ----------------------------
# 1. Atom, Bond and Spin Classes
# ----------------------------
class Atom:
    """
    Class representing an atom.

    Parameters
    ----------
    position : list or np.ndarray
        [x, y, z] position.
    species : str
        Chemical symbol.
    atom_type : str
        Atom type (can be the same as species).
    radius : float, optional
        Radius for the sphere. Default is None.
    color : str, optional
        Color used to render. Default is None.

    Attributes
    ----------
    position : np.ndarray
        [x, y, z] position.
    species : str
        Chemical symbol.
    atom_type : str
        Atom type (can be the same as species).
    radius : float
        Radius for the sphere.
    color : str
        Color used to render.
    render_info : dict
        Additional information.
    """
    def __init__(
        self,
        position: np.ndarray,
        species: str,
        atom_type: str,
        radius: float = None,
        color: str = None,
        **kwargs,
    ):
        """
        Represents an atom for visualization.

        Parameters
        ----------
        position : list or np.ndarray
            [x, y, z] position.
        species : str
            Chemical symbol.
        atom_type : str
            Atom type (can be the same as species).
        radius : float, optional
            Radius for the sphere. Default is None.
        color : str, optional
            Color used to render. Default is None.
        """
        # Read colors from a JSON file
        self.cpk_colors = load_cpk_colors()

        self.position = np.array(position)
        self.species = species
        self.atom_type = atom_type
        self.radius = radius if radius is not None else get_atomic_radius(species)
        self.color = color if color is not None else get_element_color(species, self.cpk_colors)
        self.render_info = kwargs


class Bond:
    """
    Class representing a bond rendered as a cylinder.

    Parameters
    ----------
    start : list or np.ndarray
        [x, y, z] start position.
    end : list or np.ndarray
        [x, y, z] end position.
    radius : float, optional
        Cylinder radius, by default 0.2.
    color : str, optional
        Cylinder color, by default "grey".

    Attributes
    ----------
    start : np.ndarray
        [x, y, z] start position.
    end : np.ndarray
        [x, y, z] end position.
    radius : float
        Cylinder radius.
    color : str
        Cylinder color.
    render_info : dict
        Additional information.
    """
    def __init__(
        self,
        start: np.ndarray,
        end: np.ndarray,
        radius: float = 0.2,
        color: str = "grey",
        **kwargs,
    ):
        """
        Represents a bond as a cylinder.

        Parameters
        ----------
        start : list or np.ndarray
            [x, y, z] start position.
        end : list or np.ndarray
            [x, y, z] end position.
        radius : float, optional
            Cylinder radius, by default 0.2.
        color : str, optional
            Cylinder color, by default "grey".
        """
        self.start = np.array(start)
        self.end = np.array(end)
        self.radius = radius
        self.color = color
        self.render_info = kwargs


class Spin:
    """
    Class representing a magnetic moment as an arrow.

    Parameters
    ----------
    position : list or np.ndarray
        Position of the atom.
    direction : tuple or np.ndarray, optional
        Direction of the spin (default is (0, 1, 0)).
    magnitude : float, optional
        Length scale factor (default is 0.0).
    base_radius : float, optional
        Radius of the arrow’s shaft (default is 0.1).
    tip_radius : float, optional
        Radius of the arrow’s tip (default is 0.05).
    color : str, optional
        Color for the arrow (default is "red").

    Attributes
    ----------
    position : np.ndarray
        Position of the atom.
    direction : np.ndarray
        Direction of the spin.
    magnitude : float
        Length scale factor.
    base_radius : float
        Radius of the arrow’s shaft.
    tip_radius : float
        Radius of the arrow’s tip.
    color : str
        Color for the arrow.
    render_info : dict
        Additional information.
    """
    def __init__(
        self,
        position: list | np.ndarray,
        direction: tuple | np.ndarray = (0, 1, 0),
        magnitude: float = 0.0,
        base_radius: float = 0.1,
        tip_radius: float = 0.05,
        color: str = "red",
        **kwargs,
    ):
        """
        Represents a magnetic moment as an arrow.

        Parameters
        ----------
        position : list or np.ndarray
            Position of the atom.
        direction : tuple or np.ndarray, optional
            Direction of the spin.
        magnitude : float, optional
            Length scale factor.
        base_radius : float, optional
            Radius of the arrow’s shaft.
        tip_radius : float, optional
            Radius of the arrow’s tip.
        color : str, optional
            Color for the arrow.
        """
        self.position = np.array(position)
        self.direction = np.array(direction)
        self.magnitude = magnitude
        self.base_radius = base_radius
        self.tip_radius = tip_radius
        self.color = color
        self.render_info = kwargs
