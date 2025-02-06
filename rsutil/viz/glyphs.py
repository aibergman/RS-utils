import numpy as np
from helpers import get_element_color, get_atomic_radius, load_cpk_colors


# ----------------------------
# 1. Atom, Bond and Spin Classes
# ----------------------------
class Atom:
    def __init__(self, position, species, atom_type, radius=None, color=None, **kwargs):
        """
        Represents an atom for visualization.

        Args:
            position (list or array): [x, y, z] position.
            species (str): Chemical symbol.
            atom_type (str): Atom type (can be the same as species).
            radius (float): Radius for the sphere.
            color (str): Color used to render. If None, a standard color is assigned.
            kwargs: Any additional render information.
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
    def __init__(self, start, end, radius=0.2, color="grey", **kwargs):
        """
        Represents a bond (rendered as a cylinder).

        Args:
            start (list or array): [x, y, z] start position.
            end (list or array): [x, y, z] end position.
            radius (float): Cylinder radius.
            color (str): Cylinder color.
            kwargs: Additional render-related options.
        """
        self.start = np.array(start)
        self.end = np.array(end)
        self.radius = radius
        self.color = color
        self.render_info = kwargs


class Spin:
    def __init__(
        self,
        position,
        direction=(0, 1, 0),
        magnitude=0.0,
        base_radius=0.1,
        tip_radius=0.05,
        color="red",
        **kwargs,
    ):
        """
        Represents a magnetic moment as an arrow.

        Args:
            position (list or array): Position of the atom.
            direction (tuple or array): Direction of the spin.
            magnitude (float): Length scale factor.
            base_radius (float): Radius of the arrow’s shaft.
            tip_radius (float): Radius of the arrow’s tip.
            color (str): Color for the arrow.
            kwargs: Any additional information.
        """
        self.position = np.array(position)
        self.direction = np.array(direction)
        self.magnitude = magnitude
        self.base_radius = base_radius
        self.tip_radius = tip_radius
        self.color = color
        self.render_info = kwargs
