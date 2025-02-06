import numpy as np
import pyvista as pv
from structure import Structure
from glyphs import Atom, Bond, Spin  # assuming these are defined appropriately
from itertools import product

class Visualizer:
    """
    Provides interactive toggles for visualizing different aspects of a
    structure. Interactive key bindings:
      R  : Toggle cell representation (Unit/Conventional)
      M  : Toggle spins
      C  : Toggle cell extents
      B  : Toggle bonds
      2  : Toggle 2d periodicity
      3  : Toggle 3d periodicity
    Additionally, icon buttons appear along the bottom of the window.
    """
    def __init__(self, structure: Structure):
        """
        Initializes the Visualizer.
        """
        self.structure = structure
        self.plotter = pv.Plotter()

        # --- Default data selection ---
        # We assume your structure object provides two representations:
        # atoms_uc, bonds_uc, spins_uc for the unit cell and
        # atoms_cc, bonds_cc, spins_cc for the conventional cell.
        self.use_conventional = True  # True: conventional; False: unit cell
        self._change_representation()  # sets self.atoms, self.bonds, self.spins, self.cell

        # --- Quality settings ---
        self.plotter.enable_anti_aliasing()
        self.plotter.enable_depth_peeling()
        self.plotter.set_background("white")
        self.plotter.add_light(pv.Light(position=(10, 10, 10), focal_point=(0, 0, 0), intensity=0.5))

        # --- Toggle options (initial states) ---
        self.show_bonds = False
        self.show_spins = False
        self.show_extents = True
        self.periodicity = None  # Options: None, '2d', '3d'

        # --- Set up interactive key callbacks ---
        self.plotter.add_key_event("R", self.toggle_cell)
        self.plotter.add_key_event("M", self.toggle_spins)
        self.plotter.add_key_event("C", self.toggle_extents)
        self.plotter.add_key_event("B", self.toggle_bonds)
        self.plotter.add_key_event("2", lambda: self.toggle_periodicity("2d"))
        self.plotter.add_key_event("3", lambda: self.toggle_periodicity("3d"))

        # Optionally, display instructions in text.
        instructions = ("R: Toggle Cell    M: Toggle Spins    C: Toggle Extents    "
                        "B: Toggle Bonds    2: 2d Periodicity    3: 3d Periodicity")
        self.plotter.add_text(instructions, position="upper_left", font_size=10, color="black")

        
        
    def _add_icons(self):
        """
        Adds checkbox button widgets at fixed positions in the render window.
        The positions are in pixel coordinates (relative to the lower-left corner).
        Each button is associated with a callback.
        """
        # The following positions (x, y) and sizes (in pixels) can be adjusted as needed.
        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_cell,
            value=self.use_conventional,
            position=(20, 20),  # Button for "R": Toggle Cell representation
            size=30,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("R", position=(20, 55), font_size=10, color="black")  # Label above the checkbox

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_spins,
            value=self.show_spins,
            position=(70, 20),  # Button for "M": Toggle Spins
            size=30,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("M", position=(70, 55), font_size=10, color="black")

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_extents,
            value=self.show_extents,
            position=(120, 20),  # Button for "C": Toggle Cell Extents
            size=30,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("C", position=(120, 55), font_size=10, color="black")

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_bonds,
            value=self.show_bonds,
            position=(170, 20),  # Button for "B": Toggle Bonds
            size=30,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("B", position=(170, 55), font_size=10, color="black")

        # For periodicity, we use lambda functions so that toggling one can override the other.
        self.plotter.add_checkbox_button_widget(
            callback=lambda _: self.toggle_periodicity("2d"),
            value=(self.periodicity == "2d"),
            position=(220, 20),  # Button for "2D": Toggle 2D periodicity
            size=30,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("2D", position=(220, 55), font_size=10, color="black")

        self.plotter.add_checkbox_button_widget(
            callback=lambda _: self.toggle_periodicity("3d"),
            value=(self.periodicity == "3d"),
            position=(270, 20),  # Button for "3D": Toggle 3D periodicity
            size=30,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("3D", position=(270, 55), font_size=10, color="black")

    def toggle_bonds(self, event=None):
        self.show_bonds = not self.show_bonds
        self.render()

    def toggle_spins(self, event=None):
        self.show_spins = not self.show_spins
        self.render()

    def toggle_cell(self, event=None):
        self.use_conventional = not self.use_conventional
        self._change_representation()  # update internal data from structure
        self.render()

    def toggle_extents(self, event=None):
        self.show_extents = not self.show_extents
        self.render()

    def toggle_periodicity(self, period="2d", event=None):
        if self.periodicity == period:
            self.periodicity = None
        else:
            self.periodicity = period
        self.render()

    def get_cell(self):
        return self.cell

    def expand_structure(self, cell, grid_dims):
        atoms_expanded = []
        for i in range(grid_dims[0]):
            for j in range(grid_dims[1]):
                for k in range(grid_dims[2]):
                    translation = np.dot(cell.T, [i, j, k])
                    for atom in self.atoms:
                        new_pos = atom.position + translation
                        atoms_expanded.append(
                            Atom(
                                position=new_pos,
                                species=atom.species,
                                atom_type=atom.atom_type,
                                radius=atom.radius,
                                color=atom.color,
                            )
                        )
        return atoms_expanded

    def render(self):
        self.plotter.clear()
        self._render_atoms()
        if self.show_bonds:
            self._render_bonds()
        if self.show_spins:
            self._render_spins()
        if self.show_extents:
            self._render_extents()
        if self.periodicity is not None:
            self._render_periodic_images()
        self.plotter.render()
        
        # --- Add icon buttons along the bottom ---
        self._add_icons()
        
        # -- Render ! --
        self.plotter.show(auto_close=False)


    def _render_atoms(self):
        for atom in self.atoms:
            sphere = pv.Sphere(radius=atom.radius, center=atom.position)
            self.plotter.add_mesh(sphere, color=atom.color, smooth_shading=True, specular=0.5,
                                  specular_power=10, lighting=True)

    def _render_bonds(self):
        for bond in self.bonds:
            vector = bond.end - bond.start
            length = np.linalg.norm(vector)
            center = (bond.start + bond.end) / 2
            direction = vector / length if length else [0, 0, 1]
            cylinder = pv.Cylinder(center=center, direction=direction,
                                   radius=bond.radius, height=length)
            self.plotter.add_mesh(cylinder, color=bond.color, smooth_shading=True)

    def _render_spins(self):
        for atom, spin in self.spins:
            if spin.magnitude > 0:
                arrow = pv.Arrow(start=atom.position,
                                 direction=spin.direction,
                                 tip_length=0.2 * spin.magnitude,
                                 tip_radius=spin.tip_radius,
                                 shaft_radius=spin.base_radius,
                                 scale=spin.magnitude)
                self.plotter.add_mesh(arrow, color=spin.color, smooth_shading=True)

    def _render_extents(self):
        cell = self.cell
        if cell is not None:
            origin = np.array([0, 0, 0])
            a, b, c = cell[0], cell[1], cell[2]
            corners = np.array([origin, a, b, c, a + b, a + c, b + c, a + b + c])
            xmin, ymin, zmin = corners.min(axis=0)
            xmax, ymax, zmax = corners.max(axis=0)
            box = pv.Box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))
            self.plotter.add_mesh(box, style="wireframe", color="black", smooth_shading=True)

    def _render_periodic_images(self):
        grid_dims = (2, 2, 1) if self.periodicity == "2d" else (2, 2, 2)
        cell = self.cell
        if cell is not None:
            atoms_expanded = self.expand_structure(cell, grid_dims)
            for atom in atoms_expanded:
                sphere = pv.Sphere(radius=atom.radius, center=atom.position)
                self.plotter.add_mesh(sphere, color=atom.color, opacity=0.5, smooth_shading=True)

    def _change_representation(self):
        """
        Updates the internal data from the structure based on whether
        the conventional or unit cell is used.
        """
        if self.use_conventional:
            self.atoms = self.structure.atoms_cc
            self.bonds = self.structure.bonds_cc
            self.spins = self.structure.spins_cc  # Expected as list of (Atom, Spin) tuples.
            self.cell = self.structure.conventional_cell
        else:
            self.atoms = self.structure.atoms_uc
            self.bonds = self.structure.bonds_uc
            self.spins = self.structure.spins_uc  # Expected as list of (Atom, Spin) tuples.
            self.cell = self.structure.unit_cell
