import numpy as np
import pyvista as pv
from .structure import Structure
from .glyphs import Atom  # assuming these are defined appropriately

class Visualizer:
    """
    Provides interactive toggles for visualizing different aspects of a
    structure. Interactive key bindings:
      C  : Toggle cell representation (Unit/Conventional)
      M  : Toggle spins
      E  : Toggle cell extents
      B  : Toggle bonds
      2  : Toggle 2d periodicity
      3  : Toggle 3d periodicity
      P  : Toggle perspective
      S  : Save screenshot
    Additionally, icon buttons appear along the bottom of the window.
    """
    def __init__(self, structure: Structure):
        """
        Initializes the Visualizer.
        """
        self.structure = structure
        self.plotter = pv.Plotter()
        # --- Quality settings ---
        self.plotter.add_light(pv.Light(position=(10, 10, 10), focal_point=(0, 0, 0), intensity=0.1))
        self.plotter.enable_ssao()
        self.plotter.enable_anti_aliasing()
        self.plotter.set_background("white")
        self.plotter.enable_depth_peeling()
        
        self.use_conventional = False  # Default: unit cell

        # --- Default data selection ---
        # We assume your structure object provides two representations:
        # atoms_uc, bonds_uc, spins_uc for the unit cell and
        # atoms_cc, bonds_cc, spins_cc for the conventional cell.
        if self.structure.unit_cell is not None:
            self.use_conventional = True  # True: conventional; False: unit cell
            
        self._change_representation()  # sets self.atoms, self.bonds, self.spins, self.cell

        # --- Toggle options (initial states) ---
        self.show_bonds = False
        self.show_spins = False
        self.show_extents = False
        self.periodicity = None  # Options: None, '2d', '3d'

        # --- Set up interactive key callbacks ---
        self.plotter.add_key_event("C", self.toggle_cell)
        self.plotter.add_key_event("M", self.toggle_spins)
        self.plotter.add_key_event("E", self.toggle_extents)
        self.plotter.add_key_event("B", self.toggle_bonds)
        self.plotter.add_key_event("2", lambda: self.toggle_periodicity("2d"))
        self.plotter.add_key_event("3", lambda: self.toggle_periodicity("3d"))
        self.plotter.add_key_event("P", self.toggle_perspective)
        self.plotter.add_key_event("S", self.save_screenshot)

        # Optionally, display instructions in text.
        instructions = ("C: Toggle Cell    M: Toggle Spins    E: Toggle Extents    "
                        "B: Toggle Bonds    2: 2d Periodicity    3: 3d Periodicity   "
                        "P  Toggle Perspective    S: Save Screenshot")
        self.plotter.add_text(instructions, position="upper_left", font_size=12, color="black")

        # --- Initialize lists for tracking geometry actors ---
        self._atom_actors = []
        self._bond_actors = []
        self._spin_actors = []
        self._periodic_actors = []
        self._extent_actor = None
        
        self._add_axis_glyph()  # Add a 3D axis glyph to the lower right corner

    def _add_icons(self):
        """
        Adds checkbox button widgets at fixed positions in the render window.
        The positions are in pixel coordinates (relative to the lower-left corner).
        Each button is associated with a callback.
        """
        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_cell,
            value=self.use_conventional,
            position=(20, 20),  # Button for "R": Toggle Cell representation
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("C", position=(30, 30), font_size=14, color="black", font="arial")

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_spins,
            value=self.show_spins,
            position=(90, 20),  # Button for "M": Toggle Spins
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("M", position=(100, 30), font_size=14, color="black", font="arial")

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_extents,
            value=self.show_extents,
            position=(160, 20),  # Button for "C": Toggle Cell Extents
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("E", position=(170, 30), font_size=14, color="black", font="arial")

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_bonds,
            value=self.show_bonds,
            position=(230, 20),  # Button for "B": Toggle Bonds
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("B", position=(240, 30), font_size=14, color="black", font="arial")

        self.plotter.add_checkbox_button_widget(
            callback=lambda _: self.toggle_periodicity("2d"),
            value=(self.periodicity == "2d"),
            position=(300, 20),  # Button for "2D": Toggle 2D periodicity
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("2D", position=(305, 30), font_size=14, color="black", font="arial")

        self.plotter.add_checkbox_button_widget(
            callback=lambda _: self.toggle_periodicity("3d"),
            value=(self.periodicity == "3d"),
            position=(370, 20),  # Button for "3D": Toggle 3D periodicity
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )

        self.plotter.add_text("3D", position=(375, 30), font_size=14, color="black", font="arial")

        self.plotter.add_checkbox_button_widget(
            callback=self.toggle_perspective,
            value=self.plotter.camera.parallel_projection,
            position=(440, 20),  # Button for "P": Toggle Perspective
            size=50,
            border_size=2,
            color_on="green",
            color_off="red"
        )
        self.plotter.add_text("P", position=(450, 30), font_size=14, color="black", font="arial")
        
    def _add_axis_glyph(self):
        """
        Adds a 3D axis glyph to the lower right corner of the render window.
        """
        #axes = pv.Axes(show_actor=True, actor_scale=0.15)
        #self.plotter.add_actor(axes.actor, reset_camera=False)
        self.plotter.add_axes(interactive=True, line_width=2, viewport=[0.8, 0.0, 1.0, 0.2])

    def _clear_geo_actors(self):
        """
        Remove only the geometry actors (atoms, bonds, spins, extents, and periodic images).
        This avoids wiping out lights, camera settings, and icon widgets.
        """
        for actor in self._bond_actors:
            self.plotter.remove_actor(actor)

        for actor in self._atom_actors:
            self.plotter.remove_actor(actor)

        for actor in self._spin_actors:
            self.plotter.remove_actor(actor)

        if self._extent_actor is not None:
            self.plotter.remove_actor(self._extent_actor)
            self._extent_actor = None

        for actor in self._periodic_actors:
            self.plotter.remove_actor(actor)

        self._bond_actors.clear()
        self._atom_actors.clear()
        self._spin_actors.clear()
        self._periodic_actors.clear()

    def save_screenshot(self, event=None):
        """
        Saves a screenshot of the current visualization.
        """
        filename = "screenshot.png"
        self.plotter.screenshot(filename)
        print(f"Screenshot saved as {filename}")

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

    def toggle_perspective(self, event=None):
        if self.plotter.camera.parallel_projection:
            self.plotter.camera.parallel_projection = False
        else:
            self.plotter.camera.parallel_projection = True
        self.render()

    def get_cell(self):
        return self.cell

    def expand_structure(self, cell, grid_dims):
        atoms_expanded = []
        i_min, i_max = -grid_dims[0] // 2, grid_dims[0] // 2
        j_min, j_max = -grid_dims[1] // 2, grid_dims[1] // 2
        if grid_dims[2] == 1:
            k_min, k_max = 0, 0
        else:
            k_min, k_max = -grid_dims[2] // 2, grid_dims[2] // 2 
        for i in range(i_min, i_max + 1):
            for j in range(j_min, j_max + 1):
                for k in range(k_min, k_max + 1):
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
        # Remove previously rendered geometry without clearing the whole scene.
        self._clear_geo_actors()

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

        # Start (or update) the interactive window.
        self.plotter.show(auto_close=False)

    def _render_atoms(self):
        points = np.array([atom.position for atom in self.atoms])
        radii = np.array([atom.radius for atom in self.atoms])
        colors = np.array([pv.Color(atom.color).float_rgb for atom in self.atoms])  # Convert hex to RGB

        # Create a PolyData object with points
        poly_data = pv.PolyData(points)

        # Add radii as point data (glyph will use this for scaling)
        poly_data["radii"] = radii

        # Create a sphere glyph
        sphere = pv.Sphere(theta_resolution=16, phi_resolution=16)

        # Use glyph filter to create spheres at each point
        glyphs = poly_data.glyph(scale="radii", geom=sphere)

        # **Reassign colors after glyph creation**
        glyphs["colors"] = np.repeat(colors, sphere.n_points, axis=0)  # Expand colors for all points

        # Add the glyphs to the plotter with proper coloring
        actor = self.plotter.add_mesh(glyphs, scalars="colors", rgb=True, smooth_shading=True)
        self._atom_actors.append(actor)

    def _render_bonds(self):
        cylinders = []
        for bond in self.bonds:
            vector = bond.end - bond.start
            length = np.linalg.norm(vector)
            center = (bond.start + bond.end) / 2
            direction = vector / length if length else [0, 0, 1]
            cylinder = pv.Cylinder(center=center, direction=direction,
                       radius=bond.radius, height=length)
            cylinders.append(cylinder)

        # Merge all cylinders into a single PolyData object
        if len(cylinders) > 1:
            merged_cylinders = cylinders[0].merge(cylinders[1:])  # Merge all into one
        else:
            merged_cylinders = cylinders[0]
        
        # Add the merged cylinders to the plotter
        actor = self.plotter.add_mesh(merged_cylinders, color="red", smooth_shading=True)
        self._bond_actors.append(actor)

    def _render_spins(self):
        for atom, spin in self.spins:
            if spin.magnitude > 0:
                arrow = pv.Arrow(start=atom.position,
                                 direction=spin.direction,
                                 tip_length=0.2 * spin.magnitude,
                                 tip_radius=spin.tip_radius,
                                 shaft_radius=spin.base_radius,
                                 scale=spin.magnitude)
                self._spin_actors.append(arrow)

        if self._spin_actors:
            merged_arrows = self._spin_actors[0].merge(self._spin_actors[1:])
            actor = self.plotter.add_mesh(merged_arrows, color=spin.color, smooth_shading=True)
            self._spin_actors = [actor]

    def _render_extents(self):
        cell = self.cell
        if cell is not None:
            origin = np.array([0, 0, 0])
            a, b, c = cell[0], cell[1], cell[2]
            corners = np.array([origin, a, a+b, b, c, a+c, a+b+c, b+c])
            edges = [ (0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 5), (2, 6),
                      (3, 7), (4, 5), (5, 6), (6, 7), (7, 4), ]
            lines = pv.PolyData()
            lines.points = np.array(corners)
            lines.lines = np.hstack([[2, edge[0], edge[1]] for edge in edges])
    
            self._extent_actor = self.plotter.add_mesh(lines, color="black", line_width=3)
            # corners = np.array([origin, a, b, c, a + b, a + c, b + c, a + b + c])
            # xmin, ymin, zmin = corners.min(axis=0)
            # xmax, ymax, zmax = corners.max(axis=0)
            # box = pv.Box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))
            # self._extent_actor = self.plotter.add_mesh(box, style="wireframe", color="black")

    def _render_periodic_images(self):
        grid_dims = (4, 4, 1) if self.periodicity == "2d" else (4, 4, 4)
        cell = self.cell
        if cell is not None:
            atoms_expanded = self.expand_structure(cell, grid_dims)
            points = np.array([atom.position for atom in atoms_expanded])
            radii = np.array([atom.radius for atom in atoms_expanded])
            colors = np.array([pv.Color(atom.color).float_rgb for atom in atoms_expanded])

            # Create a PolyData object with points
            poly_data = pv.PolyData(points)

            # Add radii as point data
            poly_data["radii"] = radii

            # Create a sphere glyph
            sphere = pv.Sphere(theta_resolution=16, phi_resolution=16)

            # Use glyph filter to create spheres at each point
            glyphs = poly_data.glyph(scale="radii", geom=sphere)

            # **Reassign colors after glyph creation**
            glyphs["colors"] = np.repeat(colors, sphere.n_points, axis=0)  # Expand colors for all points

            # Add the glyphs to the plotter with opacity
            actor = self.plotter.add_mesh(glyphs, scalars="colors", rgb=True, opacity=0.5, smooth_shading=True)
            self._periodic_actors.append(actor)

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
