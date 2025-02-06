import numpy as np
import pyvista as pv
from structure import Structure
from glyphs import Atom


# ----------------------------
# 4. Visualizer Class with Interactive Elements
# ----------------------------
class Visualizer:
    """
    Provides interactive toggles for visualizing different aspects of a
    structure, including bonds, spins, cell type, extents, and periodicity.
    
    Interactive key bindings (shown as text in the upper left):
      R  : Toggle cell representation (Unit/Conventional)
      M  : Toggle spins
      C  : Toggle cell extents
      B  : Toggle bonds
      2  : Toggle 2d periodicity
      3  : Toggle 3d periodicity
    """
    def __init__(self, structure: Structure):
        """
        Visualizer for a Structure using pyvista.
        """
        self.structure = structure
        self.plotter = pv.Plotter()

        # --- Default data selection ---
        # Here we assume that your Structure has already been prepared with:
        # atoms_uc, bonds_uc, spins_uc for the unit cell and
        # atoms_cc, bonds_cc, spins_cc for the conventional cell.
        #
        # Start with a default representation.
        self.use_conventional = True  # True: use conventional cell; False: unit cell
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
        self.periodicity = None  # Options: None, '2d', or '3d'

        # --- Set up interactive key callbacks ---
        # (Using single keys for periodicity: "2" for 2d, "3" for 3d.)
        self.plotter.add_key_event("R", self.toggle_cell)
        self.plotter.add_key_event("M", self.toggle_spins)
        self.plotter.add_key_event("C", self.toggle_extents)
        self.plotter.add_key_event("B", self.toggle_bonds)
        self.plotter.add_key_event("2", lambda: self.toggle_periodicity("2d"))
        self.plotter.add_key_event("3", lambda: self.toggle_periodicity("3d"))

        # Optionally display instructions in the render window.
        instructions = ("R: Toggle Cell (Unit/Conventional)    "
                        "M: Toggle Spins    "
                        "C: Toggle Extents    "
                        "B: Toggle Bonds    "
                        "2: 2d Periodicity    "
                        "3: 3d Periodicity")
        self.plotter.add_text(instructions, position="upper_left", font_size=10, color="black")

    def toggle_bonds(self):
        self.show_bonds = not self.show_bonds
        self.render()

    def toggle_spins(self):
        self.show_spins = not self.show_spins
        self.render()

    def toggle_cell(self):
        """
        Toggles the cell representation between unit and conventional.
        This updates the internal data (atoms, bonds, spins, cell) and re-renders.
        """
        self.use_conventional = not self.use_conventional
        self._change_representation()  # update self.atoms, etc.
        self.render()

    def toggle_extents(self):
        self.show_extents = not self.show_extents
        self.render()

    def toggle_periodicity(self, period="2d"):
        """
        Toggles periodicity.
        If the specified period is already enabled, it will be disabled.
        Otherwise, it sets the periodicity to the given period.
        Turning one on will override the other.
        
        Args:
            period (str): Either "2d" or "3d".
        """
        if self.periodicity == period:
            self.periodicity = None
        else:
            self.periodicity = period
        self.render()

    def get_cell(self):
        """
        Returns the current lattice cell matrix (unit or conventional)
        based on the active representation.
        """
        return self.cell

    def expand_structure(self, cell, grid_dims):
        """
        Expands the structure by replicating the cell.

        Args:
            cell (3x3 array-like): Lattice cell.
            grid_dims (tuple): (nx, ny, nz) number of replications.

        Returns:
            List of Atom objects for the periodic images.
        """
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
        
        self.plotter.show(auto_close=False)

    def _render_atoms(self):
        """
        Renders atoms as spheres.
        """
        for atom in self.atoms:
            sphere = pv.Sphere(radius=atom.radius, center=atom.position)
            self.plotter.add_mesh(sphere, color=atom.color, smooth_shading=True, specular=0.5,
                                  specular_power=10, lighting=True)

    def _render_bonds(self):
        """
        Renders bonds as cylinders.
        """
        # We assume get_bond_objects() returns bonds for the current representation.
        # (Alternatively, you could store bonds in self.bonds.)
        bonds = self.bonds
        for bond in bonds:
            vector = bond.end - bond.start
            length = np.linalg.norm(vector)
            center = (bond.start + bond.end) / 2
            direction = vector / length if length else [0, 0, 1]
            cylinder = pv.Cylinder(center=center, direction=direction,
                                   radius=bond.radius, height=length)
            self.plotter.add_mesh(cylinder, color=bond.color, smooth_shading=True)

    def _render_spins(self):
        """
        Renders spins as arrows.
        """
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
        """
        Renders the extents of the current cell as a wireframe box.
        """
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
        """
        Renders periodic images of the structure with reduced opacity.
        """
        grid_dims = (5, 5, 1) if self.periodicity == "2d" else (4, 4, 4)
        cell = self.cell
        if cell is not None:
            atoms_expanded = self.expand_structure(cell, grid_dims)
            for atom in atoms_expanded:
                sphere = pv.Sphere(radius=atom.radius, center=atom.position)
                self.plotter.add_mesh(sphere, color=atom.color, opacity=0.5, smooth_shading=True)

    def _change_representation(self):
        """
        Updates the visualized data based on the current representation toggle.
        If self.use_conventional is True, then the conventional cell data are used;
        otherwise the unit cell data are used.
        """
        if self.use_conventional:
            self.atoms = self.structure.atoms_cc
            self.bonds = self.structure.bonds_cc
            self.spins = self.structure.spins_cc  # Expect spins as a list of (Atom, Spin) tuples.
            self.cell = self.structure.conventional_cell
        else:
            self.atoms = self.structure.atoms_uc
            self.bonds = self.structure.bonds_uc
            self.spins = self.structure.spins_uc  # Expect spins as a list of (Atom, Spin) tuples.
            self.cell = self.structure.unit_cell
