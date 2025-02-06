from structure import Structure
from glyphs import Atom
import numpy as np
import pyvista as pv


# ----------------------------
# 4. Visualizer Class
# ----------------------------
class Visualizer:
    """
        Provides interactive toggles for visualizing different aspects of a 
        structure, including bonds, spins, cell type, extents, and periodicity.

        Attributes:
            structure (Structure): The structure to visualize.
            plotter (pv.Plotter): The pyvista plotter instance.
            show_bonds (bool): Toggle for showing bonds.
            show_spins (bool): Toggle for showing spins.
            use_conventional (bool): Toggle for using conventional cell.
            show_extents (bool): Toggle for showing extents.
            periodicity (str or None): Periodicity type ('2d', '3d', or None).

        Methods:
            toggle_bonds(): Toggles the display of bonds.
            toggle_spins(): Toggles the display of spins.
            toggle_cell(): Toggles between unit cell and conventional cell.
            toggle_extents(): Toggles the display of cell extents.
            toggle_periodicity(period="2d"): Toggles periodicity type.
            get_cell(): Returns the lattice cell matrix based on the toggle.
            expand_structure(cell, grid_dims): Expands the structure by 
                replicating the cell.
            render(): Clears and re-renders the scene with current toggles.
            _render_atoms(): Renders atoms in the structure.
            _render_bonds(): Renders bonds in the structure.
            _render_spins(): Renders spins in the structure.
            _render_extents(): Renders the extents of the unit or conventional 
                cell.
            _render_periodic_images(): Renders periodic images of the structure.
    """
    def __init__(self, structure: Structure):
        """
        Visualizer for a Structure using pyvista.

        Provides interactive toggles:
          - bonds / no bonds
          - unit cell / conventional cell
          - spins / no spins
          - extents / no extents
          - 2d-periodic / non-periodic
          - 3d-periodic / non-periodic
        """
        self.structure = structure
        self.plotter = pv.Plotter()

        # Toggle options (default settings)
        self.show_bonds = False
        self.show_spins = False
        self.use_conventional = True  # False means use unit cell
        self.show_extents = True 
        self.periodicity = None  # Options: None, '2d', '3d'

    def toggle_bonds(self):
        self.show_bonds = not self.show_bonds
        self.render()

    def toggle_spins(self):
        self.show_spins = not self.show_spins
        self.render()

    def toggle_cell(self):
        self.use_conventional = not self.use_conventional
        self.render()

    def toggle_extents(self):
        self.show_extents = not self.show_extents
        self.render()

    def toggle_periodicity(self, period="2d"):
        """
        Toggles periodicity.
        period: '2d' or '3d'; if already enabled with this type, then disables periodic images.
        """
        if self.periodicity == period:
            self.periodicity = None
        else:
            self.periodicity = period
        self.render()

    def get_cell(self):
        """
        Returns the lattice cell matrix based on the toggle.
        """
        if self.use_conventional and self.structure.conventional_cell is not None:
            return self.structure.conventional_cell
        elif self.structure.unit_cell is not None:
            return self.structure.unit_cell
        return None

    def expand_structure(self, cell, grid_dims):
        """
        Expands the structure by replicating the cell.

        Args:
            cell: 3x3 lattice matrix.
            grid_dims (tuple): (nx, ny, nz) number of replications.

        Returns:
            A list of new Atom objects corresponding to the periodic images.
        """
        atoms_expanded = []
        for i in range(grid_dims[0]):
            for j in range(grid_dims[1]):
                for k in range(grid_dims[2]):
                    translation = np.dot(cell.T, [i, j, k])
                    for atom in self.structure.atoms:
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
        """
        Clears and re-renders the scene with current toggles.
        """
        self.plotter.clear()

        # Render components based on toggles
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
        Renders atoms in the structure.
        """
        for atom in self.structure.atoms:
            sphere = pv.Sphere(radius=atom.radius, center=atom.position)
            self.plotter.add_mesh(
                sphere,
                color=atom.color,
                smooth_shading=True,
                specular=0.5,
                specular_power=10,
            )

    def _render_bonds(self):
        """
        Renders bonds in the structure.
        """
        bonds = self.structure.get_bond_objects()
        for bond in bonds:
            vector = bond.end - bond.start
            length = np.linalg.norm(vector)
            center = (bond.start + bond.end) / 2
            direction = vector / length if length else [0, 0, 1]
            cylinder = pv.Cylinder(
                center=center,
                direction=direction,
                radius=bond.radius,
                height=length,
            )
            self.plotter.add_mesh(cylinder, color=bond.color, smooth_shading=True)

    def _render_spins(self):
        """
        Renders spins in the structure.
        """
        for atom, spin in self.structure.get_atom_spin_objects():
            if spin.magnitude > 0:
                arrow = pv.Arrow(
                    start=atom.position,
                    direction=spin.direction,
                    tip_length=0.2 * spin.magnitude,
                    tip_radius=spin.tip_radius,
                    shaft_radius=spin.base_radius,
                    scale=spin.magnitude,
                )
                self.plotter.add_mesh(arrow, color=spin.color, smooth_shading=True)

    def _render_extents(self):
        """
        Renders the extents of the unit or conventional cell.
        """
        cell = self.get_cell()
        if cell is not None:
            origin = np.array([0, 0, 0])
            a, b, c = cell[0], cell[1], cell[2]
            corners = np.array([origin, a, b, c, a + b, a + c, b + c, a + b + c])
            xmin, ymin, zmin = corners.min(axis=0)
            xmax, ymax, zmax = corners.max(axis=0)
            box = pv.Box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))
            self.plotter.add_mesh(
                box, style="wireframe", color="black", smooth_shading=True
            )

    def _render_periodic_images(self):
        """
        Renders periodic images of the structure.
        """
        grid_dims = (5, 5, 1) if self.periodicity == "2d" else (4, 4, 4)
        cell = self.get_cell()
        if cell is not None:
            atoms_expanded = self.expand_structure(cell, grid_dims)
            for atom in atoms_expanded:
                sphere = pv.Sphere(radius=atom.radius, center=atom.position)
                self.plotter.add_mesh(
                    sphere, color=atom.color, opacity=0.5, smooth_shading=True
                )
