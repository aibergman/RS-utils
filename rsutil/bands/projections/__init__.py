"""
Orbital and site projection analysis module.

This module provides tools for analyzing band structure projected onto
specific orbitals, sites, or spin components. Includes gauge fixing,
degeneracy handling, and projection-based visualization.

Classes
-------
ProjectedBandCalculator
    Calculate band structure with orbital/site/spin projections.
    Extends BandStructureCalculator with projection capabilities.

Functions
---------
guess_basis_order
    Automatically detect the ordering of basis functions.

degeneracy_groups
    Identify degenerate band manifolds.

gauge_fix_evecs
    Apply unitary transformations to stabilize phases.

stabilize_degenerate_projections
    Stabilize projection results in degenerate band regions.

atom_orbital_indices
    Get indices of specific atoms/orbitals in basis.

parse_atom_selection
    Parse user-friendly atom/orbital selection strings.

get_projection_vrange
    Determine appropriate value range for projection plots.

compute_quantity_from_sexp
    Parse and compute custom projection quantities (advanced).

Typical Usage
-------------
Basic projection analysis:

    >>> from rsutil.bands.projections import ProjectedBandCalculator
    >>> from rsutil.bands import load_lattice, load_hamiltonian
    >>> lattice = load_lattice('lattice.nml')
    >>> hamiltonian = load_hamiltonian('rs2txt.dat')
    >>> calc = ProjectedBandCalculator(lattice, hamiltonian)
    >>> kpts, klabs = calc.setup_kpath(npoints=100)[:2]
    >>> energies, projections = calc.calculate_projected_bands(kpts)
    >>> calc.plot_projected_bands(kdist, energies, projections, klabs)

Spin-resolved projections:

    >>> energies, spin_proj = calc.calculate_projected_bands(
    ...     kpts, projection_mode='spin'
    ... )
    
Specific atom/orbital selection:

    >>> # Plot only Fe d-orbital contributions
    >>> calc.calculate_projected_bands(
    ...     kpts, atom_selection='Fe', orbital_selection='d'
    ... )
"""

from .projections import ProjectedBandCalculator
from .projection_utils import (
    guess_basis_order,
    degeneracy_groups,
    gauge_fix_evecs,
    stabilize_degenerate_projections,
    atom_orbital_indices,
    parse_atom_selection,
    get_projection_vrange,
    compute_quantity_from_sexp,
)

__all__ = [
    'ProjectedBandCalculator',
    'guess_basis_order',
    'degeneracy_groups',
    'gauge_fix_evecs',
    'stabilize_degenerate_projections',
    'atom_orbital_indices',
    'parse_atom_selection',
    'get_projection_vrange',
    'compute_quantity_from_sexp',
]
