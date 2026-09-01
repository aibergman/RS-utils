import os
import sys
import types
sys.path.insert(0, os.path.abspath('..'))

# Some macOS environments combine an older Matplotlib with a newer NumPy.  A
# broken optional 3-D plotting import should not prevent the textual API
# reference from being built.  The real plotting dependency is still required
# when users execute the visualization functions.
mplot3d = types.ModuleType('mpl_toolkits.mplot3d')
mplot3d.Axes3D = type('Axes3D', (), {})
sys.modules['mpl_toolkits.mplot3d'] = mplot3d

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RS-utils'
copyright = '2025, Anders Bergman'
author = 'Anders Bergman'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'myst_parser']
extensions.append('sphinx.ext.autodoc.typehints')

# Keep the API reference buildable in a lean documentation environment.  These
# dependencies are optional or platform-specific at import time, so autodoc
# can mock them while retaining public signatures and docstrings.
autodoc_mock_imports = [
    'f90nml', 'pymatgen', 'spglib', 'seekpath', 'matplotlib', 'scipy',
    'ase', 'pyvista', 'yaml',
]


templates_path = ['_templates']
exclude_patterns = [
    '_build',
    'dev',
    'HT_BASELINE_INVENTORY.md',
    'HT_FOUNDATION_BLUEPRINT.md',
    'HT_REMEDIATION_BLUEPRINT.md',
    'notes',
    'Thumbs.db',
    '.DS_Store',
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,  # Expand the navigation by default
}
