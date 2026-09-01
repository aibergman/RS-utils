import os
import sys
sys.path.insert(0, os.path.abspath('..'))

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

extensions = ['sphinx.ext.autodoc', 'myst_parser']
extensions.append('sphinx.ext.autodoc.typehints')


templates_path = ['_templates']
exclude_patterns = [
    '_build',
    'dev',
    'HT_BASELINE_INVENTORY.md',
    'HT_COMPATIBILITY_CONTRACT.md',
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
