# Sphinx configuration for rslmto band package
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'rslmto-band'
author = 'rslmto'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Use ReadTheDocs theme when available
try:
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
except Exception:
    html_theme = 'alabaster'

html_static_path = ['_static']
