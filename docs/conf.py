# Configuration file for the Sphinx documentation builder.
#
# For the full list of options see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Make the bobaT package importable for autodoc (repo root is one level up).
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "BoBa-T"
author = "Sarah Groves"
copyright = "2026, Sarah Groves"

# The full version, including alpha/beta/rc tags.
release = "0.1.4"
version = "0.1"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",       # pull docstrings from the source
    "sphinx.ext.autosummary",   # generate summary tables
    "sphinx.ext.napoleon",      # Google / NumPy style docstrings
    "sphinx.ext.viewcode",      # add "[source]" links
    "sphinx.ext.intersphinx",   # link to other projects' docs
    "sphinx.ext.mathjax",       # render math
    "nbsphinx",                 # render Jupyter notebooks as pages
    "nbsphinx_link",            # link to notebooks that live outside docs/
]

# Modules that are painful to install (compiled / C-extension dependencies).
# autodoc imports the package to read signatures, so we mock these instead of
# requiring them on the docs-build machine.
autodoc_mock_imports = [
    "graph_tool",
    "umap",
    "magic",
    "leidenalg",
    "numba",
]

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autosummary_generate = True

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False

# Notebooks already contain saved output; do not re-execute them at build time
# (they depend on graph-tool and example data that are not present in CI).
nbsphinx_execute = "never"
nbsphinx_allow_errors = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints", "Thumbs.db", ".DS_Store"]

master_doc = "index"

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}

# -- HTML output -------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "titles_only": False,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "BoBa-T documentation"
