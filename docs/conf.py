# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'astro-sedpy'
copyright = '2014-2022, Benjamin Johnson and Contributors'
author = 'Benjamin Johnson'

# The full version, including alpha/beta/rc tags
release = '0.3'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    "sphinx.ext.intersphinx",
    "myst_nb",
    'sphinx.ext.napoleon',
    'numpydoc'
]

myst_enable_extensions = ["dollarmath", "colon_fence"]

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autodoc_mock_imports = []
# The name of the Pygments (syntax highlighting) style to use.
#pygments_style = 'sphinx'
# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_title = "sedpy"
html_theme = 'sphinx_book_theme'
html_copy_source = True
html_show_sourcelink = True
html_theme_options = {"path_to_docs": "docs",
                      "repository_url": "https://github.com/bd-j/sedpy",
                      "repository_branch": "main",
                      "use_repository_button": True,
                      "use_edit_page_button": True,
                      "use_issues_button": True,
                      "use_download_button": True}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
#html_css_files = ['css/custom.css']
#html_logo = "_static/logo_name.png"
#html_favicon = "_static/favicon.png"

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
html_last_updated_fmt = ''
