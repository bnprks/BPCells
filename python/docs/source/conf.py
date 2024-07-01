import inspect
import os
import sys
import warnings

import bpcells

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'bpcells'
copyright = '2023, Ben Parks'
author = 'Ben Parks'

# Determine what version of docs this is for the version history dropdown
# version = bpcells.version.version
# if ".dev" in version:
#     version = "dev"
# else:
#     version = f"v{version}"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.linkcode",
    "sphinx.ext.napoleon",
    "myst_nb",
]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_context = {
    "default_mode": "light"
}

html_sidebars = {
    "**": ["search-field", "sidebar-nav-bs"],
}

# -- Options for PyData Theme -------------------------------------------------
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#references
html_theme_options = {
    "github_url": "https://github.com/bnprks/BPCells",
    # "switcher": {
    #     # "json_url": "https://bnprks.github.io/BPCells/python/_static/switcher.json",
    #     "json_url": "_static/switcher.json",
    #     "version_match": version,
    # },
    "navbar_end": ["theme-switcher", "navbar-icon-links"] #, "version-switcher"]
}


# -- Options for autodoc ----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#configuration

autosummary_generate = True

# -- Options for myst_nb ----------------------------------------------------
# https://myst-nb.readthedocs.io/en/latest/configuration.html
nb_custom_formats = {
    ".md": ["jupytext.reads", {"fmt": "mystnb"}],
}

# -- Options for myst_parser ----------------------------------------------------
# https://myst-parser.readthedocs.io/en/latest/configuration.html
myst_enable_extensions = ["colon_fence"]

# -- Options for intersphinx ----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration
intersphinx_mapping = {
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/', None),
    'pandas': ('http://pandas.pydata.org/docs/', None),
}

# based on numpy/pandas doc/source/conf.py
def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
    """
    if domain != "py":
        return None

    modname = info["module"]
    fullname = info["fullname"]

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split("."):
        try:
            with warnings.catch_warnings():
                # Accessing deprecated objects will generate noisy warnings
                warnings.simplefilter("ignore", FutureWarning)
                obj = getattr(obj, part)
        except AttributeError:
            return None

    try:
        fn = inspect.getsourcefile(inspect.unwrap(obj))
    except TypeError:
        try:  # property
            fn = inspect.getsourcefile(inspect.unwrap(obj.fget))
        except (AttributeError, TypeError):
            fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except TypeError:
        try:  # property
            source, lineno = inspect.getsourcelines(obj.fget)
        except (AttributeError, TypeError):
            lineno = None
    except OSError:
        lineno = None

    if lineno:
        linespec = f"#L{lineno}-L{lineno + len(source) - 1}"
    else:
        linespec = ""

    fn = os.path.relpath(fn, start=os.path.dirname(bpcells.__file__))

    return f"https://github.com/bnprks/bpcells/blob/main/python/src/bpcells/{fn}{linespec}"

# See: https://www.sphinx-doc.org/en/master/development/html_themes/index.html#defining-custom-template-functions
# To Do: Might be able to get the top header section nav working better  
# def setup(app):
#     app.connect("")