import os
import sys
from datetime import datetime
import re

sys.path.insert(0, os.path.abspath('../..'))
from gtranslate import __author__, __version__, __title__, __maintainer__, __url__

# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = __title__
copyright = f'{datetime.now().year}, {__maintainer__}'
author = __author__

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

github_url = __url__


# The full version, including alpha/beta/rc tags
release = __version__
version = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxarg.ext', 'sphinx.ext.napoleon', 'sphinx.ext.autodoc', 'linuxdoc.rstFlatTable',
              'myst_parser', 'sphinx_sitemap', 'nbsphinx','matplotlib.sphinxext.plot_directive']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

html_extra_path = ['robots.txt']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'canonical_url': '',
    # 'analytics_id': 'UA-84847737-2',  # Provided by Google in your dashboard
    'logo_only': True,
    'version_selector': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_context = {
    'display_github': True,
    'github_user': 'cmc-aau',
    'github_repo': 'gTranslate',
    'github_version': 'master/docs/src/'
     }

html_css_files = [
    'theme_overrides.css',
]

html_js_files = [
    'custom.js'
]

html_logo = '_static/gTranslate_logo_red.svg'
html_favicon = '_static/gtranslate_favicon.png'

# Sitemap settings
html_baseurl = 'https://cmc-aau.github.io/gTranslate/'
sitemap_url_scheme = "{link}"


# Write dynamic install instructions with correct version
install_block_path = os.path.join(os.path.dirname(__file__), 'includes', 'install_block.rst')
manually_alias_models = os.path.join(os.path.dirname(__file__), 'includes', 'manually_alias_models.rst')
css_override_path = os.path.join(os.path.dirname(__file__), '_static', 'theme_overrides.css')
os.makedirs(os.path.dirname(install_block_path), exist_ok=True)

with open(install_block_path, 'w') as f:
    f.write(f"""\
.. code-block:: bash

    # using mamba (recommended)
    mamba create -n gtranslate-{release} -c conda-forge -c bioconda gtdb-gtranslate={release}

    # or using conda
    conda create -n gtranslate-{release} -c conda-forge -c bioconda gtdb-gtranslate={release}

    # activate the environment
    conda activate gtranslate-{release}
""")

with open(manually_alias_models, 'w') as f:
    f.write(f"""\
.. code-block:: bash

    # Activate the gTranslate conda environment
    conda activate gtranslate-{release}

    # Set the environment variable to the directory containing the gTranslate models
    conda env config vars set GTRANSLATE_MODEL_PATH="/path/to/unarchived/gtranslate/models/"

    # Re-activate the environment for the variable to take effect
    conda deactivate
    conda activate gtranslate-{release}

The directory must contain the following files (see :ref:`installing#gtranslate-models`):

* ``ada_multi_class.pkl.gz``
* ``dt_multi_class.pkl.gz``
* ``knn_multi_class.pkl.gz``
* ``mlp_multi_class.pkl.gz``
* ``xgb_multi_class.pkl.gz``
* ``label_encoder.pkl.gz``
""")

if os.path.exists(css_override_path):
    # Remove the version block injected by a previous build before injecting it again,
    # otherwise each local build appends another copy to this tracked file.
    css_marker = '\n/* Dynamically injected layout fix to display version when logo_only is True */'
    with open(css_override_path) as f:
        css_content = f.read().split(css_marker)[0]
    with open(css_override_path, 'w') as f:
        f.write(css_content)
        f.write(f"""\
{css_marker}
.wy-side-nav-search > a::after {{
    content: "v{release}";
    display: block;
    text-align: center;
    margin-top: 10px;
    font-size: 90%;
    font-weight: bold;
    color: #bf2000;
}}
""")


# Replace the version placeholder in all pages, including code blocks, with the version from gtranslate/__init__.py
VERSION_PLACEHOLDER = '{{ gtranslate_version }}'


def replace_version_placeholder(app, docname, source):
    source[0] = source[0].replace(VERSION_PLACEHOLDER, release)


def setup(app):
    app.connect('source-read', replace_version_placeholder)
