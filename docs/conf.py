# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
import atexit
import re
import tempfile

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#sys.path.append(os.path.abspath('.'))

# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# -- Project information -----------------------------------------------------

project = 'KAPPA'
copyright = '2021, Lorenzo Campoli'
author = 'Lorenzo Campoli'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinxcontrib.spelling',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.restbuilder',
    'sphinx.ext.extlinks',
    'sphinx.ext.napoleon',
    'sphinx_tabs.tabs',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx' # None

intersphinx_mapping = {
     'python': ('https://docs.python.org/3', None),
     }

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'KAPPAdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

 # -- Spelling ---

if any("spell" in arg for arg in sys.argv):
    names_file = tempfile.NamedTemporaryFile(mode='w', prefix="coverage_names_", suffix=".txt")
    with open("../CONTRIBUTORS.txt") as contributors:
        names = set(re.split(r"[^\w']", contributors.read()))
        names = [n for n in names if len(n) >= 2 and n[0].isupper()]
        names_file.write("\n".join(names))
        names_file.flush()
    atexit.register(os.remove, names_file.name)

    spelling_word_list_filename = ['dict.txt', names_file.name]
    spelling_show_suggestions = False

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'KAPPA.tex', 'KAPPA Documentation',
     'Lorenzo Campoli', 'manual'),
]

extlinks = {
     # :github:`123` becomes a link to the GitHub issue, with text "issue 123".
     'github': ('https://github.com/lkampoli/kappa/issues/%s', 'issue '),
}

# Regexes for URLs that linkcheck should skip.
linkcheck_ignore = [
    # We have lots of links to GitHub, and they start refusing to serve them to linkcheck,
    # so don't bother checking them.
    r"https://github.com/lkampoli/kappa/(issues|pull)/\d+",
    # When publishing a new version, the docs will refer to the version before
    # the docs have been published.  So don't check those links.
    r"https://ml4sts.readthedocs.io/en/kappa-{}$".format(release),
]

# https://github.com/executablebooks/sphinx-tabs/pull/54
sphinx_tabs_valid_builders = ['linkcheck']

# When auto-doc'ing a class, only write the class' docstring into the class docs,
# don't automatically include the __init__ docstring.
autoclass_content = "class"

prerelease = bool(max(release).isalpha())

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'kappa', 'KAPPA Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'KAPPA', 'KAPPA Documentation',
     author, 'KAPPA', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True
