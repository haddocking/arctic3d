# -*- coding: utf-8 -*-
"""Sphinx config file."""
# from __future__ import unicode_literals

import os
import sys
from pathlib import Path

sys.path.append(Path(".."))

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.ifconfig",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "myst_parser",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

master_doc = "index"
project = "arctic3d"
year = "2022"
author = "BonvinLab"
copyright = "{0}, {1}".format(year, author)
version = release = "0.0.0"

todo_include_todos = True
pygments_style = "trac"
templates_path = ["."]
extlinks = {
    "issue": ("https://github.com/haddocking/arctic3d/issues/%s", "#"),
    "pr": ("https://github.com/haddocking/arctic3d/pull/%s", "PR #"),
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get("READTHEDOCS", None) == "True"

linkcheck_ignore = [r"https://codecov.io/*"]

# html_logo = 'img/taurenmd_logo_black.png'
html_use_smartypants = True
html_last_updated_fmt = "%b %d, %Y"
html_split_index = False
html_sidebars = {
    "**": ["searchbox.html", "globaltoc.html", "sourcelink.html"],
}
html_short_title = "%s-%s" % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

# myst options
myst_heading_anchors = 3

autosummary_generate = True
