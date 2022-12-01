# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import os
import sys
import pathlib
sys.path.insert(0,os.path.abspath('.'))
sys.path.insert(0,os.path.abspath('../..'))
sys.path.insert(0,os.path.abspath('../carbatpy'))

#sys.path.insert(0,os.path.abspath(r'C:\Users\atakan\sciebo\Python\carbatpy\carbatpy'))
#sys.path.insert(0,os.path.abspath(r'C:\Users\atakan\sciebo\Python\carbatpy'))
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
project = 'carbatpy'
copyright = '2022, B. Atakan'
author = 'B. Atakan'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon', 
              'sphinx.ext.doctest',
              'sphinx.ext.githubpages',
              'sphinx.ext.duration']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

numpydoc_show_classmembers = False
autodoc_mock_imports = ['ctREFPROP', 'CoolProp']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'nature'  # 'alabaster' 'classic'  
# html_theme_options = {
#     "rightsidebar": "true",
#     "relbarbgcolor": "black"
#     }
#html_static_path = ['_static']
