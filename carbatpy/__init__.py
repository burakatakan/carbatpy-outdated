# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:35:21 2022

@author: atakan
"""

import os
import sys

directories =("src", "src\\run_scripts", "src\\analysis", "data", "tests")
for d_name in directories:
    fpath = os.path.join(os.path.dirname(__file__), d_name)
    sys.path.append(fpath)
print(sys.path)