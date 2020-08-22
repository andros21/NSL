#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
------------------------------------------------------
        .___      _____             .__   __
      __| _/_____/ ____\____   __ __|  |_/  |_
     / __ |/ __ \   __\\__  \ |  |  \  |\   __\
    / /_/ \  ___/|  |   / __ \|  |  /  |_|  |
    \____ |\___  >__|  (____  /____/|____/__|
         \/    \/           \/

  Features:

    * load custom jupyter css      => [jupyter.css]
    * load custom matplotlib style => [theme.mplstyle]
    * change wd path               => [../out]

------------------------------------------------------
"""

import os
import numpy as np
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt
from IPython.core.display import display, HTML


with open('../../custom.css', 'r') as fs:
    css = fs.read()

html = f'<style>{css}</style>'
display(HTML(html))

plt.style.use("../../theme.mplstyle")

outPath = "../out/"
os.chdir(outPath)

# Default colors and vars for exercise 
df = dict()
yy = dict()
states = ["solid", "liquid", "gas"]
vv = ['P', 'V']
cc = cycle("rb")
