# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 15:35:15 2023

Use self-made modules from wam_ipe_plotter
"""

from wam_ipe_plotter import save_fig
import sys

filenames = sys.argv[1:]

for filename in filenames:
    save_fig(filename)

