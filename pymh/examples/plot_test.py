# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 22:36:51 2016

@author: filippo
"""
import sys
import numpy as np

# pymhfolder = '/w04d2/bfilippo/pymh'
pymhfolder = '/Users/filippo/work/pymh'
sys.path.insert(0, pymhfolder)

from pymh.vis.vis import \
    plot

# %%
ax = plot('output', clip=1e3)
