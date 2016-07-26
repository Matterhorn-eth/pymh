# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 22:36:51 2016

@author: filippo
"""
import sys
import numpy as np

sys.dont_write_bytecode = True

# pymhfolder = '/w04d2/bfilippo/pymh'
pymhfolder = '/Users/filippo/work/pymh'
sys.path.insert(0, pymhfolder)

import pymh
#from pymh.vis.vis import \
#    plot
from pymh.io.output import \
    ShotGather, VolumeBoundary

# reload(pymh.io.output)
# reload(pymh.core.classes)

# %%
# ax = plot('output', clip=1e3)

shot = ShotGather('output.su')

vb = VolumeBoundary('injection_sxx_x_2500_y_0_z_10_volume_boundary', nt=11000)

