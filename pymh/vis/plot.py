# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 11:11:57 2016

@author: bfilippo
"""

import os
import numpy as np
from pymh.param.parameters import \
    GridParam
from pymh.utils.segyread_new import \
    SEGYFile
import matplotlib.pyplot as plt
from matplotlib import gridspec

# %%
def plot(SmallGrid=GridParam,
             inprefix='model_full', inpath=os.curdir,
             outprefix='model_ibc', outpath=os.curdir,
             ext='.su', ispadding=False, isSU=True):
    """ Create smaller model from full model """

    data = SEGYFile('/'.join([inpath, inprefix + ext]), isSU=isSU)

    data = data[:]

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(data.T, aspect='auto')

    plt.show()

    return 1
