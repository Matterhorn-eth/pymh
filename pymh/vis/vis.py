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
def plot(fn,
             inprefix='model_full', inpath=os.curdir,
             outprefix='model_ibc', outpath=os.curdir,
             ext='.su', isSU=True,
             clip=1e4, aspect='auto', cmap='seismic'):
    """ Create smaller model from full model """

    data = SEGYFile('/'.join([inpath, fn + ext]), isSU=isSU)

    data = data[:]

    plotopts = {
        'vmin': -clip,
        'vmax': +clip,
        'aspect': aspect,
        'cmap': cmap
    }

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    plt.tick_params(
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        left='off',        # ticks along the top edge are off
        right='off'        # ticks along the top edge are off
    )
    im = ax.imshow(data.T, **plotopts)
    ax.xaxis.set_label_text('Horizontal location [m]')
    ax.yaxis.set_label_text('Time [s]')
    ax.set_title('Shot gather')
    # plt.xticks(np.linspace(2e3,4e3,3))

    # Colorbar
    cbax = plt.colorbar(im, orientation='vertical', shrink=0.9)
    cbax.set_ticks(range(-int(clip), int(clip+1), int(clip/5)))
    cbax.set_label('Amplitude')

    plt.show()

    return ax
