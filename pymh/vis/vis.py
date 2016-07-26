# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 11:11:57 2016

@author: bfilippo
"""

import os
import numpy as np
# from pymh.param.parameters import \
#     GridParam
from pymh.utils.segyread_new import \
    SEGYFile
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import gridspec


# %%
def plot(fn,
         inprefix='model_full', inpath=os.curdir,
         outprefix='model_ibc', outpath=os.curdir,
         ext='.su', isSU=True, endian='Little',
         clip=1e4, aspect='auto', cmap='seismic',
         interpolation='bicubic', title=None,
         colorbar=True, style='gather'):

    """ Create a plot """

    if isinstance(fn, basestring):
        data = SEGYFile('/'.join([inpath, fn + ext]), isSU=isSU, endian=endian)
        data = data[:]
    elif isinstance(fn, np.ndarray):
        data = fn
    else:
        raise TypeError('LOL!')

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')

    # Remove the ugly ticks
    plt.tick_params(
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        left='off',        # ticks along the top edge are off
        right='off'        # ticks along the top edge are off
    )

    # Set style-specific attributes
    if style == 'shotgather':
        clipmin = -clip
        clipmax = +clip
        ax.xaxis.set_label_text('Horizontal location [m]')
        ax.yaxis.set_label_text('Time [s]')
        if title is None:
            ax.set_title('Shot gather')
        else:
            ax.set_title(title)
        cb_label = 'Amplitude'
    elif style == 'slice':
        clipmin = -clip
        clipmax = +clip
        ax.xaxis.set_label_text('Horizontal location [m]')
        ax.yaxis.set_label_text('Depth [m]')
        if title is None:
            ax.set_title('Slice')
        else:
            ax.set_title(title)
        cb_label = 'Amplitude'
    elif style == 'volumeboundary':
        clipmin = -clip
        clipmax = +clip
        ax.xaxis.set_label_text('Horizontal location [m]')
        ax.yaxis.set_label_text('Time [s]')
        if title is None:
            ax.set_title('Volume boundary')
        else:
            ax.set_title(title)
        cb_label = 'Amplitude'
    elif style == 'velocity':
        clipmin = data.min()
        clipmax = data.max()
        ax.xaxis.set_label_text('Horizontal location [m]')
        ax.yaxis.set_label_text('Depth [m]')
        if title is None:
            ax.set_title('Velocity')
        else:
            ax.set_title(title)
        cb_label = 'Amplitude [m/s]'
    elif style == 'density':
        clipmin = data.min()
        clipmax = data.max()
        ax.xaxis.set_label_text('Horizontal location [m]')
        ax.yaxis.set_label_text('Depth [m]')
        if title is None:
            ax.set_title('Density')
        else:
            ax.set_title(title)
        cb_label = 'Amplitude [kg/m^3]'
    # plt.xticks(np.linspace(2e3,4e3,3))

    # Plot options
    plotopts = {
        'vmin': clipmin,
        'vmax': clipmax,
        'aspect': aspect,
        'cmap': cmap,
        'interpolation': interpolation
    }

    # Create image
    im = ax.imshow(data.T, **plotopts)

    # Colorbar
    if colorbar:
        cbax = plt.colorbar(im, orientation='vertical', shrink=0.9)
        cbax.set_ticks(range(int(clipmin+1), int(clipmax+1), int((clipmax-clipmin)/5)))
        cbax.set_label(cb_label)

    plt.show()

    return ax


def animate(full, ibc, diff,
            clip=1e3, aspect='auto', cmap='seismic', blit=False,
            interval=10, interpolation='bilinear', title=None,
            colorbar=True, style='gather'):

    """ Create a plot """
    pause = False

#    def onClick(event):
#        global pause
#        pause ^= True

    def isnap():
        t_max = full.shape[2]
        dt = 1
        t = 1
        while t < t_max:
            if not pause:
                t = t + dt
            yield t

    def snapshots(i, lol, full, ibc, diff):
        # return ibc[:, :, i]
        if not pause:
            im0.set_array(full[:, :, i].T)
            im1.set_array(ibc[:, :, i].T)
            im2.set_array(diff[:, :, i].T)
            # ax = fig.colorbar(im)
            title.set_text('{0:4}'.format(i))
            return im0, im1, im2,

    lol = 1
    fig = plt.figure()

    # lim = [-clip, clip]

    # Plot options
    plotopts = {
        'vmin': -clip,
        'vmax': +clip,
        'aspect': aspect,
        'cmap': cmap,
        'interpolation': interpolation,
        'animated': True
    }

    gs = gridspec.GridSpec(3, 1)
    ax0 = fig.add_subplot(gs[0, 0])
    title = plt.title('{0:4}'.format(0))
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[2, 0])
    i = 0
    im0 = ax0.imshow(full[:, :, i].T, **plotopts)
    im1 = ax1.imshow(ibc[:, :, i].T, **plotopts)
    im2 = ax2.imshow(diff[:, :, i].T, **plotopts)

#    fig.canvas.mpl_connect('button_press_event', onClick)
    ani = animation.FuncAnimation(fig, snapshots, isnap,
                                  fargs=(lol, full, ibc, diff),
                                  blit=blit, interval=interval, repeat=True)
    plt.show()

    return ani
