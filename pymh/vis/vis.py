# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 11:11:57 2016

@author: bfilippo
"""
from __future__ import absolute_import
import os
import numpy as np
# from pymh.param.parameters import \
#     GridParam
from pymh.utils.segyread_new import \
    SEGYFile
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import gridspec
from mpldatacursor import datacursor

__all__ = ['plot', 'animate']


# %%
def plot(fn,
         inprefix='model_full', inpath=os.curdir,
         outprefix='model_ibc', outpath=os.curdir,
         ext='.su', isSU=True, endian='Little',
         clip=1e3, aspect='auto', cmap='seismic',
         interpolation='bicubic', title=None,
         colorbar=True, style='shot_gather',
		 extent=[0, 1, 1, 0],
		 figsize=(10, 8)):

    """ Create a plot """

    if isinstance(fn, basestring):
        data = SEGYFile('/'.join([inpath, fn + ext]), isSU=isSU, endian=endian)
        data = data[:]
    elif isinstance(fn, np.ndarray):
        data = fn
    else:
        raise TypeError('LOL!')

    fig = plt.figure(figsize=figsize, facecolor='w', edgecolor='k')
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
    if style == 'shot_gather':
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
    elif style == 'sub_volume_boundary':
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
        #ax.xaxis.set_label_text('Horizontal location [m]')
        ax.xaxis.set_label_text('Inline [m]')
        #ax.yaxis.set_label_text('Crossline [m]')
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
        'interpolation': interpolation,
		'extent': extent
    }

    # Create image
    im = ax.imshow(data.T, picker=True, **plotopts)

    # Colorbar
    if colorbar:
        cbax = plt.colorbar(im, orientation='vertical', shrink=0.9)
        cbax.set_ticks(range(int(clipmin+1), int(clipmax+1), int((clipmax-clipmin)/5)))
        cbax.set_label(cb_label)


    datacursor(display='single')

    fig.canvas.mpl_connect('pick_event', onclick)
	#fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

    return ax

def onclick(event):
    #if event.xdata != None and event.ydata != None:
    #print(event.xdata, event.ydata)
    print((event.mouseevent.xdata))

def animate(full, ibc, diff,
            clip=1e3, aspect='auto', cmap='seismic', blit=False,
            interval=10, interpolation='bilinear', title=None,
            colorbar=True, style='shot_gather', extent=[0, 1, 1, 0],
			figsize=(10, 8)):

    """ Create a plot """
    pause = False

#    def onClick(event):
#        global pause
#        pause ^= True

    def isnap():
        t_max = full.shape[2]-1
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
#            ax0.title.set_text('Full (snapshot #{0:4})'.format(i))
            return im0, im1, im2,

    lol = 1
    fig = plt.figure(figsize=figsize, facecolor='w', edgecolor='k')

    # lim = [-clip, clip]

    # Plot options
    plotopts = {
        'vmin': -clip,
        'vmax': +clip,
        'aspect': aspect,
        'cmap': cmap,
        'interpolation': interpolation,
        'animated': True,
        'extent': extent
    }

    gs = gridspec.GridSpec(3, 1)
    ax0 = fig.add_subplot(gs[0, 0])
#    title0 = ax0.title.set_text('Full (snapshot #{0:4})'.format(0))
    # Remove the ugly ticks
    plt.tick_params(
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        left='off',        # ticks along the top edge are off
        right='off'        # ticks along the top edge are off
    )
    ax1 = fig.add_subplot(gs[1, 0])
#    title1 = ax1.title.set_text('IBC')
    # Remove the ugly ticks
    plt.tick_params(
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        left='off',        # ticks along the top edge are off
        right='off'        # ticks along the top edge are off
    )
    ax2 = fig.add_subplot(gs[2, 0])
#    title2 = ax2.title.set_text('Full - IBC')
    # Remove the ugly ticks
    plt.tick_params(
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        left='off',        # ticks along the top edge are off
        right='off'        # ticks along the top edge are off
    )
    i = 0
    im0 = ax0.imshow(full[:, :, i].T, **plotopts)
    im1 = ax1.imshow(ibc[:, :, i].T, **plotopts)
    im2 = ax2.imshow(diff[:, :, i].T, **plotopts)

#    ax0.xaxis.set_label_text('Horizontal location [m]')
    ax0.yaxis.set_label_text('Depth [m]')
#    ax1.xaxis.set_label_text('Horizontal location [m]')
    ax0.get_xaxis().set_ticks([])
    ax1.yaxis.set_label_text('Depth [m]')
    ax1.get_xaxis().set_ticks([])
    ax2.xaxis.set_label_text('Horizontal location [m]')
    ax2.yaxis.set_label_text('Depth [m]')

    # Annotate subplots
    ax0.annotate('Full', xy=(200, 500), xycoords='data',
                size=15,
                bbox=dict(boxstyle='round,pad=.3', fc='white', ec='black')
                )
    ax1.annotate('IBC', xy=(200, 500), xycoords='data',
                size=15,
                bbox=dict(boxstyle='round,pad=.3', fc='white', ec='black')
                )
    ax2.annotate('Full - IBC', xy=(200, 500), xycoords='data',
                size=15,
                bbox=dict(boxstyle='round,pad=.3', fc='white', ec='black')
                )

    fig.tight_layout()

#    fig.canvas.mpl_connect('button_press_event', onClick)
    ani = animation.FuncAnimation(fig, snapshots, isnap,
                                  fargs=(lol, full, ibc, diff),
                                  blit=blit, interval=interval, repeat=False)

    metadata = dict(title='Movie Test', artist='EEG',
                    comment='Testing!')
    writer = animation.writers['ffmpeg'](fps=10, bitrate=8000, metadata=metadata)
    ani.save('demo_8000k_200c.avi', writer=writer, dpi=200) #, extra_args=['-b:v', '4000k', '-bufsize', '4000k'])
#    with animation.MovieWriter.saving('myfile.mp4'):
#        animation.MovieWriter.grab_frame()


    plt.show()

    return ani
