# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 11:11:57 2016

@author: bfilippo
"""

from __future__ import absolute_import
import os
import io
import numpy as np
import struct
from pymh.param.parameters import \
    GridParam, OutputParam, DirParam
import pymh.vis.vis as vis
from pymh.utils.segyread_new import \
    SEGYFile

__all__ = ['ShotGather', 'Slice', 'SubVolumeBoundary', 'loadsnap', 'loadsnap3']


# %%
class BaseOutput(object):
    """ Base class for `Matterhorn` outputs. """

    type = None

    # def __init__(self):
    def plot(self, **kwargs):
        vis.plot(self.data, style=self.type, **kwargs)


class ShotGather(BaseOutput):
    """ Class for describing shot gathers in `Matterhorn`.

    """

    type = 'shot_gather'

    def __init__(self, fn, endian='Little', isSU=True, ext='',
                 inpath=os.curdir):
        data = SEGYFile('/'.join([inpath, fn + ext]), isSU=isSU, endian=endian)
        self.data = data[:]
        del(data)


class Slice(BaseOutput):
    """ Class for describing slices in `Matterhorn`.

    """

    type = 'slice'

    def __init__(self, fn, endian='Little', isSU=True, ext='',
                 inpath=os.curdir):
        data = SEGYFile('/'.join([inpath, fn + ext]), isSU=isSU, endian=endian)
        self.data = data[:]
        del(data)

class SubVolumeBoundary(BaseOutput):
    """ Class for describing sub volume boundaries in `Matterhorn`.

    """

    type = 'sub_volume_boundary'

    def __init__(self, fn, nt=1000, endian='Little', isSU=True, ext='',
                 inpath=os.curdir):
        self.nt = nt
        fid = io.open('/'.join([inpath, fn + ext]), mode='rb')
        header = fid.read(48)
        headerstruct = '<12i'
        self.header = struct.unpack(headerstruct, header)
        data = fid.read(self.header[10]*nt*4)
        datastruct = '<%df' % (self.header[10]*nt)
        self.data = np.reshape(np.array(struct.unpack(datastruct, data), dtype=np.float32),
                               [nt, self.header[10]]).T
        fid.close()


# %%
def loadsnap(FullGrid=GridParam, FullOutput=OutputParam, Dir=DirParam,
             inprefix='output', inpath=os.curdir,
             ext='.su', ispadding=False, isSU=True, nsnaps=None):

    """ Load full snapshots """

    start = FullOutput.parameters['start_timestep'][0]
    step = FullOutput.parameters['timestep_increment'][0]
    end = FullOutput.parameters['end_timestep'][0]
    if not nsnaps:
        nsnaps = range(start, end, step)

    # Full
    ncells_full = np.array(FullGrid.parameters['number_of_cells'], dtype=np.int32)

    # Initialize arrays
    full = np.zeros((ncells_full[0], ncells_full[2], len(nsnaps)))

    for i, isnap in enumerate(nsnaps):
        full[:, :, i] = SEGYFile('/'.join([Dir.parameters['ref'], inprefix + '_{:0>8}'.format(isnap) + ext]), isSU=True, verbose=False, endian='Little')[:]

    return full

# %%
def loadsnap3(FullGrid=GridParam, SmallGrid=GridParam,
             SmallOutput=OutputParam, Dir=DirParam,
             inprefix='output', inpath=os.curdir, small_type='ibc',
             ext='.su', ispadding=False, isSU=True, nsnaps=None):

    """ Load full and IBC/Injection snapshots """

    start = SmallOutput.parameters['start_timestep'][0]
    step = SmallOutput.parameters['timestep_increment'][0]
    end = SmallOutput.parameters['end_timestep'][0]
    if not nsnaps:
        nsnaps = range(start, end, step)

    # Full
    ncells_full = np.array(FullGrid.parameters['number_of_cells'], dtype=np.int32)
    cell_size_full = np.array(FullGrid.parameters['cell_size'], dtype=np.float32)

    # Small
    ncells_small = np.array(SmallGrid.parameters['number_of_cells'], dtype=np.int32)
    origin_small = np.array(SmallGrid.parameters['origin'], dtype=np.float32)

    iorigin_small = np.array(origin_small/cell_size_full, dtype=np.int32)
    iend_small = np.array(iorigin_small + ncells_small, dtype=np.int32)

    # Initialize arrays
    full = np.zeros((ncells_full[0], ncells_full[2], len(nsnaps)))
    small = np.zeros(full.shape)
    diff = np.zeros(full.shape)

    for i, isnap in enumerate(nsnaps):
        full[:, :, i] = SEGYFile('/'.join([Dir.parameters['ref'], inprefix + '_{:0>8}'.format(isnap) + ext]), isSU=True, verbose=False, endian='Little')[:]
        small[iorigin_small[0]:iend_small[0], iorigin_small[2]:iend_small[2], i] = SEGYFile('/'.join([Dir.parameters[small_type], inprefix + '_{:0>8}'.format(isnap) + ext]), isSU=True, verbose=False, endian='Little')[:]

    diff = full - small

    return full, small, diff