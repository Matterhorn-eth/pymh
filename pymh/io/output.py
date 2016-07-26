# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 11:11:57 2016

@author: bfilippo
"""

import os
import io
import numpy as np
import struct
from pymh.param.parameters import \
    GridParam, OutputParam, DirParam
#from pymh.core.classes import \
#    ShotGather #, Slice, VolumeBoundary
from pymh.utils.segyread_new import \
    SEGYFile


# %%
class BaseOutput(object):
    """ Base class for `Matterhorn` outputs. """

    type = None

    # def __init__(self):
    a = 1


class ShotGather(BaseOutput):
    """ Class for describing shot gathers in `Matterhorn`.

    """

    type = 'shotgather'

    def __init__(self, fn, endian='Little', isSU=True, ext='',
                 inpath=os.curdir):
        data = SEGYFile('/'.join([inpath, fn + ext]), isSU=isSU, endian=endian)
        self.data = data[:]
        del(data)


class VolumeBoundary(BaseOutput):
    """ Class for describing volume boundaries in `Matterhorn`.

    """

    type = 'volumeboundary'

    def __init__(self, fn, nt=1000, endian='Little', isSU=True, ext='',
                 inpath=os.curdir):
        self.nt = nt
        fid = io.open('injection_sxx_x_2500_y_0_z_10_volume_boundary', mode='rb')
        # fid = io.open('/'.join([inpath, fn + ext]), mode='rb')
        header = fid.read(48)
        headerstruct = '<12i'
        self.header = struct.unpack(headerstruct, header)
        data = fid.read(header[10]*nt*4)
        datastruct = '<%df' % (header[10]*nt)
        self.data = np.reshape(np.array(struct.unpack(datastruct, data), dtype=np.float32), [header[10], nt])
        fid.close()


# %%
def loadsnap(FullGrid=GridParam, SmallGrid=GridParam,
             SmallOutput=OutputParam, Dir=DirParam,
             inprefix='output', inpath=os.curdir,
             ext='.su', ispadding=False, isSU=True):

    """ Create smaller model from full model """

    start = SmallOutput.parameters['start_timestep'][0]
    step = SmallOutput.parameters['timestep_increment'][0]
    end = SmallOutput.parameters['end_timestep'][0]
    nsnaps = range(start, end ,step)

    # Full
    ncells_full = np.array(FullGrid.parameters['number_of_cells'], dtype=np.int32)
    cell_size_full = np.array(FullGrid.parameters['cell_size'], dtype=np.float32)

    # IBC
    ncells_ibc = np.array(SmallGrid.parameters['number_of_cells'], dtype=np.int32)
    origin_ibc = np.array(SmallGrid.parameters['origin'], dtype=np.float32)

    iorigin_ibc = np.array(origin_ibc/cell_size_full, dtype=np.int32)
    iend_ibc = np.array(iorigin_ibc + ncells_ibc, dtype=np.int32)

    # Initialize arrays
    full = np.zeros((ncells_full[0], ncells_full[2], len(nsnaps)))
    ibc = np.zeros(full.shape)
    diff = np.zeros(full.shape)

    for i, isnap in enumerate(nsnaps):
        full[:, :, i] = SEGYFile('/'.join([Dir.parameters['ref'], inprefix + '_{:0>8}'.format(isnap) + ext]), isSU=True, verbose=False, endian='Little')[:]
        ibc[iorigin_ibc[0]:iend_ibc[0], iorigin_ibc[2]:iend_ibc[2], i] = SEGYFile('/'.join([Dir.parameters['ibc'], inprefix + '_{:0>8}'.format(isnap) + ext]), isSU=True, verbose=False, endian='Little')[:]

    diff = full - ibc

    return full, ibc, diff