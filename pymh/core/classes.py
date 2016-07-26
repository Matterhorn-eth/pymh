# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 11:11:57 2016

@author: bfilippo
"""

import os
import io
import numpy as np
import struct
# from pymh.param.parameters import \
#     GridParam
from pymh.utils.segyread_new import \
    SEGYFile
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import gridspec


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
        fid = io.open('/'.join([inpath, fn + ext]), mode='r')
        header = fid.read(48)
        headerstruct = '<12i'
        self.header = struct.unpack(headerstruct, header)
        data = fid.read(header[10]*nt*4)
        datastruct = '<%df' % (header[10]*nt)
        self.data = np.reshape(np.array(struct.unpack(datastruct, data), dtype=np.float32), [header[10], nt])
        fid.close()

