# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 11:11:57 2016

@author: bfilippo
"""

import os
import numpy as np
from pymh.param.parameters import \
    GridParam
from pymh.utils.segyread_new import \
    SEGYFile


# %%
def truncate(SmallGrid=GridParam,
             inprefix='model_full', inpath=os.curdir,
             outprefix='model_ibc', outpath=os.curdir,
             ext='.su', ispadding=False, isSU=True):
    """ Create smaller model from full model """

    for attr in ('vp', 'rho'):
        model = SEGYFile('/'.join([inpath, inprefix + '_' + attr + ext]), isSU=isSU)

        # the 1 assumes that the FD scheme is 2nd order accurate
        if ispadding:
            pad1 = 1
            pad2 = 2 * pad1
        else:
            pad1 = pad2 = 0

        origin = np.array(SmallGrid.parameters['origin'], dtype=np.float32)
        ncells = np.array(SmallGrid.parameters['number_of_cells'], dtype=np.int32)
        cell_size = np.array(SmallGrid.parameters['cell_size'], dtype=np.float32)

        iorigin = np.array(origin/cell_size - pad1, dtype=np.int32)
        iend = np.array(iorigin + ncells + pad2, dtype=np.int32)

        model_full = model[:]
        model_small = model_full[iorigin[0]:iend[0], iorigin[2]:iend[2]]

        trhead_small = [model.trhead[i] for i in range(iorigin[0], iend[0])]
        for i in range(len(trhead_small)):
            trhead_small[i]['ns'] = ncells[2] + pad2

        model.writeSU('/'.join([outpath, outprefix + '_' + attr + ext]),
                      model_small,
                      trhead=trhead_small)

    return iorigin, iend
