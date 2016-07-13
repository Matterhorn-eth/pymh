# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:32:12 2016

@author: bfilippo
"""
# %%
import sys
import subprocess
from segyread_new import SEGYFile
import matplotlib.pyplot as plt
import os

sys.dont_write_bytecode = True

# pymhfolder = '/w04d2/bfilippo/pymh'
pymhfolder = '/Users/filippo/work/pymh'
sys.path.insert(0, pymhfolder)

import pymh
from pymh.param.parameters import \
    GridParam, DecompositionParam, TimeParam, ModelParam, SimulationParam, \
    InputParam, OutputParam, BCParam, IBCParam, DirParam, Locations
from pymh.sim.simulations import \
    BasicSim
from pymh.utils.utils import \
    isPinRectangle


reload(pymh.param.parameters)

if __name__ == '__main__':

# %% ----------------------------------------------------------------------
# Create directory structure
# -------------------------------------------------------------------------
    Dir = DirParam()
    Dir.makedirs()

# %% ----------------------------------------------------------------------
# Initialize parameters
# -------------------------------------------------------------------------

    # Common
    Time = TimeParam()
    BC = BCParam('pml')

    # Full
    FullGrid = GridParam()
    FullDecomposition = DecompositionParam()
    FullModel = ModelParam(filename_prefix=['model_full'])
    FullSimulation = SimulationParam(freesurface=[True])
    FullSimulation.parameters.update(BC.parameters)
    FullInput = InputParam()
    FullOutput = []
    FullOutput.append(OutputParam('shot_gather',
                              receiver_origin=[0, 0, 40],
                              timestep_increment=[1]))
    FullOutput.append(OutputParam('slice'))

    # IBC
    IBC = IBCParam('freesurface')

    # Injection

# %% ----------------------------------------------------------------------
# Locations
# -------------------------------------------------------------------------

    SinjLocations = Locations()
    SinjLocations.rectangle(origin=IBC.extraparameters['inj']['origin'],
                            number_of_cells=IBC.extraparameters['inj']['ncells'],
                            cell_size=IBC.extraparameters['inj']['cell_size'])
    print SinjLocations.locations
    sinjlocationsfilename = 'sinj_locations.txt'
    SinjLocations.write(sinjlocationsfilename)

    SrecLocations = Locations()
    SrecLocations.rectangle(origin=IBC.extraparameters['rec']['origin'],
                            number_of_cells=IBC.extraparameters['rec']['ncells'],
                            cell_size=IBC.extraparameters['rec']['cell_size'])
    print SrecLocations.locations
    sreclocationsfilename = 'srec_locations.txt'
    SrecLocations.write(sreclocationsfilename)

# %% ----------------------------------------------------------------------
# Is source outside or inside?
# -------------------------------------------------------------------------

    corners = []
    for loc in SrecLocations.locations:
        if loc[3]//6 == 1: corners.append(loc[:3])
    corners = [corners[i] for i in [0, 1, 3, 2]]
    print corners

    IBC.extraparameters['source_inside'] = \
        isPinRectangle(corners, FullInput.parameters['location'])




# %% ----------------------------------------------------------------------
# Create files
# -------------------------------------------------------------------------

    if IBC.extraparameters['source_inside'] =
    FullOutput.append(OutputParam('shot_gather',
                              receiver_origin=[0, 0, 40],
                              timestep_increment=[1]))
    fullfilename = 'full.txt'

    Full = BasicSim(FullGrid,
                    FullDecomposition,
                    Time,
                    FullModel,
                    FullSimulation,
                    FullInput,
                    FullOutput)

    Full.create(fullfilename)

# %%
    # subprocess.call(['matterhorn', fullfilename])

# %%

    # shot = SEGYFile('shot.su', isSU=True)

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # im = ax.imshow(shot.readTraces().T, aspect='auto')
    # plt.show()
