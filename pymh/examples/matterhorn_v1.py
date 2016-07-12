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

pymhfolder = '/w04d2/bfilippo/pymh'
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

    # -------------------------------------------------------------------------
    # Create directory structure
    # -------------------------------------------------------------------------
    Dir = DirParam()
    Dir.makedirs()

    # -------------------------------------------------------------------------
    # Full
    # -------------------------------------------------------------------------

    # Initialize parameters
    Grid = GridParam()
    Decomposition = DecompositionParam()
    Time = TimeParam()
    Model = ModelParam(filename_prefix=['model_full'])
    Simulation = SimulationParam()
    BC = BCParam('pml')
    Simulation.parameters.update(BC.parameters)
    Input = []
    Input.append(InputParam())
    Output = []
    Output.append(OutputParam('shot_gather',
                              receiver_origin=[0, 0, 40],
                              timestep_increment=[1]))
    Output.append(OutputParam('slice'))

    fullfilename = 'full.txt'

    Full = BasicSim(Grid,
                    Decomposition,
                    Time,
                    Model,
                    Simulation,
                    Input,
                    Output)

    Full.create(fullfilename)

# %% Locations
    IBC = IBCParam('freesurface')
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

# %%
    corners = []
    for loc in SrecLocations.locations:
        if loc[3] in range(6,10): corners.append(loc[:3])
    corners = corners[:2] + corners[3:1:-1]
    print corners
    
# %%    

    a = isPinRectangle(corners, [500, 0, 1400])
    print a

# %%
    # subprocess.call(['matterhorn', fullfilename])

# %%

    # shot = SEGYFile('shot.su', isSU=True)

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # im = ax.imshow(shot.readTraces().T, aspect='auto')
    # plt.show()
