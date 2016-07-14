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
# pymhfolder = '/Users/filippo/work/pymh'
sys.path.insert(0, pymhfolder)

import pymh
from pymh.param.parameters import \
    GridParam, DecompositionParam, TimeParam, ModelParam, SimulationParam, \
    InputParam, OutputParam, BCParam, IBCParam, InjectionParam, \
    DirParam, Locations
from pymh.sim.simulations import \
    BasicSim
from pymh.utils.utils import \
    isPinRectangle

reload(pymh.param.parameters)

# Abbreviations
# fn = filename
# bd = boundary


# %%
if __name__ == '__main__':

    # ----------------------------------------------------------------------
    # Create directory structure
    # -------------------------------------------------------------------------
    Dir = DirParam()
    Dir.makedirs()

    # Add cluster

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
                                  receiver_origin=[0, 0, 40]))
    FullOutput.append(OutputParam('slice',
                                  timestep_increment=[50]))

    # IBC
    IBC = IBCParam('rigidboundary')
    IBCGrid = GridParam(origin=IBC.extraparameters['inj']['origin'],
                        number_of_cells=IBC.extraparameters['inj']['ncells'])
    IBCDecomposition = DecompositionParam(number_of_cells_per_node_x=[IBC.extraparameters['inj']['ncells'][0]],
                                          number_of_cells_per_node_y=[IBC.extraparameters['inj']['ncells'][1]],
                                          number_of_cells_per_node_z=[IBC.extraparameters['inj']['ncells'][2]])
    IBCModel = ModelParam(filename_prefix=['model_ibc'])
    IBCSimulation = SimulationParam()
    IBCInjection = InjectionParam()
    IBCOutput = []
    IBCOutput.append(OutputParam('slice',
                                 timestep_increment=[50]))

    # Injection
    Injection = InjectionParam('mps')
    InjectionGrid = GridParam(origin=Injection.extraparameters['origin'],
                        number_of_cells=Injection.extraparameters['ncells'])
    InjectionDecomposition = DecompositionParam(number_of_cells_per_node_x=[Injection.extraparameters['ncells'][0]],
                                                number_of_cells_per_node_y=[Injection.extraparameters['ncells'][1]],
                                              number_of_cells_per_node_z=[Injection.extraparameters['ncells'][2]])
    InjectionModel = ModelParam(filename_prefix=['model_inj'])
    InjectionOutput = []
    InjectionOutput.append(OutputParam('slice',
                                 timestep_increment=[50]))
    # Green's functions
    GFInput = InputParam(wavelet=['delta'],
                         spread=['single_point_staggered_velocity'])
    GFOutput = []
    GFOutput.append(OutputParam('sub_volume_boundary',
                                stagger_on_sub_volume=[True],
                                boundary_thickness=[1]))

    # %% ----------------------------------------------------------------------
    # Locations
    # -------------------------------------------------------------------------

    SinjLocations = Locations()
    SinjLocations.rectangle(origin=IBC.extraparameters['inj']['origin'],
                            number_of_cells=IBC.extraparameters['inj']['ncells'],
                            cell_size=IBC.extraparameters['inj']['cell_size'])
    # print SinjLocations.locations
    sinjlocations_fn = 'sinj_locations.txt'
    SinjLocations.write(sinjlocations_fn)

    SrecLocations = Locations()
    SrecLocations.rectangle(origin=IBC.extraparameters['rec']['origin'],
                            number_of_cells=IBC.extraparameters['rec']['ncells'],
                            cell_size=IBC.extraparameters['rec']['cell_size'])
    # print SrecLocations.locations
    sreclocations_fn = 'srec_locations.txt'
    SrecLocations.write(sreclocations_fn)

    # %% ----------------------------------------------------------------------
    # Is source outside or inside?
    # -------------------------------------------------------------------------

    corners = []
    for loc in SrecLocations.locations:
        if loc[3]//6 == 1:
            corners.append(loc[:3])
    corners = [corners[i] for i in [0, 1, 3, 2]]
    # print corners

    IBC.extraparameters['source_inside'] = \
        isPinRectangle(corners, FullInput.parameters['location'])

    # %% ----------------------------------------------------------------------
    # Extrapolation?
    # -------------------------------------------------------------------------
    IBC.extraparameters['extrap'] = True

    # %% ----------------------------------------------------------------------
    # Create files
    # -------------------------------------------------------------------------

    # Full
    if not IBC.extraparameters['source_inside']:
        injection_sxx_fn = 'injection_sxx_x_{}_y_{}_z_{}'.format(*FullInput.parameters['location'])
        FullOutput.append(OutputParam('sub_volume_boundary',
                                      receiver_locations=[sinjlocations_fn],
                                      filename_prefix=[injection_sxx_fn],
                                      boundary_thickness=[1],
                                      attribute=['S00XX']))

        injection_vn_fn = 'injection_vn_x_{}_y_{}_z_{}'.format(*FullInput.parameters['location'])
        FullOutput.append(OutputParam('sub_volume_boundary',
                                      receiver_locations=[sinjlocations_fn],
                                      filename_prefix=[injection_vn_fn],
                                      boundary_thickness=[1],
                                      attribute=['normal_velocity'],
                                      stagger_on_sub_volume=[True]))

        if IBC.parameters['type'][0] is 'freesurface':
            injection_sxx_staggered_fn = 'injection_sxx_staggered_x_{}_y_{}_z_{}'.format(*FullInput.parameters['location'])
            FullOutput.append(OutputParam('sub_volume_boundary',
                                          receiver_locations=[sinjlocations_fn],
                                          filename_prefix=[injection_sxx_staggered_fn],
                                          boundary_thickness=[1],
                                          attribute=['S00XX'],
                                          stagger_on_sub_volume=[True]))

    full_input_fn = 'full.txt'

    FullSim = BasicSim(FullGrid,
                       FullDecomposition,
                       Time,
                       FullModel,
                       FullSimulation,
                       FullInput,
                       FullOutput)

    FullSim.create(full_input_fn, path=Dir.parameters['ref'])

    # IBC
    if IBC.extraparameters['extrap']:
        ibc_extrap_sxx_fn = 'ebc_extrap_sxx'
        IBCOutput.append(OutputParam('sub_volume_boundary',
                                     receiver_locations=[sreclocations_fn],
                                     filename_prefix=[ibc_extrap_sxx_fn],
                                     boundary_thickness=[1],
                                     attribute=['S00XX']))

        ibc_extrap_vn_fn = 'ebc_extrap_vn'
        IBCOutput.append(OutputParam('sub_volume_boundary',
                                     receiver_locations=[sreclocations_fn],
                                     filename_prefix=[ibc_extrap_vn_fn],
                                     boundary_thickness=[1],
                                     attribute=['normal_velocity'],
                                     stagger_on_sub_volume=[True]))

    ibc_input_fn = 'ibc.txt'

    if IBC.extraparameters['source_inside']:
        InpOrInj = FullInput
    else:
        InpOrInj = IBCInjection
        if IBC.parameters['type'][0] is 'freesurface':
            InpOrInj.parameters.pop('injection_filelist_di')
        else:
            InpOrInj.parameters.pop('injection_filelist_mono')

    IBCSim = BasicSim(IBCGrid,
                      IBCDecomposition,
                      Time,
                      IBCModel,
                      IBCSimulation,
                      IBC,
                      InpOrInj,
                      IBCOutput)

    IBCSim.create(ibc_input_fn, path=Dir.parameters['ibc'])

    # Injection
    if IBC.extraparameters['extrap']:
        inj_extrap_sxx_fn = 'inj_extrap_sxx'
        InjectionOutput.append(OutputParam('sub_volume_boundary',
                                           receiver_locations=[sreclocations_fn],
                                           filename_prefix=[inj_extrap_sxx_fn],
                                           boundary_thickness=[1],
                                           attribute=['S00XX']))

        inj_extrap_vn_fn = 'inj_extrap_vn'
        InjectionOutput.append(OutputParam('sub_volume_boundary',
                                           receiver_locations=[sreclocations_fn],
                                           filename_prefix=[inj_extrap_vn_fn],
                                           boundary_thickness=[1],
                                           attribute=['normal_velocity'],
                                           stagger_on_sub_volume=[True]))

    injection_input_fn = 'injection.txt'

    InjectionSim = BasicSim(InjectionGrid,
                            InjectionDecomposition,
                            Time,
                            InjectionModel,
                            FullSimulation,
                            Injection,
                            InjectionOutput)

    InjectionSim.create(injection_input_fn, path=Dir.parameters['inj'])

    # Green's functions
    if IBC.parameters['type'][0] is 'freesurface':
        GFOutput[0].parameters['attribute'] = ['S00XX']
    else:
        GFOutput[0].parameters['attribute'] = ['normal_velocity']

    GFOutput[0].parameters['receiver_locations'] = [sinjlocations_fn]

    if IBC.extraparameters['extrap']:
        GFOutput.append(OutputParam('shot_gather',
                                    receiver_origin=[0, 0, 20],
                                    receiver_increment=[8, 0, 0],
                                    number_of_receivers=[125]))

    facedict = {0: ('x_di', 'x_source'),
                1: ('x_di', 'x_source'),
                2: ('z_di', 'z_source'),
                3: ('z_di', 'z_source')
                }

    for i in SrecLocations.locations:
        face = i[3] % 6
        for (j, k) in [('mono', 'isotropic_stress_source'), facedict[face]]:

            gf_input_fn = 'input_{}_x_{}_y_{}_z_{}.txt'.format(j, *i[:3])

            GFInput.parameters['type'] = [k]
            GFInput.parameters['location'] = list(i[:3])

            gf_fn = 'GF_{}_x_{}_y_{}_z_{}'.format(j, *i[:3])
            gf_extrap_fn = 'GF_extrap_{}_x_{}_y_{}_z_{}'.format(j, *i[:3])
            gf_extrap_fn2 = 'extrap_GF_{}_x_{}_y_{}_z_{}'.format(j, *i[:3])

            GFOutput[0].parameters['filename_prefix'] = [gf_fn]
            GFOutput[1].parameters['filename_prefix'] = [gf_extrap_fn]

            GFSim = BasicSim(FullGrid,
                             FullDecomposition,
                             Time,
                             FullModel,
                             FullSimulation,
                             GFInput,
                             GFOutput)

            GFSim.create(gf_input_fn, path=Dir.parameters['ibc_gf'])
            # diff_string = 'diff /w04d2/bfilippo/pymh/pymh/examples/IBC_simulation/GF_files/' + gf_extrap_fn + '.su /w04d2/bfilippo/matterhorn_filippo/tests/EBC/EBC_injection_v3/freesurface/EBC_simulation/GF_files/' + gf_extrap_fn2 + '.su'
            # diff_string = 'diff /w04d2/bfilippo/pymh/pymh/examples/IBC_simulation/GF_files/' + gf_fn + '_volume_boundary /w04d2/bfilippo/matterhorn_filippo/tests/EBC/EBC_injection_v3/freesurface/EBC_simulation/GF_files/' + gf_fn + '_volume_boundary'
            # diff_string = ['diff', 'IBC_simulation/GF_files/' + gf_fn + '_volume_boundary', 'IBC_simulation/GF_files/v3/EBC_simulation/GF_files/' + gf_fn + '_volume_boundary']
            # print diff_string
            # os.system(diff_string)
            # subprocess.call(diff_string)
# %%

#    for i in SrecLocations.locations:
        #subprocess.call(['diff', ])

# %%

    # shot = SEGYFile('shot.su', isSU=True)

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # im = ax.imshow(shot.readTraces().T, aspect='auto')
    # plt.show()
