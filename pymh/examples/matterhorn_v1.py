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
import copy

sys.dont_write_bytecode = True

pymhfolder = '/w04d2/bfilippo/pymh'
# pymhfolder = '/Users/filippo/work/pymh'
sys.path.insert(0, pymhfolder)

import pymh
from pymh.param.parameters import \
    GridParam, DecompositionParam, TimeParam, ModelParam, SimulationParam, \
    InputParam, OutputParam, BCParam, IBCParam, InjectionParam, \
    DirParam, LocationsParam
from pymh.sim.simulations import \
    BasicSim
from pymh.utils.utils import \
    isPinRectangle

reload(pymh.param.parameters)

# Abbreviations
# fn = filename
# bd = boundary
# vb = volume boundary


# %%
if __name__ == '__main__':

    # ----------------------------------------------------------------------
    # Create directory structure
    # -------------------------------------------------------------------------
    Dir = DirParam()
    Dir.makedirs()

    # Add cluster directory

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

    Locations = {}
    locations_fn = {}
    for surf in ('inj', 'rec'):
        Locations[surf] = LocationsParam()
        Locations[surf].rectangle(origin=IBC.extraparameters[surf]['origin'],
                                  number_of_cells=IBC.extraparameters[surf]['ncells'],
                                  cell_size=IBC.extraparameters[surf]['cell_size'])
        # print Locations[surf].locations
        locations_fn[surf] = 's' + surf + '_locations.txt'
        Locations[surf].write(locations_fn[surf])

    # Temporary, eventually I need to fix this
    # It is only needed when using gffu to reorder injection boundary data
    Locations_for_inj_util = copy.deepcopy(Locations['rec'])
    Locations_for_inj_util.locations[1:] = []
    locations_inj_fn = 's' + 'rec' + '_locations_injection.txt'
    Locations_for_inj_util.write(locations_inj_fn, path=Dir.parameters['ref'])

    # %% ----------------------------------------------------------------------
    # Is source outside or inside?
    # -------------------------------------------------------------------------

    # This also needs to be improved a lot
    # Right now, it only check if the source is inside srec
    corners = []
    for loc in Locations['rec'].locations:
        if loc[3]//6 == 1:
            corners.append(loc[:3])
    corners = [corners[i] for i in [0, 1, 3, 2]]
    # print corners

    IBC.extraparameters['source_inside'] = \
        isPinRectangle(corners, FullInput.parameters['location'])

    # %% ----------------------------------------------------------------------
    # Extrapolation?
    # -------------------------------------------------------------------------

    # Does it really have to be included in IBC?
    IBC.extraparameters['extrap'] = True

    # %% ----------------------------------------------------------------------
    # Create files
    # -------------------------------------------------------------------------

    # Full
    if not IBC.extraparameters['source_inside']:
        injection_mono_fn = 'injection_sxx_x_{}_y_{}_z_{}'.format(*FullInput.parameters['location'])
        FullOutput.append(OutputParam('sub_volume_boundary',
                                      receiver_locations=[locations_fn['inj']],
                                      filename_prefix=[injection_mono_fn],
                                      boundary_thickness=[1],
                                      attribute=['S00XX']))

        injection_di_fn = 'injection_vn_x_{}_y_{}_z_{}'.format(*FullInput.parameters['location'])
        FullOutput.append(OutputParam('sub_volume_boundary',
                                      receiver_locations=[locations_fn['inj']],
                                      filename_prefix=[injection_di_fn],
                                      boundary_thickness=[1],
                                      attribute=['normal_velocity'],
                                      stagger_on_sub_volume=[True]))

        if IBC.parameters['type'][0] is 'freesurface':
            injection_mono_staggered_fn = 'injection_sxx_staggered_x_{}_y_{}_z_{}'.format(*FullInput.parameters['location'])
            FullOutput.append(OutputParam('sub_volume_boundary',
                                          receiver_locations=[locations_fn['inj']],
                                          filename_prefix=[injection_mono_staggered_fn],
                                          boundary_thickness=[1],
                                          attribute=['S00XX'],
                                          stagger_on_sub_volume=[True]))

    injection_fn = {}
    injection_fn['mono'] = [injection_mono_fn]
    injection_fn['di'] = [injection_di_fn]

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
        ibc_extrap_mono_fn = 'ebc_extrap_sxx'
        IBCOutput.append(OutputParam('sub_volume_boundary',
                                     receiver_locations=[locations_fn['rec']],
                                     filename_prefix=[ibc_extrap_mono_fn],
                                     boundary_thickness=[1],
                                     attribute=['S00XX']))

        ibc_extrap_di_fn = 'ebc_extrap_vn'
        IBCOutput.append(OutputParam('sub_volume_boundary',
                                     receiver_locations=[locations_fn['rec']],
                                     filename_prefix=[ibc_extrap_di_fn],
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
        inj_extrap_mono_fn = 'inj_extrap_sxx'
        InjectionOutput.append(OutputParam('sub_volume_boundary',
                                           receiver_locations=[locations_fn['rec']],
                                           filename_prefix=[inj_extrap_mono_fn],
                                           boundary_thickness=[1],
                                           attribute=['S00XX']))

        inj_extrap_di_fn = 'inj_extrap_vn'
        InjectionOutput.append(OutputParam('sub_volume_boundary',
                                           receiver_locations=[locations_fn['rec']],
                                           filename_prefix=[inj_extrap_di_fn],
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

    GFOutput[0].parameters['receiver_locations'] = [locations_fn['inj']]

    if IBC.extraparameters['extrap']:
        GFOutput.append(OutputParam('shot_gather',
                                    receiver_origin=[0, 0, 20],
                                    receiver_increment=[8, 0, 0],
                                    number_of_receivers=[125]))

    facedict = {0: (('mono', 'isotropic_stress_source', 'mono'), ('x_di', 'x_source', 'di')),
                1: (('mono', 'isotropic_stress_source', 'mono'), ('x_di', 'x_source', 'di')),
                2: (('mono', 'isotropic_stress_source', 'mono'), ('z_di', 'z_source', 'di')),
                3: (('mono', 'isotropic_stress_source', 'mono'), ('z_di', 'z_source', 'di'))
                }

    gf_input_fn = {}
    gf_fn = {}
    gf_extrap_fn = {}
    for i in ('mono', 'di'):
        gf_input_fn[i] = []
        gf_fn[i] = []
        gf_extrap_fn[i] = []

    for (i, loc) in enumerate(Locations['rec'].locations):
        face = loc[3] % 6
        for (j, k, l) in facedict[face]:

            gf_input_fn[l].append('input_{}_x_{}_y_{}_z_{}.txt'.format(j, *loc[:3]))

            GFInput.parameters['type'] = [k]
            GFInput.parameters['location'] = list(loc[:3])

            gf_fn[l].append('GF_{}_x_{}_y_{}_z_{}'.format(j, *loc[:3]))
            gf_extrap_fn[l].append('GF_extrap_{}_x_{}_y_{}_z_{}'.format(j, *loc[:3]))
            # gf_extrap_fn2 = 'extrap_GF_{}_x_{}_y_{}_z_{}'.format(j, *loc[:3])

            GFOutput[0].parameters['filename_prefix'] = [gf_fn[l][i]]
            GFOutput[1].parameters['filename_prefix'] = [gf_extrap_fn[l][i]]

            GFSim = BasicSim(FullGrid,
                             FullDecomposition,
                             Time,
                             FullModel,
                             FullSimulation,
                             GFInput,
                             GFOutput)

            GFSim.create(gf_input_fn[l][i], path=Dir.parameters['ibc_gf'])

# %%
    # GF_mono_volume_boundary_list.txt
    # GF_di_volume_boundary_list.txt
    gf_vb_fn = IBCSim.volume_boundary(prefix='GF', list_fn=gf_fn,
                   path=Dir.parameters['ibc_gf'])

# %%
    # GF_mono_utility.txt
    # GF_di_utility.txt
    gf_util_fn = IBCSim.utility(prefix='GF', nt=Time.parameters['number_of_timesteps'][0],
                   sinj='sinj_locations.txt',
                   srec='srec_locations.txt',
                   path=Dir.parameters['ibc_gf'])

# %%
    # injection_mono_volume_boundary_list.txt
    # injection_di_volume_boundary_list.txt
    inj_vb_fn = IBCSim.volume_boundary(prefix='injection', list_fn=injection_fn,
                   path=Dir.parameters['ref'])

# %%
    # injection_mono_utility.txt
    # injection_di_utility.txt
    InjectionSim.utility(prefix='injection', nt=Time.parameters['number_of_timesteps'][0],
               sinj='sinj_locations.txt',
               srec='srec_locations_injection.txt',
               path=Dir.parameters['ref'])

# %%

#    for i in SrecLocations.locations:
        #subprocess.call(['diff', ])

            # diff_string = 'diff /w04d2/bfilippo/pymh/pymh/examples/IBC_simulation/GF_files/' + gf_extrap_fn + '.su /w04d2/bfilippo/matterhorn_filippo/tests/EBC/EBC_injection_v3/freesurface/EBC_simulation/GF_files/' + gf_extrap_fn2 + '.su'
            # diff_string = 'diff /w04d2/bfilippo/pymh/pymh/examples/IBC_simulation/GF_files/' + gf_fn + '_volume_boundary /w04d2/bfilippo/matterhorn_filippo/tests/EBC/EBC_injection_v3/freesurface/EBC_simulation/GF_files/' + gf_fn + '_volume_boundary'
            # diff_string = ['diff', 'IBC_simulation/GF_files/' + gf_fn + '_volume_boundary', 'IBC_simulation/GF_files/v3/EBC_simulation/GF_files/' + gf_fn + '_volume_boundary']
            # print diff_string
            # os.system(diff_string)
            # subprocess.call(diff_string)
# %%

    # shot = SEGYFile('shot.su', isSU=True)

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # im = ax.imshow(shot.readTraces().T, aspect='auto')
    # plt.show()
