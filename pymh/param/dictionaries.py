# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 17:48:12 2016

@author: bfilippo
"""

# %% Dictionaries

filedict = {
    'sinj': ['sinj_locations.txt'],
    'srec': ['srec_locations.txt']
    }

dirdict = {
    'base': ['REF_simulation'],
    'ref': ['REF_simulation'],
    'inj': ['INJ_simulation'],
    'ibc': ['IBC_simulation'],
    'ibc_gf': ['IBC_simulation/GF_files']
    }

userdict = {
    'user': ['bfilippo'],
    'platform': ['aug']
    }

griddict = {
    'origin': [0, 0, 0],
    'number_of_cells': [251, 1, 201],
    'cell_size': [4, 4, 4]
    }

decompositiondict = {
    'number_of_nodes': [1, 1, 1],
    'number_of_cells_per_node_x': [griddict['number_of_cells'][0]],
    'number_of_cells_per_node_y': [griddict['number_of_cells'][1]],
    'number_of_cells_per_node_z': [griddict['number_of_cells'][2]]
    }

timedict = {
    'timestep_size': [8e-4],
    'number_of_timesteps': [3200]
    }

modeldict = {
    'type': ['acoustic'],
    'format': ['su'],
    }

simulationdict = {
    'type': ['acoustic_2d'],
    'order_time_operator': [2],
    'spatial_operator': ['taylor_2'],
    'freesurface': [False]
    }

pmlbcdict = {
    'bc_type': ['pml'],
    'pml_width_in_gridpoints': [10],
    'pml_power': [4.0],
    'pml_frequency': [25.0],
    'pml_damping_vel': [2000.0],
    }

rigidboundarybcdict = {
    'bc_type': ['bc_specified_per_edge'],
    'bc_x_min': ['rigidboundary'],
    'bc_x_max': ['rigidboundary'],
    'bc_z_min': ['rigidboundary'],
    'bc_z_max': ['rigidboundary']
    }

freesurfacebcdict = {
    'bc_type': ['bc_specified_per_edge'],
    'bc_x_min': ['freesurface'],
    'bc_x_max': ['freesurface'],
    'bc_z_min': ['freesurface'],
    'bc_z_max': ['freesurface']
    }

ibcdict = {
    'type': ['rigidboundary'],
    'sinj_locations': ['sinj_locations.txt'],
    'srec_locations': ['srec_locations.txt'],
    'greensfunction_filelist_mono': ['GF_mono_file_list'],
    'greensfunction_filelist_di': ['GF_di_file_list']
    }

ibcextradict = {
    'extrap': [False],
    'inj': {},
    'rec': {},
    'ngpts': [3],
    }

ibcextradict['inj'] = {
    'origin': [400, 0, 300],
    'ncells': [101, 1, 75],
    'cell_size': griddict['cell_size']
    }

ibcextradict['rec'] = {
    'origin': [
        ibcextradict['inj']['origin'][i] +
        ibcextradict['ngpts'][0]*griddict['cell_size'][i] for i in range(3)
        ],
    'ncells': [
        ibcextradict['inj']['ncells'][i] -
        2*ibcextradict['ngpts'][0] for i in range(3)
        ],
    'cell_size': griddict['cell_size']
    }

injectiondict = {
    'method': ['mps'],
    'injection_locations': ['sinj_locations.txt'],
    'injection_filelist_mono': ['injection_mono_file_list'],
    'injection_filelist_di': ['injection_di_file_list']
    }

injectionextradict = {
    'ngpts': [10],
    'origin': [],
    'ncells': []
    }

injectionextradict['origin'] = [
    ibcextradict['inj']['origin'][i] -
    (injectionextradict['ngpts'][0])*griddict['cell_size'][i] for i in range(3)
    ]

injectionextradict['ncells'] = [
    ibcextradict['inj']['ncells'][i] +
    2*injectionextradict['ngpts'][0] for i in range(3)
    ]

inputdict = {
    'type': ['isotropic_stress_source'],
    'wavelet': ['ricker'],
    'location': [200, 0, 250],
    'spread': ['trilinear'],
    'scale_factor': [1.0],
    'shift': [1.0],
    'central_frequency': [25.0]
    }

outputdict = {
    'attribute': ['S00XX'],
    'timestep_increment': [1],
    'start_timestep': [0],
    'end_timestep': [timedict['number_of_timesteps'][0] - 1],
    'filename_prefix': ['output']
    }

# %%

fullsimdict = {
    0: 'Grid',
    1: 'Decomposition',
    2: 'Time',
    3: 'Model',
    4: 'Simulation',
    5: 'Input',
    6: 'Output'}

fullsimtuple = tuple(fullsimdict.values())
