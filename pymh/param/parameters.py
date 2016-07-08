import os

# __all__ = ['DomainBase', 'RectangularDomain', 'DomainBC', 'Neumann',
# 'Dirichlet', 'PML']

# Mapping between dimension and the key labels for Cartesian domain
_cart_keys = {2: [(0, 'x'), (1, 'z')],
              3: [(0, 'x'), (1, 'y'), (2, 'z')]
              }

# %% Dictionaries

filedict = {
    'sinj': ['sinj_locations.txt'],
    'srec': ['srec_locations.txt']
    }

dirdict = {
    'ref': ['REF_simulation'],
    'inj': ['INJ_simulation'],
    'ebc': ['EBC_simulation'],
    'ebc_gf': ['EBC_simulation/GF_files']
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
    'freesurface': [True]
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
    'inj': {},
    'rec': {},
    'ngpts': [3],
    }

ibcdict['inj'] = {
    'origin': [400, 0, 300],
    'ncells': [101, 1, 75]
    }

ibcdict['rec'] = {
    'origin': [
        ibcdict['inj']['origin'][i] +
        ibcdict['ngpts'][0]*griddict['cell_size'][i] for i in range(3)
        ],
    'ncells': [
        ibcdict['inj']['ncells'][i] -
        2*ibcdict['ngpts'][0] for i in range(3)
        ],
    }

injectiondict = {
    'ngpts': [10],
    'origin': [],
    'ncells': [],
    }

injectiondict['origin'] = [
    ibcdict['inj']['origin'][i] -
    (injectiondict['ngpts'][0])*griddict['cell_size'][i] for i in range(3)
    ]

injectiondict['ncells'] = [
    ibcdict['inj']['ncells'][i] +
    2*injectiondict['ngpts'][0] for i in range(3)
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
    'timestep_increment': [50],
    'start_timestep': [0],
    'end_timestep': [timedict['number_of_timesteps'][0] - 1],
    'filename_prefix': ['output']
    }


# %%
class BaseParam(object):
    """ Base class for `Matterhorn` parameters. """

    type = None

    # def __init__(self):
    # self.dim = len(configs)

    def write(self, filename):
        """ Write section to file """

        filename.write('<%s>\n' % self.type)
        for key in self.parameters:
            # filename.writelines('{} {}\n'.format(k,v)
            # for k, v in self.parameters.items())
            values = ' '.join(str(self.parameters[key][k])
                for k in range(len(self.parameters[key])))
            filename.write('{} {}\n'.format(key, values))
        filename.write('</%s>\n\n' % self.type)


class GridParam(BaseParam):
    """ Class for describing grid parameters in `Matterhorn`.

    """

    type = 'grid'

    def __init__(self, **kwargs):
        self.parameters = griddict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class DecompositionParam(BaseParam):
    """ Class for describing decomposition parameters in `Matterhorn`.

    """

    type = 'decomposition'

    def __init__(self, **kwargs):
        self.parameters = decompositiondict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class TimeParam(BaseParam):
    """ Class for describing time parameters in `Matterhorn`.

    """

    type = 'time'

    def __init__(self, **kwargs):
        self.parameters = timedict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class ModelParam(BaseParam):
    """ Class for describing model parameters in `Matterhorn`.

    """

    type = 'model'

    def __init__(self, **kwargs):
        self.parameters = modeldict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)


class SimulationParam(BaseParam):
    """ Class for describing simulation parameters in `Matterhorn`.

    """

    type = 'simulation'

    def __init__(self, **kwargs):
        self.parameters = simulationdict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class BCParam(BaseParam):
    """ Class for describing boundary conditions parameters in `Matterhorn`.

    """

    type = 'bc'

    def __init__(self, bc, **kwargs):
        self.bc = bc
        if self.bc == 'pml':
            self.parameters = pmlbcdict.copy()
        elif self.bc == 'rigidboundary':
            self.parameters = rigidboundarybcdict.copy()
        elif self.bc == 'freesurface':
            self.parameters = freesurfacebcdict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)


class IBCParam(BaseParam):
    """ Class for describing IBC parameters in `Matterhorn`.

    """

    type = 'ibc'

    def __init__(self, **kwargs):
        self.parameters = ibcdict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class InjectionParam(BaseParam):
    """ Class for describing injection parameters in `Matterhorn`.

    """

    type = 'injection'

    def __init__(self, **kwargs):
        self.parameters = injectiondict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class InputParam(BaseParam):
    """ Class for describing input parameters in `Matterhorn`.

    """

    type = 'input'

    def __init__(self, **kwargs):
        self.parameters = inputdict.copy()
        for key in kwargs:
            self.parameters.update(kwargs)

class OutputParam(BaseParam):
    """ Class for describing output parameters in `Matterhorn`.

    """
    
    type = 'output'

    def __init__(self, output, **kwargs):
        self.parameters = outputdict.copy()
        self.output = output
        self.parameters['type'] = [self.output]
        if self.output is 'shot_gather':
            self.parameters['format'] = ['su']
            self.parameters['receiver_origin'] = [0, 0, 20]
            self.parameters['receiver_increment'] = [4, 0, 0]
            self.parameters['number_of_receivers'] = [250]
        elif self.output == 'slice':
            self.parameters['format'] = ['su']
            self.parameters['axis'] = ['y']
            self.parameters['slice_index'] = [0]
        elif self.output == 'sub_volume_boundary':
            self.parameters['receiver_locations'] = ['srec_locations.txt']
            self.parameters['boundary_thickness'] = [1]
        for key in kwargs:
            self.parameters.update(kwargs)

# %%
class DirParam(object):
    """
    Class for creating the directory tree for
    a `Matterhorn` simulation.
    """

    type = 'dir'

    def __init__(self, *configs):
        self.parameters = dirdict.copy()

        for d in self.parameters:
            self.parameters[d] = \
                '/'.join([os.getcwd(), str(self.parameters[d][0])])

    def makedirs(self):
        """ Create directories """

        for key in self.parameters:
            try:
                os.makedirs(self.parameters[key])
            except OSError:
                print 'Directory %s already exists\n' % self.parameters[key]

# %%
        # Setup access in the parameters dict by both letter and dimension number
        #for (i,k) in _cart_keys[self.dim]:

            #config = configs[i]

            #if len(config) == 4:
            #    L, R, lbc, rbc = config # min, max, left BC, right BC
            #    unit = None
            #elif len(config) == 5:
            #    L, R, lbc, rbc, unit = config # min, max, left BC, right BC, unit

            #param = Bunch(lbound=float(L), rbound=float(R), lbc=lbc, rbc=rbc, unit=unit)

            #param.length = float(R) - float(L)


            # access the dimension data by index, key, or shortcut
            #self.parameters[i] = param # d.dim[-1]
            #self.parameters[k] = param # d.dim['z']
            #self.__setattr__(k, param) # d.z

#    def collect(self, prop):
#        """ Collects a ndim-tuple of property `prop`."""
#        return tuple([self.parameters[i][prop] for i in xrange(self.dim)])
#
#    def get_lengths(self):
#        """Returns the physical size of the domain as a ndim-tuple."""
#        return self.collect('length')


        