# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:14:19 2016

@author: bfilippo
"""

from __future__ import absolute_import
import os
from pymh.param.dictionaries import *
from pymh.param.parameters import OutputParam

__all__ = ['BaseSim', 'BasicSim', 'FullSim', 'IBCSim']


# %%
class BaseSim(object):
    """ Base class for `Matterhorn` simulations. """

    type = None

    def __init__(self,
                 simtuple,
                 Param,
                 extrap,
                 source_inside,
                 ibc_type):

        self.parameters = {}
        self.tuple = simtuple
        for i in range(len(self.tuple)):
            self.parameters[self.tuple[i]] = Param[i]
        self.extrap = extrap
        self.source_inside = source_inside
        self.ibc_type = ibc_type

    def create(self, filename, path=os.curdir):
        """ Create simulation file """

        self.file = open('/'.join([path, filename]), 'w')

        for param in self.tuple:
            try:
                for i in self.parameters[param]:
                    i.write(self.file)
            except:
                self.parameters[param].write(self.file)

        self.file.close()

    def utility(self, dim='2d', prefix='LOL', nt=(),
                sinj='sinj_locations.txt',
                srec='srec_locations.txt',
                path=os.curdir):
        """ Create utility files """

        fn = {}
        for md in ('mono', 'di'):
            fn[md] = prefix + '_' + md + '_utility.txt'
            fid = open('/'.join([path, fn[md]]), 'w')
            fid.write(dim + '\n')
            fid.write('input_file_list_name {}\n'.format(prefix + '_' + md + '_volume_boundary_list.txt'))
            fid.write('output_file_name_prefix {}\n'.format(prefix + '_' + md))
            fid.write('sinj_locations {}\n'.format(sinj))
            fid.write('srec_locations {}\n'.format(srec))
            fid.write('number_of_timesteps {}\n'.format(nt))
            fid.close()

        return fn

    def volume_boundary(self, list_fn={}, prefix='LOL', path=os.curdir):
        """ Create volume boundary files """

        fn = {}
        for md in ('mono', 'di'):
            fn[md] = prefix + '_' + md + '_volume_boundary_list.txt'
            fid = open('/'.join([path, fn[md]]), 'w')
            for i in list_fn[md]:
                fid.write('{}_volume_boundary\n'.format(i))
            fid.close()

        return fn


class BasicSim(BaseSim):
    """ Class for describing basic simulations in `Matterhorn`.

    """

    type = 'basic'

    def __init__(self, *args):
        self.parameters = args
#        for key in kwargs:
#            self.parameters.update(kwargs)


class FullSim(BaseSim):
    """ Class for describing full simulations in `Matterhorn`.

    """

    type = 'full'

    def __init__(self,
                 Param,
                 extrap=False,
                 source_inside=True,
                 ibc_type='freesurface',
                 locations_fn='sinj_locations.txt',
                 *args):

        super(FullSim, self).__init__(fullsimtuple,
                                      Param,
                                      extrap,
                                      source_inside,
                                      ibc_type)

        if not self.source_inside:
            self.injection_mono_fn = 'injection_sxx_x_{}_y_{}_z_{}'.format(*self.parameters['Input'].parameters['location'])
            self.parameters['Output'].append(OutputParam('sub_volume_boundary',
                                          receiver_locations=[locations_fn],
                                          filename_prefix=[self.injection_mono_fn],
                                          boundary_thickness=[1],
                                          attribute=['S00XX'],
                                          end_timestep=[self.parameters['Time'].parameters['number_of_timesteps'][0] - 1]))

            self.injection_di_fn = 'injection_vn_x_{}_y_{}_z_{}'.format(*self.parameters['Input'].parameters['location'])
            self.parameters['Output'].append(OutputParam('sub_volume_boundary',
                                          receiver_locations=[locations_fn],
                                          filename_prefix=[self.injection_di_fn],
                                          boundary_thickness=[1],
                                          attribute=['normal_velocity'],
                                          stagger_on_sub_volume=[True],
                                          end_timestep=[self.parameters['Time'].parameters['number_of_timesteps'][0] - 1]))

            if self.ibc_type is 'freesurface':
                self.injection_mono_staggered_fn = 'injection_sxx_staggered_x_{}_y_{}_z_{}'.format(*self.parameters['Input'].parameters['location'])
                self.parameters['Output'].append(OutputParam('sub_volume_boundary',
                                              receiver_locations=[locations_fn],
                                              filename_prefix=[self.injection_mono_staggered_fn],
                                              boundary_thickness=[1],
                                              attribute=['S00XX'],
                                              stagger_on_sub_volume=[True],
                                              end_timestep=[self.parameters['Time'].parameters['number_of_timesteps'][0] - 1]))

class IBCSim(BaseSim):
    """ Class for describing IBC simulations in `Matterhorn`.

    """

    type = 'IBC'

    def __init__(self,
                 Param,
                 extrap=False,
                 source_inside=True,
                 ibc_type='freesurface',
                 locations_fn='srec_locations.txt',
                 *args):

        super(IBCSim, self).__init__(ibcsimtuple,
                                     Param,
                                     extrap,
                                     source_inside,
                                     ibc_type)

        if not self.extrap:
            self.ibc_extrap_mono_fn = 'ibc_extrap_sxx'
            self.parameters['Output'].append(OutputParam('sub_volume_boundary',
                                          receiver_locations=[locations_fn],
                                          filename_prefix=[self.ibc_extrap_mono_fn],
                                          boundary_thickness=[1],
                                          attribute=['S00XX'],
                                          end_timestep=[self.parameters['Time'].parameters['number_of_timesteps'][0] - 1]))

            self.ibc_extrap_di_fn = 'ibc_extrap_vn'
            self.parameters['Output'].append(OutputParam('sub_volume_boundary',
                                          receiver_locations=[locations_fn],
                                          filename_prefix=[self.ibc_extrap_di_fn],
                                          boundary_thickness=[1],
                                          attribute=['normal_velocity'],
                                          stagger_on_sub_volume=[True],
                                          end_timestep=[self.parameters['Time'].parameters['number_of_timesteps'][0] - 1]))