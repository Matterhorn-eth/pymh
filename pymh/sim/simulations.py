# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:14:19 2016

@author: bfilippo
"""

from __future__ import absolute_import
import os

__all__ = ['BaseSim', 'BasicSim']


# %%
class BaseSim(object):
    """ Base class for `Matterhorn` simulations. """

    type = None

    # def __init__(self):
    # self.dim = len(configs)

    def create(self, filename, path=os.curdir):
        """ Create simulation file """

        self.file = open('/'.join([path, filename]), 'w')

        for param in self.parameters:
            try:
                for i in param:
                    i.write(self.file)
            except:
                param.write(self.file)

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
