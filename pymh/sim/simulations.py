# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:14:19 2016

@author: bfilippo
"""
import os


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


class BasicSim(BaseSim):
    """ Class for describing basic simulations in `Matterhorn`.

    """

    type = 'basic'

    def __init__(self, *args):
        self.parameters = args
#        for key in kwargs:
#            self.parameters.update(kwargs)
