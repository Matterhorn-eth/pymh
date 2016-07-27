# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:14:19 2016

@author: bfilippo
"""

from __future__ import absolute_import

__all__ = ['isPinRectangle']


# %%
def isPinRectangle(r, P):
    """
        r: A list of four points, each has a x- and a y- coordinate
        P: A point
    """

    # print P[0], P[1]
    areaRectangle = 0.5*abs(
        #                 y_A      y_C      x_D      x_B
                        (r[0][2]-r[2][2])*(r[3][0]-r[1][0])
        #                  y_B     y_D       x_A     x_C
                      + (r[1][2]-r[3][2])*(r[0][0]-r[2][0])
                    )

    ABP = 0.5*abs(
             r[0][0]*(r[1][2]-P[2])
            +r[1][0]*(P[2]-r[0][2])
            +P[0]*(r[0][2]-r[1][2])
          )
    BCP = 0.5*abs(
             r[1][0]*(r[2][2]-P[2])
            +r[2][0]*(P[2]-r[1][2])
            +P[0]*(r[1][2]-r[2][2])
          )
    CDP = 0.5*abs(
             r[2][0]*(r[3][2]-P[2])
            +r[3][0]*(P[2]-r[2][2])
            +P[0]*(r[2][2]-r[3][2])
          )
    DAP = 0.5*abs(
             r[3][0]*(r[0][2]-P[2])
            +r[0][0]*(P[2]-r[3][2])
            +P[0]*(r[3][2]-r[0][2])
          )

    # print 'area: %s\n' % areaRectangle
    # print 'sum: %s %s %s %s\n' % (ABP, BCP, CDP, DAP)
    return areaRectangle == (ABP+BCP+CDP+DAP)