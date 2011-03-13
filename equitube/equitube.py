# -*- coding: utf-8 -*-

#########################################################################
# Copyright (C) 2011 by John Harris <harrisj@mnstate.edu>               #
# This program is free software; you can redistribute it and#or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation; either version 2 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program; if not, write to the                         #
# Free Software Foundation, Inc.,                                       #
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
#########################################################################

import matplotlib.pyplot as plt
import numpy as np
import math
import random
import optparse


class Equitube:
    def __init__(self,argv):
        
        usage = "usage: %prog [options]"
        parser = optparse.OptionParser(usage=usage)
        variables, arguments = self._parseOptions(argv, parser)

        self._density = variables.density

    def _parseOptions(self, argv, parser):

        density_help_list = [
            "This option allows you to specify the nanotube ",
            "density for the system. This is used to approximate ",
            "the number of tubes created during the simulation."
            ]
        parser.add_option('--density', '-d', default=0.5,
            help=''.join(density_help_list))
        area_help_list = [
            "This option allows you to specity the size of the ",
            "field on which the nanotubes are oriented."
            ]
        parser.add_option('--area', '-a', default=10,
            help=''.join(area_help_list))
       
        return parser.parse_args()


class Tube:
    def __init__(self, argv):
        self._m = None
        self._l = None
        self._b = None
        self._xcm = None
        self._xmin = None
        self._xmax = None
        self._ycm = None
        self._ymin = None
        self._ymax = None
        self._theta = None

    def create_line():
        self._m = random.randint(-20,20)
        self._l = random.uniform(5,12)
        self._xcm, self._ycm = random.uniform(10,15)
        self._theta = np.arctan(math.sqrt(self._m*self._m))
        self._xmin = self._xcm - self._l*np.cos(self._theta)
        self._xmax = self._xcm + self._l*np.cos(self._theta)
        self._ymin = self._ycm - self._l*np.sin(self._theta)
        self._ymax = self._ycm + self._l*np.sin(self._theta)
        self._b = self._ycm - self._m*self._xcm 

        return self._m, self._l, self._theta, self._xmin, self._xcm, self._xmax, self._ymin, self._ycm, self._ymax


    
class EquitubeException(Exception):
    def __init__(self, message, *args):
        super(EquitubeException, self).__init__(args)
        self._message = message

    def GetMessage(self):
        return self._message


'''
def intersect(x_min, x_max, b1, b2, m1, m2):
    
    # This function evaluates 2 node to see if they intersect within
    # the node range x_min and x_max both need to come from the same
    # node, it doesn't which node.
 
    x_int = (b2-b1)/(m1-m2)
    if x_min <= x_int and x_int <= x_max:
        return True
    return False
'''


'''
test
fig = plt.figure()
ax = fig.add_subplot(111)
x, y = np.random.rand(2, 2)
x1,y1 = np.random.rand(2, 2)
line = ax.plot(x, y, 'bs-')
line1 = ax.plot(x1, y1, 'bs-', picker=5)

plt.draw()
plt.show()
'''
