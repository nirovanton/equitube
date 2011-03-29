# -*- coding: utf-8 -*-

########################################################################
# Copyright (C) 2011 by John Harris <harrisj@mnstate.edu>              #
#                                                                      #
# This program is free software; you can redistribute it and#or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 2 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program; if not, write to the                        #
# Free Software Foundation, Inc.,                                      #
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.            #
########################################################################

from equitube.tube import Tube
from equitube.field import Field
import sys
import time

try:
    import matplotlib.pyplot as plt
except:
    print >> sys.stderr, "Could not find matplotlib.  Perhaps, install matplotlib?"

import numpy as np

class Plot:
    
    def __init__(self, length = False):
        
        self._length = length

    def plotField(self, tube_dict):
        """A function to plot the tubes within the system.

        Generates a 2D square field of length = width, in this field
        it generates a periodic boundry condition and plots individual
        line segments for each tube object.
        """
        #plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axis([0,self._length,0,self._length])
        for key in tube_dict.keys():
            P = tube_dict[key].getParams()['P']
            Q = tube_dict[key].getParams()['Q']
            lines, = ax.plot([P[0],Q[0]],[P[1],Q[1]], '-o')
        plt.show()
        #plt.draw()
        #fig.canvas.draw()
        return None

