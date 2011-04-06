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

import numpy as np
import random

class Tube:
    def __init__(self):
        """ Initializes a clean tube
        
        the length parameter specifies the dimensions
        of the field, It is passed to the class to ensure
        that the randomly generated tubes have a center pt
        inside the bounds of the field.
        
        regardless of slope point P is always on left side.
        P(x,y) 0----------------0 Q(x,y)
        """
        self._m = None
        self._l = None
        self._b = None
        self._theta = None
        self._neighbors = {}
        self._cm = [] #[Xcm,Ycm]
        self._P = [] #[Xp,Yp]
        self._Q = [] #[Xq,Yq]

    def createLine(self,center,theta):
        """ Generates a Line segment

        This function generates the initial line segment, It also
        modifies the variables that are based on slope and center
        point after each iteration.
        """
        if self._l == None:
            self._l = random.uniform(4,8)
        self._theta = theta
        self._m = np.tan(theta)
        self._b = center[1] - self._m*center[0]
        self._neighbors = {}
        self._cm = center

        xmax = self._cm[0] + abs(self._l*np.cos(theta))
        xmin = self._cm[0] - abs(self._l*np.cos(theta))
        ymax = self._cm[1] + abs(self._l*np.sin(theta))
        ymin = self._cm[1] - abs(self._l*np.sin(theta))
        if self._m <= 0:
            self._P = [xmin,ymax]
            self._Q = [xmax,ymin]
        else:
            self._P = [xmin,ymin]
            self._Q = [xmax,ymax]

        return None

    def addNeighbors(self, key, angle, x_intercept):
        """ Function for adding neighbors
        
        The Field class keeps track of the intersecting
        tubes via a dict'd list. {key:[angle,x_intecept]}
        """
        self._neighbors[key] = list()
        self._neighbors[key] = [angle,x_intercept]
     
        return None
    
    def getParams(self):
        """ A function that returns requested parameters.

        This function takes all of the individual tube
        attributes, and returns them as a dict.
        """
        params = {
            'm':self._m,
            'l':self._l,
            'b':self._b,
            'cm':self._cm,
            'theta':self._theta,
            'P':self._P,
            'Q':self._Q,
            'neighbors':self._neighbors
            }

        return params

