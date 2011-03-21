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
import math
import random

class Tube:
    def __init__(self,length = False):
        """ Initializes a clean tube
        
        the length parameter specifies the dimmentions
        of the field, It is passed to the class to ensure
        that the randomly generated tubes have a center pt
        inside the bounds of the field.
        """
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
        self._length = length
        self._neighbours = {}

    def createLine(self):
        """ Generating a random Line segment
        """
        self._m = random.randint(-50,50)
        self._l = random.randint(3,10)
        self._xcm = random.uniform(0,self._length)
        self._ycm = random.uniform(0,self._length)
        self._theta = np.arctan(math.sqrt(self._m*self._m))
        self._xmin = self._xcm - self._l*np.cos(self._theta)
        self._xmax = self._xcm + self._l*np.cos(self._theta)
        self._ymin = self._ycm - self._l*np.sin(self._theta)
        self._ymax = self._ycm + self._l*np.sin(self._theta)
        self._b = self._ycm - self._m*self._xcm

        return None

    def addNeighbours(self, key, angle, x_intercept):
        """ Function for adding neighbors
        
        The Field class keeps track of the tubes in the field with
        a dict, key->Tube. This function populates another dict,
        (key of intercepting tube) -> (x coord of intercept)
        y = mx+b is then used to find y coord of intercept
        """
        self._neighbours[key] = dict()
        self._neighbours[key][angle] = x_intercept 
     
        return None
    
    def getParams(self):
        """ A function that returns requested parameters.
        """
        params = {'m':self._m,'l':self._l,'b':self._b,'xcm':self._xcm,'ycm':self._ycm,'xmin':self._xmin,'xmax':self._xmax,'ymin':self._ymin,'ymax':self._ymax,'theta':self._theta,'neighbours':self._neighbours}

        return params

