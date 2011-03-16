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
        self._length = lenth

    def create_line():
        self._m = random.randint(-20,20)
        self._l = random.randint(3,10)
        self._xcm, self._ycm = random.uniform(0,self.length)
        self._theta = np.arctan(math.sqrt(self._m*self._m))
        self._xmin = self._xcm - self._l*np.cos(self._theta)
        self._xmax = self._xcm + self._l*np.cos(self._theta)
        self._ymin = self._ycm - self._l*np.sin(self._theta)
        self._ymax = self._ycm + self._l*np.sin(self._theta)
        self._b = self._ycm - self._m*self._xcm
        return self._m, self._l, self._theta, self._xmin, self._xcm, self._xmax, self._ymin, self._ycm, self._ymax
