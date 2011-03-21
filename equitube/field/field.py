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

class Field:
    """An object for creating the field that holds the tube network.
    
    The Field class only keeps track of the different instances of
    the tube objects. It also uses this list to gain access to each
    tube object and calculate the energies of the entire field.
    """

    def __init__(self,length):
        """Initializes the field
        """
        self._tubes = {}
        self._length = length

    def addTubes(self, number = 1):
        """This function adds a tube to the system 

        This function will initialize another tube object and will add it
        to the list of active tubes.
        """
        
        for x in range(number):
            self._tubes[x] = Tube(self._length)
            self._tubes[x].createLine()
        
        return None

    def getTubes(self):
        """Returns the generated tube array
        """
        return self._tubes
   
    def getTube(self,index):
        """Return the tube assocaited with the passed index
        """
        return self.tubes[index]
 
    def getVanderPotential(self):
        """Calculates the Van der Waals Potential for the system

        This function uses the attributes stored within each tube object
        to calculate the Van der Waals energies as a function of theta.
        for the entire system. It will have the form 1/sin(theta)
        """
        return None

    def getSpringPotential(self):
        """Calculates the Spring Potential energy stored in the substate.
        
        This function will consider small changes in center and end
        point position and use it to calculate the spring potential energy.
        of the entire system.
        """
        return None

    def calculateIntercepts(self):
        """Calculate the intercept of all the tubes.
        
        This function takes the collection of tubes and determines where
        they intercept as well as the angle of intercept. all intercept
        angles are calculated to the right side of each intercept.
        """
        
        for dex1 in range(len(self._tubes)):
            tube1 = self._tubes[dex1].getParams()
            for dex2 in range(len(self._tubes)):
                tube2 = self._tubes[dex2].getParams()
                if tube1['m'] == tube2['m']:
                    continue
                x = (tube2['b']-tube1['b'])/(tube1['m']-tube2['m'])
   
                # Verify that the tubes intersect.
                #TODO Make this more elegant.
                if x >= tube1['xmin'] and x <= tube1['xmax'] \
                and x >= tube2['xmin'] and x <= tube2['xmax'] :
                    
                    # Both >= 0
                    if tube1['m'] >= 0 and tube2['m'] >= 0:
                        if tube1['m'] > tube2['m']:
                            angle = tube1['theta']-tube2['theta']
                        else:
                            angle = tube2['theta']-tube1['theta']
                    
                    # Both < 0
                    if tube1['m'] < 0 and tube2['m'] < 0:
                        if tube1['m'] < tube2['m']:
                            angle = tube1['theta']-tube2['theta']
                        else:
                            angle = tube2['theta']-tube1['theta']
                    
                    # One < 0 and One >= 0
                    if tube1['m'] >= 0 and tube2['m'] < 0:
                        angle = tube1['theta']+tube2['theta']
                    if tube1['m'] < 0 and tube2['m'] >= 0:
                        angle = tube1['theta']+tube2['theta']
                    
                    self._tubes[dex1].addNeighbours(dex2,angle,x)
                    continue
                else:
                    continue
        
        return None 
