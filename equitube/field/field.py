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
import numpy as np

class Field:
    """An object for creating the field that holds the tube network.
    
    The Field class only keeps track of the different instances of
    the tube objects. It also uses this list to gain access to each
    tube object and calculate the energies of the entire field.
    """

    def __init__(self,length):
        """Initializes the field
        """
        self._startP = {}
        self._startQ = {}
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
            params = self._tubes[x].getParams()
            self._startP[x] = list()
            self._startP[x] = params['P'] 
            self._startQ[x] = list()
            self._startQ[x] = params['Q']
        return None

    def getTubes(self):
        """Returns the generated tube array
        """
        return self._tubes
   
    def getTube(self,index):
        """Return the tube assocaited with the passed index
        """
        return self._tubes[index]
 
    def getVanderPotential(self):
        """Calculates the Van der Waals Potential for the system

        The intersections for each tube apply a torque that acts to 
        rotate the tube about the intersection point. This function
        calculates the effective torque and pivot point using the
        intersectons. It also approximates the total potential energy
        of the torques across the system. The function returns the total
        Van der Waals potential for the system as well as the information
        needed to rotate the tubes about the calculated pivot points.
        """
        for tube_id in self._tubes.keys():
            Fp = 0.0
            Fq = 0.0
            P = self._tubes[tube_id].getParams()['P'] #P[x,y]
            Q = self._tubes[tube_id].getParams()['Q'] #Q[x,y]
            slope = self._tubes[tube_id].getParams()['m']
            neighbour_dict = self._tubes[tube_id].getParams()['neighbours']
            
            for index in neighbour_dict.keys():
                xint = neighbour_dict[index][1]
                yint = slope*xint + self._tubes[tube_id].getParams()['b']
                phi = neighbour_dict[index][0]
                torque = 1/np.sin(phi)-1
                slope2 = self._tubes[index].getParams()['m']
                
                if slope >= 0:
                    Rp = np.sqrt((xint-P[0])**2+(yint-P[1])**2)
                    Rq = np.sqrt((Q[0]-xint)**2+(Q[1]-yint)**2)
                else:
                    Rp = np.sqrt((xint-P[0])**2+(P[1]-yint)**2)
                    Rq = np.sqrt((Q[0]-xint)**2+(yint-Q[1])**2)
                
                if slope > slope2:
                    if phi > np.pi/2:
                        Fp += -1*torque/Rp
                        Fq += torque/Rq
                        continue
                    Fp += torque/Rp
                    Fq += -1*torque/Rq
                    continue
                elif slope < slope2:
                    if phi > np.pi/2:
                        Fp += torque/Rp
                        Fq += -1*torque/Rq
                        continue
                    Fp += -1*torque/Rp
                    Fq += torque/Rq
                    continue

        #return potential,pivot_dict
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
                if x >= tube1['P'][0] and x <= tube1['Q'][0] \
                and x >= tube2['P'][0] and x <= tube2['Q'][0] :
                    if tube1['m'] > tube2['m']:
                        angle = tube1['theta'] - tube2['theta']
                    else:
                        angle = tube2['theta'] - tube1['theta']
                    self._tubes[dex1].addNeighbours(dex2,angle,x)
        return None

