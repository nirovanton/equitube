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
import random
import time

class Field:
    """An object for creating the field that holds the tube network.
    
    The Field class only keeps track of the different instances of
    the tube objects. It also houses several methods to calculate
    intercepts, point forces, traversal paths, and rotate the tubes.
    """

    def __init__(self,length,k,dm):
        """Initializes the field
        """
        self._dm = dm # +/- change in slope per iteration
        self._k = k # spring constant
        self._startP = {}
        self._startQ = {}
        self._tubes = {}
        self._paths = {}
        self._length = length

    def addTubes(self, number = 1):
        """This function adds a tube to the system 

        This function will initialize another tube object and will add it
        to the list of active tubes. The non-random nature of the computer 
        generated random #s is evident here as many tubes share a common 
        approximate length and position. the time.sleep statements are my 
        feeble attempt to increase the randomness of the generated lines.
        """
       
        for x in range(number):
            time.sleep(random.uniform(.001,0))
            self._tubes[x] = Tube()
            m = random.uniform(-7,7)
            time.sleep(random.uniform(.001,0))
            cm = [random.uniform(0,self._length-.5),random.uniform(0,self._length-.5)]
            self._tubes[x].createLine(m,cm)
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
   
    def getPointForces(self):
        """Calculates the Point forces for each tube in the system.

        The intersections for each tube apply a torque that acts to 
        rotate the tube about the intersection point. This function
        calculates the effective torque and pivot point using the
        intersectons. It also factors in the displacement from the
        original position and uses it to calculate the spring force acting
        on the ends of each tube. This is factored into the torque
        calculations as well.
        """
        
        force_dict = {}
        for tube_id in self._tubes.keys():
            force_dict[tube_id] = list()
            Fp = 0.0
            Fq = 0.0
            P = self._tubes[tube_id].getParams()['P'] #P[x,y]
            Q = self._tubes[tube_id].getParams()['Q'] #Q[x,y]
            slope = self._tubes[tube_id].getParams()['m']
            
            p_dist = np.sqrt((P[0]-self._startP[tube_id][0])**2 \
                           + (P[1]-self._startP[tube_id][1])**2)
            q_dist = np.sqrt((Q[0]-self._startQ[tube_id][0])**2 \
                           + (Q[1]-self._startQ[tube_id][1])**2)

            if self._startP[tube_id][1] > P[1]:
                Fp += self._k*p_dist
            else:
                Fp += -1*self._k*p_dist

            if self._startQ[tube_id][1] > Q[1]:
                Fq += self._k*q_dist
            else:
                Fq += -1*self._k*q_dist
                        
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
                    elif phi < np.pi/2:    
                        Fp += torque/Rp
                        Fq += -1*torque/Rq
                elif slope < slope2:
                    if phi > np.pi/2:
                        Fp += torque/Rp
                        Fq += -1*torque/Rq
                    elif phi < np.pi/2:
                        Fp += -1*torque/Rp
                        Fq += torque/Rq
                
            force_dict[tube_id] = [Fp,Fq]
        return force_dict

    def rotateTubes(self, force_dict):
        """This function rotates the tubes about a calculated point.
        
        Each intersect for a given tube generates a torque about the
        intersection point. The getVanderPotential() fuction calculates
        the total torque on a tube as well as the relative pivot point.
        This function takes the direction of rotation and the pivot
        and uses them to rotate the tobe by a small amount.

        l = 1/2 length of tube
        d = distance from pivot to center
        """
        # TODO Ponder a more graceful method.. 'cuz damn.
        for tube_id in self._tubes.keys():
            l = self._tubes[tube_id].getParams()['l']
            m = self._tubes[tube_id].getParams()['m']
            cm = self._tubes[tube_id].getParams()['cm']
            theta = abs(self._tubes[tube_id].getParams()['theta'])
            Fp = force_dict[tube_id][0]
            Fq = force_dict[tube_id][1]
            #print tube_id , Fp,Fq, m
            if Fp == 0 and Fq == 0:
                continue
            
            if m > 0:
                # M+ P+ Q+
                if Fp > 0 and Fq > 0:
                    if Fp < Fq:
                        d = ((Fq+Fp)/(Fq-Fp))*l
                        m += self._dm
                        x_piv = cm[0] - d*np.cos(theta)
                        y_piv = cm[1] - d*np.sin(theta)
                        xcm = x_piv + d*np.cos(np.arctan(m))
                        ycm = y_piv + d*np.sin(np.arctan(m))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue    
                    elif Fp > Fq:
                        d = ((Fp+Fq)/(Fp-Fq))*l
                        m += -1*self._dm
                        x_piv = cm[0] + d*np.cos(theta)
                        y_piv = cm[1] + d*np.sin(theta)
                        xcm = x_piv - d*np.cos(np.arctan(m))
                        ycm = y_piv - d*np.sin(np.arctan(m))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                # M+ P- Q-
                elif Fp < 0 and Fq < 0:
                    if Fp < Fq:
                        d = ((abs(Fq)+abs(Fp))/(abs(Fq)-abs(Fp)))*l
                        m += self._dm
                        x_piv = cm[0] + d*np.cos(theta)
                        y_piv = cm[1] + d*np.sin(theta)
                        xcm = x_piv - d*np.cos(np.arctan(m))
                        ycm = y_piv - d*np.sin(np.arctan(m))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                    elif Fp > Fq:
                        d = ((abs(Fp)+abs(Fq))/(abs(Fp)-abs(Fq)))*l
                        m += -1*self._dm
                        x_piv = cm[0] - d*np.cos(theta)
                        y_piv = cm[1] - d*np.sin(theta)
                        xcm = x_piv + d*np.cos(np.arctan(m))
                        ycm = y_piv + d*np.sin(np.arctan(m))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                # M+ P+/- Q+/-
                elif abs(Fp) < abs(Fq):
                    d = ((abs(Fq)-abs(Fp))/(abs(Fq)+abs(Fp)))*l
                    if Fp > 0:
                        m += -1*self._dm
                    elif Fp < 0:
                        m += self._dm
                    x_piv = cm[0] - d*np.cos(theta)
                    y_piv = cm[1] - d*np.sin(theta)
                    xcm = x_piv + d*np.cos(np.arctan(m))
                    ycm = y_piv + d*np.sin(np.arctan(m))
                    self._tubes[tube_id].createLine(m,[xcm,ycm])
                    continue
                elif abs(Fp) > abs(Fq):
                    d = ((abs(Fp)-abs(Fq))/(abs(Fp)+abs(Fq)))*l
                    if Fp > 0:
                        m += -1*self._dm
                    elif Fp < 0:
                        m += self._dm
                    x_piv = cm[0] + d*np.cos(theta)
                    y_piv = cm[1] + d*np.sin(theta)
                    xcm = x_piv - d*np.cos(np.arctan(m))
                    ycm = y_piv - d*np.sin(np.arctan(m))
                    self._tubes[tube_id].createLine(m,[xcm,ycm])
                    continue
            
            elif m < 0:
                # M- P+ Q+
                if Fp > 0 and Fq > 0:
                    if Fp < Fq:
                        d = ((Fq+Fp)/(Fq-Fp))*l
                        m += self._dm
                        x_piv = cm[0] - d*np.cos(theta)
                        y_piv = cm[1] + d*np.sin(theta)
                        xcm = x_piv + d*np.cos(np.arctan(abs(m)))
                        ycm = y_piv - d*np.sin(np.arctan(abs(m)))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                    elif Fp > Fq:
                        d = ((Fp+Fq)/(Fp-Fq))*l
                        m += -1*self._dm
                        x_piv = cm[0] + d*np.cos(theta)
                        y_piv = cm[1] - d*np.sin(theta)
                        xcm = x_piv - d*np.cos(np.arctan(abs(m)))
                        ycm = y_piv + d*np.sin(np.arctan(abs(m)))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                # M- P- Q-
                if Fp < 0 and Fq < 0:
                    if Fp < Fq:
                        d = ((abs(Fq)+abs(Fp))/(abs(Fq)-abs(Fp)))*l
                        m += self._dm
                        x_piv = cm[0] + d*np.cos(theta)
                        y_piv = cm[1] - d*np.sin(theta)
                        xcm = x_piv - d*np.cos(np.arctan(abs(m)))
                        ycm = y_piv + d*np.sin(np.arctan(abs(m)))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                    elif Fp > Fq:
                        d = ((abs(Fp)+abs(Fq))/(abs(Fp)-abs(Fq)))*l
                        m += -1*self._dm
                        x_piv = cm[0] - d*np.cos(theta)
                        y_piv = cm[1] + d*np.sin(theta)
                        xcm = x_piv + d*np.cos(np.arctan(abs(m)))
                        ycm = y_piv - d*np.sin(np.arctan(abs(m)))
                        self._tubes[tube_id].createLine(m,[xcm,ycm])
                        continue
                # M- P+/- Q+/-
                elif abs(Fp) < abs(Fq):
                    d = ((abs(Fq)-abs(Fp))/(abs(Fq)+abs(Fp)))*l
                    if Fp > 0: 
                        m += -1*self._dm
                    elif Fp < 0:
                        m += self._dm
                    x_piv = cm[0] - d*np.cos(theta)
                    y_piv = cm[1] + d*np.sin(theta)
                    xcm = x_piv + d*np.cos(np.arctan(abs(m)))
                    ycm = y_piv - d*np.sin(np.arctan(abs(m)))
                    self._tubes[tube_id].createLine(m,[xcm,ycm])
                    continue
                elif abs(Fp) > abs(Fq):
                    d = ((abs(Fp)-abs(Fq))/(abs(Fp)+abs(Fq)))*l
                    if Fp > 0: 
                        m += -1*self._dm
                    elif Fp < 0:
                        m += self._dm
                    x_piv = cm[0] + d*np.cos(theta)
                    y_piv = cm[1] - d*np.sin(theta)
                    xcm = x_piv - d*np.cos(np.arctan(abs(m)))
                    ycm = y_piv + d*np.sin(np.arctan(abs(m)))
                    self._tubes[tube_id].createLine(m,[xcm,ycm])
                    continue
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
    
    def traverseNeighbours(self, index, prev_tubes):
        """Calculates the total number of paths from left to right.

        All tubes with ( P[0] <= 0 ) intersect the left side of
        the field. These are the tubes initially passed to this
        function. This function then gathers neighbours and pass
        them recursively. prev_tubes is a list of tube ID #'s. 
        It keeps track of the previous tubes along a given path of
        recursion. It prevents backtracking.
        """
        prev_tubes.append(index)
        print prev_tubes

        neighbours = self._tubes[index].getParams()['neighbours']
        neighbour_key_list = neighbours.keys()
        count = 0
        if self._tubes[index].getParams()['Q'][0] >= self._length:
            count += 1
            print prev_tubes
        for key in neighbour_key_list:
            if prev_tubes.count(key) > 0:
                continue
            count += self.traverseNeighbours(key,prev_tubes)
        return count
