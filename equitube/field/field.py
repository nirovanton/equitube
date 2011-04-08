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

    def __init__(self,length,spring_const,
                 vander_const,radius,increment):
        """Initializes the field
        """
        self._k = spring_const
        self._V = vander_const
        self._a = radius*10**-9
        self._inc = increment

        self._startP = {}
        self._startQ = {}
        self._tubes = {}
        self._length = length

    def addTubes(self, number = 1):
        """This function adds a tube to the system.

        This function will initialize another tube object and will add it
        to the list of active tubes. The non-random nature of the computer 
        generated random #s is evident here as many tubes share a common 
        approximate length and position. the time.sleep statements are my 
        feeble attempt to increase the randomness of the generated lines.
        """
        for x in range(number):
            self._tubes[x] = Tube()
            theta  = random.uniform(0,2*np.pi)
            time.sleep(random.uniform(random.uniform(.0001,0),0))
            cm = [random.uniform(0,self._length-.5),random.uniform(0,self._length-.5)]
            self._tubes[x].createLine(cm,theta)
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
   
    def calculateIntercepts(self):
        """Calculate the intercept of all the tubes.
        
        This function takes the collection of tubes and determines where
        they intercept as well as the angle of intercept. all intercept
        angles are calculated to the right side of each intercept.

        The nested if statements that convert a theta to a phi are
        important steps for the ease of future calculations of the 
        point forces used in relaxational rotations.
        """
        for dex1 in range(len(self._tubes)):
            tube1 = self._tubes[dex1].getParams()
            for dex2 in range(len(self._tubes)):
                tube2 = self._tubes[dex2].getParams()
                if tube1['m'] == tube2['m']:
                    continue
                x = (tube2['b']-tube1['b'])/(tube1['m']-tube2['m'])
                if x >= tube1['P'][0] and x <= tube1['Q'][0] \
                and x >= tube2['P'][0] and x <= tube2['Q'][0]:
                    angle = abs(tube1['theta']-tube2['theta'])
                    self._tubes[dex1].addNeighbors(dex2,angle,x)
        return None
    

    def potential(self):
        """ Calculates the energy of the network.

        This function considers the endpoint displacement
        and the angle of interaction to generate an
        approximate numerical value for the energy.
        """
        energy = 0.0
        for tube_id in self._tubes.keys():
            params = self._tubes[tube_id].getParams()
           
            # spring potential  ~ kr^2
            P_dist_sqrd = (params['P'][0]-self._startP[tube_id][0])**2 + \
                          (params['P'][1]-self._startP[tube_id][1])**2
            Q_dist_sqrd = (params['Q'][0]-self._startQ[tube_id][0])**2 + \
                          (params['Q'][1]-self._startQ[tube_id][1])**2
            energy += (P_dist_sqrd + Q_dist_sqrd) * self._k
            
            # vanderwaals potental ~ aV/sin(theta)
            for key in params['neighbors'].keys():
                theta = params['neighbors'][key][0]
                energy += (self._a*self._V)/abs(np.sin(theta))

        return energy
           
    
    def relax(self):
        """ A Bruteforce method for relaxing the network.

        This offsets each variable (x,y,theta) by an incremental amount
        it then checks to see if the resulting change reduced the energy
        of the system. If a +dx and a -dx both resulted in an increase in
        energy the program considers the incremental value to be too great,
        it halves the increment and continues to the next variable.

        increment_values is a dict() that stores the individual dx dy dtheta
        increments for each tube.
        """
        relaxed = False

        # prepare the increment dictionary.
        increment_values = {}
        for tube_id in self._tubes.keys():
            increment_values[tube_id] = list()
            increment_values[tube_id] = [self._inc,self._inc,self._inc]

        # relax the network. 
        while not relaxed:
            self.calculateIntercepts()
            previous_potential = self.potential()
            
            for tube_id in self._tubes.keys():
                x_lock = False
                y_lock = False
                theta_lock = False
                
                dx = increment_values[tube_id][0]
                dy = increment_values[tube_id][1]
                dtheta = increment_values[tube_id][2]
                
                x,y = self._tubes[tube_id].getParams()['cm']
                theta = self._tubes[tube_id].getParams()['theta']
                                
                """ Minimize the x coordinate.
                """
                self._tubes[tube_id].createLine([x+dx,y],theta)
                self.calculateIntercepts()
                if self.potential() < previous_potential:
                    x_lock = True
                    x += dx
                if x_lock == False:
                    self._tubes[tube_id].createLine([x-dx,y],theta)
                    self.calculateIntercepts()
                    if self.potential() < previous_potential:
                        x_lock = True
                        x += -dx
                """ If increasing and decreasing the y coordinate and
                still not smaller energy: increment is too big. Cut
                incrementer in half.
                """
                if x_lock == False:
                    increment_values[tube_id][0] = dx/2

                """ Minimize the y coordinate.
                """
                self._tubes[tube_id].createLine([x,y+dy],theta)
                self.calculateIntercepts()
                if self.potential() < previous_potential:
                    y_lock = True
                    y += dy
                if y_lock == False:
                    self._tubes[tube_id].createLine([x,y-dy],theta)
                    self.calculateIntercepts()
                    if self.potential() < previous_potential:
                        y_lock = True
                        y += -dy
                """ If increasing and decreasing the y coordinate and
                still not smaller energy: increment is too big. Cut 
                incrementer in half.
                """
                if y_lock == False:
                    increment_values[tube_id][1] = dy/2

                """ Minimizing theta
                """
                self._tubes[tube_id].createLine([x,y],theta + dtheta)
                self.calculateIntercepts()
                if self.potential() < previous_potential:
                    theta_lock = True
                    theta += dtheta
                if theta_lock == False:
                    self._tubes[tube_id].createLine([x,y],theta - dtheta)
                    self.calculateIntercepts()
                    if self.potential() < previous_potential:
                        theta_lock = True
                        theta += -dtheta
                """ If increasing or decreasing the angle, and still 
                not smaller energy: Increment is too big. Cut incrementer
                in half.
                """
                if theta_lock == False:
                    increment_values[tube_id][2] = dtheta/2
                    self._tubes[tube_id].createLine([x,y],theta)

 
            
            """ Verify that all of the increments are very small. If they
            are <= 10^-20 consider the system relaxed.
            """
            find = True
            for inc_list in increment_values.values():
                if (round(inc_list[0],20) != 0) \
                or (round(inc_list[1],20) != 0) \
                or (round(inc_list[2],20) != 0):
                    find = False
                    continue
            if find:
                relaxed = True

        return None            

