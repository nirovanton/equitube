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


class Field():
    """An object for creating the field that holds the tube network.
    
    The Field class only keeps track of the different instances of
    the tube objects. It also uses this list to gain access to each
    tube object and calculate the energies of the entire field.
    """

    def __init__(self):
        """Initializes the field
        """
        self._tubes = []

    def addTubes(number = 1):
        """This function adds a tube to the system 

        This function will initialize another tube object and will add it
        to the list of active tubes.
        """

        self.tubes = [ Tube for x in range(0, number) ]
        return None 

    def getVanderPotential():
        """Calculates the Van der Waals Potential for the system

        This function uses the attributes stored within each tube object
        to calculate the Van der Waals energies as a function of theta.
        for the entire system. It will have the form 1/sin(theta)
        """
        return None

    def getSpringPotential():
        """Calculuates the Spring Potential energy stored in the substate.
        
        This function will consider small changes in center and end
        point position and use it to calculate the spring potential energy.
        of the entire system.
        """
        return None
