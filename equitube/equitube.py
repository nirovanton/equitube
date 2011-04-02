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

import sys

try:
    import matplotlib.pyplot as plt
except:
    print >> sys.stderr, "Could not load matplotlib.  Perhaps, install matplotlib?"

import optparse
import numpy as np
from field import Field
from plot import Plot
from tube import Tube

class Equitube:
    def __init__(self,argv):
        
        usage = "usage: %prog [options]"
        parser = optparse.OptionParser(usage=usage)
        variables, arguments = self._parseOptions(argv, parser)

        self._count = int(variables.count)
        self._springconst = float(variables.spring)
        self._length = int(variables.length)
        self._deltaslope = float(variables.deltaslope)

    def _parseOptions(self, argv, parser):

        count_help_list = [
            "This option allows you to specify the nanotube ",
            "count for the system. This is the number of tubes",
            "created during the simulation."
            ]
        parser.add_option('--count', '-c', default=50,
            help=''.join(count_help_list))
        length_help_list = [
            "This option allows you to specify the length of the ",
            "field on which the nanotubes are oriented, the field ",
            "is square for simplicity."
            ]
        parser.add_option('--length', '-l', default=10,
            help=''.join(length_help_list))
        spring_help_list = [
            "This option allos you to specify the strength of the ",
            "spring constant and control the elasticity of the ",
            "substrate the tubes are adhered to."
            ]
        parser.add_option('--spring', '-k', default=.001,
            help=''.join(spring_help_list))
        deltaslope_help_list = [
            "This option allows you to specify the amount of change ",
            "the slope undergoes at each iteration"
            ]
        parser.add_option('--deltaslope', '-m', default = 0.05,
            help=''.join(deltaslope_help_list))
        return parser.parse_args()
        
    def Run(self):
        try:
            field = Field(self._length, self._springconst, self._deltaslope)
            field.addTubes(self._count)
            plot = Plot(self._length)
            tubes = field.getTubes()

            end = 0
            while end < 1:
                field.calculateIntercepts()
                point_forces = field.getPointForces()
                traverses = 0
                
                for index in tubes.keys():
                    if tubes[index].getParams()['P'][0] <= 0:
                        traverses += field.traverseNeighbours(index,[])
                print traverses, len(field.getPaths())
                for path in field.getPaths():
                    print path,"::",len(path)
                field.resetPaths()
                print len(field.getPaths())


                
                #if end % 1 == 0:
                #    plot.plotField(tubes)
                field.rotateTubes(point_forces)
                end += 1
                
        except EquitubeException, e:
            raise EquitubeException(e.get_message())

class EquitubeException(Exception):
    def __init__(self, message, *args):
        super(EquitubeException, self).__init__(args)
        self._message = message

    def GetMessage(self):
        return self._message

