/*######################################################################
# Copyright (C) 2011 by John Harris <john.harris@ndsu.edu>             #
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
# Fr ee Software Foundation, Inc.,                                      #
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.            #
######################################################################*/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

// global vars.
//---------------------
int tube_count = 50;
int field_size = 10;
float vand_const = 1.0;
float spring_const = 1.0;
float radius = 1e-9;
float** tube_array = new float*[tube_count];
float** p_initial = new float*[tube_count];
float** q_initial = new float*[tube_count];
//---------------------


float calculateEnergy()
{
    float potential = 0.0;
    
    for(int k=0; k<tube_count; k++)
    {
        float m1 = tan(tube_array[k][2]);
        float b1 = tube_array[k][1]-m1*tube_array[k][0];
 
        float p1[2];
        float q1[2];
        p1[0] = tube_array[k][0] - abs(tube_array[k][3]*cos(tube_array[k][2]));
        q1[0] = tube_array[k][0] + abs(tube_array[k][3]*cos(tube_array[k][2]));
        if(m1 >= 0)
        {
            p1[1] = tube_array[k][1] - abs(tube_array[k][3]*sin(tube_array[k][2]));
            q1[1] = tube_array[k][1] + abs(tube_array[k][3]*sin(tube_array[k][2]));
        }
        else
        {
            p1[1] = tube_array[k][1] + abs(tube_array[k][3]*sin(tube_array[k][2]));
            q1[1] = tube_array[k][1] - abs(tube_array[k][3]*sin(tube_array[k][2]));
        }
            
        // spring
        potential += spring_const*(pow((p_initial[k][0]-p1[0]),2)+pow((p_initial[k][1]-p1[1]),2));
        potential += spring_const*(pow((q_initial[k][0]-q1[0]),2)+pow((q_initial[k][1]-q1[1]),2));

        // Vanderwaals
        for(int j=0; j<tube_count; j++)
        {
            float m2 = tan(tube_array[j][2]);
            if(m1 == m2)
                continue;
            float b2 = tube_array[j][1]-m2*tube_array[j][0];            
            
            float p2[2];
            float q2[2];
            p2[0] = tube_array[j][0] - abs(tube_array[j][3]*cos(tube_array[j][2]));
            q2[0] = tube_array[j][0] + abs(tube_array[j][3]*cos(tube_array[j][2]));
            if(m2 >= 0)
            {
                p2[1] = tube_array[j][1] - abs(tube_array[j][3]*sin(tube_array[j][2]));
                q2[1] = tube_array[j][1] + abs(tube_array[j][3]*sin(tube_array[j][2]));
            }   
            else
            {   
                p2[1] = tube_array[j][1] + abs(tube_array[j][3]*sin(tube_array[j][2]));
                q2[1] = tube_array[j][1] - abs(tube_array[j][3]*sin(tube_array[j][2]));
            } 

            float x_int = (b2-b1)/(m1-m2);
            if(x_int >= p1[0] && x_int >= p2[0] && x_int <= q1[0] && x_int <= q2[0])
                potential += -1*(radius*vand_const)/abs(sin(tube_array[k][2]-tube_array[j][2]));
        }
    }
    
    return potential;
}

// TODO 
// ------------
// void relaxNetwork()
// void compressNetwork()
// float calculateCurrent()
// ------------

int main ()
{
    srand(time(0));
   
    float max_angle = M_PI/2;
    float min_angle = -M_PI/2;
    float energy;

    // Create the tubes.
    for(int i = 0; i < tube_count; i++)
    {
        // Random variables
        float angle = (min_angle-max_angle)*(float)rand()/(float)(RAND_MAX)+max_angle;
        float x_cm = field_size*(float)rand()/(float)(RAND_MAX);
        float y_cm = field_size*(float)rand()/(float)(RAND_MAX);
        float length = 3-(float)rand()/(float)(RAND_MAX);
      
        // Populate the tube array.
        tube_array[i] = new float[4];
        tube_array[i][0] = x_cm;
        tube_array[i][1] = y_cm;
        tube_array[i][2] = angle;
        tube_array[i][3] = length;

        // (P) O------------O (Q) ( Px < Qx )
        // populate the initial endpoint arrays.
        p_initial[i] = new float[2];
        q_initial[i] = new float[2];
        p_initial[i][0] = x_cm - abs(length*cos(angle));
        q_initial[i][0] = x_cm + abs(length*cos(angle));
        if(tan(angle) >= 0)
        {
            p_initial[i][1] = y_cm - abs(length*sin(angle));
            q_initial[i][1] = y_cm + abs(length*sin(angle));
        }
        else
        {
            p_initial[i][1] = y_cm + abs(length*sin(angle));
            q_initial[i][1] = y_cm - abs(length*sin(angle));
        }
    }
    
    energy = calculateEnergy();
    cout << energy <<endl;
    return 0;
}

