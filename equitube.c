/*######################################################################
# Copyright (C) 2011 by John Harris <john.harrisj@ndsu.edu>            #
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
######################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mygraph.h>

//GLOBALS 
//=========================================
#define tube_count 10
int field_size = 10;
double energy;
double radius = 1.0;
double vand_const = 1.0;
double spring_const = 1.0;
double tube_array[tube_count][4];
double p_initial[tube_count][4];
double iteration_array[tube_count][9];
int done=0,pause=1,sstep=1;
//=========================================

double calculateEnergy()
{
    double potential = 0.0;
    for(int k=0; k<tube_count; k++)
    {
        double m1 = tan(tube_array[k][2]);
        double b1 = tube_array[k][1]-m1*tube_array[k][0];
  
        double p1[4];
        p1[0] = tube_array[k][0] - fabs(tube_array[k][3]*cos(tube_array[k][2]));
        p1[2] = tube_array[k][0] + fabs(tube_array[k][3]*cos(tube_array[k][2]));
        if(m1 >= 0)
        {
            p1[1] = tube_array[k][1] - fabs(tube_array[k][3]*sin(tube_array[k][2]));
            p1[3] = tube_array[k][1] + fabs(tube_array[k][3]*sin(tube_array[k][2]));
        }
        else
        {
            p1[1] = tube_array[k][1] + fabs(tube_array[k][3]*sin(tube_array[k][2]));
            p1[3] = tube_array[k][1] - fabs(tube_array[k][3]*sin(tube_array[k][2]));
        }
            
        // spring
        potential += spring_const*(pow((p_initial[k][0]-p1[0]),2)+pow((p_initial[k][1]-p1[1]),2));
        potential += spring_const*(pow((p_initial[k][2]-p1[2]),2)+pow((p_initial[k][3]-p1[3]),2));

        // Vanderwaals
        for (int j=0; j<tube_count; j++)
        {
            double m2 = tan(tube_array[j][2]);
            if(m1 == m2)
                continue;
            double b2 = tube_array[j][1]-m2*tube_array[j][0];            
            
            double p2[4];
            p2[0] = tube_array[j][0] - fabs(tube_array[j][3]*cos(tube_array[j][2]));
            p2[2] = tube_array[j][0] + fabs(tube_array[j][3]*cos(tube_array[j][2]));
            if(m2 >= 0)
            {
                p2[1] = tube_array[j][1] - fabs(tube_array[j][3]*sin(tube_array[j][2]));
                p2[3] = tube_array[j][1] + fabs(tube_array[j][3]*sin(tube_array[j][2]));
            }   
            else
            {   
                p2[1] = tube_array[j][1] + fabs(tube_array[j][3]*sin(tube_array[j][2]));
                p2[3] = tube_array[j][1] - fabs(tube_array[j][3]*sin(tube_array[j][2]));
            } 

            double x_int = (b2-b1)/(m1-m2);
            if(x_int >= p1[0] && x_int >= p2[0] && x_int <= p1[2] && x_int <= p2[2])
                potential += -1*(radius*vand_const)/fabs(sin(tube_array[k][2]-tube_array[j][2]));
        }
    }
    return potential;
}


void relaxNetwork()
{
    double starting_energy = calculateEnergy();
    for(int j = 0; j < tube_count; j++)
    {
        // ADJUST THE ITERATIONS
        // ---------------------
        if(iteration_array[j][1] > 4 || iteration_array[j][2] > 4)
        {
            double tmp = iteration_array[j][0]+.25*iteration_array[j][0];
            iteration_array[j][0] = tmp;
        }
        if(iteration_array[j][4] > 4 || iteration_array[j][5] > 4)
        {
            double tmp = iteration_array[j][3]+.25*iteration_array[j][3];
            iteration_array[j][3] = tmp;
        }
        if(iteration_array[j][7] > 4 || iteration_array[j][8] > 4)
        {
            double tmp = iteration_array[j][6]+.25*iteration_array[j][6];
            iteration_array[j][6] = tmp;
        }
    
        // MINIMIZE X
        // ----------
        int x_min = 0;
        // <incease x>
        tube_array[j][0] += iteration_array[j][0];
        if(calculateEnergy() > starting_energy)
            tube_array[j][0] += -1*iteration_array[j][0];
        else if(calculateEnergy() <= starting_energy)
        {
            starting_energy = calculateEnergy();                
            iteration_array[j][1] += 1;
            iteration_array[j][2] = 0;
            x_min = 1;
        }
        if(x_min == 0)
        {
            // <decrease x>
            tube_array[j][0] += -1*iteration_array[j][0];
            if(calculateEnergy() > starting_energy)
                tube_array[j][0] += iteration_array[j][0];
            else if(calculateEnergy() <= starting_energy)
            {
                starting_energy = calculateEnergy();
                iteration_array[j][1] = 0;
                iteration_array[j][2] += 1;
                x_min = 1;
            }
            if(x_min == 0)
            {
                double tmp = iteration_array[j][0]/2;
                iteration_array[j][0] = tmp;
                iteration_array[j][1] = 0;
                iteration_array[j][2] = 0;
            }
        }
             
        // MINIMIZE Y
        // ----------
        int y_min = 0;
        // <incease y>
        tube_array[j][1] += iteration_array[j][3];
        if(calculateEnergy() > starting_energy)
            tube_array[j][1] += -1*iteration_array[j][3];
        if(calculateEnergy() <= starting_energy)
        {
            starting_energy = calculateEnergy();
            iteration_array[j][4] += 1;
            iteration_array[j][5] = 0;
            y_min = 1;
        }
        if(y_min == 0)
        {
            // <decrease y>
            tube_array[j][1] += -1*iteration_array[j][3];
            if(calculateEnergy() > starting_energy)
                tube_array[j][1] += iteration_array[j][3];
            if(calculateEnergy() <= starting_energy)
            {
                starting_energy = calculateEnergy();
                iteration_array[j][4] = 0;
                iteration_array[j][5] += 1;
                y_min = 1;
            }
            if(y_min == 0)
            {
                double tmp = iteration_array[j][3]/2;
                iteration_array[j][3] = tmp;
                iteration_array[j][4] = 0;
                iteration_array[j][5] = 0;
            }
        }

        // MINIMIZE THETA
        // --------------
        int theta_min = 0;
        // <incease theta>
        tube_array[j][2] += iteration_array[j][6];
        if(calculateEnergy() > starting_energy)
            tube_array[j][2] += -1*iteration_array[j][6];
        if(calculateEnergy() <= starting_energy)
        {
            starting_energy = calculateEnergy();
            iteration_array[j][7] += 1;
            iteration_array[j][8] = 0;
            theta_min = 1;
        }
        if(theta_min == 0)
        {
            // <decrease theta>
            tube_array[j][2] += -1*iteration_array[j][6];
            if(calculateEnergy() > starting_energy)
                tube_array[j][2] += iteration_array[j][6];
            if(calculateEnergy() <= starting_energy)
            {
                starting_energy = calculateEnergy();
                iteration_array[j][7] = 0;
                iteration_array[j][8] += 1;
                theta_min = 1;
            }
            if(theta_min == 0)
            {
                double tmp = iteration_array[j][6]/2;
                iteration_array[j][6] = tmp;
                iteration_array[j][7] = 0;
                iteration_array[j][8] = 0;
            }
        }
    }
}

// TODO 
// ------------
// double calculateCurrent()
// ------------


void compressNetwork()
{
    for(int a = 0; a < tube_count; a++)
    {
        if(tube_array[a][0]<field_size/2)
        {
            double offset = 1-2*tube_array[a][0]/field_size;
            double tmp = tube_array[a][0] + offset;
            tube_array[a][0] = tmp;
        }
        if(tube_array[a][0]>field_size/2)
        {
            double offset = 1-2*(field_size-tube_array[a][0])/field_size;
            double tmp = tube_array[a][0] - offset;
            tube_array[a][0] = tmp;
        }
    }
}


void Initialize()
{
    srand(time(0));

    double max_angle = M_PI/2;
    double min_angle = -M_PI/2;
 
    // Create the tubes.
    for(int i = 0; i < tube_count; i++)
    {
        // Random variables
        double angle = (min_angle-max_angle)*(double)rand()/(double)(RAND_MAX)+max_angle;
        double x_cm = field_size*(double)rand()/(double)(RAND_MAX);
        double y_cm = field_size*(double)rand()/(double)(RAND_MAX);
        double length = 3-(double)rand()/(double)(RAND_MAX);

        // Initialize the iteration array.
        // 0 -> x_iter  1-> x+  2-> x-
        // 3 -> y_iter  4-> y+  5-> y-
        // 6 -> t_iter  7-> t+  8-> t-
        iteration_array[i][0] = .05;
        iteration_array[i][1] = 0;
        iteration_array[i][2] = 0;
        iteration_array[i][3] = .05;
        iteration_array[i][4] = 0;
        iteration_array[i][5] = 0;
        iteration_array[i][6] = .05;
        iteration_array[i][7] = 0;
        iteration_array[i][8] = 0;

        // Populate the tube array.
        tube_array[i][0] = x_cm;
        tube_array[i][1] = y_cm;
        tube_array[i][2] = angle;
        tube_array[i][3] = length;

        // (P) O------------O (Q) ( Px < Qx )
        // populate the initial endpoint arrays.
        p_initial[i][0] = x_cm - fabs(length*cos(angle));
        p_initial[i][2] = x_cm + fabs(length*cos(angle));
        if(tan(angle) >= 0)
        {
            p_initial[i][1] = y_cm - fabs(length*sin(angle));
            p_initial[i][3] = y_cm + fabs(length*sin(angle));
        }
        else
        {
            p_initial[i][1] = y_cm + fabs(length*sin(angle));
            p_initial[i][3] = y_cm - fabs(length*sin(angle));
        }
    }
}


void GUI()
{
    static int two=2;
    static char Name[tube_count][50];
    
    for(int i=0;i<tube_count;i++)
    {
        sprintf(Name[i],"tube %i",i);
        DefineGraphN_RxR(Name[i],p_initial[i],&two,NULL);
    }

    StartMenu("Nanotubes",1);
    DefineDouble("Energy",&energy);
    DefineGraph(curve2d_,"Tubes");
    DefineFunction("Initialize",&Initialize);
    DefineBool("Pause",&pause);
    DefineBool("Single Step",&sstep);
    DefineBool("Done",&done);
    EndMenu();
}


int main ()
{
    Initialize();
    energy = calculateEnergy();
    GUI();

    while (!done)
    {
        Events(1);
        DrawGraphs();
        if (!(pause||sstep))
        {
            sstep=0;
            energy = calculateEnergy();
            compressNetwork();
            //relaxNetwork();
        } 
        else
            sleep(1);
    }
}
