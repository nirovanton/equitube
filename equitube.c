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


//---------------------
#define tube_count 10
int field_size = 10;
double radius = 1e-9;
double vand_const = 1.0;
double spring_const = 1.0;
double tube_array[tube_count][4];
double p_initial[tube_count][4];
double energy;
//---------------------

int done=0,pause=1,sstep=1;


double calculateEnergy()
{
    double potential = 0.0;
    
    for(int k=0; k<tube_count; k++)
    {
        double m1 = tan(tube_array[k][2]);
        double b1 = tube_array[k][1]-m1*tube_array[k][0];
 
        double p1[4];
        p1[0] = tube_array[k][0] - abs(tube_array[k][3]*cos(tube_array[k][2]));
        p1[2] = tube_array[k][0] + abs(tube_array[k][3]*cos(tube_array[k][2]));
        if(m1 >= 0)
        {
            p1[1] = tube_array[k][1] - abs(tube_array[k][3]*sin(tube_array[k][2]));
            p1[3] = tube_array[k][1] + abs(tube_array[k][3]*sin(tube_array[k][2]));
        }
        else
        {
            p1[1] = tube_array[k][1] + abs(tube_array[k][3]*sin(tube_array[k][2]));
            p1[3] = tube_array[k][1] - abs(tube_array[k][3]*sin(tube_array[k][2]));
        }
            
        // spring
        potential += spring_const*(pow((p_initial[k][0]-p1[0]),2)+pow((p_initial[k][1]-p1[1]),2));
        potential += spring_const*(pow((p_initial[k][2]-p1[2]),2)+pow((p_initial[k][3]-p1[3]),2));

        // Vanderwaals
        for(int j=0; j<tube_count; j++)
        {
            double m2 = tan(tube_array[j][2]);
            if(m1 == m2)
                continue;
            double b2 = tube_array[j][1]-m2*tube_array[j][0];            
            
            double p2[4];
            p2[0] = tube_array[j][0] - abs(tube_array[j][3]*cos(tube_array[j][2]));
            p2[2] = tube_array[j][0] + abs(tube_array[j][3]*cos(tube_array[j][2]));
            if(m2 >= 0)
            {
                p2[1] = tube_array[j][1] - abs(tube_array[j][3]*sin(tube_array[j][2]));
                p2[3] = tube_array[j][1] + abs(tube_array[j][3]*sin(tube_array[j][2]));
            }   
            else
            {   
                p2[1] = tube_array[j][1] + abs(tube_array[j][3]*sin(tube_array[j][2]));
                p2[3] = tube_array[j][1] - abs(tube_array[j][3]*sin(tube_array[j][2]));
            } 

            double x_int = (b2-b1)/(m1-m2);
            if(x_int >= p1[0] && x_int >= p2[0] && x_int <= p1[2] && x_int <= p2[2])
                potential += -1*(radius*vand_const)/abs(sin(tube_array[k][2]-tube_array[j][2]));
        }
    }
    
    return potential;
}



void relaxNetwork()
{
    
    // Initialize Increment arrays
    // ---------------------------

    double x_iteration_array[tube_count][3];
    double y_iteration_array[tube_count][3];
    double theta_iteration_array[tube_count][3];
    for(int a=0; a<tube_count; a++)
    {
        // [ dx , increase_counter , decrease_counter ]
        x_iteration_array[a][0] = .05;
        x_iteration_array[a][1] = 0;
        x_iteration_array[a][2] = 0;
        y_iteration_array[a][0] = .05;
        y_iteration_array[a][1] = 0;
        y_iteration_array[a][2] = 0;
        theta_iteration_array[a][0] = .05;
        theta_iteration_array[a][1] = 0;
        theta_iteration_array[a][2] = 0;
    }
    

    // Relax the system
    // ----------------
    int relaxed = 0;
    while(relaxed == 0)
    {   
        double starting_energy = calculateEnergy();
        for(int j = 0; j < tube_count; j++)
        {
            
            // MINIMIZE X
            // ----------
            int x_min = 0;
            // <incease x>
            tube_array[j][0] += x_iteration_array[j][0];
            if(calculateEnergy() > starting_energy)
            {
                tube_array[j][0] += -1*x_iteration_array[j][0];
            }
            else if(calculateEnergy() <= starting_energy)
            {

                starting_energy = calculateEnergy();                
                x_iteration_array[j][1] += 1;
                x_iteration_array[j][2] = 0;
                x_min = 1;
            }
            if(x_min == 0)
            {
                // <increase x>
                tube_array[j][0] += -1*x_iteration_array[j][0];
                if(calculateEnergy() > starting_energy)
                    tube_array[j][0] += x_iteration_array[j][0];
                else if(calculateEnergy() <= starting_energy)
                {
                    starting_energy = calculateEnergy();
                    x_iteration_array[j][2] += 1;
                    x_iteration_array[j][1] = 0;
                    x_min = 1;
                }
                if(x_min == 0)
                {
                    double tmp = x_iteration_array[j][0]/2;
                    x_iteration_array[j][0] = tmp;
                    x_iteration_array[j][1] = 0;
                    x_iteration_array[j][2] = 0;
               }
            }
             
            // MINIMIZE Y
            // ----------
            int y_min = 0;
            // <incease y>
            tube_array[j][1] += y_iteration_array[j][0];
            if(calculateEnergy() > starting_energy)
                tube_array[j][1] += -1*y_iteration_array[j][0];
            if(calculateEnergy() <= starting_energy)
            {
                starting_energy = calculateEnergy();
                y_iteration_array[j][1] += 1;
                y_iteration_array[j][2] = 0;
                y_min = 1;
            }
            if(y_min == 0)
            {
                // <decrease y>
                tube_array[j][1] += -1*y_iteration_array[j][0];
                if(calculateEnergy() > starting_energy)
                    tube_array[j][1] += y_iteration_array[j][0];
                if(calculateEnergy() <= starting_energy)
                {
                    starting_energy = calculateEnergy();
                    y_iteration_array[j][2] += 1;
                    y_iteration_array[j][1] = 0;
                    y_min = 1;
                }
                if(y_min == 0)
                {
                    double tmp = y_iteration_array[j][0]/2;
                    y_iteration_array[j][0] = tmp;
                    y_iteration_array[j][1] = 0;
                    y_iteration_array[j][2] = 0;
                }
            }

            // MINIMIZE THETA
            // --------------
            int theta_min = 0;
            // <incease theta>
            tube_array[j][2] += theta_iteration_array[j][0];
            if(calculateEnergy() > starting_energy)
                tube_array[j][2] += -1*theta_iteration_array[j][0];
            if(calculateEnergy() <= starting_energy)
            {
                starting_energy = calculateEnergy();
                theta_iteration_array[j][1] += 1;
                theta_iteration_array[j][2] = 0;
                theta_min = 1;
            }
            if(theta_min == 0)
            {
                // <decrease theta>
                tube_array[j][2] += -1*theta_iteration_array[j][0];
                if(calculateEnergy() > starting_energy)
                    tube_array[j][2] += theta_iteration_array[j][0];
                if(calculateEnergy() <= starting_energy)
                {
                    starting_energy = calculateEnergy();
                    theta_iteration_array[j][2] += 1;
                    theta_iteration_array[j][1] = 0;
                    theta_min = 1;
                }
                if(theta_min == 0)
                {
                    double tmp = theta_iteration_array[j][0]/2;
                    theta_iteration_array[j][0] = tmp;
                    theta_iteration_array[j][1] = 0;
                    theta_iteration_array[j][2] = 0;
                }
            }
        }
    }
}

// TODO 
// ------------
// void compressNetwork()
// double calculateCurrent()
// ------------

void Initialize(){
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
      
        // Populate the tube array.
        tube_array[i][0] = x_cm;
        tube_array[i][1] = y_cm;
        tube_array[i][2] = angle;
        tube_array[i][3] = length;

        // (P) O------------O (Q) ( Px < Qx )
        // populate the initial endpoint arrays.
        p_initial[i][0] = x_cm - abs(length*cos(angle));
        p_initial[i][2] = x_cm + abs(length*cos(angle));
        if(tan(angle) >= 0)
        {
            p_initial[i][1] = y_cm - abs(length*sin(angle));
            p_initial[i][3] = y_cm + abs(length*sin(angle));
        }
        else
        {
            p_initial[i][1] = y_cm + abs(length*sin(angle));
            p_initial[i][3] = y_cm - abs(length*sin(angle));
        }
    }
}

void GUI(){
  static int two=2;
  static char Name[tube_count][50];
  for (int i=0;i<tube_count;i++){
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
  GUI();
  while (!done){
    Events(1);
    DrawGraphs();
    if (!(pause||sstep)){
      sstep=0;
      energy = calculateEnergy();
      relaxNetwork();
    } else {sleep(1);}
  }
}

