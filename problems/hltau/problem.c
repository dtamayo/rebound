/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the IAS15 integrator
 * to integrate the outer planets of the solar system. The initial 
 * conditions are taken from Applegate et al 1986. Pluto is a test
 * particle. This example is a good starting point for any long term orbit
 * integrations.
 *
 * You probably want to turn off the visualization for any serious runs.
 * Just go to the makefile and set `OPENGL=0`. 
 *
 * The example also works with the Wisdom-Holman symplectic integrator.
 * Simply change the integrator to `integrator_wh.c` in the Makefile.
 * 
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein, Shangfei Liu, Dave Spiegel
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"

#include <string.h>

const double mjup = 9.54e-4; // solar masses

const double a[9] = // in AU
{
	0.,		// placeholder for star (have to calc COM)
  	9.,   	// gap 1
  	23.4, 	// gap 2
  	30.,  	// gap 3
  	36.6, 	// gap 4
  	46.2, 	// gap 5
  	55.2, 	// gap 6
  	66.,  	// gap 7
  	77.4  	// gap 8
};

const double massfac = 8.e-2;

const double m[9] = // in solar masses                                                                               
{
	0.55,				// star
    massfac*mjup,   	// gap 1                                                                                                  
    massfac*mjup,   	// gap 2                                                                                                   
   	massfac*mjup,   	// gap 3  
    massfac*mjup,   	// gap 4  
    massfac*mjup,   	// gap 5  
    massfac*mjup,   	// gap 6  
    massfac*mjup,  		// gap 7  
    massfac*mjup,   	// gap 8                                                                                                      
};  

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
    dt  		= 0.5;				// in years.  Innermost would have P ~25 yrs for 1 solar mass star.  IAS15 is adaptive anyway
	tmax		= 1e6;		
	G		  	= 4*M_PI*M_PI;		// units of years, AU and solar masses.
	
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(200); 			// Init box with width 1000 astronomical units (max a = 80 AU to start)

	struct particle p[9]; 			// also include star
	
	// Initial conditions
	
	struct timeval tim;
	gettimeofday(&tim, NULL);
	srand(tim.tv_usec);
	
	double mx[3] = {0.,0.,0.};
	double mxdot[3] = {0.,0.,0.}; // to hold the total of m_i * x_i and m_i * v_i for all the planets, so that can set star's ini. conds s.t. com is fixed at 0
	 
	for (int i=1;i<=8;i++){			// initialize the planets first, then initialize star so that center of mass is fixed at 0
		double phi = (float)rand()/RAND_MAX*2*M_PI; // choose random azimuthal angle [0,2PI]
		
		p[i].x  = a[i]*cos(phi); 			p[i].y  = a[i]*sin(phi);	 		p[i].z  = 0.;
		mx[0] += m[i]*p[i].x;				mx[1] += m[i]*p[i].y;				mx[2] += m[i]*p[i].z;
		
		double vkep = sqrt(G*m[0]/a[i]);
		p[i].vx = -vkep*sin(phi); 			p[i].vy = vkep*cos(phi);	 		p[i].vz = 0.;
		mxdot[0] += m[i]*p[i].vx;			mxdot[1] += m[i]*p[i].vy;			mxdot[2] += m[i]*p[i].vz;
		
		p[i].ax = 0; 						p[i].ay = 0; 						p[i].az = 0;
		p[i].m  = m[i];
	}
	
	// Now initialize star such that center of mass is fixed at origin (so total mx & mxdot = 0)
	
	p[0].x = -mx[0]/m[0];					p[0].y = -mx[1]/m[0];				p[0].z = -mx[2]/m[0];
	p[0].vx = -mxdot[0]/m[0];				p[0].vy = -mxdot[1]/m[0];			p[0].vz = -mxdot[2]/m[0];
	p[0].ax = 0;							p[0].ay = 0;						p[0].az = 0;
	p[0].m = m[0];
	
	for (int i=0;i<=8;++i){	particles_add(p[i]); }
	
	printf("Sum of m_i * x_i at t=0 is (%f, %f, %f)\n", mx[0]+p[0].m*p[0].x, mx[1]+p[0].m*p[0].y, mx[2]+p[0].m*p[0].z);
	printf("Sum of m_i * v_i at t=0 is (%f, %f, %f)\n", mx[0]+p[0].m*p[0].x, mx[1]+p[0].m*p[0].y, mx[2]+p[0].m*p[0].z);
	
	tools_move_to_center_of_momentum();

	system("rm -f orbits.txt");
	
	/*int x = 2;
	char* test = "pngs/img0000.png";
	char y = (char)x;
	printf("**** %c ******\n", y);*/
}

void problem_inloop(){
}

void problem_output(){
	
	if (output_check(1000.*dt)){
		output_timing();
	}
	if (output_check(100.)){ 	// output heliocentric orbital elements every 10000 years
		output_append_orbits("orbits.txt");
	}
}

void problem_finish(){
}
