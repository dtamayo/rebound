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
#include <sys/time.h>
#include "main.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"
#include "problem.h"
#include "input.h"

void append_orbits();
const double mjup = 9.54e-4; // solar masses

double outputcadence = 500.;
double Nplanets = 5;
double a[] = // in AU
{
	0.,		// placeholder for star (have to calc COM)
  	9.,   	// gap 1
  	23., 	// gap 2
  	46., 	// gap 5
  	55.,	// gap 6
  	66.  	// gap 7
};

double E0;
double* L0;
double ethresh=0.1;

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
    
    double massfac = input_get_double(argc,argv,"mass",1); // in jup masses
    double starmass = input_get_double(argc,argv,"starmass",0.55); // in solar masses
    double sigma_e = input_get_double(argc,argv,"sigmae",0.01); 			// scale of eccentricity rayleigh distribution
    double sigma_i = input_get_double(argc,argv,"sigmai",sigma_e); 			// scale of inclination rayleigh distribution.  By default use same as sigma_e value
    double it = input_get_double(argc,argv,"it",0);							// iteration (running several realizations of same set of parameters)

    dt  		= 0.01;				// in years.  Innermost would have P ~25 yrs for 1 solar mass star.  IAS15 is adaptive anyway
	tmax		= 1e6;
	G		  	= 4*M_PI*M_PI;		// units of years, AU and solar masses.
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(400); 			// Using no boundary conditions now, so particles aren't removed beyond this dist.  This is just for graphics box size

	// Initial conditions
	
	struct timeval tim;
	gettimeofday(&tim, NULL);
	srand(tim.tv_usec);
	
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = starmass;			// This is a sub-solar mass star
	particles_add(star);

    for (int i=1;i<=Nplanets;i++){			// initialize the planets first, then initialize star so that center of mass is fixed at 0
		double Omega = (float)rand()/RAND_MAX*2*M_PI;;
		double omega = (float)rand()/RAND_MAX*2*M_PI;;
		double f = (float)rand()/RAND_MAX*2*M_PI;
		double inc = tools_rayleigh(sigma_i);
		double e = tools_rayleigh(sigma_e);
		struct particle p = tools_init_orbit3d(starmass, massfac*mjup, a[i], e, inc, Omega, omega, f);
		particles_add(p);
	}
		
    particles_assign_IDs(); // keep track of which particle is which
    boundaries_track_conservation(); // keep track of conserved quantities

	tools_move_to_center_of_momentum();
    
	output_prepare_directory();

	E0 = tools_get_total_E();
	L0 = calloc(sizeof(double),3);
	tools_get_total_L_vec(L0);

    output_append_orbits("ini.txt");
	
#ifdef LIBPNG
    system("mkdir pngs");
#endif //LIBPNG
    
    output_double("Planet masses [mjup]",massfac);
    output_double("Star's mass [msolar]", starmass);
    output_double("Eccentricity scale parameter", sigma_e);
    output_double("Inclination scale parameter", sigma_i);
    output_double("Iteration number", it);
    system("cat config.log");
}

void problem_inloop(){
}

void append_orbits(){
	struct particle com = particles[0];
	for (int i=1;i<N;i++){
		char filename[20];
		sprintf(filename, "orbit%d.txt", particle_IDs[i]);
		FILE* of = fopen(filename,"a");
		if (of==NULL){
			printf("\n\nError while opening file '%s'.\n",filename);
			return;
		}
		struct orbit o = tools_p2orbit(particles[i],com);
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
		fclose(of);
		if(1.-o.e < ethresh){
			char efile[20];
			sprintf(efile, "restart%.1e.bin", ethresh);
			output_binary(efile);
			of = fopen("restarts.txt", "a");
			if (of==NULL){
				printf("\n\nError while opening file '%s'.\n",filename);
				return;
			}
			fprintf(of, "%e\t%d\t%e\n", t, particle_IDs[i], 1.- o.e);
			ethresh /= 10.;
			fclose(of);
		}
		com = tools_get_center_of_mass(com,particles[i]);
	}


}

void problem_output(){
	if (output_check(outputcadence)){
		output_timing();
		append_orbits();
#ifdef LIBPNG
        //output_png("pngs/");
#endif
	}
}

void problem_finish(){
	char* eos = "eos.txt"; // end of simulation time
	FILE* of = fopen(eos, "w");
	if (of==NULL){
	    printf("\n\nError while opening file '%s'.\n", eos);
	    return;
	}
	fprintf(of, "%.3e", t);
	fclose(of);

	if(t > tmax-outputcadence){
		char* end = "endcondition.txt"; // end of simulation time
		of = fopen(end, "w");
		if (of==NULL){
		    printf("\n\nError while opening file '%s'.\n", end);
		    return;
		}
		fprintf(of, "1\n"); // simulation ended because we reached tmax
		fclose(of);
	}

    char* finalstate = "finalstate.txt"; // orb elements of particles at end of simulation
    of = fopen(finalstate, "w");
    if (of==NULL){
    	printf("\n\nError while opening file '%s'.\n", finalstate);
    	return;
    }
    struct particle com = particles[0];
    for (int i=1;i<N;i++){
    	struct orbit o = tools_p2orbit(particles[i],com);
    	fprintf(of,"%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,particle_IDs[i],o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
    	com = tools_get_center_of_mass(com,particles[i]);
    }
    fclose(of);

    char* cons = "conservation.txt";
    of = fopen(cons, "w");
    if (of==NULL){
     	printf("\n\nError while opening file '%s'.\n", cons);
       	return;
    }
    double dE = tools_get_total_E() + Eadj - E0; // take final energy, add the adjustments from ejections and recenterings to new COMs and subtract initial E
    fprintf(of, "%.15e\n", dE/E0);

    double* Lf = calloc(sizeof(double),3);
    double* dL = calloc(sizeof(double),3);
    tools_get_total_L_vec(Lf);
    for(int j=0;j<3;j++){
    	dL[j] = Lf[j] + Ladj[j] - L0[j];
    	fprintf(of, "%.15e\n", dL[j]/L0[j]);
    }
    fclose(of);

    output_binary("restart.bin");
}
