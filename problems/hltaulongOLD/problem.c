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

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
    
    double it = input_get_double(argc,argv,"it",0);							// iteration (running several realizations of same set of parameters)

    dt  		= 0.01;				// in years.  Innermost would have P ~25 yrs for 1 solar mass star.  IAS15 is adaptive anyway
	tmax		= 3.2e6;
	G		  	= 4*M_PI*M_PI;		// units of years, AU and solar masses.
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(200); 			// Using no boundary conditions now, so particles aren't removed beyond this dist.  This is just for graphics box size

	if (input_check_restart(argc,argv)!=1){
		printf("Must provide restart.bin file\n");
		exit(1);
	}
	// Initial conditions
	
	particles_assign_IDs(); // keep track of which particle is which

    tools_move_to_center_of_momentum();
    
    output_prepare_directory();

    FILE* of = fopen("EandL.txt","a");
    if (of==NULL){
    	printf("\n\nError while opening EandL file");
    	exit(1);
    }

    fprintf(of, "%.15e\t%.15e\n", tools_get_total_E(), tools_get_total_L());
	
#ifdef LIBPNG
    system("mkdir pngs");
#endif //LIBPNG
    
    output_double("Iteration number", it);
    system("cat config.log");
}

void problem_inloop(){
}

void append_orbits(){
	struct particle com = particles[0];
	for (int i=1;i<N;i++){
		char filename[50];
		sprintf(filename, "orbit%d.txt", particle_IDs[i]);
		FILE* of = fopen(filename,"a");
		if (of==NULL){
			printf("\n\nError while opening file '%s'.\n",filename);
			return;
		}
		struct orbit o = tools_p2orbit(particles[i],com);
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
		fclose(of);
		com = tools_get_center_of_mass(com,particles[i]);
	}
}

void problem_output(){
	if (output_check(250.)){
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

    char* finalstate = "finalstate.txt"; // orb elements of particles at end of simulation
    FILE* of2 = fopen(finalstate, "w");
    if (of2==NULL){
    	printf("\n\nError while opening file '%s'.\n", eos);
    	return;
    }
    struct particle com = particles[0];
    for (int i=1;i<N;i++){
    	struct orbit o = tools_p2orbit(particles[i],com);
    	fprintf(of2,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
    	com = tools_get_center_of_mass(com,particles[i]);
    }
    fclose(of2);
    output_binary("restart.bin");

    FILE* of3 = fopen("EandL.txt","a");
    if (of3==NULL){
       	printf("\n\nError while opening EandL file");
    }

    fprintf(of3, "%.15e\t%.15e\n", tools_get_total_E(), tools_get_total_L());
    fclose(of3);
}
