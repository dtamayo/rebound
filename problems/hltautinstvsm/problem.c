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
void check_jumps();
const double mjup = 9.54e-4; // solar masses

const int Nplanets = 3;
const double dr_thresh = 5; // if deltar changes by more than this thresh, exit (either because ae > thresh, or delta a > thresh between checks)
double a[4];
double it, delta, massfac,a3,daOvera,q,P2;

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
    
	a3 = input_get_double(argc,argv,"a3",66.);
    daOvera = input_get_double(argc,argv,"daOvera", 0.45);
    q = input_get_double(argc,argv,"q",1.);
    massfac = input_get_double(argc,argv,"mass",1); // in jup masses
    double starmass = input_get_double(argc,argv,"starmass",0.55); // in solar masses
    double sigma_e = input_get_double(argc,argv,"sigmae",1.e-6); 			// scale of eccentricity rayleigh distribution
    double sigma_i = input_get_double(argc,argv,"sigmai",sigma_e); 			// scale of inclination rayleigh distribution.  By default use same as sigma_e value
    it = input_get_double(argc,argv,"it",0);							// iteration (running several realizations of same set of parameters)
    delta = daOvera/pow((2./3.*massfac*mjup/starmass), (1./3.));

    double a2 = a3/(1 + daOvera);
    double a1 = a2/(1 + q*daOvera);

    dt  		= 0.01;				// in years.  Innermost would have P ~25 yrs for 1 solar mass star.  IAS15 is adaptive anyway
    G		  	= 4*M_PI*M_PI;		// units of years, AU and solar masses.
    P2			= sqrt(a2*a2*a2/G/starmass); // period of 2nd planet (want to refer to inner planet, like da/a_inner)
    tmax		= 1.e6*P2;

#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(200); 			// Using no boundary conditions now, so particles aren't removed beyond this dist.  This is just for graphics box size

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

	a[1] = a1;	a[2] = a2;	a[3] = a3;

	printf("a1, a2, a3 = %f, %f, %f\n M = %f delta = %f, da/a = %f", a1,a2,a3, massfac, delta, (a3 - a2)/a2);

    for (int i=1;i<=Nplanets;i++){			// initialize the planets first, then initialize star so that center of mass is fixed at 0
		double Omega = (float)rand()/RAND_MAX*2*M_PI;;
		double omega = (float)rand()/RAND_MAX*2*M_PI;;
		double f = (float)rand()/RAND_MAX*2*M_PI;;
		double inc = tools_rayleigh(sigma_i);
		double e = tools_rayleigh(sigma_e);

		struct particle p = tools_init_orbit3d(starmass, massfac*mjup, a[i], e, inc, Omega, omega, f);
		particles_add(p);
	}
		
    tools_move_to_center_of_momentum();
    output_prepare_directory();
    output_binary("ini.bin");
	
#ifdef LIBPNG
    system("mkdir pngs");
#endif //LIBPNG
}


void problem_inloop(){
}

void append_orbits(){
	struct particle com = particles[0];
	for (int i=1;i<N;i++){
		char filename[20];
		sprintf(filename, "orbit%d.txt", i);
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

void check_jumps(){
	struct particle com = particles[0];
	for (int i=1;i<N;i++){
		struct orbit o = tools_p2orbit(particles[i],com);
		if(o.a*o.e > dr_thresh){ // quit if one of the eccentricity excursions grows beyond thresh (widths are ~5AU so set thresh ~ 5)
			printf("\nPlanet %d had a*e > %f AU\n",i,dr_thresh);
			exit_simulation=1;
		}

		if(abs(o.a - a[i]) > dr_thresh){
			printf("\nPlanet %d had semimajor axis jump by more than %f AU", i, dr_thresh);
			exit_simulation=1;
		}
		a[i] = o.a; // update a for next check

		com = tools_get_center_of_mass(com,particles[i]);
	}

	return;
}

void problem_output(){
	if (output_check(1000.*dt)){
		check_jumps();
#ifdef LIBPNG
        //output_png("pngs/");
#endif
	}
	if (output_check(tmax/1.e3)){
		append_orbits();
		output_timing();
	}
}

void problem_finish(){
	char filename[40];
	sprintf(filename, "../eos_a3_%0.2f_da_%0.2f.txt", a3, daOvera);
	FILE* of = fopen(filename, "a");
	if (of==NULL){
	    printf("\n\nError while opening file '%s'.\n", filename);
	    return;
	}
	fprintf(of, "%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", massfac, delta, q, it, t);
	fclose(of);
}
