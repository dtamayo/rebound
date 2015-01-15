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

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double* tau_i;  /**< Inclination damping timescale in years for all particles */

void problem_migration_forces();

double starmass = 0.55;
const int Nplanets = 3;
const double dr_thresh = 5; // if deltar changes by more than this thresh, exit (either because ae > thresh, or delta a > thresh between checks)
double a[4] = {0.,0.,0.,0.};
double deltacheck = 1000.; // number of years between checks for jumps and updating of masses
const double mearth = 3.e-6; // mass of Earth in solar masses
double mthresh;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt  		= 0.1;				// in years.  Innermost would have P ~25 yrs for 1 solar mass star.  IAS15 is adaptive anyway
	tmax		= 1e8;
	G		  	= 4*M_PI*M_PI;		// units of years, AU and solar masses.


	double it = input_get_double(argc,argv,"it",0);
	mthresh = input_get_double(argc,argv,"mthresh", 300.); // cut off rise in mass at this mass (in earth masses)
	mthresh *= mearth; // get it into solar masses
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL

	init_boxwidth(1000); 			// Init box with width 1000 astronomical units (max a = 80 AU to start)

	if (input_check_restart(argc,argv)!=1){
		printf("Must provide restart.bin file\n");
		exit(1);
	}
    
	tau_a = calloc(sizeof(double),N);
	tau_e = calloc(sizeof(double),N);
	tau_i = calloc(sizeof(double),N);
	problem_additional_forces = problem_migration_forces; 	//Set function pointer to add dissipative forces.

	for (int i=1;i<=Nplanets;i++){			// initialize the planets first, then initialize star so that center of mass is fixed at 0
		tau_e[i] = 1.e5;
		tau_i[i] = 1.e5;
		if(i==1){ // only set the inner planet to have an outward migration
			tau_a[i] = -1.5e7;
		}
		else{
			tau_a[i] = 0;
		}
	}

	tools_move_to_center_of_momentum();
    
    output_prepare_directory();

    output_append_orbits("ini.txt");
	
#ifdef LIBPNG
    system("mkdir pngs");
#endif //LIBPNG
    
    output_double("Iteration number", it);
    printf("N = %d\n", N);
}

void problem_migration_forces(){
	struct particle com = particles[0]; // calculate migration forces with respect to center of mass;
	for(int i=1;i<N;i++){
		if (tau_e[i]!=0||tau_a[i]!=0){
			struct particle* p = &(particles[i]);
			const double dvx = p->vx-com.vx;
			const double dvy = p->vy-com.vy;
			const double dvz = p->vz-com.vz;

			if (tau_a[i]!=0){ 	// Migration
				p->ax -=  dvx/(2.*tau_a[i]);
				p->ay -=  dvy/(2.*tau_a[i]);
				p->az -=  dvz/(2.*tau_a[i]);
			}
			if (tau_e[i]!=0 || tau_i[i]!=0){ 	// Need h and e vectors if there is either ecc or inc damping
				const double mu = G*(com.m + p->m);
				const double dx = p->x-com.x;
				const double dy = p->y-com.y;
				const double dz = p->z-com.z;

				const double hx = dy*dvz - dz*dvy;
				const double hy = dz*dvx - dx*dvz;
				const double hz = dx*dvy - dy*dvx;
				const double h = sqrt ( hx*hx + hy*hy + hz*hz );
				const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
				const double r = sqrt ( dx*dx + dy*dy + dz*dz );
				const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
				const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
				const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
				const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
				const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity

				if (tau_e[i]!=0){					// Eccentricity damping
					const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
					const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
					const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/tau_e[i]/1.5;
					p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
					p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
					p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
				}
				if (tau_i[i]!=0){	// Inclination damping
					p->az += -2*dvz/tau_i[i];
					const double prefac = (hx*hx + hy*hy)/h/h/tau_i[i];
					p->ax += prefac*dvx;
					p->ay += prefac*dvy;
					p->az += prefac*dvz;
				}
			}

		}
		com = tools_get_center_of_mass(com,particles[i]);
	}
}

void problem_inloop(){
}

void append_orbits(char *filename){
	FILE* of = fopen(filename,"a");
	if (of==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}
	struct particle com = particles[0];
	for (int i=1;i<N;i++){
		struct orbit o = tools_p2orbit(particles[i],com);
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f, particles[i].m);
		if(o.a > 1e4 || o.a < 0){
			printf("\nPlanet %d became unbound\n", i);
			char* ofname = "../ej.txt";
			FILE* of2 = fopen(ofname,"a");
			if (of2==NULL){
				printf("\n\nError while opening file '%s'.\n",filename);
				return;
			}
			fprintf(of2,"%d\n",i);
			fclose(of2);
			exit_simulation=1;
		} // quit if one of the semimajor axes jumps beyond 10000 AU or becomes hyperbolic

		com = tools_get_center_of_mass(com,particles[i]);
	}
	fclose(of);
}

void check_jumps(){
	struct particle com = particles[0];
	double deltaM;
	for (int i=1;i<N;i++){
		struct orbit o = tools_p2orbit(particles[i],com);

		if(o.a*o.e > dr_thresh){ // quit if one of the eccentricity excursions grows beyond thresh (widths are ~5AU so set thresh ~ 5)
			printf("\nPlanet %d had a*e > %f AU\n",i,dr_thresh);
			exit_simulation=1;
		}

		if(a[i] == 0.){	a[i] = o.a;}
		if(abs(o.a - a[i]) > dr_thresh){
			printf("\nPlanet %d had semimajor axis jump by more than %f AU", i, dr_thresh);
			exit_simulation=1;
		}
		a[i] = o.a; // update a for next check

		if(i == N-1){ // calculate the libration timescale for the outer guy since it will have the longest orbital period and thus longest timescale
			if(particles[i].m < mthresh){
				double tlib = 0.078*pow(particles[i].m/starmass, -2./3.)*o.P;
				deltaM = particles[i].m * deltacheck / 10. / tlib;
			}
			else{
				deltaM = 0.;
			}

		}

		com = tools_get_center_of_mass(com,particles[i]);
	}

	for (int i=1;i<N;i++){
		particles[i].m += deltaM;
	}

	tools_move_to_center_of_momentum();

	return;
}

void problem_output(){
	if (output_check(250.)){
	    output_timing();
		append_orbits("orbits.txt");
#ifdef LIBPNG
        //output_png("pngs/");
#endif
	}
	if (output_check(deltacheck)){
		check_jumps();
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
    output_binary("restart.bin");
}
