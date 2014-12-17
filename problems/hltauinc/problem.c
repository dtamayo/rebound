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

double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
void problem_migration_forces();
void append_orbits();
const double mjup = 9.54e-4; // solar masses

const double a[6] = // in AU
{
	0.,		// placeholder for star (have to calc COM)
  	9.,   	// gap 1
  	23.4, 	// gap 2
  	50.7, 	// gap 5 average of 2 gaps in big one
  	66.,  	// gap 7
  	77.4  	// gap 8
};

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
    
    double massfac = input_get_double(argc,argv,"mass",1); // in jup masses
    double starmass = input_get_double(argc,argv,"starmass",0.55); // in solar masses
    double taue = input_get_double(argc,argv,"taue",10000);
    double k = input_get_double(argc,argv,"k",100);
    double sigma_e = input_get_double(argc,argv,"sigmae",0.01); 			// scale of eccentricity rayleigh distribution
    double sigma_i = input_get_double(argc,argv,"sigmai",0.01); 			// scale of inclination rayleigh distribution

    dt  		= 0.1;				// in years.  Innermost would have P ~25 yrs for 1 solar mass star.  IAS15 is adaptive anyway
	tmax		= 1e6;		
	G		  	= 4*M_PI*M_PI;		// units of years, AU and solar masses.
	
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(400); 			// Init box with width 1000 astronomical units (max a = 80 AU to start)

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

    for (int i=1;i<=5;i++){			// initialize the planets first, then initialize star so that center of mass is fixed at 0
		double Omega = (float)rand()/RAND_MAX*2*M_PI;;
		double omega = (float)rand()/RAND_MAX*2*M_PI;;
		double f = (float)rand()/RAND_MAX*2*M_PI;;
		double inc = tools_rayleigh(sigma_i);
		double e = tools_rayleigh(sigma_e);

		struct particle p = tools_init_orbit3d(starmass, massfac*mjup, a[i], e, inc, Omega, omega, f);
		particles_add(p);
	}
		
    tau_a = calloc(sizeof(double),N);
	tau_e = calloc(sizeof(double),N);
    
    for(int j=1;j<=5;++j){
        tau_a[j] = taue*k;
        tau_e[j] = taue;
    }
    
	problem_additional_forces = problem_migration_forces; 	//Set function pointer to add dissipative forces.
    
	tools_move_to_center_of_momentum();
    
    output_prepare_directory();

    output_append_orbits("ini.txt");
	
#ifdef LIBPNG
    system("mkdir pngs");
#endif //LIBPNG
    
    output_double("semimajor axis decay timescale [yrs]", taue*k);
    output_double("eccentricity decay timescale [yrs]", taue);
    output_double("K ratio between tau_a and tau_e", k);
    output_double("Planet masses [mjup]",massfac);
    output_double("Star's mass [msolar]", starmass);
    output_double("Eccentricity scale parameter", sigma_e);
    output_double("Inclination scale parameter", sigma_i);
    system("cat config.log");
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
			if (tau_e[i]!=0){ 	// Eccentricity damping
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
				const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
				p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
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
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
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
			char* eos = "eos.txt"; // end of simulation time
			FILE* of3 = fopen(eos, "w");
			if (of3==NULL){
			    printf("\n\nError while opening file '%s'.\n", eos);
			    return;
			}
			fprintf(of3, "%.3e", t);
			fflush(stdout);
			fclose(of3);
			exit_simulation=1; } // quit if one of the semimajor axes jumps beyond 10000 AU or becomes hyperbolic
		com = tools_get_center_of_mass(com,particles[i]);
	}
	fclose(of);
}
void problem_output(){
	if (output_check(1000*dt)){
		output_timing();
		append_orbits("orbits.txt");
		tools_move_to_center_of_momentum();
#ifdef LIBPNG
        //output_png("pngs/");
#endif
	}
}

void problem_finish(){
    char* finalstate = "finalstate.txt"; // xyz of particles at end of simulation
    output_ascii(finalstate);
    output_binary("restart.bin");
}
