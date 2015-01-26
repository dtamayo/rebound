/**
 * @file 	libias15.c
 * @brief 	Declarations and functions for shared library access to the IAS15 integrator.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Dave Spiegel <dave@ias.edu>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Dave Spiegel
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
#include <string.h>
#include "particle.h"
#include "integrator.h"
#include "gravity.h"
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "tools.h"

// Default values of parameters and constants
double dt 	= 0.01;	
double t 	= 0;
double tmax	= 0;
double G 	= 1;
double softening = 0;
extern int Nmax;	

extern int N3allocated;
extern double dt_last_success;
extern double* at;
extern double* x0;
extern double* v0;
extern double* a0;
extern double* csx;
extern double* csv;
extern double* g[7];
extern double* b[7];
extern double* e[7];
extern double* br[7];
extern double* er[7];

// pointers for damping timescales
double *tau_a, *tau_e, *tau_i;

// Function pointer to additional forces
void (*problem_additional_forces) () = NULL;

// Particle structure
struct particle* particles = NULL;

// Particle getter/setter methods.
void setp(struct particle* _p){
	free(particles);
	particles = malloc(sizeof(struct particle)*N);
	memcpy(particles,_p,sizeof(struct particle)*N);
}
struct particle particle_get(int i){
	return particles[i];
}
struct particle* particles_get(){
	return particles;
}
void set_additional_forces(void (* _cb)()){
	problem_additional_forces = _cb;
}

// Integrate for 1 step
void step(){ 
	if (N<=0){
		fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
		return;
	}
	integrator_part1();
	gravity_calculate_acceleration();
	if (problem_additional_forces) problem_additional_forces();

	integrator_part2();
}

// Integrate for 1 step
void reset(){ 
	dt 	= 0.01;	
	t 	= 0;
	tmax	= 0;
	G 	= 1;
	softening = 0;
	N = 0;
	Nmax = 0;
	N_active = -1;
	N_megno = 0;
	free(particles);
	particles = NULL;
	N3allocated = 0;
	dt_last_success = 0;
	for (int l=0;l<7;++l) {
		free(g[l]);
		g[l] = NULL;
		free(b[l]);
		b[l] = NULL;
		free(e[l]);
		e[l] = NULL;
		free(br[l]);
		br[l] = NULL;
		free(er[l]);
		er[l] = NULL;
	}
	free(at);
	at =  NULL;
	free(x0);
	x0 =  NULL;
	free(v0);
	v0 =  NULL;
	free(a0);
	a0 =  NULL;
	free(csx);
	csx=  NULL;
	free(csv);
	csv=  NULL;
	struct timeval tim;
	gettimeofday(&tim, NULL);
	srand ( tim.tv_usec + getpid());
	free(tau_a);
	free(tau_e);
	free(tau_i);
}

// Integrate until t=_tmax
void integrate(double _tmax){
	tmax = _tmax;
	double dt_last_done = dt;
	while(t<tmax){
		if (N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return;
		}
		step();
		if (t+dt>=tmax){
			dt = tmax-t;
		}else{
			dt_last_done = dt;
		}
	}
	dt = dt_last_done;
}
	 
// The following code is needed to leave the original REBOUND files unchanged. 
// Currently the libias15 library only supports an infinitely large box.
// Infinite box size
int root_nx = 1;
int root_ny = 1;
int root_nz = 1;
double boxsize = 0;
double boxsize_x = 0;
double boxsize_y = 0;
double boxsize_z = 0;

int boundaries_particle_is_in_box(struct particle p){
	return 1;
}

// No ghost boxes for now.
int nghostx = 0;	
int nghosty = 0;
int nghostz = 0;

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = 0;
	gb.shifty = 0;
	gb.shiftz = 0;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}

void additional_forces(){
	struct particle com = particles[0]; // calculate migration forces with respect to center of mass;
	for(int i=1;i<N;i++){
		if (tau_e[i] != 0. || tau_a[i] != 0.){
			struct particle* p = &(particles[i]);
			const double dvx = p->vx-com.vx;
			const double dvy = p->vy-com.vy;
			const double dvz = p->vz-com.vz;

			if (tau_a[i] != 0.){
				p->ax -=  dvx/(2.*tau_a[i]);
				p->ay -=  dvy/(2.*tau_a[i]);
				p->az -=  dvz/(2.*tau_a[i]);
			}
			if (tau_e[i] != 0. || tau_i[i] != 0.){ 	// Need h and e vectors if there is either ecc or inc damping
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

				if (tau_e[i] != 0.){					// Eccentricity damping
					const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
					const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
					const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/tau_e[i]/1.5;
					p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
					p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
					p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
				}
				if (tau_i[i] != 0.){	// Inclination damping
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

void init_damping_forces(){
	set_additional_forces(additional_forces);
	tau_a = calloc(sizeof(double),N);
	tau_e = calloc(sizeof(double),N);
	tau_i = calloc(sizeof(double),N);
	return;
}

void add_migration(double *_tau_a)
{
	for(int i=0; i<N; ++i){
		tau_a[i] = _tau_a[i];
	}
}

void add_e_damping(double *_tau_e)
{
	for(int i=0; i<N; ++i){
		tau_e[i] = _tau_e[i];
	}
}

void add_i_damping(double *_tau_i)
{
	for(int i=0; i<N; ++i){
		tau_i[i] = _tau_i[i];
	}
}
