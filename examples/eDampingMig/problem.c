/**
 * Velocity dependent drag force
 *
 * This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "libreboundxf.h"
#include "tools.h"

void heartbeat(struct reb_simulation* const r);

double tmax = 40.;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->dt 			= 1e-4;		// initial timestep.
	r->integrator	= REB_INTEGRATOR_WHFAST;

	// Setup callback function for velocity dependent forces.
	r->additional_forces 	= reboundxf_forces;
	r->force_is_velocity_dependent = 1;
	// Setup callback function for outputs.
	r->heartbeat	= heartbeat;
	//r->usleep		= 1;		// Slow down integration (for visualization only)

	double tmax = 1.e4;
	double e0 = 0.1;
	double ainner = 1.;
	double aouter = 10.;

	struct reb_particle p; 
	p.m  	= 1.;	
	reb_add(r, p); 

	struct reb_particle p1 = reb_tools_init_orbit2d(r->G, 1., 1.e-6, ainner, e0, 0., 0.);	
	struct reb_particle p2 = reb_tools_init_orbit2d(r->G, 1., 1.e-6, aouter, e0, 0., 0.);	
	reb_add(r,p1);
	reb_add(r,p2);

	reb_move_to_com(r);
	
	// Delete previous output
	system("rm -v r.txt");	

	// Do the integration
	reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* const r){
	// Output some information to the screen every 100th timestep
	if(reb_output_check(r, 100.*r->dt)){
		reb_output_timing(r, tmax);
	}
	// Output the particle position to a file every timestep.
	const struct reb_particle* const particles = r->particles;
	FILE* f = fopen("r.txt","a");
	fprintf(f,"%e\t%e\t%e\n",r->t,particles[0].x, particles[1].vx);
	fclose(f);
}
