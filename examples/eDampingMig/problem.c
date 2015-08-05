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

double tmax = 1.e5;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->dt 			= 0.012;		// initial timestep.
	r->integrator	= REB_INTEGRATOR_WHFAST;

	// Setup callback function for velocity dependent forces.
	//r->post_timestep_modifications 	= reboundxf_modify_elements;
	r->force_is_velocity_dependent = 1;
	r-> G = 4*M_PI*M_PI;
	// Setup callback function for outputs.
	r->heartbeat	= heartbeat;
	//r->usleep		= 1;		// Slow down integration (for visualization only)

	struct reb_particle p; 
	p.m  	= 1.;	
	reb_add(r, p); 

	struct reb_particle p1 = reb_tools_init_orbit2d(r->G, 1., 1.e-8, 1.0, 0.0, 0., 0.);	
	struct reb_particle p2 = reb_tools_init_orbit2d(r->G, 1., 1.e-5, pow(2.1,(2./3.)), 0., 0., 0.);	
	reb_add(r,p1);
	reb_add(r,p2);

	double tau_e[3] = {0.,1.e5,100.};
	double tau_a[3] = {0., 1.e10, 1.e7};
	double tau_po[3] = {0.,1.e5,100.};

	/*reboundxf_set_migration(tau_a, 3);
	reboundxf_set_e_damping_p(1.);
	reboundxf_set_e_damping(tau_e,3);
	reboundxf_set_peri_precession(tau_po,3);
	*/
	reb_move_to_com(r);

	// Delete previous output
	system("rm -v r.txt");	

	// Do the integration
	reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* const r){
	// Output some information to the screen every 100th timestep
	if(reb_output_check(r, 100.*r->dt)){
		struct reb_orbit o1 = reb_tools_p2orbit(r->G, r->particles[1], r->particles[0]);
		struct reb_orbit o2 = reb_tools_p2orbit(r->G, r->particles[2], r->particles[0]);
		printf("%f\t%f\t%f\n", r->t, o1.a, o2.a);
		//reb_output_timing(r, tmax);
	}
	// Output the particle position to a file every timestep.
	/*const struct reb_particle* const particles = r->particles;
	FILE* f = fopen("r.txt","a");
	fprintf(f,"%e\t%e\t%e\n",r->t,particles[0].x, particles[1].vx);
	fclose(f);*/
}
