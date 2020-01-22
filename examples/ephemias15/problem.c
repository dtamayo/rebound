/**
 * Highly eccentric orbits
 *
 * This example uses the IAS15 integrator to simulate
 * a very eccentric planetary orbit. The integrator
 * automatically adjusts the timestep so that the pericentre passages
 * resolved with high accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double tmax;
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->G            = 1;        // Gravitational constant
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;

    const double mu = 1.e-3;
    struct reb_particle star = {0}; 
    star.m  = 1.-mu;
    struct reb_particle planet; 
    planet.m  = mu;

    star.x = -mu;
    star.vy = -mu;
    planet.x = 1-mu;
    planet.vy = 1-mu;
    star.r=1./200.;
    planet.r=1./2000.;

    struct reb_particle tp = {0};
    tp.x = 1.5;
    tp.vy = 0.4; // 1/sqrt(1.5); // for circular orbit

    reb_add(r, star); 
    reb_add(r, planet); 
    reb_add(r, tp); 
    r->N_active=2;    
    r->ri_ias15.epsilon_global=0;
    // initial timestep
    tmax            = 10000.;

    reb_integrate(r, tmax);    
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r,tmax/10000.)){        // outputs to the screen
        reb_output_timing(r, tmax);
    }
    // Output the time and the current timestep. Plot it to see how IAS15 automatically reduces the timestep at pericentre. 
    FILE* of = fopen("timestep.txt","a"); 
    fprintf(of,"%e\t%e\t\n",r->t/tmax,r->dt/tmax);
    fclose(of);
}

