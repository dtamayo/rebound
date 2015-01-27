/**
 * @file 	tools.h
 * @brief 	Tools for creating distributions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#ifndef TOOLS_H
#define TOOLS_H
#include "particle.h"
/**
 * Struct representing a Keplerian orbit.
 */
struct orbit {
	double a;
	double r;	// Radial distance from central object
	double h;	// Angular momentum
	double P;	// Orbital period
	double l;
	double e;
	double inc;
	double Omega; 	// longitude of ascending node
	double omega; 	// argument of perihelion
	double f; 	// true anomaly
};


/**
 * Calculates a random variable in a given range.
 * @param min Minimum value.
 * @param max Maximum value.
 */
double tools_uniform(double min, double max);

/**
 * Calculates a random variable drawn form a powerlaw distribution.
 * @param min Minimum value.
 * @param max Maximum value.
 * @param slop Slope of powerlaw distribution.
 */
double tools_powerlaw(double min, double max, double slope);

/**
 * Calculate a random number with normal distribution.
 * Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
 * @param variance Variance of normal distribution.
 * @return Random number with normal distribution (mean 0). 
 */
double tools_normal(double variance);

/**
 * Calculates a random variable drawn form a Rayleigh distribution.  Calculated as described on Rayleigh distribution wikipedia page
 * @param sigma Scale parameter.
 */
double tools_rayleigh(double sigma);

/**
 * This function sets up a Plummer sphere.
 * @param _N Number of particles in the plummer sphere.
 * @param M Total mass of the cluster.
 * @param R Characteristic radius of the cluster.
 */

void tools_init_plummer(int _N, double M, double R);

/**
 * Initialize a particle on an orbit in the xy plane.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param omega Pericenter of the particle.
 * @param f true anomaly of the particle.
 */
struct particle tools_init_orbit2d(double M, double m, double a, double e, double omega, double f);

/**
 * Initialize a particle on a 3D orbit.  See Fig. 2.13 of Murray & Dermott Solar System Dynamics for diagram.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param i inclination of the particle to the reference plane.
 * @param Omega Longitude of the ascending node of the particle.
 * @param omega argument of pericenter of the particle.
 * @param f true anomaly of the particle.
 */

struct particle tools_init_orbit3d(double M, double m, double a, double e, double i, double Omega, double omega, double f);

/**
 * This function calculated orbital elements for a given particle. 
 * @param p Particle for which the orbit is calculated.
 * @param star Star or central object particle
 * @return Orbital parameters. 
 */
struct orbit tools_p2orbit(struct particle p, struct particle star);

/**
 * Returns the energy associated with a particular particle index, i.e., kinetic energy plus interaction energy with all the other particles
 * @param i index of the particle for which the energy is calculated
 */
double tools_particle_E(int i);

/**
 * This function returns p/a, a conserved quantity in the general 3-body problem.  This value determines whether the bodies can undergo close encounters or not (see Gladman 1993 Eq. 12)
 * @param G Gravitational constant.
 * @param p1 first particle
 * @param p2 second particle
 * @param p3 third particle
 */
double tools_p_over_a(double G, struct particle p1, struct particle p2, struct particle p3);

/**
 * This functions populates L with the center of mass's angular momentum
 * @param L pointer to a vector of 3 doubles (to be populated by function)
 */
void tools_com_L_vec(double* L);

/**
 * This function returns the center of mass's kinetic energy
 */
double tools_com_ke();

/**
 * This function returns the system's total energy
 */
double tools_get_total_E();

/**
 * This function populates the L pointer with the particle with index i's angular momentum components
 * @param i index of the particle for which to calculate the angular momentum
 * @param L Pointer to an array of 3 doubles (to be populated by function)
 */
void tools_particle_L_vec(int i, double* L);

/**
 * This function populates the L pointer with the system's angular momentum components
 * @param L Pointer to an array of 3 doubles (to be populated by function)
 */
void tools_get_total_L_vec(double* L);

/**
 * This function returns the system's total angular momentum
 */
double tools_get_total_L();

/**
 * Move to center of momentum and center of mass frame.
 */
void tools_move_to_center_of_momentum();

/**
 * Returns the center of mass of particle p1 and p2.
 */
struct particle tools_get_center_of_mass(struct particle p1, struct particle p2);

#endif 	// TOOLS_H
