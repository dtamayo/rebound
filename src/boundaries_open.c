/**
 * @file 	boundaries.c
 * @brief 	Implementation of open boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * This file implements open boundary conditions. Every particle that leaves 
 * the box is removed from the simulation 
 * 
 * 
 * @section LICENSE
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"
#include "particle.h"
#include "tools.h"

int nghostx = 0;
int nghosty = 0;
int nghostz = 0;

int track_IDs = 0;
int track_conservation = 0;

void boundaries_check(){
	for (int i=0;i<N;i++){
		int removep = 0;
		if(particles[i].x>boxsize_x/2.){
			removep = 1;
		}
		if(particles[i].x<-boxsize_x/2.){
			removep = 1;
		}
		if(particles[i].y>boxsize_y/2.){
			removep = 1;
		}
		if(particles[i].y<-boxsize_y/2.){
			removep = 1;
		}
		if(particles[i].z>boxsize_z/2.){
			removep = 1;
		}
		if(particles[i].z<-boxsize_z/2.){
			removep = 1;
		}
		if (removep==1){
			if (N==1){
				printf("Last particle removed. Exiting.\n");
				exit(0);
			}
#ifndef TREE
			if(track_conservation == 1){
				FILE* of = fopen("ejections.txt","a");
				if (of==NULL){
				   	printf("\n\nError while opening ejections file");
				   	exit(1);
				}

				fprintf(of, "%d\n",particle_IDs[i]);
				fclose(of);

				of = fopen("ej_cons.txt","a");
				if (of==NULL){
				   	printf("\n\nError while opening the ejection conservation file");
				   	exit(1);
				}

				fprintf(of, "***Particle %d was ejected***\n", i);
				double Eprev = tools_get_total_E();
				double Eejected = tools_particle_E(i);

				double* Lprev = calloc(sizeof(double), 3);
				double* Lejected = calloc(sizeof(double), 3);
				tools_get_total_L_vec(Lprev);
				tools_particle_L_vec(i, Lejected);

				particles[i] = particles[N-1];
				if(track_IDs == 1){
					particle_IDs[i] = particle_IDs[N-1];
				}
				i--;
				N--;

				double Eafter = tools_get_total_E();
				fprintf(of, "Frac. E of ej. particle = %.3e\n", Eejected/Eprev);
				fprintf(of, "delta E = %.15e\n", Eafter + Eejected - Eprev);

				double* Lafter = calloc(sizeof(double), 3);
				tools_get_total_L_vec(Lafter);
				fprintf(of, "Frac. L of ej. particle = %.3e\n", sqrt(Lejected[0]*Lejected[0] + Lejected[1]*Lejected[1] + Lejected[2]*Lejected[2])/sqrt(Lprev[0]*Lprev[0] + Lprev[1]*Lprev[1] + Lprev[2]*Lprev[2]));
				fprintf(of, "delta (Lx, Ly, Lz) = (%.15e, %.15e, %.15e)\n\n", Lafter[0] + Lejected[0] - Lprev[0], Lafter[1] + Lejected[1] - Lprev[1], Lafter[2] + Lejected[2] - Lprev[2]);

				double Ecom = tools_com_ke();
				double* Lcom = calloc(sizeof(double), 3);
				tools_com_L_vec(Lcom);

				tools_move_to_center_of_momentum();

				double Efinal =  tools_get_total_E();
				double* Lfinal = calloc(sizeof(double), 3);
				tools_get_total_L_vec(Lfinal);

				fprintf(of, "Frac. E in COM = %.3e\n", Ecom/Eafter);
				fprintf(of, "delta E = %.15e\n", Efinal + Ecom - Eafter);
				fprintf(of, "Frac. L in COM = %.3e\n", sqrt(Lcom[0]*Lcom[0] + Lcom[1]*Lcom[1] + Lcom[2]*Lcom[2])/sqrt(Lafter[0]*Lafter[0] + Lafter[1]*Lafter[1] + Lafter[2]*Lafter[2]));
				fprintf(of, "delta (Lx, Ly, Lz) = (%.15e, %.15e, %.15e)\n\n", Lfinal[0] + Lcom[0] - Lafter[0], Lfinal[1] + Lcom[1] - Lafter[1], Lfinal[2] + Lcom[2] - Lafter[2]);

				Eadj += Eejected + Ecom; // update the energy adjustment from ejected particles
				for(int j=0;j<3;j++){
					Ladj[j] += Lejected[j] + Lcom[j];
				}
			}
			else{
				particles[i] = particles[N-1];
				i--;
				N--;
			}

			if(N==2){
				exit_simulation = 1;
				char* end = "endcondition.txt"; // end of simulation time
				FILE* of = fopen(end, "w");
				if (of==NULL){
				    printf("\n\nError while opening file '%s'.\n", end);
				    return;
				}
				fprintf(of, "2\n"); // simulation ended because went down to 2 planets
				fclose(of);
			}

			if(N==3){
				double poveracrit = 1. + pow(3, 4./3.)*particles[1].m*particles[2].m/pow(particles[0].m, 2./3.)/pow(particles[1].m+particles[2].m, 4./3.);

				if(tools_p_over_a(G, particles[0], particles[1], particles[2]) > poveracrit){
					exit_simulation = 1;
					char* end = "endcondition.txt"; // end of simulation time
					FILE* of = fopen(end, "w");
					if (of==NULL){
					    printf("\n\nError while opening file '%s'.\n", end);
					    return;
					}
					fprintf(of, "3\n"); // simulation ended because went down to 3 planets with H and E above critical value for Hill stability (Gladman 1993)
					fclose(of);
				}
			}
#endif
		}
	}
}

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = boxsize_x*(double)i;
	gb.shifty = boxsize_y*(double)j;
	gb.shiftz = boxsize_z*(double)k;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}

/**
 * Checks if a given particle is within the computational domain.
 * @param p Particle to be checked.
 * @return Return value is 1 if particle is inside the box and 0 otherwise.
 */
int boundaries_particle_is_in_box(struct particle p){
	if(p.x>boxsize_x/2.){
		return 0;
	}
	if(p.x<-boxsize_x/2.){
		return 0;
	}
	if(p.y>boxsize_y/2.){
		return 0;
	}
	if(p.y<-boxsize_y/2.){
		return 0;
	}
	if(p.z>boxsize_z/2.){
		return 0;
	}
	if(p.z<-boxsize_z/2.){
		return 0;
	}
	return 1;
}

void boundaries_track_IDs(){
	track_IDs = 1;
}

void boundaries_track_conservation(){
	track_conservation = 1;
	Eadj = 0.; // initialize adjustment values of E and L for when particles leave box
	Ladj = calloc(sizeof(double),3);
}
