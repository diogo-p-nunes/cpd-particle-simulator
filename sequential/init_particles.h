/***********************************************************************************************
 *                                                                                              *
 *                                          Grupo: 18 *
 *                                                                                              *
 *                                   Beatriz Marques , 80809 * Carlos  Carvalho,
 *81395                                    * Diogo   Nunes   , 85184 *
 *                	              							       *
 *		    Copyright (c) 2019 Beatriz, Carlos e Diogo. All rights
 *reserved.           *
 *                                                                                              *
 ***********************************************************************************************/

#ifndef INIT_PARTICLES_H
#define INIT_PARTICLES_H

#include <stdlib.h>

#define RND0_1 ((double)random() / ((long long)1 << 31))
#define G 6.67408e-11
#define EPSLON 0.0005

typedef struct particle_t {
  double x, y, vx, vy, m;
  double fx, fy;
} particle_t;

void init_particles(long seed, long ncside, long long n_part, particle_t *par);

#endif
