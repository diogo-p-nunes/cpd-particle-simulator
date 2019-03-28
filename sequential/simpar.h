/***********************************************************************************************
 *                                                                                              *
 *                                          Grupo: 18                                           *
 *                                                                                              *
 *                                   Beatriz Marques , 80809                                    *
 *                                   Carlos  Carvalho, 81395                                    *
 *                                   Diogo   Nunes   , 85184                                    *
 *                	              							                                    *
 *		    Copyright (c) 2019 Beatriz, Carlos e Diogo. All rights reserved.                    *
 *                                                                                              *
 ***********************************************************************************************/

#ifndef SIMPAR_H
#define SIMPAR_H

#include "init_particles.h"

#define PI (3.141592653589793)

/*
All type definitions
*/

typedef struct {
    // if cell has npar == 0, then it does not have any par inside
    double x, y, m;
    long npar;
} cell_t;

void print_particles(long long n_part, particle_t *par);

void print_cells(long ncside, cell_t **cells);

void calc_and_print_overall_cm(long long n_part, particle_t *par);

double euclidean_distance(double x1, double x2, double y1, double y2);

void free_memory(int ncside, cell_t **cells, particle_t *par);

void update_force(cell_t *cell, particle_t *par);

void calc_all_particle_force(long ncside, cell_t **cells, long long n_part, particle_t *par);

void update_vel(double acc_x, double acc_y, particle_t *particle);

void update_pos(double acc_x, double acc_y, particle_t *particle);

void calc_all_particle_new_values(long ncside, long long n_part, particle_t *par);

int calc_cell_number(double pos, double interval, long ncside);

void calc_all_cells_cm(long ncside, cell_t **cells, long long n_part, particle_t *par);

int wrap_around(int index, long min, long max);

void init_cells_matrix(long ncside, cell_t **cells);

void create_cells_matrix(long ncside, cell_t **cells);

void init_particle_force(particle_t *par, long long n_part);

#endif
