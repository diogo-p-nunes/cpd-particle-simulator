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

#ifndef SIMPAR_MPI_H
#define SIMPAR_MPI_H

#include "init_particles.h"

#define PI (3.141592653589793)

/*
All type definitions
*/

typedef struct {
    // if cell has npar == 0, then it does not have any par inside
    double x, y, m;
    int cx, cy;
    long npar;

} cell_t;

void calc_and_print_overall_cm(long long n_part, particle_t *par, double **result);

void free_memory(int ncside, cell_t **cells, particle_t *par);

void update_force(cell_t *cell, particle_t *par);

void calc_all_particle_force(long ncside, cell_t **cells, long long n_part, int *cx, int *cy,
            particle_t *par, int cols, int rows, cell_t **id_received_cells_map,int *dims, int*counter);

void update_vel(double acc_x, double acc_y, particle_t *particle);

void update_pos(double acc_x, double acc_y, particle_t *particle);

void calc_all_particle_new_values(long ncside, long long n_part, particle_t *par);

int calc_cell_number(double pos, double interval, long ncside);

void calc_all_cells_cm(cell_t **cells, long long n_part, particle_t *par, int cols, int rows,
        long ncside, double interval);

void init_cells_matrix(int cols, int rows, cell_t **cells, int xmin, int ymin);

void create_cells_matrix(int cols, int rows, cell_t **cells);

void init_particle_force(particle_t *par, long long num_par, int id, particle_t **containsParticleZero);


// new funcs

void get_processor_c(int id, int dim, int size, long ncside, int *c);

void init_processors_particles(particle_t **arr, int *processors_particles_size, int p);

void generate_sending_data(long long n_part, particle_t* par, long ncside, particle_t **processors_particles,
        int* processors_particles_sizes, int*dims);

void add_particle_to_processor_array(particle_t **array, particle_t *par, int id, int* processors_particles_sizes);

void distribute_processors_particles(particle_t **processors_particles, int* processors_particles_sizes, int p);

void print_processors_particles(particle_t **processors_particles, int *processors_particles_sizes, int p);

void print_particle(particle_t *par);

void calculate_c(int *cx, int *cy, int id, int ncside, int*dims);


#endif
