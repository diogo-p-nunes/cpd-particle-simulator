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
    long npar;

} cell_t;

void calc_and_print_overall_cm(long long n_part, particle_t *par);

void free_memory(int ncside, cell_t **cells, particle_t *par);

void update_force(cell_t *cell, particle_t *par);

void calc_all_particle_force(long ncside, cell_t **cells, long long n_part, particle_t *par);

void update_vel(double acc_x, double acc_y, particle_t *particle);

void update_pos(double acc_x, double acc_y, particle_t *particle);

void calc_all_particle_new_values(long ncside, long long n_part, particle_t *par);

int calc_cell_number(double pos, double interval, long ncside);

void calc_all_cells_cm(long ncside, cell_t **cells, long long n_part, particle_t *par);

void init_cells_matrix(long ncside, cell_t **cells);

void create_cells_matrix(long ncside, cell_t **cells);

void init_particle_force(particle_t *par, long long n_part);


// new funcs

void get_processor_c(int id, int dim, int size, long ncside, int *c);

void init_processors_particles(particle_t **arr, int *processors_particles_size, int p);

void generate_sending_data(long long n_part, particle_t* par, long ncside, particle_t **processors_particles,
        int **id_map, int* processors_particles_sizes);

void add_particle_to_processor_array(particle_t **array, particle_t *par, int id, int* processors_particles_sizes);

void distribute_processors_particles(particle_t **processors_particles, int* processors_particles_sizes, int p);

void init_id_map(int **id_map, long ncside);

void populate_id_map(int **id_map, int dims[], int sizes[], long ncside, int p);

void print_id_map(int **id_map,long ncside);

void print_processors_particles(particle_t **processors_particles, int *processors_particles_sizes, int p);

void print_particle(particle_t *par);


#endif
