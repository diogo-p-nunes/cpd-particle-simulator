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

#include "simpar-mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stddef.h>


#define euclidean(x1,x2,y1,y2)       ((sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2))))
#define wrap_around(index, min, max) (index < min ? max : (index > max ) ? min : index)



void get_processor_c(int id, int dim, int size, long ncside, int *c) {
    int coord = id % dim, cmin = coord * size, cmax;

    if(coord == dim-1) cmax = ncside-1;
    else cmax = (coord+1) * size -1;
    c[0] = cmin; c[1] = cmax;
}


void init_processors_particles(particle_t ***arr, int *processors_particles_sizes, int p) {
    int i=0;
    for(i; i<p; i++) {
        arr[i] = NULL;
        processors_particles_sizes[i] = 0;
    }
}

void generate_sending_data(long long n_part, particle_t* par, long ncside, particle_t ***processors_particles,
        int **id_map, int* processors_particles_sizes) {

    // go through each particle and determine to which processor it belongs to
    // add that particle to the dynamic array of particles of that processor

    int cellx, celly, id;
    double interval = 1.0 / ncside;

    for(int i=0; i<n_part; i++) {
        cellx = calc_cell_number(par[i].x, interval, ncside);
        celly = calc_cell_number(par[i].y, interval, ncside);
        id = id_map[cellx][celly];
        add_particle_to_processor_array(&processors_particles[id], &par[i], id, processors_particles_sizes);
    }
}

void add_particle_to_processor_array(particle_t ***array, particle_t *par, int id, int* processors_particles_sizes) {

    if(processors_particles_sizes[id] == 0) {
        (*array) = (particle_t**) malloc(sizeof(particle_t*));
        (*array)[processors_particles_sizes[id]] = par;

        processors_particles_sizes[id] += 1;
    }
    else {
        (*array) = (particle_t**) realloc((*array), (processors_particles_sizes[id]+1) * sizeof(particle_t*));
        (*array)[processors_particles_sizes[id]] = par;

        processors_particles_sizes[id] += 1;
    }

}


void distribute_processors_particles(particle_t ***processors_particles, int dims[], int sizes[], long ncside) {

}


void init_id_map(int **id_map, long ncside) {
    for(int i=0; i<ncside; i++) {
        id_map[i] = (int*) malloc(ncside* sizeof(int));
    }
}

void populate_id_map(int **id_map, int dims[], int sizes[], long ncside, int p) {
    int *cx, *cy;
    for(int k=0; k<p; k++){

        cx = (int*) malloc(2*sizeof(int));
        cy = (int*) malloc(2*sizeof(int));

        get_processor_c(k, dims[0], sizes[0], ncside, cx);
        get_processor_c(k, dims[1], sizes[1], ncside, cy);
        //printf("id:%d, cx=[%d,%d], cy=[%d,%d]\n", k, cx[0], cx[1], cy[0], cy[1]);

        for(int i=cx[0]; i<=cx[1]; i++) {
            for(int j=cy[0]; j<=cy[1]; j++) {
                id_map[i][j] = k;
            }
        }
    }
}


void print_id_map(int **id_map,long ncside) {
    for(int i=0; i<ncside; i++) {
        for(int j=0; j<ncside; j++) {
            printf("[%d,%d] -> %d\n", i,j, id_map[i][j]);
            fflush(stdout);
        }
    }
}


void print_processors_particles(particle_t ***processors_particles, int *processors_particles_sizes, int p) {
    for(int i=0; i<p; i++) {
        printf("id:%d - ", i);
        fflush(stdout);
        particle_t **array = processors_particles[i];

        for(int j=0; j<processors_particles_sizes[i]; j++) {
            print_particle(array[j]);
        }
        printf("\n");
    }
}


void print_particle(particle_t *par) {
    printf("(%.2f,%.2f) ", par->x, par->y);
    fflush(stdout);
}


/*****************
 * CELL FUNCTIONS
 *****************/

int calc_cell_number(double pos, double interval, long ncside) {
    return labs(((int) floor(pos / interval)) % ncside);
}

void init_cells_matrix(long ncside, cell_t **cells) {
    // default value to 0
    int i, j;
    cell_t *cell;

    for (i = 0; i < ncside; i++) {
        for (j = 0; j < ncside; j++) {
            cell = &cells[i][j];
            cell->x = 0.0;
            cell->y = 0.0;
            cell->m = 0.0;
            cell->npar = 0;
        }
    }
}

void create_cells_matrix(long ncside, cell_t **cells) {
    // access cell (x, y) => cells[x][y]
    int i;
    for (i = 0; i < ncside; i++) {
        cells[i] = malloc(ncside * sizeof(cell_t));
    }
}


/*
void build_mpi_cell_type(cell_t *cell, MPI_Datatype *mpi_cell_type, MPI_Datatype *mpi_particle_type) {

    // Create a MPI type for struct cell_t
    // Build a derived datatype consisting of three doubles, a long and one mpi_particle_type
    const int nitems=5;

    // Specify the number of elements of each type
    int blocklengths[5] = {1,1,1,1,1};

    // First specify the types
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_LONG, (*mpi_particle_type)};

    MPI_Aint displacements[5];
    MPI_Aint addresses[6];

    // Calculate the displacements of the members relative to cell
    MPI_Get_address(cell, &addresses[0]);
    MPI_Get_address(&(cell->x), &addresses[1]);
    MPI_Get_address(&(cell->y), &addresses[2]);
    MPI_Get_address(&(cell->m), &addresses[3]);
    MPI_Get_address(&(cell->npar), &addresses[4]);
    MPI_Get_address(&(cell->par), &addresses[5]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];
    displacements[2] = addresses[3] - addresses[0];
    displacements[3] = addresses[4] - addresses[0];
    displacements[4] = addresses[5] - addresses[0];

    // Create the derived type
    MPI_Type_create_struct(nitems, blocklengths, displacements, types, mpi_cell_type);

    // Commit it so that it can be used
    MPI_Type_commit(mpi_cell_type);
}


void build_mpi_particle_type(particle_t *par, MPI_Datatype *mpi_particle_type) {

    // Create a MPI type for struct particle_t
    // Build a derived datatype consisting of three doubles, a long and one mpi_particle_type
    const int nitems=8;

    // Specify the number of elements of each type
    int blocklengths[8] = {1,1,1,1,1,1,1,1};

    // TODO: E aqui que esta a ciganice. Estou a definir um tipo novo baseado nele mesmo.
    // TODO: Nao sei como dar a volta a isto, se calhar esta mal a forma de pensar.
    // First specify the types
    MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, (*mpi_particle_type)};

    MPI_Aint displacements[8];
    MPI_Aint addresses[9];

    // Calculate the displacements of the members relative to cell
    MPI_Get_address(par, &addresses[0]);
    MPI_Get_address(&(par->x), &addresses[1]);
    MPI_Get_address(&(par->y), &addresses[2]);
    MPI_Get_address(&(par->vx), &addresses[3]);
    MPI_Get_address(&(par->vy), &addresses[4]);
    MPI_Get_address(&(par->m), &addresses[5]);
    MPI_Get_address(&(par->fx), &addresses[6]);
    MPI_Get_address(&(par->fy), &addresses[7]);
    MPI_Get_address(&(par->next), &addresses[8]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];
    displacements[2] = addresses[3] - addresses[0];
    displacements[3] = addresses[4] - addresses[0];
    displacements[4] = addresses[5] - addresses[0];
    displacements[5] = addresses[6] - addresses[0];
    displacements[6] = addresses[7] - addresses[0];
    displacements[7] = addresses[8] - addresses[0];

    // Create the derived type
    MPI_Type_create_struct(nitems, blocklengths, displacements, types, mpi_particle_type);

    // Commit it so that it can be used
    MPI_Type_commit(mpi_particle_type);
}
*/

/*******
 * MAIN
 *******/

int main(int argc, char *argv[]) {

    int id, p;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // init values
    long nseed = strtol(argv[1], NULL, 10);
    long ncside = strtol(argv[2], NULL, 10);
    long long n_part = strtol(argv[3], NULL, 10);
    long n_tsteps = strtol(argv[4], NULL, 10);

    // information regarding the whole system
    particle_t *par;
    cell_t **cells;

    /*
     * First we have to initiate all particles in non-parallel
     * mode so that we have the same base to start with. Only process 0 will do this.
     * It then has to distribute these particles among the other processes.
     */

    // Get best dimensions
    int dims[2] = {0,0};

    // x = dims[0], y = dims[1]
    MPI_Dims_create(p, 2, dims);

    // x = sizes[0], y = sizes[1]
    int sizes[2] = {floor(ncside/dims[0]), floor(ncside/dims[1])};

    // Define the processor id for the coordinates
    int **id_map = (int**) malloc(sizeof(int*)*ncside);
    init_id_map(id_map, ncside);
    populate_id_map(id_map, dims, sizes, ncside, p);

    if(id == 0) {

        // init particles locally
        par = malloc(n_part * sizeof(particle_t));
        init_particles(nseed, ncside, n_part, par);

        // init array with particles that belong to each processor id
        particle_t ***processors_particles = (particle_t***) malloc(p * sizeof(particle_t**));

        // keep in memory the number of particles in each processor id array above
        int *processors_particles_sizes = (int*) malloc(p * sizeof(int));

        // in initialize array and sizes side array
        init_processors_particles(processors_particles, processors_particles_sizes, p);

        // populate the processors particles array with the corresponding data of each processor
        generate_sending_data(n_part, par, ncside, processors_particles, id_map, processors_particles_sizes);

        print_processors_particles(processors_particles, processors_particles_sizes, p);

        // distribute particles across all processors
        //distribute_processors_particles(processors_particles, dims, sizes, ncside);
    }

    // All processes must wait until process 0 has finished initializing the system
    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * From here on the computation is done identically at each process.
     */

    // TODO: Acho que nao da para nativamente enviar estruturas no MPI, portanto temos de definir
    // TODO: nos mesmos as estruturas. Eu tentei desta forma, mas estou a fazer umas ciganices que dao
    // TODO: SEGMENTATION FAULT
    // build MPI data types for transfer

    /*
    MPI_Datatype mpi_particle_type;
    build_mpi_particle_type(&par[0], &mpi_particle_type);

    MPI_Datatype mpi_cell_type;
    build_mpi_cell_type(&cells[0][0], &mpi_cell_type, &mpi_particle_type);

     */

    MPI_Finalize();
    return EXIT_SUCCESS;
}
