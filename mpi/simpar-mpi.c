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


/*********
 * HELPERS
 *********/

void print_id_map(int **id_map,long ncside) {
    for(int i=0; i<ncside; i++) {
        for(int j=0; j<ncside; j++) {
            printf("[%d,%d] -> %d\n", i,j, id_map[i][j]);
            fflush(stdout);
        }
    }
}


void print_processors_particles(particle_t **processors_particles, int *processors_particles_sizes, int p) {
    for(int i=0; i<p; i++) {
        printf("id:%d - ", i);
        fflush(stdout);
        particle_t *array = processors_particles[i];

        for(int j=0; j<processors_particles_sizes[i]; j++) {
            print_particle(&array[j]);
        }
        printf("\n");
    }
}


void print_particle(particle_t *par) {
    printf("(%.2f,%.2f) ", par->x, par->y);
    fflush(stdout);
}



/**************************
 * DATA GENERATION/DISTRIB
 **************************/

void get_processor_c(int id, int dim, int size, long ncside, int *c) {
    int coord = id % dim, cmin = coord * size, cmax;

    if(coord == dim-1) cmax = ncside-1;
    else cmax = (coord+1) * size -1;
    c[0] = cmin; c[1] = cmax;
}

void init_processors_particles(particle_t **arr, int *processors_particles_sizes, int p) {
    int i=0;
    for(i; i<p; i++) {
        arr[i] = NULL;
        processors_particles_sizes[i] = 0;
    }
}

void generate_sending_data(long long n_part, particle_t* par, long ncside, particle_t **processors_particles,
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

void add_particle_to_processor_array(particle_t **array, particle_t *par, int id, int* processors_particles_sizes) {

    if(processors_particles_sizes[id] == 0) {
        (*array) = (particle_t*) malloc(sizeof(particle_t));
        (*array)[processors_particles_sizes[id]] = (*par);

        processors_particles_sizes[id] += 1;
    }
    else {
        (*array) = (particle_t*) realloc((*array), (processors_particles_sizes[id]+1) * sizeof(particle_t));
        (*array)[processors_particles_sizes[id]] = (*par);

        processors_particles_sizes[id] += 1;
    }

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

        for(int i=cx[0]; i<=cx[1]; i++) {
            for(int j=cy[0]; j<=cy[1]; j++) {
                id_map[i][j] = k;
            }
        }
    }
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


void build_mpi_particle_type(MPI_Datatype *mpi_particle_type) {

    // Data type we want to duplicate into MPI
    particle_t *par;

    // Create a MPI type for struct particle_t
    // Build a derived datatype consisting of three doubles, a long and one mpi_particle_type
    const int nitems=7;

    // Specify the number of elements of each type
    int blocklengths[7] = {1,1,1,1,1,1,1};

    // First specify the types
    MPI_Datatype types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,};

    MPI_Aint displacements[7];
    MPI_Aint addresses[8];

    // Calculate the displacements of the members relative to cell
    MPI_Get_address(par, &addresses[0]);
    MPI_Get_address(&(par->x), &addresses[1]);
    MPI_Get_address(&(par->y), &addresses[2]);
    MPI_Get_address(&(par->vx), &addresses[3]);
    MPI_Get_address(&(par->vy), &addresses[4]);
    MPI_Get_address(&(par->m), &addresses[5]);
    MPI_Get_address(&(par->fx), &addresses[6]);
    MPI_Get_address(&(par->fy), &addresses[7]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];
    displacements[2] = addresses[3] - addresses[0];
    displacements[3] = addresses[4] - addresses[0];
    displacements[4] = addresses[5] - addresses[0];
    displacements[5] = addresses[6] - addresses[0];
    displacements[6] = addresses[7] - addresses[0];

    // Create the derived type
    MPI_Type_create_struct(nitems, blocklengths, displacements, types, mpi_particle_type);

    // Commit it so that it can be used
    MPI_Type_commit(mpi_particle_type);
}


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

    // information regarding the problem
    particle_t *par;

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

    // Define the mapping of processor id for processor coordinates
    int **id_map = (int**) malloc(sizeof(int*)*ncside);
    init_id_map(id_map, ncside);
    populate_id_map(id_map, dims, sizes, ncside, p);

    // Define the number of particles for each processor (local)
    int num_par = 0;

    // Define data sending TAGS for validation: IMPORTANT
    int NUM_PAR_TAG = 0;
    int PAR_TAG = 1;

    // build MPI data types for transfer
    MPI_Datatype mpi_particle_type;
    build_mpi_particle_type(&mpi_particle_type);

    if(id == 0) {

        // init particles locally
        par = malloc(n_part * sizeof(particle_t));
        init_particles(nseed, ncside, n_part, par);

        // init array with particles that belong to each processor id
        particle_t **processors_particles = (particle_t**) malloc(p * sizeof(particle_t*));

        // keep in memory the number of particles in each processor id particle array
        int *processors_particles_sizes = (int*) malloc(p * sizeof(int));

        // in initialize array and sizes side array
        init_processors_particles(processors_particles, processors_particles_sizes, p);

        // populate the processors particles array with the corresponding data of each processor
        generate_sending_data(n_part, par, ncside, processors_particles, id_map, processors_particles_sizes);

        //print_processors_particles(processors_particles, processors_particles_sizes, p);

        // Distribute particles across all processors
        for(int k=1; k<p; k++) {
            // First each process needs to know how many particles it is going to receive in the array
            MPI_Send(&processors_particles_sizes[k], 1, MPI_INT, k, NUM_PAR_TAG, MPI_COMM_WORLD);

            // Now we can send the particles to the processors
            MPI_Send(processors_particles[k], processors_particles_sizes[k], mpi_particle_type, k, PAR_TAG, MPI_COMM_WORLD);
        }

        // Now processor id=0 can set its own particles (no longer knows about all particles)
        num_par = processors_particles_sizes[0];
        par = processors_particles[0];

    }
    else {

        // Receive respective particles from processor id=0 who initialized the system
        // First each process needs to know how many particles it is going to receive in the array
        MPI_Recv(&num_par, 1, MPI_INT, 0, NUM_PAR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // allocate memory for incoming number of particles
        par = (particle_t*) malloc(num_par * sizeof(particle_t));

        // all other processors receive the distributed data, respectively
        MPI_Recv(par, num_par, mpi_particle_type, 0, PAR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // print processor received particles
    printf("id:%d - num_par=%d\n", id, num_par);
    for(int i=0; i<num_par; i++) {
        printf("id:%d (%.2f, %.2f)\n", id, par[i].x, par[i].y);
    }



    /*
     * From here on the computation is done identically at each process.
     */
    cell_t **cells;






    MPI_Finalize();
    return EXIT_SUCCESS;
}
