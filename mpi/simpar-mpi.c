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

#define euclidean(x1,x2,y1,y2)       ((sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2))))
#define wrap_around(index, min, max) (index < min ? max : (index > max ) ? min : index)

/**********
 * HELPERS
 **********/

void print_particles(particle_t *par) {
    /*
     *  Helper function to visualize each particle's position and mass
     *  print = (x, y) - m
     * */
    printf("[Particles]\n");
    if(par == NULL) {
        printf("\t-\n");
    }
    else {
        while(par->next != NULL) {
            printf("\tp=(%.2f, %.2f) | f=(%.2f, %.2f)\n", par->x, par->y, par->fx, par->fy);
            par = par->next;
        }
        printf("\tp=(%.2f, %.2f) | f=(%.2f, %.2f)\n", par->x, par->y, par->fx, par->fy);
    }
}

void print_cells(long ncside, cell_t **cells) {
    /*
     *  Helper function to visualize each cells's position and mass
     *  print = (x, y) - m
     * */
    int i, j;
    for (i = 0; i < ncside; i++) {
        for (j = 0; j < ncside; j++) {
            printf("\n[Cell (x=%d, y=%d)]  ", i, j);
            printf("(%.2f, %.2f) - %.2f, %ld\n", cells[i][j].x, cells[i][j].y, cells[i][j].m, cells[i][j].npar);
            print_particles(cells[i][j].par);
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
            cell->par = NULL;
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

void allocate_particles_to_cells(long ncside, long long n_part, cell_t **cells, particle_t *par) {
    // go through all particles and allocate each to the cell
    // it belongs to
    long long k;
    int cx, cy;
    double interval = 1.0 / ncside;

    for (k = 0; k < n_part; k++) {
        cx = calc_cell_number(par[k].x, interval, ncside);
        cy = calc_cell_number(par[k].y, interval, ncside);
        add_particle_to_cell_linked_list(&cells[cx][cy], &par[k]);
    }
}

void add_particle_to_cell_linked_list(cell_t *cell, particle_t *par) {
    /*
     * Add new particle to the beginning of
     * the linked list of particles of the given cell
     */

    particle_t *head = cell->par;
    cell->par = par;
    par->next = head;
}

/*******
 * MAIN
 *******/

int main(int argc, char *argv[]) {

    int id = 0, p;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // init values
    long nseed = strtol(argv[1], NULL, 10);
    long ncside = strtol(argv[2], NULL, 10);
    long long n_part = strtol(argv[3], NULL, 10);
    long n_tsteps = strtol(argv[4], NULL, 10);
    particle_t *par;
    cell_t **cells;

    /*
     * First we have to initiate all particles in non-parallel
     * mode so that we have the same base to start with. Only process 0 will do this.
     * We then have to distribute these particles among the other processes.
     */
    if(id == 0) {
        printf("==== Particles init allocation by process 0 ========\n");
        // init particles locally
        par = malloc(n_part * sizeof(particle_t));
        init_particles(nseed, ncside, n_part, par);

        // init cell matrix locally
        cells = malloc(ncside * sizeof(cell_t *));
        create_cells_matrix(ncside, cells);
        init_cells_matrix(ncside, cells);
        allocate_particles_to_cells(ncside, n_part, cells, par);

        print_cells(ncside, cells);
    }

    // TODO: Insert WAIT function here because process 0 has to finish first

    /* From here on the computation is done identically at each process.
     * First we have to distribute the particles (that only process 0 contains,
     * across the network.
     * Only then can each process start the computation.
     */

    printf("==== Particles received by process %d ========\n", id);

    /*
    cell_t *column; // receiving variable
    int cols_per_proc = ncside / p;

    MPI_Scatter(cells, cols_per_proc, cell_t , column,
                cols_per_proc, cell_t, 0, MPI_COMM_WORLD);


    print_cells();

    */

    MPI_Finalize();

    return EXIT_SUCCESS;
}
