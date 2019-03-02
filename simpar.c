#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init_particles.c"
#include "types.h"


void print_particles(long long n_part, particle_t *par) {
	/*
	Helper function to visualize each particle's position and mass
	print = (x, y) - m
	*/
    printf("[Particles]\n");
    int i;
    for(i = 0; i < n_part; i++) {
        printf("\t(%f, %f) - %f\n", par[i].x, par[i].y, par[i].m);
    }
}

void print_cells(long ncside, cell_t **cells) {
	/*
	Helper function to visualize each cells's position and mass
	print = (x, y) - m
	*/
    printf("[Cells]\n");
    int i, j;
    for(i = 0; i < ncside; i++) {
		for(j = 0; j < ncside; j++) {
        	printf("\t(%f, %f) - %f\n", cells[i][j].x, cells[i][j].y, cells[i][j].m);
		}
    }
}

void calc_all_cells_cm(long ncside, cell_t **cells, long long n_part, particle_t *par) {
	/*
	Determine center of mass of all cells.
	For each particle, determine to which cell it belongs,
	add to total mass, total number of particles and total position

	CM = (for all i SUM((xi, yi) * mi)) / SUM(for all i, mi)
	*/
	int i, j, cellx, celly;
	double interval = 1.0 / ncside;
    for(i = 0; i < n_part; i++) {
		cellx = ((int) floor(par[i].x / interval)) % ncside ;
		celly = ((int) floor(par[i].y / interval)) % ncside ;
		//printf("[Part] (%f, %f) => [Cell] (%d, %d)\n", par[i].x, par[i].y, cellx, celly);

		cells[cellx][celly].x += par[i].x * par[i].m;
		cells[cellx][celly].y += par[i].y * par[i].m;
		cells[cellx][celly].m += par[i].m;
		cells[cellx][celly].npar++;
    }

	// after all total masses and positions have been determined
	for(i = 0; i < ncside; i++) {
		for(j = 0; j < ncside; j++) {
			if(cells[i][j].npar != 0) {
				cells[i][j].x = cells[i][j].x / cells[i][j].m ;
				cells[i][j].y = cells[i][j].y / cells[i][j].m;
			}
		}
	}
}

void init_cells_matrix(long ncside, cell_t **cells) {
	// access cell (x, y) => cells[x][y]
	int i, j;
	for(i=0; i < ncside; i++) {
		cells[i] = malloc(ncside * sizeof(cell_t));
	}

	// default value to 0
	for(i=0; i < ncside; i++) {
		for(j=0; j < ncside; j++) {
			cells[i][j].x = 0;
			cells[i][j].y = 0;
			cells[i][j].x = 0;
			cells[i][j].npar = 0;
		}
	}
}

int main(int argc, char *argv[]) {

    // receive exactly 4 arguments (first is file name)
    if(argc != 5) return EXIT_FAILURE;

    // init values
    long nseed = strtol(argv[1], NULL, 10);
    long ncside = strtol(argv[2], NULL, 10);
    long long n_part = strtol(argv[3], NULL, 10);
    particle_t *par = malloc(n_part * sizeof(particle_t));
    init_particles(nseed, ncside, n_part, par);

	// init cell matrix
	cell_t **cells = malloc(ncside * sizeof(cell_t*));
	init_cells_matrix(ncside, cells);

	// determine center of mass of all cells
	calc_all_cells_cm(ncside, cells, n_part, par);

	// TODO: compute the gravitational force applied to each particle

	// TODO: calculate the new velocity and then the new position of each particle

	print_particles(n_part, par);
	print_cells(ncside, cells);

	// free memory
	free(par);
	int i;
	for(i=0; i < ncside; i++) {
		free(cells[i]);
	}
	free(cells);

    return EXIT_SUCCESS;
}