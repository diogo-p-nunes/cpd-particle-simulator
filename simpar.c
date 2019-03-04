/***********************************************************************************************
*                                                                                              *
*                                          Grupo: 18                                           *
*                                                                                              *
*                                   Beatriz Marques , 80809                                    *
*                                   Carlos  Carvalho, 81395                                    *
*                                   Diogo   Nunes   , 85184                                    *
*                	              							       *
*		    Copyright (c) 2019 Beatriz, Carlos e Diogo. All rights reserved.           *
*                                                                                              *
***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "init_particles.h"
#include "simpar.h"

void print_particles(long long n_part, particle_t *par) {
    /*
     *  Helper function to visualize each particle's position and mass
     *  print = (x, y) - m
     * */
    printf("[Particles]\n");
    int i;
    for(i = 0; i < n_part; i++) {
        printf("\tp=(%f, %f) | f=(%f, %f)\n", par[i].x, par[i].y, par[i].fx, par[i].fy);
    }
}

void print_cells(long ncside, cell_t **cells) {
    /*
     *  Helper function to visualize each cells's position and mass
     *  print = (x, y) - m
     * */
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
     *  Determine center of mass of all cells.
     *  For each particle, determine to which cell it belongs,
     *  add to total mass, total number of particles and total position
     *
     *  CM = (for all i SUM((xi, yi) * mi)) / SUM(for all i, mi)
     * */
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
				cells[i][j].x = cells[i][j].x / cells[i][j].m;
				cells[i][j].y = cells[i][j].y / cells[i][j].m;
			}
        }
    }
}

double euclidean_distance(double x1, double x2, double y1, double y2){
    /*
     *  Compute Euclidean distance between two points (x1,y1) and (x2,y2)
     *
     * */
    double dx = (x1-x2);
    double dy = (y1-y2);
    return sqrt(dx*dx + dy*dy);
}

void update_force(int i, int cellx, int celly, cell_t **cells, particle_t *par){
    /*
     * Update the gravitacional force of particle i
     * applied by the center of mass from cell (cellx, celly).
     *
     * */
    double dist;
	double fdirx, fdiry, magnitude, norm;
	
    if(cells[cellx][celly].npar != 0){
        dist = euclidean_distance(cells[cellx][celly].x, par[i].x, cells[cellx][celly].y, par[i].y);
		if(dist >= EPSLON) {
			magnitude = (G * par[i].m * cells[cellx][celly].m) / (dist * dist);
		}
		else {
			magnitude = 0;
		}

		// determine force direction and norm
		fdirx = cells[cellx][celly].x - par[i].x; 
		fdiry = cells[cellx][celly].y - par[i].y;
		norm = sqrt(fdirx*fdirx + fdiry*fdiry);

		// normalize direction vector
		fdirx = fdirx / norm; 
		fdiry = fdiry / norm;

		// update particle's force vector
		par[i].fx += (fdirx * magnitude); 
		par[i].fy += (fdiry * magnitude);
    }
}

void calc_all_particle_force(long ncside, cell_t **cells, long long n_part, particle_t *par){
    /*
     *  Determine gravitational force for each particle.
     *  A particle is influenced by gravity from all the center of masses
     *  of the adjacenct cells and the center of mass of the cell that belongs to.
     *
     *  GF = G * ma * mb / dab^2
     * */
    int i, cellx, celly, cx, cy;
    double interval = 1.0 / ncside;
	
    for(i = 0; i < n_part; i++) {
        cellx = ((int) floor(par[i].x / interval)) % ncside;
		celly = ((int) floor(par[i].y / interval)) % ncside;

		// TODO: This whole part can be refactored into a few loops

		// actual cell of the particle
		update_force(i, cellx, celly, cells, par);

		// right cell
		cy = celly+1;
		if(cy == ncside)
		cy = 0;
		update_force(i, cellx, cy, cells, par);

		// left cell
		cy = celly-1;
		if(cy == -1)
			cy = ncside-1;
		update_force(i, cellx, cy, cells, par);

		// up cell
		cx = cellx-1;
		if(cx == -1)
			cx = ncside-1;
		update_force(i, cx, celly, cells, par);

		// diagonal up-left
		cy = celly-1;
		if(cy == -1)
			cy = ncside-1;
		update_force(i, cx, cy, cells, par);

		// diagonal up-right
		cy = celly+1;
		if(cy == ncside)
			cy = 0;
		update_force(i, cx, cy, cells, par);

		// down cell
		cx = cellx+1;
		if(cx == ncside)
			cx = 0;
		update_force(i, cx, celly, cells, par);

		// diagonal down-right
		cy = celly+1;
		if(cy == ncside)
			cy = 0;
		update_force(i, cx, cy, cells, par);

		// diagonal down_left
		cy = celly-1;
		if(cy == -1)
			cy = ncside-1;
		update_force(i, cx, cy, cells, par);
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
	    cells[i][j].m = 0;
	    cells[i][j].npar = 0;
        }
    }
}

void free_memory(int ncside, cell_t **cells, particle_t *par){
    free(par);
    int i;
    for(i=0; i < ncside; i++)
        free(cells[i]);
    free(cells);
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
	
    // compute the gravitational force applied to each particle
    calc_all_particle_force(ncside, cells, n_part, par);

    // TODO: calculate the new velocity and then the new position of each particle

    print_particles(n_part, par);
    print_cells(ncside, cells);

    // free memory
    free_memory(ncside, cells, par);

    return EXIT_SUCCESS;
}
