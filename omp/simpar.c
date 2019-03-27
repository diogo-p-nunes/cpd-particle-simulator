/***********************************************************************************************
 *                                                                                              *
 *                                          Grupo: 18 *
 *                                                                                              *
 *                                   Beatriz Marques , 80809 * Carlos  Carvalho,
 *81395                                    * Diogo   Nunes   , 85184 *
 *                	              							                                   *
 *		    Copyright (c) 2019 Beatriz, Carlos e Diogo. All rights
 *reserved.                   *
 *                                                                                              *
 ***********************************************************************************************/

#include "simpar.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**********
 * HELPERS
 **********/

void print_particles(long long n_part, particle_t *par) {
  /*
   *  Helper function to visualize each particle's position and mass
   *  print = (x, y) - m
   * */
  printf("[Particles]\n");
  int i;
  for (i = 0; i < n_part; i++) {
    printf("\tp=(%f, %f) | f=(%f, %f)\n", par[i].x, par[i].y, par[i].fx,
           par[i].fy);
  }
}

void print_cells(long ncside, cell_t **cells) {
  /*
   *  Helper function to visualize each cells's position and mass
   *  print = (x, y) - m
   * */
  printf("[Cells]\n");
  int i, j;
  for (i = 0; i < ncside; i++) {
    for (j = 0; j < ncside; j++) {
      printf("\t(%f, %f) - %f\n", cells[i][j].x, cells[i][j].y, cells[i][j].m);
    }
  }
}

/****************
 * AUX FUNCTIONS
 ****************/

void calc_and_print_overall_cm(long long n_part, particle_t *par) {
  int i;
  double cmx = 0, cmy = 0, cmm = 0;

  #pragma omp parallel for reduction (+:cmx,cmy,cmm)
  for (i = 0; i < n_part; i++) {
    cmx += par[i].x * par[i].m;
    cmy += par[i].y * par[i].m;
    cmm += par[i].m;
  }

  cmx = cmx / cmm;
  cmy = cmy / cmm;

  printf("%.2f %.2f\n", cmx, cmy);
}

double euclidean_distance(double x1, double x2, double y1, double y2) {
  /*
   *  Compute Euclidean distance between two points (x1,y1) and (x2,y2)
   * */
  double dx = (x1 - x2);
  double dy = (y1 - y2);
  return sqrt(dx * dx + dy * dy);
}

void free_memory(int ncside, cell_t **cells, particle_t *par) {
  free(par);
  int i;

  #pragma omp parallel for
  for (i = 0; i < ncside; i++)
    free(cells[i]);
  free(cells);
}

int wrap_around(int index, long min, long max) {
  if (index < min)
    return (max);
  if (index > (max))
    return min;
  else
    return index;
}

/*********************
 * PARTICLE FUNCTIONS
 *********************/

void update_force(cell_t *cell, particle_t *particle) {
  /*
   * Update the gravitacional force of particle i
   * applied by the center of mass from cell (cellx, celly).
   *
   * */
  double dist, fx, fy, magnitude, norm;

  if (cell->npar != 0) {
    dist = euclidean_distance(cell->x, particle->x, cell->y, particle->y);

    if (dist >= EPSLON) {
      magnitude = (G * particle->m * cell->m) / (dist * dist);
    } else {
      magnitude = 0;
    }

    // determine force direction and norm
    fx = cell->x - particle->x;
    fy = cell->y - particle->y;
    norm = sqrt(fx * fx + fy * fy);

    if (norm != 0) {
      // normalize direction vector
      fx = fx / norm;
      fy = fy / norm;
    }

    // update particle's force vector
    particle->fx += (fx * magnitude);
    particle->fy += (fy * magnitude);
  }
}

void calc_all_particle_force(long ncside, cell_t **cells, long long n_part,
                             particle_t *par) {
  /*
   *  Determine gravitational force for each particle.
   *  A particle is influenced by gravity from all the center of masses
   *  of the adjacenct cells and the center of mass of the cell that belongs to.
   *
   *  GF = G * ma * mb / dab^2
   * */
  int i, cellx, celly, cx, cy;
  double interval = 1.0 / ncside;

  #pragma omp parallel for
  for (i = 0; i < n_part; i++) {
    cellx = calc_cell_number(par[i].x, interval, ncside);
    celly = calc_cell_number(par[i].y, interval, ncside);

    int k, j;
    // start on the top left cell in relation to this one
    cx = wrap_around(cellx - 1, 0, ncside - 1);
    cy = wrap_around(celly + 1, 0, ncside - 1);
    for (k = 0; k < 3; k++) {
      for (j = 0; j < 3; j++) {
        // update force
        update_force(&cells[cx][cy], &par[i]);

        // move to next cell on the right
        cx = wrap_around(cx + 1, 0, ncside - 1);
      }
      // move to the row below, left-most cell
      cy = wrap_around(cy - 1, 0, ncside - 1);
      cx = wrap_around(cellx - 1, 0, ncside - 1);
    }
  }
}

void update_vel(double acc_x, double acc_y, particle_t *particle) {
  /*
   * Update the velocity vector of particle i
   * */
  particle->vx += acc_x;
  particle->vy += acc_y;
}

void update_pos(double acc_x, double acc_y, particle_t *particle) {
  /*
   * Update the velocity and then the position vector of particle i
   * */

  // update velocity before updating position
  update_vel(acc_x, acc_y, particle);

  double x_new, y_new;
  x_new = particle->x + particle->vx + (1 / 2) * acc_x;
  y_new = particle->y + particle->vy + (1 / 2) * acc_y;

  // wrap around the positions
  // IMPORTANT - Ranges of the positions -> x = [0, 1] and y = [0, 1]
  if (x_new < 0)
    x_new = 1 - fabs(x_new);
  else if (x_new > 1)
    x_new = x_new - 1;

  if (y_new < 0)
    y_new = 1 - fabs(y_new);
  else if (y_new > 1)
    y_new = y_new - 1;

  particle->x = x_new;
  particle->y = y_new;
}

void calc_all_particle_new_values(long ncside, long long n_part,
                                  particle_t *par) {
  /*
   *  Calculate the new velocity and then the new position of each particle.
   * */
  int i;
  double acc_x, acc_y;

  #pragma omp parallel for
  for (i = 0; i < n_part; i++) {
    // acceleration of the particle
    acc_x = par[i].fx / par[i].m;
    acc_y = par[i].fy / par[i].m;

    // actual cell of the particle
    update_pos(acc_x, acc_y, &par[i]);
  }
}

void init_particle_force(particle_t *par, long long n_part) {
  int i;

  #pragma omp parallel for
  for (i = 0; i < n_part; i++) {
    par[i].fx = 0;
    par[i].fy = 0;
  }
}

/*****************
 * CELL FUNCTIONS
 *****************/

int calc_cell_number(double pos, double interval, long ncside) {
  return labs(((int)floor(pos / interval)) % ncside);
}

void calc_all_cells_cm(long ncside, cell_t **cells, long long n_part,
                       particle_t *par) {
  /*
   *  Determine center of mass of all cells.
   *  For each particle, determine to which cell it belongs,
   *  add to total mass, total number of particles and total position
   *
   *  CM = (for all i SUM((xi, yi) * mi)) / SUM(for all i, mi)
   * */
  int i, j, cellx, celly;
  double interval = 1.0 / ncside;

  #pragma omp parallel for reduction (+:&cells)
  for (i = 0; i < n_part; i++) {
    cellx = calc_cell_number(par[i].x, interval, ncside);
    celly = calc_cell_number(par[i].y, interval, ncside);

    cells[cellx][celly].x += par[i].x * par[i].m;
    cells[cellx][celly].y += par[i].y * par[i].m;
    cells[cellx][celly].m += par[i].m;
    cells[cellx][celly].npar++;
  }
  // after all total masses and positions have been determined
  for (i = 0; i < ncside; i++) {
    for (j = 0; j < ncside; j++) {
      if (cells[i][j].npar != 0) {
        cells[i][j].x = cells[i][j].x / cells[i][j].m;
        cells[i][j].y = cells[i][j].y / cells[i][j].m;
      }
    }
  }
}

void init_cells_matrix(long ncside, cell_t **cells) {
  // default value to 0
  int i, j;

  #pragma omp parallel for
  for (i = 0; i < ncside; i++) {
    for (j = 0; j < ncside; j++) {
      cells[i][j].x = 0;
      cells[i][j].y = 0;
      cells[i][j].m = 0;
      cells[i][j].npar = 0;
    }
  }
}

void create_cells_matrix(long ncside, cell_t **cells) {
  // access cell (x, y) => cells[x][y]
  int i;

  #pragma omp parallel for
  for (i = 0; i < ncside; i++) {
    cells[i] = malloc(ncside * sizeof(cell_t));
  }
}

/*******
 * MAIN
 *******/

int main(int argc, char *argv[]) {

  // receive exactly 4 arguments (first is file name)
  if (argc != 5)
    return EXIT_FAILURE;

  long nseed, ncside, n_tsteps;
  long long n_part;
  particle_t *par;
  cell_t **cells;

  // init values
  nseed = strtol(argv[1], NULL, 10);
  ncside = strtol(argv[2], NULL, 10);
  n_part = strtol(argv[3], NULL, 10);
  n_tsteps = strtol(argv[4], NULL, 10);

  // init particles
  par = malloc(n_part * sizeof(particle_t));
  init_particles(nseed, ncside, n_part, par);

  // init cell matrix
  cells = malloc(ncside * sizeof(cell_t *));
  create_cells_matrix(ncside, cells);
  init_cells_matrix(ncside, cells);

  int i;
  for (i = 0; i < n_tsteps; i++) {
    // determine center of mass of all cells
    calc_all_cells_cm(ncside, cells, n_part, par);

    // compute the gravitational force applied to each particle
    calc_all_particle_force(ncside, cells, n_part, par);

    // print_particles(n_part, par);
    calc_all_particle_new_values(ncside, n_part, par);

    // init cells and particles aplied forces for next timestep
    init_cells_matrix(ncside, cells);
    init_particle_force(par, n_part);
  }

  // Print the desired outputs
  printf("%.2f %.2f\n", par[0].x, par[0].y);
  calc_and_print_overall_cm(n_part, par);

  // free memory
  free_memory(ncside, cells, par);

  return EXIT_SUCCESS;
}
