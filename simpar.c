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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI (3.141592653589793)

#include "init_particles.h"
#include "simpar.h"

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

void calc_and_print_overall_cm(long long n_part, particle_t *par) {
  int i;
  double cmx = 0, cmy = 0, cmm = 0;
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

void update_force(int i, int cellx, int celly, cell_t **cells,
                  particle_t *par) {
  /*
   * Update the gravitacional force of particle i
   * applied by the center of mass from cell (cellx, celly).
   *
   * */
  double dist, fx, fy, magnitude, norm;

  if (cells[cellx][celly].npar != 0) {
    dist = euclidean_distance(cells[cellx][celly].x, par[i].x,
                              cells[cellx][celly].y, par[i].y);

    if (dist >= EPSLON) {
      magnitude = (G * par[i].m * cells[cellx][celly].m) / (dist * dist);
    } else {
      magnitude = 0;
    }

    // determine force direction and norm
    fx = cells[cellx][celly].x - par[i].x;
    fy = cells[cellx][celly].y - par[i].y;
    norm = sqrt(fx * fx + fy * fy);

    if (norm != 0) {
      // normalize direction vector
      fx = fx / norm;
      fy = fy / norm;
    }

    // update particle's force vector
    par[i].fx += (fx * magnitude);
    par[i].fy += (fy * magnitude);
  }
}

int validate_cell(int index, long ncside) {
  if (index < 0)
    return (ncside - 1);
  if (index > (ncside - 1))
    return 0;
  else
    return index;
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

  for (i = 0; i < n_part; i++) {
    cellx = calc_cell_number(par[i].x, interval, ncside);
    celly = calc_cell_number(par[i].y, interval, ncside);

    int k, j;
    // start on the top left cell in relation to this one
    cx = validate_cell(cellx - 1, ncside);
    cy = validate_cell(celly + 1, ncside);
    for (k = 0; k < 3; k++) {
      for (j = 0; j < 3; j++) {
        // update force
        update_force(i, cx, cy, cells, par);

        // move to next cell on the right
        cx = validate_cell(cx + 1, ncside);
      }
      // move to the row below, left-most cell
      cy = validate_cell(cy - 1, ncside);
      cx = validate_cell(cellx - 1, ncside);
    }
  }
}

void update_vel(int i, double acc_x, double acc_y, particle_t *par) {
  /*
   * Update the velocity vector of particle i
   * */
  par[i].vx += acc_x;
  par[i].vy += acc_y;
}

void update_pos(int i, long ncside, double acc_x, double acc_y,
                particle_t *par) {
  /*
   * Update the velocity and then the position vector of particle i
   * */

  // update velocity before updating position
  update_vel(i, acc_x, acc_y, par);

  double x_new, y_new;

  x_new = par[i].x + par[i].vx + (1 / 2) * acc_x;
  y_new = par[i].y + par[i].vy + (1 / 2) * acc_y;

  // wrap around the positions
  // IMPORTANT - Ranges of the positions -> x = [0, 1] and y = [0, 1]
  if (x_new < 0)
    x_new = 1 - fabs(x_new);
  if (y_new < 0)
    y_new = 1 - fabs(y_new);

  if (x_new > 1)
    x_new = x_new - 1;
  if (y_new > 1)
    y_new = y_new - 1;

  par[i].x = x_new;
  par[i].y = y_new;
}

void calc_all_particle_new_values(long ncside, long long n_part,
                                  particle_t *par) {
  /*
   *  Calculate the new velocity and then the new position of each particle.
   * */
  int i;
  double acc_x, acc_y;

  for (i = 0; i < n_part; i++) {
    // acceleration of the particle
    acc_x = par[i].fx / par[i].m;
    acc_y = par[i].fy / par[i].m;

    // actual cell of the particle
    update_pos(i, ncside, acc_x, acc_y, par);
  }
}

void init_cells_matrix(long ncside, cell_t **cells) {
  // access cell (x, y) => cells[x][y]
  int i, j;
  for (i = 0; i < ncside; i++) {
    cells[i] = malloc(ncside * sizeof(cell_t));
  }

  // default value to 0
  for (i = 0; i < ncside; i++) {
    for (j = 0; j < ncside; j++) {
      cells[i][j].x = 0;
      cells[i][j].y = 0;
      cells[i][j].m = 0;
      cells[i][j].npar = 0;
    }
  }
}

void init_cells_and_forces_on_par(cell_t **cells, particle_t *par,
                                  long long n_part, long ncside) {
  // default value to 0
  int i, j;
  for (i = 0; i < ncside; i++) {
    for (j = 0; j < ncside; j++) {
      cells[i][j].x = 0;
      cells[i][j].y = 0;
      cells[i][j].m = 0;
      cells[i][j].npar = 0;
    }
  }

  // dont forget to zero down particles applied force for the next timestep
  for (i = 0; i < n_part; i++) {
    par[i].fx = 0;
    par[i].fy = 0;
  }
}

void free_memory(int ncside, cell_t **cells, particle_t *par) {
  free(par);
  int i;
  for (i = 0; i < ncside; i++)
    free(cells[i]);
  free(cells);
}

int main(int argc, char *argv[]) {

  // receive exactly 4 arguments (first is file name)
  if (argc != 5)
    return EXIT_FAILURE;

  // init values
  long nseed = strtol(argv[1], NULL, 10);
  long ncside = strtol(argv[2], NULL, 10);
  long long n_part = strtol(argv[3], NULL, 10);
  long n_tsteps = strtol(argv[4], NULL, 10);

  int i;
  particle_t *par = malloc(n_part * sizeof(particle_t));
  init_particles(nseed, ncside, n_part, par);

  // init cell matrix
  cell_t **cells = malloc(ncside * sizeof(cell_t *));
  init_cells_matrix(ncside, cells);

  for (i = 0; i < n_tsteps; i++) {
    // determine center of mass of all cells
    calc_all_cells_cm(ncside, cells, n_part, par);

    // compute the gravitational force applied to each particle
    calc_all_particle_force(ncside, cells, n_part, par);

    // print_particles(n_part, par);
    calc_all_particle_new_values(ncside, n_part, par);

    // init cells and particles aplied forces for next timestep
    init_cells_and_forces_on_par(cells, par, n_part, ncside);
  }

  // Print the desired outputs
  printf("%.2f %.2f\n", par[0].x, par[0].y);
  calc_and_print_overall_cm(n_part, par);

  // print_particles(n_part, par);
  // print_cells(ncside, cells);

  // free memory
  free_memory(ncside, cells, par);

  return EXIT_SUCCESS;
}
