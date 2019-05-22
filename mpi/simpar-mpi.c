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
#define in_process(cx,cy,cols,rows) (cx >= 0 && cx <= cols && cy >= 0 && cy <= rows)
#define get_id(dim0, dim1, x, y, ncside) ((dim0*x)/ncside + ((dim1*y)/ncside)*dim0)



/**********************************
 * DATA GENERATION AND DISTRIBUTION
 **********************************/

void calculate_c(int *cx, int *cy, int id, int ncside, int*dims){

    /*
     *  Calculates range of X and Y coordinates of given process with ID = id
     *
     */

    int x, y, i, last, y_aux, x_aux;
    int entered = 0;

    for(x = 0; x < ncside; x++) {
        for (y = 0; y < ncside; y++) {
            if(id == get_id(dims[0],dims[1],x,y,ncside)) {
                cx[0] = x;
                cy[0] = y;

                last = y;

                // fix x
                for(i = 0; i < ncside; i++){
                    y_aux = y+i+1;
                    if(id != get_id(dims[0],dims[1],x,y_aux,ncside)) cy[1] = last;
                    else last = y_aux;
                }
                last = x;

                // fix y
                for(i = 0; i < ncside; i++){
                    x_aux = x+i+1;
                    if(id != get_id(dims[0],dims[1],x_aux,y,ncside)) cx[1] = last;
                    else last = x_aux;
                }

                entered = 1;
            }
            if(entered) break;
        }
        if(entered) break;
    }
}

void init_processors_particles(particle_t **arr, int *processors_particles_sizes, int p) {
    for(int i = 0; i<p; i++) {
        arr[i] = NULL;
        processors_particles_sizes[i] = 0;
    }
}

void generate_sending_data(long long n_part, particle_t* par, long ncside, particle_t **processors_particles,
        int* processors_particles_sizes, int*dims) {

    // go through each particle and determine to which processor it belongs to
    // add that particle to the dynamic array of particles of that processor

    int cellx, celly, id;
    double interval = 1.0 / ncside;

    for(int i=0; i<n_part; i++) {
        cellx = calc_cell_number(par[i].x, interval, ncside);
        celly = calc_cell_number(par[i].y, interval, ncside);
        id = get_id(dims[0], dims[1], cellx, celly, ncside);
        add_particle_to_processor_array(&processors_particles[id], &par[i], id, processors_particles_sizes);
    }
}

void add_particle_to_processor_array(particle_t **array, particle_t *par, int id, int* processors_particles_sizes) {

    if(processors_particles_sizes[id] == 0)
        (*array) = (particle_t *) malloc(sizeof(particle_t));

    else
        (*array) = (particle_t *) realloc((*array), (processors_particles_sizes[id] + 1) * sizeof(particle_t));

    (*array)[processors_particles_sizes[id]] = (*par);
    processors_particles_sizes[id] += 1;

}

void build_contiguous_array(cell_t *send_cells, cell_t **cells, int tsize, int i, int j, int rows, int id) {

    // we only have this problem when (i=0,j=1) or (i=1,j=0)
    int entry = j==0 ? 0 : rows-1;

    for(int k=0; k<tsize; k++) {
        send_cells[k] = cells[k][entry];
    }
}

void send_and_receive_cells(cell_t **id_received_cells_map, int *counter, cell_t **cells, int*dims,
        int *cx, int *cy, int ncside, int cols, int rows, MPI_Datatype mpi_cell_type, int id, int p) {

    int tx, ty, i, j, k;
    int tid, tsize;
    int CELLS_TAG = 666;

    for(i = 0; i < p; i++) counter[i] = 0;

    // SENDING AND RECEIVING ALL DATA CELLS NECESSARY BY OTHERS AND ME
    for(i=0; i<2; i++) {
        for(j=0; j<2; j++) {
            for(k=0; k<2; k++) {
                // i -> cx, j -> cy

                // determine direction im sending data
                int xdir = i == 0 ? -1 : 1;
                int ydir = j == 0 ? -1 : 1;
                tsize = 1;

                // if not diagonal
                if(k == 1) {
                    xdir = (i+j) == 1 ? 0 : xdir;
                    ydir = (i+j) == 0 || (i+j) == 2 ? 0 : ydir;
                    tsize = (i+j) == 1 ? cols : rows;
                }

                // target processor info
                tx = wrap_around(cx[i]+xdir, 0, ncside-1);
                ty = wrap_around(cy[j]+ydir, 0, ncside-1);
                tid = get_id(dims[0],dims[1],tx,ty,ncside);

                // data to send
                cell_t *send_cells;
                if(((i == 0 && j == 1) || (i == 1 && j == 0)) && (k==1)) {
                    send_cells = (cell_t*) malloc(tsize * sizeof(cell_t));
                    build_contiguous_array(send_cells, cells, tsize, i, j, rows, id);
                }
                else {
                    send_cells = &cells[i*(cols-1)][0];
                }

                if(tsize == 1) {
                    int entry = j==0 ? 0 : rows-1;
                    send_cells = &cells[i*(cols-1)][entry];
                }

                // Async - Send - and Wait
                MPI_Send(send_cells, tsize, mpi_cell_type, tid, CELLS_TAG, MPI_COMM_WORLD);

                //printf("id:%d sent - id:%d - tsize:%d\n", id, tid, tsize);
                //print_cells(id, cx, cy, &send_cells, 1, tsize);
                //flush(stdout);
            }
        }
    }


    // RECEIVING ALL DATA CELLS NECESSARY FOR ME
    for(i=0; i<2; i++) {
        for(j=0; j<2; j++) {
            for(k=0; k<2; k++) {
                // i -> cx, j -> cy

                // determine direction im receiving data
                int xdir = i == 0 ? -1 : 1;
                int ydir = j == 0 ? -1 : 1;
                tsize = 1;

                // if not diagonal
                if(k == 1) {
                    xdir = (i+j) == 1 ? 0 : xdir;
                    ydir = (i+j) == 0 || (i+j) == 2 ? 0 : ydir;
                    tsize = (i+j) == 1 ? cols : rows;
                }

                // target processor info
                tx = wrap_around(cx[i]+xdir, 0, ncside-1);
                ty = wrap_around(cy[j]+ydir, 0, ncside-1);
                tid = get_id(dims[0],dims[1],tx,ty,ncside);

                if(counter[tid] > 0){

                    cell_t* tmp_aux = (cell_t*) realloc(id_received_cells_map[tid], (counter[tid] + tsize) * sizeof(cell_t));
                    cell_t* tmp_receive = (cell_t*) malloc(tsize * sizeof(cell_t));

                    id_received_cells_map[tid] = tmp_aux;

                    MPI_Recv(tmp_receive, tsize, mpi_cell_type, tid, CELLS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    int c = 0;
                    for(int m = counter[tid]; m < counter[tid]+tsize; m++)
                        id_received_cells_map[tid][m] = tmp_receive[c++];

                    counter[tid] += tsize;
                }
                else {
                    id_received_cells_map[tid] = (cell_t*) malloc(tsize * sizeof(cell_t));
                    counter[tid] = tsize;
                    MPI_Recv(id_received_cells_map[tid], tsize, mpi_cell_type, tid, CELLS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                //printf("id:%d received - id:%d\n", id, tid);
                //fflush(stdout);
                //print_cells(id, cx, cy, &id_received_cells_map[tid], 1, tsize);
            }
        }
    }


}

void send_and_receive_moved_particles(particle_t **par, long long *num_par, long ncside, int *dims, int id, int p,
        MPI_Datatype mpi_particle_type) {

    //printf("id:%d Sending particles ...\n", id);
    // automatically set to zero
    long long *counter = (long long*) calloc(p, sizeof(long long));
    particle_t **sending_part = (particle_t**) malloc(p * sizeof(particle_t*));

    int x, y, tid;
    double interval = 1.0 / ncside;

    //printf("INICIO id:%d - num_par:%lld\n", id, (*num_par));

    // determine what and how many particles to send to each processor
    for(int i=0; i<(*num_par); i++) {
        // real cell coordinates
        //printf("interval=%f ncside=%ld\n", interval, ncside);
        x = calc_cell_number((*par)[i].x, interval, ncside);
        y = calc_cell_number((*par)[i].y, interval, ncside);
        tid = get_id(dims[0], dims[1], x, y, ncside);

        //printf("id=%d p=(%f,%f) c=(%d,%d) tid=%d\n", id, (*par)[i].x, (*par)[i].y, x, y, tid);

        if(counter[tid] == 0) {
            sending_part[tid] = (particle_t*) malloc(sizeof(particle_t));
        }
        else {
            sending_part[tid] = (particle_t*) realloc(sending_part[tid], (counter[tid] + 1) * sizeof(particle_t));
        }
        counter[tid] += 1;
        sending_part[tid][counter[tid]-1] = (*par)[i];

    }

    (*num_par) = counter[id];
    (*par) = (particle_t*) malloc((*num_par) * sizeof(particle_t));
    (*par) = sending_part[id];

    // send the particles
    int PART_TAG = 999;
    for(int i=0; i<p; i++) {
        if(i==id) continue;

        int new_num_par = 0;

        // First each process needs to know how many particles it is going to receive in the array
        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Irecv(&new_num_par, 1, MPI_INT, i, PART_TAG, MPI_COMM_WORLD, &request);
        MPI_Send(&counter[i], 1, MPI_INT, i, PART_TAG, MPI_COMM_WORLD);
        MPI_Wait(&request, MPI_STATUS_IGNORE);

        request = MPI_REQUEST_NULL;
        long long old_size;
        particle_t* tmp_receive;

        if(new_num_par > 0) {
            //printf("[RECEIVING] %d -> id=%d new_num_par=%lld\n", i, id, new_num_par);
            fflush(stdout);

            old_size = (*num_par);
            (*num_par) = (*num_par) + new_num_par;

            particle_t* tmp_aux = (particle_t*) realloc((*par), (*num_par) * sizeof(particle_t));
            tmp_receive = (particle_t*) malloc(new_num_par * sizeof(particle_t));
            (*par) = tmp_aux;
            MPI_Irecv(tmp_receive, new_num_par, mpi_particle_type, i, PART_TAG, MPI_COMM_WORLD, &request);
        }

        if(counter[i] > 0) {
            MPI_Send(sending_part[i], counter[i], mpi_particle_type, i, PART_TAG, MPI_COMM_WORLD);
            //printf("id=%d -> %d new_num_par=%lld\n", id, i, counter[i]);
        }

        if(new_num_par > 0) {
            MPI_Wait(&request, MPI_STATUS_IGNORE);

            int c = 0;
            for(int m = old_size; m < (*num_par); m++) (*par)[m] = tmp_receive[c++];
        }

    }


    //printf("FIM id:%d - num_par:%lld\n", id, (*num_par));
}



/*****************
 * CELL FUNCTIONS
 *****************/

int calc_cell_number(double pos, double interval, long ncside) {
    return labs(((int) floor(pos / interval)) % ncside);
}

int calc_processor_cell_number(double pos, double interval, long ncside, int max) {
    return labs((((int) floor(pos / interval)) % ncside) % max);
}

void init_cells_matrix(int cols, int rows, cell_t **cells, int xmin, int ymin) {
    int i, j;
    cell_t *cell;

    for (i = 0; i < cols; i++) {
        for (j = 0; j < rows; j++) {
            cell = &cells[i][j];
            cell->x = 0.0;
            cell->y = 0.0;
            cell->m = 0.0;
            cell->npar = 0;

            cell->cx = i+xmin;
            cell->cy = j+ymin;
        }
    }
}

void create_cells_matrix(int cols, int rows, cell_t **cells) {
    // access cell (x, y) => cells[x][y]
    for (int i = 0; i < cols; i++)
        cells[i] = malloc(rows * sizeof(cell_t));
}

void build_mpi_particle_type(MPI_Datatype *mpi_particle_type) {

    // Data type we want to duplicate into MPI
    particle_t *par;

    // Create a MPI type for struct particle_t
    // Build a derived datatype consisting of three doubles, a long and one mpi_particle_type
    const int nitems=8;

    // Specify the number of elements of each type
    int blocklengths[8] = {1,1,1,1,1,1,1,1};

    // First specify the types
    MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG_LONG};

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
    MPI_Get_address(&(par->tag), &addresses[8]);

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

void build_mpi_cell_type(MPI_Datatype *mpi_cell_type) {

    // Data type we want to duplicate into MPI
    cell_t *cell;

    // Create a MPI type for struct particle_t
    // Build a derived datatype consisting of three doubles and a long
    const int nitems=6;

    // Specify the number of elements of each type
    int blocklengths[6] = {1,1,1,1,1,1};

    // First specify the types
    MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_LONG};

    MPI_Aint displacements[6];
    MPI_Aint addresses[7];

    // Calculate the displacements of the members relative to cell
    MPI_Get_address(cell, &addresses[0]);
    MPI_Get_address(&(cell->x), &addresses[1]);
    MPI_Get_address(&(cell->y), &addresses[2]);
    MPI_Get_address(&(cell->m), &addresses[3]);
    MPI_Get_address(&(cell->cx), &addresses[4]);
    MPI_Get_address(&(cell->cy), &addresses[5]);
    MPI_Get_address(&(cell->npar), &addresses[6]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];
    displacements[2] = addresses[3] - addresses[0];
    displacements[3] = addresses[4] - addresses[0];
    displacements[4] = addresses[5] - addresses[0];
    displacements[5] = addresses[6] - addresses[0];


    // Create the derived type
    MPI_Type_create_struct(nitems, blocklengths, displacements, types, mpi_cell_type);

    // Commit it so that it can be used
    MPI_Type_commit(mpi_cell_type);
}





/*****************************
 * PROBLEM SPECIFIC FUNCTIONS
 *****************************/

void calc_all_cells_cm(cell_t **cells, long long n_part, particle_t *par, int cols, int rows,
        long ncside, double interval) {
    /*
     *  Determine center of mass of all cells of given processor.
     *  For each particle, determine to which cell it belongs,
     *  add to total mass, total number of particles and total position
     *
     *  CM = (for all i SUM((xi, yi) * mi)) / SUM(for all i, mi)
     * */

    int i, j, x, y;
    cell_t *cell;

    for (i = 0; i < n_part; i++) {
        x = calc_processor_cell_number(par[i].x, interval, ncside, cols);
        y = calc_processor_cell_number(par[i].y, interval, ncside, rows);
        cell = &cells[x][y];

        cell->x += par[i].x * par[i].m;
        cell->y += par[i].y * par[i].m;
        cell->m += par[i].m;
        cell->npar++;
    }

    // after all total masses and positions have been determined
    for (i = 0; i < cols; i++) {
        for (j = 0; j < rows; j++) {
            cell = &cells[i][j];
            if (cell->npar != 0) {
                cell->x = cell->x / cell->m;
                cell->y = cell->y / cell->m;
            }
        }
    }
}

void update_force(cell_t *cell, particle_t *particle) {
    /*
     * Update the gravitacional force of particle i
     * applied by the center of mass from cell (cellx, celly).
     *
     * */
    double dist, fx, fy, magnitude, norm;

    dist = euclidean(cell->x, particle->x, cell->y, particle->y);

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

void calc_all_particle_force(long ncside, cell_t **cells, long long n_part, int *cx, int *cy,particle_t *par,
        int cols, int rows, cell_t **id_received_cells_map, int *dims, int*counter) {
    /*
     *  Determine gravitational force for each particle.
     *  A particle is influenced by gravity from all the center of masses
     *  of the adjacenct cells and the center of mass of the cell that belongs to.
     *
     *  GF = G * ma * mb / dab^2
     */
    int i, x, y, c_x, c_y, index = 0;
    cell_t *c;

    double interval = 1.0 / ncside;

    for (i = 0; i < n_part; i++) {
        x = calc_cell_number(par[i].x, interval, ncside);
        y = calc_cell_number(par[i].y, interval, ncside);

        int k, j, l;

        // start on the top left cell in relation to this one
        c_x = wrap_around(x - 1, 0, ncside - 1); // real coordinates
        c_y = wrap_around(y + 1, 0, ncside - 1); // real coordinates

        for (k = 0; k < 3; k++) {
            for (j = 0; j < 3; j++) {

                // check if coordinates are in process or not
                if(!in_process(c_x-cx[0], c_y-cy[0], cols-1, rows-1)){
                    int id = get_id(dims[0],dims[1],c_x,c_y,ncside);

                    for(l = 0; l < counter[id]; l++)
                        if(id_received_cells_map[id][l].cx == c_x && id_received_cells_map[id][l].cy == c_y) {
                            c = &id_received_cells_map[id][l];
                            break;
                        }

                    update_force(c, &par[i]);
                }
                else {
                    update_force(&cells[c_x-cx[0]][c_y-cy[0]], &par[i]); // convert to process coordinates
                }

                // move to next cell on the right
                c_x = wrap_around(c_x + 1, 0, ncside - 1);
            }
            // save the last ones

            // move to the row below, left-most cell
            c_y = wrap_around(c_y - 1, 0, ncside - 1);
            c_x = wrap_around(x - 1, 0, ncside - 1);
        }
    }
}

void calc_all_particle_new_values(long ncside, long long n_part, particle_t *par) {
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
        update_pos(acc_x, acc_y, &par[i]);
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

void init_particle_force(particle_t *par, long long num_par, int id, particle_t **containsParticleZero) {
    for(int i=0; i<num_par; i++) {
        //printf("id:%d i:%d\n", id, i);
        //fflush(stdout);
        par[i].fx = 0;
        par[i].fy = 0;
        if(par[i].tag == 0) {
            (*containsParticleZero) = &par[i];
        }
    }
}

void calc_and_print_overall_cm(long long n_part, particle_t *par, double **result) {
    int i;
    double cmx = 0, cmy = 0, cmm = 0;

    for (i = 0; i < n_part; i++) {
        cmx += par[i].x * par[i].m;
        cmy += par[i].y * par[i].m;
        cmm += par[i].m;
    }

    (*result)[0] = cmx;
    (*result)[1] = cmy;
    (*result)[2] = cmm;
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

    // Define the number of particles for each processor (local)
    long long num_par = 0;

    // Define data sending TAGS for validation: IMPORTANT
    int NUM_PAR_TAG = 0;
    int PAR_TAG = 1;

    // build MPI data types for transfer
    MPI_Datatype mpi_particle_type, mpi_cell_type;
    build_mpi_particle_type(&mpi_particle_type);
    build_mpi_cell_type(&mpi_cell_type);

    long long max_generating_part = 10 * 1000000;
    long long processed_n_part = 0;
    long long processing_n_part = n_part < max_generating_part ? n_part : max_generating_part;

    if(id != 0) {
        // Receive generated particles in separate iterations to avoid memory overflow
        while (processed_n_part < n_part) {
            MPI_Request request = MPI_REQUEST_NULL;
            int new_num_par = 0;

            // Receive respective particles from processor id=0 who initialized the system
            // First each process needs to know how many particles it is going to receive in the array
            MPI_Irecv(&new_num_par, 1, MPI_INT, 0, NUM_PAR_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            request = MPI_REQUEST_NULL;

            //printf("id:%d receiveing new_num_par = %d\n", id, new_num_par);
            //fflush(stdout);

            // allocate memory for incoming number of particles
            if(num_par == 0) {
                par = (particle_t*) malloc(new_num_par * sizeof(particle_t));
            }
            else {
                par = (particle_t*) realloc(par, (num_par+new_num_par) * sizeof(particle_t));
            }

            // all other processors receive the distributed data, respectively
            MPI_Irecv(&par[num_par], new_num_par, mpi_particle_type, 0, PAR_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);

            //printf("id:%d done received new_num_par = %d\n", id, new_num_par);
            //fflush(stdout);

            // update number of particles of processor
            num_par += new_num_par;

            // update for the next iteration of particles
            processed_n_part += processing_n_part;
            processing_n_part = (n_part - processed_n_part) < max_generating_part ? (n_part - processed_n_part) : max_generating_part;
        }
    }
    if(id == 0) {
        // init particles locally
        particle_t *aux_par = malloc(n_part * sizeof(particle_t));
        init_particles(nseed, ncside, n_part, aux_par);

        // Generate sending data in separate loops to avoid memory overflow
        while (processed_n_part < n_part) {

            printf("processed:%lld, processing now:%lld\n", processed_n_part, processing_n_part);
            fflush(stdout);

            // init array with particles that belong to each processor id
            particle_t **processors_particles  = (particle_t**) malloc(p * sizeof(particle_t*));

            // keep in memory the number of particles in each processor id particle array
            // initialize array and sizes side array
            int *processors_particles_sizes = (int*) malloc(p * sizeof(int));
            init_processors_particles(processors_particles, processors_particles_sizes, p);

            //printf("processors particles init done\n");
            //fflush(stdout);

            // populate the processors particles array with the corresponding data of each processor
            generate_sending_data(processing_n_part, &aux_par[processed_n_part], ncside, processors_particles, processors_particles_sizes, dims);

            //printf("generate data done\n");
            //fflush(stdout);

            // Distribute particles across all processors
            for(int i=1; i<p; i++) {
                //printf("Sending processor %d num_par = %d\n", i, processors_particles_sizes[i]);
                // First each process needs to know how many particles it is going to receive in the array
                MPI_Send(&processors_particles_sizes[i], 1, MPI_INT, i, NUM_PAR_TAG, MPI_COMM_WORLD);

                // Now we can send the particles to the processors
                MPI_Send(processors_particles[i], processors_particles_sizes[i], mpi_particle_type, i, PAR_TAG, MPI_COMM_WORLD);
            }

            //printf("sent done\n");
            //fflush(stdout);

            // Now processor id=0 can set its own particles (no longer knows about all particles)
            int new_num_par = processors_particles_sizes[0];
            if(num_par == 0) {
                par = (particle_t*) malloc(new_num_par * sizeof(particle_t));
            }
            else {
                par = (particle_t*) realloc(par, (num_par+new_num_par) * sizeof(particle_t));
            }

            // populate id=0 particle array
            int c = 0;
            for(int m=num_par; m < (num_par + new_num_par); m++) {
                par[m] = processors_particles[0][c++];
            }

            // update number of particles of id=0
            num_par += new_num_par;

            // free the used memory
            for(int i=1; i<p; i++) free(processors_particles[i]);
            free(processors_particles_sizes);

            //printf("populate 0 done\n");
            //fflush(stdout);

            // update for the next iteration of particles
            processed_n_part += processing_n_part;
            processing_n_part = (n_part - processed_n_part) < max_generating_part ? (n_part - processed_n_part) : max_generating_part;
        }
    }

    //printf("id:%d Going on\n", id);
    //fflush(stdout);
    //MPI_Barrier(MPI_COMM_WORLD);

    /* From here on the computation is done identically at each process. */

    // Determine which cells belong to me
    int *cx = (int*) malloc(2*sizeof(int));
    int *cy = (int*) malloc(2*sizeof(int));
    calculate_c(cx, cy, id, ncside, dims);
    int cols = cx[1] - cx[0] +1;
    int rows = cy[1] - cy[0] +1;

    // init cells matrix
    cell_t **cells;
    cells = (cell_t**) malloc(cols * sizeof(cell_t*));
    create_cells_matrix(cols, rows, cells);
    init_cells_matrix(cols, rows, cells, cx[0], cy[0]);

    double interval = 1.0 / ncside;

    int* counter = (int*) calloc(p, sizeof(int));
    particle_t *containsParticleZero;

    for (int i = 0; i < n_tsteps; i++) {
        // determine center of mass of all cells
        calc_all_cells_cm(cells, num_par, par, cols, rows, ncside, interval);

        //MPI_Barrier(MPI_COMM_WORLD);

        //printf("id:%d Sending cells\n", id);
        //fflush(stdout);
        // First we have to receive all the cells that we will need from the surrounding cells
        cell_t **id_received_cells_map = (cell_t**) malloc(p * sizeof(cell_t**));
        send_and_receive_cells(id_received_cells_map, counter, cells, dims, cx, cy, ncside, cols, rows, mpi_cell_type, id, p);
        //printf("id:%d Received cells\n", id);
        //fflush(stdout);

        //MPI_Barrier(MPI_COMM_WORLD);

        calc_all_particle_force(ncside, cells, num_par, cx, cy, par, cols, rows, id_received_cells_map,dims, counter);
        calc_all_particle_new_values(ncside, num_par, par);

        MPI_Barrier(MPI_COMM_WORLD);

        // determine if any particle has moved to another processors designated cells
        // if yes, send it
        send_and_receive_moved_particles(&par, &num_par, ncside, dims, id, p, mpi_particle_type);

        MPI_Barrier(MPI_COMM_WORLD);

        // flag to know which processor has particle zero
        containsParticleZero = NULL;

        // init cells and particles aplied forces for next timestep
        init_particle_force(par, num_par, id, &containsParticleZero);
        init_cells_matrix(cols, rows, cells, cx[0], cy[0]);

        free(id_received_cells_map);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Print the desired outputs
    if(containsParticleZero != NULL) {
        printf("%.2f %.2f\n", containsParticleZero->x, containsParticleZero->y);
        fflush(stdout);
    }

    // calc overall CM
    int PROC_CM = 909;
    double *result = (double*) malloc(3 * sizeof(result));
    calc_and_print_overall_cm(num_par, par, &result);

    if(id != 0) {
        // First each process needs to know how many particles it is going to receive in the array
        MPI_Send(result, 3, MPI_DOUBLE, 0, PROC_CM, MPI_COMM_WORLD);
    }
    else {
        // Receive and calc overall CM
        double cmx = result[0];
        double cmy = result[1];
        double cmm = result[2];

        for(int i=1; i<p; i++) {
            MPI_Recv(result, 3, MPI_DOUBLE, i, PROC_CM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            cmx += result[0];
            cmy += result[1];
            cmm += result[2];

            result[0] = 0;
            result[1] = 0;
            result[2] = 0;
        }

        printf("%.2f %.2f\n", cmx/cmm, cmy/cmm);
        fflush(stdout);
    }

    free(par);
    free(cells);
    free(counter);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
