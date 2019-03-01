#include <stdio.h>
#include <stdlib.h>
#include "init_particles.c"

int main(int argc, char *argv[]) {

    // receive exactly 4 arguments (first is file name)
    if(argc != 5) return EXIT_FAILURE;

    // init values
    particle_t *par;
    nseed = (long) argv[1];
    ncside = (long) argv[2];
    n_part = (long long) argv[3];

    init_particles(nseed, ncside, n_part, par);
    

    /*
        ADD CODE HERE
    */

   
    return EXIT_SUCCESS;
}