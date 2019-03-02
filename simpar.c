#include <stdio.h>
#include <stdlib.h>
#include "init_particles.c"


void print_particles(long long n_part, particle_t *par) {
    printf("[Particles]\n");
    int i;
    for(i = 0; i < n_part; i++) {
        printf("\t(%f, %f) - %f\n", par[i].x, par[i].y, par[i].m);
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

    print_particles(n_part, par);

    return EXIT_SUCCESS;
}