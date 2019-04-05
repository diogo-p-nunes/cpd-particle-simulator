################################################################################################
#                                                                                              *
#                                          Grupo: 18                                           *
#                                                                                              *
#                                   Beatriz Marques , 80809                                    *
#                                   Carlos  Carvalho, 81395                                    *
#                                   Diogo   Nunes   , 85184                                    *
#                							         	       *
#		    Copyright (c) 2019 Beatriz, Carlos e Diogo. All rights reserved.           *
#                                                                                              *
################################################################################################


# Makefile
# Parallel and Distributed Computation (2019)

#							-std=gnu99  to allow c++ comments ("//")
CFLAGS = -fopenmp -g -pedantic -std=gnu99
CC = gcc-8

all: clean simpar-omp

simpar-omp: simpar-omp.o init_particles.o
	$(CC) $(CFLAGS) -o simpar-omp simpar-omp.o init_particles.o -lm

simpar-omp.o: simpar-omp.c simpar-omp.h init_particles.h
	$(CC) $(CFLAGS) -c simpar-omp.c

init_particles.o: init_particles.c init_particles.h
	$(CC) $(CFLAGS) -c init_particles.c

clean:
	rm -rf *.o simpar-omp *.inc *.rc *.mod.c *.tab.c *.dSYM
