################################################################################################
#                                                                                              *
#                                          Grupo: 18                                           *
#                                                                                              *
#                                   Beatriz Marques , 80809                                    *
#                                   Carlos  Carvalho, 81395                                    *
#                                   Diogo   Nunes   , 85184                                    *
#                							         	       								   *
#		    Copyright (c) 2019 Beatriz, Carlos e Diogo. All rights reserved.           		   *
#                                                                                              *
################################################################################################


# Makefile
# Parallel and Distributed Computation (2019)

#							-std=gnu99  to allow c++ comments ("//")

CFLAGS = -g -Wall -pedantic -std=gnu99
CC = mpicc

all: clean simpar-mpi

simpar-mpi: simpar-mpi.o init_particles.o
	$(CC) $(CFLAGS) -o simpar-mpi simpar-mpi.o init_particles.o -lm

simpar-mpi.o: simpar-mpi.c simpar-mpi.h init_particles.h
	$(CC) $(CFLAGS) -c simpar-mpi.c

init_particles.o: init_particles.c init_particles.h
	$(CC) $(CFLAGS) -c init_particles.c

clean:
	rm -f *.o simpar *.dSYM
