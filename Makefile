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
CFLAGS = -g -Wall -pedantic -std=gnu99
CC = gcc

all: simpar

simpar: simpar.o init_particles.o
	$(CC) $(CFLAGS) -o simpar simpar.o init_particles.o

simpar.o: simpar.c simpar.h init_particles.h
	$(CC) $(CFLAGS) -c simpar.c

init_particles.o: init_particles.c init_particles.h
	$(CC) $(CFLAGS) -c init_particles.c

zip: clean
	zip cpd_18.zip simpar.c simpar.h init_particles.c init_particles.h Makefile README.md

clean:
	rm -f *.o simpar
