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

#ifndef SIMPAR_H
#define SIMPAR_H

/*
All type definitions
*/

typedef struct cell_t{
    // if cell has npar == 0, then it does not have any par inside
    double x, y, m;
    int npar;
} cell_t;

#endif
