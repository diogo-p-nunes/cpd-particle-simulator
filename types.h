/*
All type definitions
*/

typedef struct cell_t {
    // if cell has npar == 0, then it does not have any par inside
	double x, y, m;
    int npar;
} cell_t;