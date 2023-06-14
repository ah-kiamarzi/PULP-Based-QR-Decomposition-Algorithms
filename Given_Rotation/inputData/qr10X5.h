#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__


#define N_ROW 10 // row size
#define N_COL 5 // column size

void factorization(float Qt[][N_ROW], float Rt[][N_COL]);
void givens_rotation(float x, float y, float *c, float *s);

#endif
