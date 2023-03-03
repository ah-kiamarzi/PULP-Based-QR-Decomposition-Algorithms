#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__

#define N_CHANNELS 4 // row size
#define EV_WINDOWS_SIZE 100 // column size


void factorization(float Qt[][N_CHANNELS], float Rt[][EV_WINDOWS_SIZE]);
void givens_rotation(float x, float y, float *c, float *s);

#endif 
