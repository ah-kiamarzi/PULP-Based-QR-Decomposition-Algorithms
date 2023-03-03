#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__

#define N_CHANNELS 100 // row size
#define EV_WINDOWS_SIZE 4 // column size


void qr_gramSmidt(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE], float input[][N_CHANNELS]);
float norm(float *v, int row, int column);

#endif 
