#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__

#define N_CHANNELS 16 // row size
#define EV_WINDOWS_SIZE 100 // column size

void qr_household(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE]);
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O);

#endif 
