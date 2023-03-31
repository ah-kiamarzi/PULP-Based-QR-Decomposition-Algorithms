#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__

//#define N_CHANNELS 4 // row size
//#define EV_WINDOWS_SIZE 100 // column size

#define N_CHANNELS 5 // row size
#define EV_WINDOWS_SIZE 10 // column size



/*void factorization(float Qt[][N_CHANNELS], float Rt[][EV_WINDOWS_SIZE]);
void givens_rotation(float x, float y, float *c, float *s);*/
//void qr_household(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE]);
void qr_household(double Q[][N_CHANNELS], double R[][EV_WINDOWS_SIZE]);
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O);

#endif 
