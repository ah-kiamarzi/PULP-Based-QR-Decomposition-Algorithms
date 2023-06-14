#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__

#define N_ROW 10 // row size
#define N_COL 5 // column size

void qr_household(float Q[][N_ROW], float R[][N_COL]);
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O);

#endif
