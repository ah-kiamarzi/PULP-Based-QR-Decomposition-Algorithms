#include <stdio.h>
#include <math.h>

#ifndef __QR__
#define __QR__


//There is no need to change the N_ROW and N_COL. it should be the same as number of rows and number of columns of original matrix. 
#define N_ROW 10 // row size
#define N_COL 5 // column size

void qr_gramSmidt(float Q[][N_ROW], float R[][N_COL], float input[][N_ROW]);


#endif
