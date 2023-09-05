#include "pmsis.h"
#include "util.h"
#include "inputData/qr40X40.h"
#include "inputData/data40X40.inc"
#include "math.h"

#define MAX(a,b) ((a) > (b))? (a): (b)

PI_L1 float aR[MAX(N_COL,N_ROW)];
PI_L1 float bR[MAX(N_COL,N_ROW)];
PI_L1 float aQ[MAX(N_COL,N_ROW)];
PI_L1 float bQ[MAX(N_COL,N_ROW)];


//Q is transposed
PI_L1 float Q[N_ROW][N_ROW];
PI_L1 float R[N_ROW][N_COL];

PI_L1 float one = 1.0f;
PI_L1 float zero = 0.0f;
PI_L1 float two = 2.0f;

PI_L1 float c=0;
PI_L1 float s=0;
PI_L1 float xi;
PI_L1 float xj;

PI_L1 int ii = 0;

float res[N_ROW][N_COL];

int numBarr = 0;

definePrefArr


void cluster_main();

void pe_entry(void *arg){
	cluster_main();
}

void cluster_entry(void *arg){
    pi_cl_team_fork((NUM_CORES), pe_entry, 0);
}

static int test_entry(){

    struct pi_device cluster_dev;
    struct pi_cluster_conf cl_conf;
    struct pi_cluster_task cl_task;

    pi_cluster_conf_init(&cl_conf);
    pi_open_from_conf(&cluster_dev, &cl_conf);
    if (pi_cluster_open(&cluster_dev)){
        return -1;
    }

    pi_cluster_send_task_to_cl(&cluster_dev, pi_cluster_task(&cl_task, cluster_entry, NULL));

    pi_cluster_close(&cluster_dev);

    return 0;
}

static void test_kickoff(void *arg){
    int ret = test_entry();
    pmsis_exit(ret);
}

int main(){
 
 
    return pmsis_kickoff((void *)test_kickoff);
}

void cluster_main(){


    pi_perf_conf(
    (1<<PI_PERF_CYCLES) | 
    (1<<PI_PERF_INSTR)  | 
    (1<<PI_PERF_ACTIVE_CYCLES) | 
    (1<<PI_PERF_LD_EXT) | 
    (1<<PI_PERF_TCDM_CONT) | 
    (1<<PI_PERF_LD_STALL) | 
    (1<<PI_PERF_IMISS) |
    (1<<0x11) | 
    (1<<0x12) |
    (1<<0x13) | 
    (1<<0x14));
   
    pi_cl_team_barrier();

    if (pi_core_id()==0){
		initialMatrixEYE(&Q[0][0],N_ROW,N_ROW);
		initialMatrixMatrix(&R[0][0],N_ROW,N_COL,&input[0][0]);
    }


	pi_cl_team_barrier();
		
	pi_perf_reset();
	pi_perf_start();

	factorization(Q, R);

	pi_perf_stop();
	printPerfCounters
	pi_cl_team_barrier();



	// for(ii = 0; ii < 1; ii++){
	// 	if (pi_core_id()==0){
	// 		initialMatrixEYE(&Q[0][0],N_ROW,N_ROW);
	// 		initialMatrixMatrix(&R[0][0],N_ROW,N_COL,&input[0][0]);
	// 	}		
	// 	pi_cl_team_barrier();
	// 	if(ii == 0){
	// 		pi_cl_team_barrier();
				
	// 		pi_perf_reset();
	// 		pi_perf_start();

	// 		factorization(Q, R);

	// 		pi_perf_stop();
	// 		printPerfCounters
	// 		pi_cl_team_barrier();

	// 	}else{
	// 		pi_cl_team_barrier();

	// 		factorization(Q, R);

	// 		pi_cl_team_barrier();
	// 	}

	// }

	if(pi_core_id() == 0){
		printf("NUMBARR = %d\n",numBarr);
	}

    pi_cl_team_barrier();
    
	#ifdef DEBUG
 	if(pi_core_id()==0){
		
		printMatrix(&Q[0][0],N_ROW,N_ROW,"Q",1);
		printMatrix(&R[0][0],N_ROW,N_COL,"R",0);

		for (int i = 0; i < N_ROW; i++) {
			for (int j = 0; j < N_COL; j++) {
				res[i][j] = 0;
				for (int k = 0; k < N_ROW; k++) {
					res[i][j] += Q[k][i] * R[k][j];
				}
			}
		}
		printMatrix(&res[0][0],N_ROW,N_COL,"Res",0);
	}

	pi_cl_team_barrier();
	#endif
    
	return;
}

inline float Sqrt(float x) {
        float res;
        asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
        return res;
}


void inline givens_rotation(float x, float y, float *c, float *s){
    float r;
    
    if ((fabs(x)) >= (fabs(y))){ 
        r = x / y;
        *s = one/sqrtf(one + r*r);
        *c = (*s)*r;
	}
    else{
        r = y / x;
        *c = one /sqrtf(one + r*r);
        *s = (*c)*r;
	}
    return;
}

void factorization(float Q[][N_ROW], float R[][N_COL]){
	int core_id = pi_core_id();
	
	int i;

	
	#if NUM_CORES > 1
	int blockSize_column = N_COL/NUM_CORES;
	int start_column = core_id*blockSize_column;

	int blockSize_row = N_ROW/ NUM_CORES;
	int start_row = core_id*blockSize_row; 

	if(core_id==(NUM_CORES - 1)){
		blockSize_column = N_COL - (NUM_CORES - 1) * blockSize_column;
		blockSize_row = N_ROW - (NUM_CORES - 1) * blockSize_row;
	}
	int end_row = start_row + blockSize_row;
	int end_column = start_column + blockSize_column;
	int start = 0;

	#endif

	for(int k=0; k<N_COL; k++){
		#if NUM_CORES > 1

		if(start_column >= k){
			start = start_column;
		}else if(end_column >= k){
			start = k;
		}else{
			start = end_column+1;
		}
		#endif

		for(int j=k+1; j<N_ROW; j++){ 

			if(R[j][k]!=0){
				if(core_id == 0){
					xi=R[k][k];
					xj=R[j][k];
					givens_rotation(xi, xj, &c, &s);
				}

				#if NUM_CORES > 1
				for(i = start; (i<end_column); i++){
					float Rki0 = R[k][i];
					float Rji0 = R[j][i];
					aR[i] = Rki0;
					bR[i] = Rji0;
				}

				for(i = start_row;(i+1<end_row); i+=2){
					float a0 = Q[k][i];
					float a1 = Q[k][i+1];
					aQ[i] = a0;
					aQ[i+1] = a1;
				}
				if((i < end_row)){
					float a0 = Q[k][i];
					aQ[i] = a0;
				}

				for(i = start_row; (i+1<end_row); i+=2){
					float b0 = Q[j][i];
					float b1 = Q[j][i+1];						
					bQ[i] = b0;
					bQ[i+1] = b1;
				}
				if((i < end_row)){
					float b0 = Q[j][i];
					bQ[i] = b0;
				}


				pi_cl_team_barrier();
				BarrierCounter

				for(i = start; (i<end_column); i++){
					float ai0 = aR[i]; 
					float bi0 = bR[i];
					R[k][i]     = (c * ai0) + (s * bi0);
					R[j][i]     = (c * bi0) - (s * ai0);
				}

				#else

				for(i=k; i+1<N_COL; i+=2){

					float Rki0 = R[k][i]; 
					float Rji0 = R[j][i]; 

					float Rki1 = R[k][i+1]; 
					float Rji1 = R[j][i+1]; 

					R[k][i]       = (c * Rki0) + (s * Rji0);
					R[j][i]       = (c * Rji0) - (s * Rki0);
					R[k][i+1]     = (c * Rki1) + (s * Rji1);
					R[j][i+1]     = (c * Rji1) - (s * Rki1);
				}
				if(i < N_COL){
					float Rki0 = R[k][i]; 
					float Rji0 = R[j][i]; 
					
					R[k][i]       = (c * Rki0) + (s * Rji0);
					R[j][i]       = (c * Rji0) - (s * Rki0);

				}

				#endif
				
				
				#if NUM_CORES > 1

				// for(i = start_row; (i+1<end_row); i+=2){
				// 	float a0 = aQ[i];float b0 = bQ[i];
				// 	float a1 = aQ[i+1];float b1 = bQ[i+1];
				// 	Q[k][i] = (c * a0)+(s * b0);
				// 	Q[k][i+1] = (c * a1)+(s * b1);
				// 	Q[j][i] = (c * b0)-(s * a0);
				// 	Q[j][i+1] = (c * b1)-(s * a1);
				// }
				// if((i < end_row)){
				// 	float a0 = aQ[i];float b0 = bQ[i];
				// 	Q[k][i] = (c * a0)+(s * b0);
				// 	Q[j][i] = (c * b0)-(s * a0);
				// }


				for(i = start_row; (i<end_row); i++){
					float a0 = aQ[i];float b0 = bQ[i];
					Q[k][i] = (c * a0)+(s * b0);
					Q[j][i] = (c * b0)-(s * a0);
				}
				
				#else

				for(i=0; (i+1)<N_ROW; i+=2){

					float Qki0 = Q[k][i]  ;float Qji0 = Q[j][i];
					float Qki1 = Q[k][i+1];float Qji1 = Q[j][i+1];
					
					Q[k][i]   = (c * Qki0)+(s * Qji0);
					Q[j][i] = (c * Qji0)-(s * Qki0);

					Q[k][i+1] = (c * Qki1)+(s * Qji1);
					Q[j][i+1] = (c * Qji1)-(s * Qki1);
				}if(i < N_ROW){
					float Qki0 = Q[k][i]  ;float Qji0 = Q[j][i];

					Q[k][i]   = (c * Qki0)+(s * Qji0);
					Q[j][i] = (c * Qji0)-(s * Qki0);					
				}

				#endif
			}
		}
    }
    return;
}

