#include "pmsis.h"
#include "util.h" 
#include "inputData/T/qr40X40.h"
#include "inputData/T/data40X40.inc"

#include <math.h>

//Q should be reversed
PI_L1 float Q[N_ROW][N_ROW]; 
PI_L1 float Q_temp[N_ROW][N_ROW];

//R should be reversed => 
PI_L1 float R[N_COL][N_ROW];
PI_L1 float R_temp[N_COL][N_ROW];




PI_L1 float v[N_ROW];
PI_L1 float b[N_COL];
PI_L1 float v_temp[N_ROW][N_ROW];

float res[N_ROW][N_COL];

int numBarr = 0;


definePrefArr



void cluster_main();
void qr_household(float Q[][N_ROW], float R[][N_ROW]);
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O, int core_id);

void pe_entry(void *arg){cluster_main();}

void cluster_entry(void *arg){pi_cl_team_fork((NUM_CORES), pe_entry, 0);}

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
	pmsis_exit(ret);}

int main(){
	return pmsis_kickoff((void *)test_kickoff);}


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
  
	if (pi_core_id()==0){ 

		initialMatrixEYE(&Q[0][0],N_ROW,N_ROW);
		initialMatrixMatrix(&R[0][0],N_ROW,N_COL,&input[0][0]);
		initialMatrixConst(&R_temp[0][0],N_ROW,N_COL,0);
	}
    
	pi_cl_team_barrier();
    
	pi_perf_reset();
	pi_perf_start();
    
    qr_household(Q, R);

	pi_perf_stop();
	
	printPerfCounters

	if(pi_core_id() == 0){
		printf("NUMBARR = %d\n",numBarr);
	}
	
	pi_cl_team_barrier();
	#ifdef DEBUG
 	if(pi_core_id()==0){
		
		printMatrix(&Q[0][0],N_ROW,N_ROW,"Q",1);
		printMatrixU(&R[0][0],N_ROW,N_COL,"R",1);

		for (int i = 0; i < N_ROW; i++) {
			for (int j = 0; j < N_COL; j++) {
				res[i][j] = 0;
				for (int k = 0; k < j+1; k++) {
					res[i][j] += Q[k][i] * R[j][k];
				}
			}
		}
		printMatrix(&res[0][0],N_ROW,N_COL,"Res",0);
	}

	pi_cl_team_barrier();
	#endif
}


inline float Sqrt(float x) {
        float res;
        asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
        return res;
}

void house(float *v,float *b, float *x, int n){
	float sigma = 0;
	float x0 = x[0];
	if(n == 1){
		sigma = 0;
	}else{
		for (int i = 1; i < n; i++){
			sigma = sigma + x[i]*x[i];
		}
	}
	v[0] = one;
	for (int i = 1; i < n; i++){
		v[i] = x[i];
	}
	if(sigma == 0){
		*b = 0;
	}else{
		float mu = Sqrt(x0*x0 + sigma);
		if(x0 <= 0){
			v[0] = x0 - mu;
		}else{
			v[0] = -sigma/(x0+mu);
		}
	}
	float v0 = v[0]; 
	float v0p2 = v0*v0;
	*b = two*v0p2/(sigma+v0p2);
	v[0] = 1;
	for (int i = 1; i < n; i++){
		v[i] = v[i]/v0;
	}
}

void house_opt(float *v,float *b, float *x, int n){
	float sigma = 0;
	float x0 = x[0];
	if(n == 1){
		sigma = 0;
	}else{
		int i;
		for (i = 1; i+1 < n; i+=2){
			float x0 = x[i];
			float x1 = x[i+1];
			float x0p2 = x0*x0;
			float x1p2 = x1*x1;
			sigma = sigma + x0p2+x1p2;
		}
		if(i < n){
			float x0 = x[i];
			float x0p2 = x0*x0;
			sigma = sigma + x0p2;
		}
	}
	v[0] = one;
	int i;
	for (i = 1; i+1 < n; i+=2){
		float x0 = x[i];
		float x1 = x[i+1];
		v[i] = x0;
		v[i+1] = x1;
	}
	if(i < n){
		float x0 = x[i];
		v[i] = x0;
	}

	if(sigma == 0){
		*b = 0;
	}else{
		float mu = Sqrt(x0*x0 + sigma);
		if(x0 <= 0){
			v[0] = x0 - mu;
		}else{
			v[0] = -sigma/(x0+mu);
		}
	}
	float v0 = v[0]; 
	float v0p2 = v0*v0;
	*b = two*v0p2/(sigma+v0p2);
	v[0] = one;
	for (int i = 1; i < n; i++){
		v[i] = v[i]/v0;
	}
}

void qr_household(float Q[][N_ROW], float R[][N_ROW]){
	int core_id = pi_core_id();
	int l;
	float bj;
	#if NUM_CORES > 1
		int blockSize_ROW = N_ROW/NUM_CORES;
		int start_ROW = core_id*blockSize_ROW;

		if(core_id==(NUM_CORES - 1)){
			blockSize_ROW = N_ROW - (NUM_CORES - 1)* blockSize_ROW;}
		int end_ROW = start_ROW+blockSize_ROW;

		int blockSize_COL = N_COL/NUM_CORES;
		int start_COL = core_id*blockSize_COL;

		if(core_id==(NUM_CORES - 1)){
			blockSize_COL = N_COL - (NUM_CORES - 1)* blockSize_COL;}
		int end_COL = start_ROW+blockSize_COL;
		int start = 0;

	#endif

	for (int j = 0; j < N_COL; j++){
		if(core_id == 0){
			house_opt(&v[0],&b[j], &R[j][j], N_ROW-j);
		}
		#if NUM_CORES > 1
			pi_cl_team_barrier();
			BarrierCounter
		#endif
		bj = b[j];
		#if NUM_CORES > 1
			for (int i = start_ROW; ((i < N_ROW-j) && (i < end_ROW)); i++){
				float vi = v[i];
				for (int k = 0; k < N_ROW-j; k++){
					if(i == k){
						v_temp[i][k] =	one - bj * vi*vi;
					}else{
						v_temp[i][k] = -bj * vi*v[k];
					}
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for (int i = 0; i < N_ROW-j; i++){
				float vi = v[i];
				for (int k = 0; k < N_ROW-j; k++){
					if(i == k){
						v_temp[i][k] =	one - bj * vi*vi;
					}else{
						v_temp[i][k] = -bj * vi*v[k];
					}
				}
			}
		#endif

		#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}
			for(int i = start; i < end_ROW; i++){
				for (int k = j; k < N_COL; k++){
					float temp = 0;

					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i-j][l-j]*R[k][l];
					}
					R_temp[k][i] = temp;
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j; i < N_ROW; i++){
				for (int k = j; k < N_COL; k++){
					float temp = 0;
					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i-j][l-j]*R[k][l];
					}
					R_temp[k][i] = temp;
				}
			}
		#endif

		#if NUM_CORES > 1
			if(start_COL >= j){
				start = start_COL;
			}else if(end_COL >= j){
				start = j;
			}else{
				start = end_COL+1;
			}

			for(int i = start; i < end_COL; i++){
				int k;
				for (k = j; k+1 < N_ROW; k+=2){
					float R_temp_ik0 = R_temp[i][k];
					float R_temp_ik1 = R_temp[i][k+1];
					R[i][k] = R_temp_ik0;
					R[i][k+1] = R_temp_ik1;
				}
				if(k < N_ROW){
					float R_temp_ik0 = R_temp[i][k];
					R[i][k] = R_temp_ik0;
				}
			}

		#else
			for(int i = j; i < N_COL; i++){
				for (int k = j; k < N_ROW; k++){
					R[i][k] = R_temp[i][k];
				}
			}
		#endif

		if(j < N_ROW){
		#if NUM_CORES > 1
			if(start_ROW >= j+1){
				start = start_ROW;
			}else if(end_ROW >= j+1){
				start = j+1;
			}else{
				start = end_ROW+1;
			}

			int i;
			for(i = start; i+1 < end_ROW; i+=2){
				float v0 = v[1+i-j-1]; 
				float v1 = v[1+i+1-j-1]; 
				R[j][i] = v0;
				R[j][i+1] = v1;
			}
			if(i < end_ROW){
				float v0 = v[1+i-j-1]; 
				R[j][i] = v0;
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j+1; i < N_ROW; i++){
				R[j][i] = v[1+i-j-1];
			}

		#endif
		}
	}


	for (int j = N_COL-1; j >= 0; j--){
		v[j] = one;
		#if NUM_CORES > 1
			if(start_ROW >= j+1){
				start = start_ROW;
			}else if(end_ROW >= j+1){
				start = j+1;
			}else{
				start = end_ROW+1;
			}

			for(int i = start; i < end_ROW; i++){
				v[i] = R[j][i];
			}
		#else
			for(int i = j+1; i < N_ROW; i++){
				v[i] = R[j][i];
			}
		#endif


		#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}

			float bj = b[j];
			for (int i = start; i < end_ROW; i++){
				float vi = v[i];
				for (int k = j; k < N_ROW; k++){
					if(i == k){
						v_temp[i][k] = one - bj * vi*vi;
					}else{
						v_temp[i][k] =   - bj * vi*v[k];
					}
				}
			}
		#else
			float bj = b[j];
			for (int i = j; i < N_ROW; i++){
				float vi = v[i];
				for (int k = j; k < N_ROW; k++){
					if(i == k){
						v_temp[i][k] = one - bj * vi*vi;
					}else{
						v_temp[i][k] =   - bj * vi*v[k];
					}
				}
			}
		#endif

		#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}
			for(int i = start; i < end_ROW; i++){
				for (int k = j; k < N_ROW; k++){
					float temp = 0;
					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i][l]*Q[k][l];
					}
					Q_temp[k][i] = temp;
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j; i < N_ROW; i++){
				for (int k = j; k < N_ROW; k++){
					float temp = 0;
					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i][l]*Q[k][l];
					}
					Q_temp[k][i] = temp;
				}
			}
		#endif
	
		#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}		
			for(int i = start; i < end_ROW; i++){
				for (int k = j; k < N_ROW; k++){
					Q[i][k] = Q_temp[i][k];
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j; i < N_ROW; i++){
				for (int k = j; k < N_ROW; k++){
					Q[i][k] = Q_temp[i][k];
				}
			}		
		#endif
	}
	return;
}

