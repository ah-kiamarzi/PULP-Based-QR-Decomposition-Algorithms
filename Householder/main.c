#include "pmsis.h"
#include "util.h" 
#include "inputData/qr40X40.h"
#include "inputData/data40X40.inc"
#include <math.h>


PI_L1 float Q[N_ROW][N_ROW]; 
PI_L1 float R[N_ROW][N_COL];
PI_L1 float v[N_ROW];
PI_L1 float sel_col[N_ROW];
PI_L1 float Q_temp[N_ROW][N_ROW];
PI_L1 float R_temp[N_ROW][N_COL];
PI_L1 float H[N_ROW][N_ROW];
PI_L1 float n;
PI_L1 int i2;
PI_L1 float zero = 0;
PI_L1 float one = 1;
PI_L1 float two = 2;
PI_L1 float temp;
PI_L1 static float buffer[NUM_CORES];
PI_L1 static int idx[NUM_CORES];

float res[N_ROW][N_COL];


definePrefArr



void cluster_main();
void qr_household(float Q[][N_ROW], float R[][N_COL]);
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
		// for(int i=0; i<N_ROW; i++){
		// 	for(int j=0; j<N_ROW; j++){
		// 		if(i == j){
		// 			Q[i][j]=1;
		// 		}else{
		// 			Q[i][j]=0;
		// 		}
		// 	}
		// }


		// for(int i=0; i<N_ROW; i++){
		// 	for(int j=0; j<N_COL; j++){
		// 		R[i][j]=input[i][j];            
		// 	}
		// }
		initialMatrixEYE(&Q[0][0],N_ROW,N_ROW);
		initialMatrixMatrix(&R[0][0],N_ROW,N_COL,&input[0][0]);
	}
    
	pi_cl_team_barrier();
    
	pi_perf_reset();
	pi_perf_start();
    
    qr_household(Q, R);

	pi_perf_stop();
	
	printPerfCounters

	pi_cl_team_barrier();
	#ifdef DEBUG
 	if(pi_core_id()==0){
		
		printMatrix(&Q[0][0],N_ROW,N_ROW,"Q",0);
		printMatrix(&R[0][0],N_ROW,N_COL,"R",0);

		for (int i = 0; i < N_ROW; i++) {
			for (int j = 0; j < N_COL; j++) {
				res[i][j] = 0;
				for (int k = 0; k < N_ROW; k++) {
					res[i][j] += Q[i][k] * R[k][j];
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


//float norm(float *v, int row){
void norm(float *v, int row,float *temp){
	i2=0;
	n = 0.0f;
	int j;


	#if NUM_CORES > 1
	int blockSize_row = row/NUM_CORES;
	int start_row = pi_core_id()*blockSize_row;

	if(pi_core_id()==(NUM_CORES - 1)){
		blockSize_row = row - (NUM_CORES - 1)* blockSize_row;}

	buffer[pi_core_id()]=0;
	idx[pi_core_id()]=0;
	
	for(j = start_row; (j<row) && (j<start_row + blockSize_row); j++){
		buffer[pi_core_id()] = buffer[pi_core_id()] + v[idx[pi_core_id()]+start_row]*v[idx[pi_core_id()]+start_row];
		idx[pi_core_id()]+=1;
	}
		
	pi_cl_team_barrier();
			
	if(pi_core_id()==0){
		for(j=0; j<NUM_CORES; j++){
			n += buffer[j];
		}
		*temp = Sqrt(n);
	}
	pi_cl_team_barrier();
	#else

	for(j=0; j<row;j++){
		n = n + v[j]*v[j];
	}
		*temp = Sqrt(n);
	#endif

	// return sqrtf(n);

}



#pragma GCC push_options
#pragma GCC optimize ("-fivopts")
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O,int core_id) {


    int i = M; // loop counter for M
    int j = N; // loop counter for N
    int k = O; // loop counter for O

    if(M<=0 || N<=0 || O<=0) return;

    for (k = core_id; k < O/2; k += NUM_CORES) {
        for (i = 0; i < M/2; i++) {

            float sum00 = 0;
            float sum01 = 0;
            float sum10 = 0;
            float sum11 = 0;

#ifdef UNROLL_INNER_LOOP
            for (j = 0; j < N/2; j++) {
                float AVal0 = pSrcA[i * 2 * N + (j*2)];
                float AVal1 = pSrcA[i * 2 * N + N + (j*2)];
                float BVal0 = pSrcB[(j*2) * O + (k * 2)];
                float BVal1 = pSrcB[(j*2) * O + (k * 2 + 1)];

                float AVal2 = pSrcA[i * 2 * N + (j*2+1)];
                float AVal3 = pSrcA[i * 2 * N + N + (j*2+1)];
                float BVal2 = pSrcB[(j*2+1) * O + (k * 2)];
                float BVal3 = pSrcB[(j*2+1) * O + (k * 2 + 1)];

                sum00 = sum00 + AVal0 * BVal0 + AVal2 * BVal2;
                sum01 = sum01 + AVal0 * BVal1 + AVal2 * BVal3;
                sum10 = sum10 + AVal1 * BVal0 + AVal3 * BVal2;
                sum11 = sum11 + AVal1 * BVal1 + AVal3 * BVal3;	        	
            }
	    
#else	   
            for (j = 0; j < N; j++) {
                float AVal0 = pSrcA[i * 2 * N + (j)];
                float AVal1 = pSrcA[i * 2 * N + N + (j)];

                float BVal0 = pSrcB[j * O + (k * 2)];
                float BVal1 = pSrcB[j * O + (k * 2 + 1)];

                sum00 = sum00 + (float) AVal0 * (float) BVal0;
                sum01 = sum01 + (float) AVal0 * (float) BVal1;
                sum10 = sum10 + (float) AVal1 * (float) BVal0;
                sum11 = sum11 + (float) AVal1 * (float) BVal1;
            }
#endif	    

            pDstC[(i * 2) * O + k * 2] = sum00;
            pDstC[(i * 2) * O + k * 2 + 1] = sum01;
	    	pDstC[(i * 2 + 1) * O + k * 2] = sum10;
            pDstC[(i * 2 + 1) * O + k * 2 + 1] = sum11;
        } // i 
    } // k
    // clean up code
    i = i * 2;
#ifdef UNROLL_INNER_LOOP    
    j = j * 2;
#endif
    k = k * 2;
    // check if every index is nicely finished
	if (i == M && j == N && k >= O) {

    } else {
		if (core_id >= O/2){
			i = 0;
		}
        uint32_t iEnd = i;
        uint32_t jEnd = j;
        uint32_t kEnd = k >= O ? O : k;

        // clean up for j
        if (jEnd != N) {
            for (i = 0; i < iEnd; i++) {
                for (k = 0; k < kEnd; k += NUM_CORES) {
                    float sum = 0;
                    for (j = jEnd; j < N; j++) {
                        sum += sum + pSrcA[i * N + j] * pSrcB[j * O + k];
                    }
                    pDstC[i * O + k] += sum;
                }
            }
        }

        // clean up for i
        if (iEnd != M) {
            for (k = core_id; k < kEnd; k += NUM_CORES) {
                for (i = iEnd; i < M; i++) {
                    float sum = 0;
                    for (j = 0; j < N; j++) {
                        sum = sum + pSrcA[i * N + j] * pSrcB[j * O + k];
                    }
                    pDstC[i * O + k] = sum;
                }
            }
        }

        // clean up for k
        for (k = kEnd; k < O; k += NUM_CORES) {	
            for (i = 0; i < M; i++) {
                float sum = 0;
                for (j = 0; j < N; j++) {
                    sum = sum + pSrcA[i * N + j] * pSrcB[j * O + k];
                }
                pDstC[i * O + k] = sum;
            }
        }
    }

    pi_cl_team_barrier(0);

}

#pragma GCC pop_options

void qr_household(float Q[][N_ROW], float R[][N_COL]){
	int core_id = pi_core_id();
	#if NUM_CORES > 1
	int blockSize_ROW = N_ROW/NUM_CORES;
	int start_ROW = core_id*blockSize_ROW;

	if(core_id==(NUM_CORES - 1)){
		blockSize_ROW = N_ROW - (NUM_CORES - 1)* blockSize_ROW;}
	#endif

	for(int i=0;i<N_COL;i++){
    	#if NUM_CORES > 1

		for(int k = start_ROW; k < (start_ROW + blockSize_ROW); k++){
			for(int j=0; j<N_ROW; j++){
				if(k == j){
					H[k][j]=one;
				}else{
					H[k][j]=zero;
				}
			}
		}

		#else

		for(int k=0; k<N_ROW; k++){
			for(int j=0; j<N_ROW; j++){
				if(k == j){
					H[k][j]=one;
				}else{
					H[k][j]=zero;
				}
			}
		}

		#endif


    	#if NUM_CORES > 1

		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			if(j<i){
				sel_col[j] = zero;
			}else{
				sel_col[j] = R[j][i];
			}
			// float Rji0 = R[j][i];			
			// if(j<i){
			// 	sel_col[j] = zero;
			// }else if((j+1) <(start_ROW + blockSize_ROW)){
			// 	float Rj0, Rj1;
			// 	Rj0 = R[j][i];
			// 	Rj1 = R[j+1][i];
			// 	sel_col[j] = Rj0;
			// 	sel_col[j+1] = Rj1;
			// 	j++;
			// }else{
			// 	sel_col[j] = Rji0;
			// }
		}
		pi_cl_team_barrier();
		#else

		for(int j=0;j<N_ROW;j++){
			if(j<i){
				sel_col[j] = zero;
			}else{
				sel_col[j] = R[j][i];
			}
			// float Rji0 = R[j][i];
			// if(j<i){
			// 	sel_col[j] = zero;
			// }else if((j+1) < N_ROW){
			// 	float Rj0, Rj1;
			// 	Rj0 = R[j][i];
			// 	Rj1 = R[j+1][i];
			// 	sel_col[j] = Rj0;
			// 	sel_col[j+1] = Rj1;
			// 	j++;
			// }else{
			// 	sel_col[j] = Rji0;
			// }
		}

		#endif

		pi_cl_team_barrier();
		
		norm(&sel_col[0],N_ROW, &temp);
		
		pi_cl_team_barrier();

    	#if NUM_CORES > 1

		float sel_coli;
		if(i < N_ROW){
			sel_coli = sel_col[i];
		}
		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			if(j == i){
				if(sel_coli >= 0){
					v[j] = sel_coli + temp;
				}else{
					v[j] = sel_coli - temp;
				}
			}else{
				v[j] = sel_col[j];
			}
		}

		pi_cl_team_barrier();

		#else

		float sel_coli;
		if(i < N_ROW){
			sel_coli = sel_col[i];
		}
		for(int j=0;j<N_ROW;j++){
			if(j == i){
				if(sel_coli >= 0){
					v[j] = sel_coli + temp;
				}else{
					v[j] = sel_coli - temp;
				}
			}else{
				v[j] = sel_col[j];
			}
		}


		#endif

		pi_cl_team_barrier();
		if(core_id == 0){
			temp = zero;
		}
		pi_cl_team_barrier();

    	#if NUM_CORES > 1
        buffer[core_id]=0;
        idx[core_id]=0;


		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){

			buffer[core_id] = buffer[core_id] + v[idx[core_id]+start_ROW]*v[idx[core_id]+start_ROW];
			idx[core_id] = idx[core_id] + 1;
		}

		pi_cl_team_barrier();
        if(core_id==0){
            for(int j=0; j<NUM_CORES; j++){
                temp += buffer[j];
            }
		}
		pi_cl_team_barrier();

		#else

		for(int j=0; (j+1)<N_ROW; j+=2){
			float vj0 = v[j];float vj1 = v[j+1];
			float sumTemp = vj0*vj0;
			sumTemp = sumTemp + vj1*vj1;
			temp = temp + sumTemp;
		}
			
		if(N_ROW%2){
			float t0 = v[N_ROW-1];
			float sumTemp = t0*t0;
			temp = temp + sumTemp;
		}


		#endif

    	#if NUM_CORES > 1

		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			int k;
			float vj = v[j];
			float vjx2 = vj*two;
			if(j >= i){
				for(k = i; k+1 < N_ROW; k+=2){
					float Hjk0 ,Hjk1;
					float vk0 ,vk1;
					vk0 = v[k];
					Hjk0 = H[j][k];
					vk1 = v[k+1];
					Hjk1 = H[j][k+1];
					H[j][k] = Hjk0 - vk0*vjx2/temp;
					H[j][k+1] = Hjk1 - vk1*vjx2/temp;
				}
				if(k < N_ROW){
					float Hjk0;
					float vk0;
					vk0 = v[k];
					Hjk0 = H[j][k];
					H[j][k] = Hjk0 - vk0*vjx2/temp;					
				}
			}
		}
		pi_cl_team_barrier();
		#else

		for (int j = i; j < N_ROW; j++){
			int k;
			float vj = v[j];
			float vjx2 = vj*two;
			for(k = i; k+1 < N_ROW; k+=2){
				float Hjk0 ,Hjk1;
				float vk0 ,vk1;
				vk0 = v[k];
				Hjk0 = H[j][k];
				vk1 = v[k+1];
				Hjk1 = H[j][k+1];
				H[j][k] = Hjk0 - vk0*vjx2/temp;
				H[j][k+1] = Hjk1 - vk1*vjx2/temp;
				
			}
			if(k < N_ROW){
				float Hjk0;
				float vk0;
				vk0 = v[k];
				Hjk0 = H[j][k];
				H[j][k] = Hjk0 - vk0*vjx2/temp;
			}
		}



		#endif

		pi_cl_team_barrier();

		matMul(&H[0][0],&R[0][0],&R_temp[0][0],N_ROW,N_ROW,N_COL,core_id);
		matMul(&Q[0][0],&H[0][0],&Q_temp[0][0],N_ROW,N_ROW,N_ROW,core_id);

		pi_cl_team_barrier();

		// if (pi_core_id() == 0)
		// {

		// 	for (int i1 = 0; i1 < N_ROW; i1++) {
		// 		for (int j1 = 0; j1 < N_COL; j1++) {
		// 			R_temp[i1][j1] = 0;
		// 			for (int k1 = 0; k1 < N_ROW; k1++) {
		// 				R_temp[i1][j1] += H[i1][k1] * R[k1][j1];
		// 			}
		// 			//printf("%f ", R_temp[i][j]);
		// 		}
		// 		//printf("\n");
		// 	}


		// 	for (int i1 = 0; i1 < N_ROW; i1++) {
		// 		for (int j1 = 0; j1 < N_ROW; j1++) {
		// 			Q_temp[i1][j1] = 0;
		// 			for (int k1 = 0; k1 < N_ROW; k1++) {
		// 				Q_temp[i1][j1] += Q[i1][k1] * H[k1][j1];
		// 			}
		// 			//printf("%f ", R_temp[i][j]);
		// 		}
		// 		//printf("\n");
		// 	}			
		// }

		// pi_cl_team_barrier();

		#if NUM_CORES > 1

		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			int k;
			for(k = 0; k+1 < N_COL; k+=2){
				float R_temp0, R_temp1;
				R_temp0 = R_temp[j][k];
				R_temp1 = R_temp[j][k+1];
				R[j][k] = R_temp0;
				R[j][k+1] = R_temp1;
			}
			if(k < N_COL){
				float R_temp0 = R_temp[j][k];
				R[j][k] = R_temp0;
			}
		}

		#else

		for(int j = 0;j < N_ROW; j++){//??????
			int k;
			for(k = 0; k+1 < N_COL; k+=2){
				float R_temp0, R_temp1;
				R_temp0 = R_temp[j][k];
				R_temp1 = R_temp[j][k+1];
				R[j][k] = R_temp0;
				R[j][k+1] = R_temp1;
			}
			if(k < N_COL){
				float R_temp0 = R_temp[j][k];
				R[j][k] = R_temp0;
			}			
		}

		#endif
		
		#if NUM_CORES > 1

		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			int k;
			for(k = 0; (k+1) < N_ROW; k+=2){
				float Q_temp0, Q_temp1;
				Q_temp0 = Q_temp[j][k];
				Q_temp1 = Q_temp[j][k+1];
				Q[j][k] = Q_temp0;
				Q[j][k+1] = Q_temp1;
			}
			if(k < N_ROW){
				float Q_temp0 = Q_temp[j][k];
				Q[j][k] = Q_temp0;
			}
		}
		
		
		pi_cl_team_barrier();

		#else

		for(int j = 0;j < N_ROW; j++){//??????
			int k;
			for(k = 0; (k+1) < N_ROW; k+=2){
				float Q_temp0, Q_temp1;
				Q_temp0 = Q_temp[j][k];
				Q_temp1 = Q_temp[j][k+1];
				Q[j][k] = Q_temp0;
				Q[j][k+1] = Q_temp1;
			}
			if(k < N_ROW){
				float Q_temp0 = Q_temp[j][k];
				Q[j][k] = Q_temp0;
			}
		}

		#endif

	}
	return;
}

