#include "pmsis.h"
#include "stats.h" 
#include "inputData/qr40X40.h"
#include "inputData/data40X40.inc"
#include <math.h>

#define STATS 
#define STACK_SIZE 2048
#define MOUNT 1
#define UNMOUNT 0
#define CID 0



// PI_L1 float Q[N_ROW][N_ROW]; 
// PI_L1 float R[N_ROW][N_COL];
// PI_L1 float v[N_ROW];
// PI_L1 float sel_col[N_ROW];
// PI_L1 float Q_temp[N_ROW][N_ROW];
// PI_L1 float R_temp[N_ROW][N_COL];
// PI_L1 float H[N_ROW][N_ROW];
// PI_L1 float n;
// PI_L1 int i2;
// PI_L1 float zero = 0;
// PI_L1 float one = 1;
// PI_L1 float two = 2;
// PI_L1 float temp;
// PI_L1 static float buffer[NUM_CORES];
// PI_L1 static int idx[NUM_CORES];

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


int cyclesArr[NUM_CORES];
int instrArr[NUM_CORES];
int activeArr[NUM_CORES];
int ldEXTArr[NUM_CORES];
int TCDMArr[NUM_CORES];
int ldstallArr[NUM_CORES];
int imissArr[NUM_CORES];


void cluster_main();
//void qr_household(float Q[][N_ROW], float R[][N_COL]);
void qr_household(float Q[][N_ROW], float R[][N_COL]);
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O);

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

	unsigned long cycles = 0;
	unsigned long instr = 0;
	unsigned long active = 0;
	unsigned long ldext = 0;
	unsigned long tcdmcont = 0;
	unsigned long ldstall = 0;
	unsigned long imiss = 0;
	unsigned long apu_cont = 0;
	unsigned long apu_dep = 0;
	unsigned long apu_type = 0;
	unsigned long apu_wb = 0;


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
		for(int i=0; i<N_ROW; i++){
			for(int j=0; j<N_ROW; j++){
				if(i == j){
					Q[i][j]=1;
				}else{
					Q[i][j]=0;
				}
			}
		}


		for(int i=0; i<N_ROW; i++){
			for(int j=0; j<N_COL; j++){
				R[i][j]=input[i][j];            
			}
		}
	}
    
	pi_cl_team_barrier();
    
	pi_perf_reset();
	pi_perf_start();
    
    qr_household(Q, R);

	pi_perf_stop();
	
	cycles   = pi_perf_read (PI_PERF_CYCLES);
	instr    = pi_perf_read (PI_PERF_INSTR);
	active   = pi_perf_read (PI_PERF_ACTIVE_CYCLES);
	ldext    = pi_perf_read (PI_PERF_LD_EXT);
	tcdmcont = pi_perf_read (PI_PERF_TCDM_CONT);
	ldstall  = pi_perf_read (PI_PERF_LD_STALL);
	imiss    = pi_perf_read (PI_PERF_IMISS);
	
	//apu_cont = pi_perf_read (0x12);
	//apu_dep  = pi_perf_read (0x13);
	//apu_type = pi_perf_read __SPRREAD (0x791);
	//apu_wb   = pi_perf_read (0x14);
   
    
	cyclesArr[pi_core_id()] = cycles/REPEAT;
	instrArr[pi_core_id()] = instr/REPEAT;
	activeArr[pi_core_id()] = active/REPEAT;
	ldEXTArr[pi_core_id()] = ldext/REPEAT;
	TCDMArr[pi_core_id()] = tcdmcont/REPEAT;
	ldstallArr[pi_core_id()] = ldstall/REPEAT;
	imissArr[pi_core_id()] = imiss/REPEAT;

	if(pi_core_id() == 0){
		long int AVGCycles = 0;
		long int AVGInstr = 0;
		long int AVGActive = 0;
		long int AVGldEXT = 0;
		long int AVGTCDM = 0;
		long int AVGLdstall = 0;
		long int AVGImiss = 0;
		for (int i = 0; i < NUM_CORES; i++){
			AVGCycles += cyclesArr[i];
			AVGInstr += instrArr[i];
			AVGActive += activeArr[i];
			AVGldEXT += ldEXTArr[i];
			AVGTCDM += TCDMArr[i];
			AVGLdstall += ldstallArr[i];
			AVGImiss += imissArr[i];
		}
		printf("AVGCycles = %lu\n",AVGCycles/NUM_CORES);
		printf("AVGInstr = %lu\n",AVGInstr/NUM_CORES);
		printf("AVGActive = %lu\n",AVGActive/NUM_CORES);
		printf("AVGldEXT = %lu\n",AVGldEXT/NUM_CORES);
		printf("AVGTCDM = %lu\n",AVGTCDM/NUM_CORES);
		printf("AVGLdstall = %lu\n",AVGLdstall/NUM_CORES);
		printf("AVGImiss = %lu\n",AVGImiss/NUM_CORES);


		printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",AVGCycles/NUM_CORES,AVGInstr/NUM_CORES,AVGActive/NUM_CORES,AVGldEXT/NUM_CORES
		,AVGTCDM/NUM_CORES,AVGLdstall/NUM_CORES,AVGImiss/NUM_CORES);
	}

	pi_cl_team_barrier();

	int id = pi_core_id();
	printf("[%d] cycles = %lu\n", id, cycles/REPEAT);
	printf("[%d] instr = %lu\n", id, instr/REPEAT);
	printf("[%d] active cycles = %lu\n", id, active/REPEAT);
	printf("[%d] ext load = %lu\n", id, ldext/REPEAT);
	printf("[%d] TCDM cont = %lu\n", id, tcdmcont/REPEAT);
	printf("[%d] ld stall = %lu\n", id, ldstall/REPEAT);
	printf("[%d] imiss = %lu\n", id, imiss/REPEAT);
	//printf("[%d] apu_cont = %lu\n", id, apu_cont/REPEAT);
	//printf("[%d] apu_dep = %lu\n", id, apu_dep/REPEAT);
	//printf("[%d] apu_type = %lu\n", id, apu_type/REPEAT);
	//printf("[%d] apu_wb = %lu\n", id, apu_wb/REPEAT);
    
	pi_cl_team_barrier();
    
	// if(pi_core_id()==0){
	// 	printf("\n\nQ = \n");
	// 	for(int i=0; i<N_ROW; i++){
	// 		for(int j=0; j<N_ROW; j++)
	// 			printf("%.20f ", Q[i][j]);
	// 		printf("\n");
	// 	}

	// 	printf("\n\nR = \n");
	// 	for(int i=0; i<N_ROW; i++){
	// 		for(int j=0; j<N_COL; j++)
	// 			printf("%.20f ", R[i][j]);
	// 		printf("\n");
	// 	}
	// 	printf("\n\nRes = \n");
	// 	float res[N_ROW][N_COL];
	// 	for (int i = 0; i < N_ROW; i++) {
	// 		for (int j = 0; j < N_COL; j++) {
	// 			res[i][j] = 0;
	// 			for (int k = 0; k < N_ROW; k++) {
	// 				res[i][j] += Q[i][k] * R[k][j];
	// 			}
	// 			printf("%.20f ", res[i][j]);
	// 		}
	// 		printf("\n");
	// 	}
	// }
	pi_cl_team_barrier();
}


inline float Sqrt(float x) {
        float res;
        asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
        return res;
}


//float norm(float *v, int row){
float norm(float *v, int row){
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
	}
	pi_cl_team_barrier();
	#else


	// for(j=0; j<row;j++){
	// 	n = n + v[i2]*v[i2];
	// 	i2++;
	// }
	for(j=0; j<row;j++){
		n = n + v[j]*v[j];
	}
	#endif

	// return Sqrt(n);
	return sqrtf(n);

}



#pragma GCC push_options
#pragma GCC optimize ("-fivopts")
void __attribute__ ((noinline)) matMul(float * __restrict__ pSrcA, float * __restrict__ pSrcB, float * __restrict__ pDstC, int M, int N, int O) {


    int i = M; // loop counter for M
    int j = N; // loop counter for N
    int k = O; // loop counter for O

    int core_id = pi_core_id();
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
	for(int i=0;i<N_COL;i++){
		


    	#if NUM_CORES > 1
        int blockSize_ROW = N_ROW/NUM_CORES;
        int start_ROW = pi_core_id()*blockSize_ROW;

        if(pi_core_id()==(NUM_CORES - 1)){
            blockSize_ROW = N_ROW - (NUM_CORES - 1)* blockSize_ROW;}
		#endif



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


		// pi_cl_team_barrier();
		// if(pi_core_id() == 0){
		// 	printf("H\n");
		// 	for(int j = 0; j < N_ROW; j++){
		// 		for(int k = 0; k < N_ROW; k++){
		// 			// printf("%.7f\t",sel_col[j]);
		// 			printf("%e\t",H[j][k]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }
		// pi_cl_team_barrier();



		#endif


    	#if NUM_CORES > 1

		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			// if(j<i){
			// 	sel_col[j] = zero;
			// }else{
			// 	sel_col[j] = R[j][i];
			// }
			float Rji0 = R[j][i];			
			if(j<i){
				sel_col[j] = zero;
			}else if((j+1) <(start_ROW + blockSize_ROW)){
				float Rj0, Rj1;
				Rj0 = R[j][i];
				Rj1 = R[j+1][i];
				sel_col[j] = Rj0;
				sel_col[j+1] = Rj1;
				j++;
			}else{
				sel_col[j] = Rji0;
			}
		}
		pi_cl_team_barrier();
		#else

		for(int j=0;j<N_ROW;j++){
			// if(j<i){
			// 	sel_col[j] = zero;
			// }else{
			// 	sel_col[j] = R[j][i];
			// }
			float Rji0 = R[j][i];
			if(j<i){
				sel_col[j] = zero;
			}else if((j+1) < N_ROW){
				float Rj0, Rj1;
				Rj0 = R[j][i];
				Rj1 = R[j+1][i];
				sel_col[j] = Rj0;
				sel_col[j+1] = Rj1;
				j++;
			}else{
				sel_col[j] = Rji0;
			}
		}

		#endif

		// pi_cl_team_barrier();
		// if(pi_core_id() == 0){
		// 	printf("sel_col\n");
		// 	for(int j = 0; j < N_ROW; j++){
		// 		printf("%.10f\t",sel_col[j]);
		// 		// printf("%e\t",sel_col[j]);
		// 	}
		// 	printf("\n");
		// }
		pi_cl_team_barrier();

		temp = norm(&sel_col[0],N_ROW);

		pi_cl_team_barrier();
		// if(pi_core_id() == 0){
		// 	printf("norm = %.10f\n",temp);
		// 	// printf("norm = %e\n",temp);
		// }
		// pi_cl_team_barrier();


    	#if NUM_CORES > 1

		// for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
		// 	if(j == i){
		// 		if(sel_col[j] >= 0){
		// 			v[j] = sel_col[j] + temp;
		// 		}else{
		// 			v[j] = sel_col[j] - temp;
		// 		}
		// 	}else{
		// 		v[j] = sel_col[j];
		// 	}
		// }

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
			// }else if((j+1) == i){
			// 	float sel_col0;
			// 	sel_col0 = sel_col[j];
			// 	v[j] = sel_col0;
			}else if(j+1 < (start_ROW + blockSize_ROW)){
				float sel_col0, sel_col1;
				sel_col0 = sel_col[j];
				sel_col1 = sel_col[j+1];
				v[j] = sel_col0;
				v[j+1] = sel_col1;
				j++;
			}else{
				float sel_col0;
				sel_col0 = sel_col[j];
				v[j] = sel_col0;				
			}
		}

		pi_cl_team_barrier();

		#else

		// for(int j=0;j<N_ROW;j++){
		// 	if(j == i){
		// 		if(sel_col[j] >= 0){
		// 			v[j] = sel_col[j] + temp;
		// 		}else{
		// 			v[j] = sel_col[j] - temp;
		// 		}
		// 	}else{
		// 		v[j] = sel_col[j];
		// 	}
		// }


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
			}else if(j+1 < N_ROW){
				float sel_col0, sel_col1;
				sel_col0 = sel_col[j];
				sel_col1 = sel_col[j+1];
				v[j] = sel_col0;
				v[j+1] = sel_col1;
				j++;
			}
			else{
				float sel_col0;
				sel_col0 = sel_col[j];
				v[j] = sel_col0;
			}
		}


		#endif

		// pi_cl_team_barrier();
		// if(pi_core_id() == 0){
		// 	printf("v\n");
		// 	for(int j = 0; j < N_ROW; j++){
		// 		printf("%.20f\t",v[j]);
		// 		// printf("%e\t",v[j]);
		// 	}
		// 	printf("\n");
		// }
		// pi_cl_team_barrier();



		pi_cl_team_barrier();

		temp = zero;

		pi_cl_team_barrier();

    	#if NUM_CORES > 1
        buffer[pi_core_id()]=0;
        idx[pi_core_id()]=0;


		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){

			buffer[pi_core_id()] = buffer[pi_core_id()] + v[idx[pi_core_id()]+start_ROW]*v[idx[pi_core_id()]+start_ROW];
			idx[pi_core_id()] = idx[pi_core_id()] + 1;
		}

		pi_cl_team_barrier();
        if(pi_core_id()==0){
            for(int j=0; j<NUM_CORES; j++){
                temp += buffer[j];
            }
		}
		pi_cl_team_barrier();

		#else

		// for (int j = 0; j < N_ROW; j++){
		// 	temp = temp + v[j]*v[j];
		// }

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

		// pi_cl_team_barrier();
		// if(pi_core_id() == 0){
		// 	printf("temp = %.10f\n",temp);
		// 	// printf("temp = %e\n",temp);
		// }
		// pi_cl_team_barrier();

    	#if NUM_CORES > 1

		// for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
		// 	if(j >= i){
		// 		for(int k = i; k < N_ROW; k++){
		// 			H[j][k] = H[j][k] - two*v[k]*v[j]/temp;
		// 		}
		// 	}
		// }

		for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
			int k;
			float vj = v[j];
			float vjx2 = vj*two;
			if(j >= i){
				for(k = i; k+1 < N_ROW; k+=2){
					// printf("id = %d\tstart_ROW = %d\tblockSize_ROW = %d\ti = %d\t j = %d\t k = %d\n",pi_core_id(),start_ROW,blockSize_ROW,i,j,k);	
					float Hjk0 ,Hjk1;
					float vk0 ,vk1;
					vk0 = v[k];
					Hjk0 = H[j][k];
					vk1 = v[k+1];
					Hjk1 = H[j][k+1];
					// H[j][k] = Hjk0 - vk0*vjx2;
					// H[j][k+1] = Hjk1 - vk1*vjx2;
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

		// for (int j = i; j < N_ROW; j++){
		// 	for(int k = i; k < N_ROW; k++){
		// 		H[j][k] = H[j][k] - (two*v[k]*v[j])/temp;
		// 	}
		// }

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
				// if(k == 0){
				// 	printf("K is 0\n");
				// }
				
			}
			if(k < N_ROW){
				float Hjk0;
				float vk0;
				vk0 = v[k];
				Hjk0 = H[j][k];
				H[j][k] = Hjk0 - vk0*vjx2/temp;
				// if(k == 0){
				// 	printf("K is 0\n");
				// }
			}
		}



		#endif

		// pi_cl_team_barrier();
		// if(pi_core_id() == 0){
		// 	printf("H\n");
		// 	for(int j = 0; j < N_ROW; j++){
		// 		for(int k = 0; k < N_ROW; k++){
		// 			printf("%.10f\t",H[j][k]);
		// 			// printf("%e\t",H[j][k]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }
		pi_cl_team_barrier();

		matMul(&H[0][0],&R[0][0],&R_temp[0][0],N_ROW,N_ROW,N_COL);
		matMul(&Q[0][0],&H[0][0],&Q_temp[0][0],N_ROW,N_ROW,N_ROW);

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

		// for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
		// 	for(int k = 0; k < N_COL; k++){
		// 		R[j][k] = R_temp[j][k];
		// 	}
		// }

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

		// for(int j = 0;j < N_ROW; j++){//??????
		// 	for(int k = 0; k < N_COL; k++){
		// 		R[j][k] = R_temp[j][k];
		// 	}
		// }

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

		// for(int j = start_ROW; j < (start_ROW + blockSize_ROW); j++){
		// 	for(k = 0; k < N_ROW; k++){
		// 		Q[j][k] = Q_temp[j][k];
		// 	}
		// }

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

		// for(int j = 0;j < N_ROW; j++){//??????
		// 	for(k = 0; k < N_ROW; k++){
		// 		Q[j][k] = Q_temp[j][k];
		// 	}
		// }


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

