#include "pmsis.h"
#include "stats.h" 
#include "qr.h"
#include "data.inc"
#include <math.h>

#define STATS 
#define STACK_SIZE 2048
#define MOUNT 1
#define UNMOUNT 0
#define CID 0



PI_L1 float Q[N_CHANNELS][N_CHANNELS]; 
PI_L1 float R[N_CHANNELS][EV_WINDOWS_SIZE];
PI_L1 float v[N_CHANNELS];
PI_L1 float sel_col[N_CHANNELS];
PI_L1 float Q_temp[N_CHANNELS][N_CHANNELS];
PI_L1 float R_temp[N_CHANNELS][EV_WINDOWS_SIZE];
PI_L1 float H[N_CHANNELS][N_CHANNELS];
PI_L1 float input_temp[N_CHANNELS][EV_WINDOWS_SIZE];
PI_L1 float n;
PI_L1 int i2;
PI_L1 float zero = 0;
PI_L1 float one = 1;
PI_L1 float two = 2;
PI_L1 float temp;
PI_L1 static float buffer[NUM_CORES];
PI_L1 static int idx[NUM_CORES];



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

void cluster_main();
void qr_household(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE]);
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
		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<N_CHANNELS; j++){
				if(i == j){
					Q[i][j]=1;
				}else{
					Q[i][j]=0;
				}
			}
		}


		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<EV_WINDOWS_SIZE; j++){
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
	
	apu_cont = pi_perf_read (0x12);
	apu_dep  = pi_perf_read (0x13);
	//apu_type = pi_perf_read __SPRREAD (0x791);
	apu_wb   = pi_perf_read (0x14);
   
    
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
    
	/*if(pi_core_id()==0){
	printf("\n\nQ = \n");
		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<N_CHANNELS; j++)
				printf("%f ", Q[i][j]);
			printf("\n");
		}

	printf("\n\nR = \n");
		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<EV_WINDOWS_SIZE; j++)
				printf("%f ", R[i][j]);
			printf("\n");
		}
	}*/

	pi_cl_team_barrier();
}


inline float Sqrt(float x) {
        float res;
        asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
        return res;
}


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
	
	//for(j = start_row; (j<row) && (j<start_row + blockSize_row); j++){
	for(j = start_row; (j<row) && (j<start_row + blockSize_row); j++){
		//printf("id = %d\tblockSize_row = %d\tstart_row = %d\n",pi_core_id(),blockSize_row,start_row);

		buffer[pi_core_id()] = buffer[pi_core_id()] + v[idx[pi_core_id()]+start_row]*v[idx[pi_core_id()]+start_row];
		//idx[pi_core_id()]+=row;
inline float Sqrt(float x) {
        float res;
        asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
        return res;
}
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


	for(j=0; j<row;j++){
		n = n + v[i2]*v[i2];
		i2++;
	}
	
	#endif

	return sqrt(n);
	
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
        uint32_t iEnd = i;
        uint32_t jEnd = j;
        uint32_t kEnd = k >= O ? O : k;

        // clean up for j
        if (jEnd != N) {
            for (i = 0; i < iEnd; i++) {
                for (k = 0; k < kEnd; k += NUM_CORES) {
                    int32_t sum = 0;
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
                    int32_t sum = 0;
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
                int32_t sum = 0;
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



void qr_household(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE]){
	int j;

	for(int i=0;i<EV_WINDOWS_SIZE;i++){
		


    	#if NUM_CORES > 1
        int blockSize_NC = N_CHANNELS/NUM_CORES;
        int start_NC = pi_core_id()*blockSize_NC;

        if(pi_core_id()==(NUM_CORES - 1)){
            blockSize_NC = N_CHANNELS - (NUM_CORES - 1)* blockSize_NC;}
		#endif



    	#if NUM_CORES > 1

		for(int k = start_NC; k < (start_NC + blockSize_NC); k++){
			for(int j=0; j<N_CHANNELS; j++){
				if(k == j){
					H[k][j]=one;
				}else{
					H[k][j]=zero;
				}
			}
		}

		#else

		for(int k=0; k<N_CHANNELS; k++){
			for(int j=0; j<N_CHANNELS; j++){
				if(k == j){
					H[k][j]=one;
				}else{
					H[k][j]=zero;
				}
			}
		}

		#endif


    	#if NUM_CORES > 1

		for(int j = start_NC; j < (start_NC + blockSize_NC); j++){
			if(j<i){
				sel_col[j] = zero;
			}else{
				sel_col[j] = R[j][i];
			}
		}

		#else

		for(j=0;j<N_CHANNELS;j++){
			if(j<i){
				sel_col[j] = zero;
			}else{
				sel_col[j] = R[j][i];
			}
		}

		#endif

		pi_cl_team_barrier();
		temp = norm(&sel_col[0],N_CHANNELS);
		pi_cl_team_barrier();

    	#if NUM_CORES > 1

		for(int j = start_NC; j < (start_NC + blockSize_NC); j++){
			if(j == i){
				if(sel_col[0] >= 0){
					v[j] = sel_col[j] + temp;
				}else{
					v[j] = sel_col[j] - temp;
				}
			}else{
				v[j] = sel_col[j];
			}
		}
		pi_cl_team_barrier();

		#else

		for(j=0;j<N_CHANNELS;j++){
			if(j == i){
				if(sel_col[0] >= 0){
					v[j] = sel_col[j] + temp;
				}else{
					v[j] = sel_col[j] - temp;
				}
			}else{
				v[j] = sel_col[j];
			}
		}

		#endif

		temp = zero;

		pi_cl_team_barrier();

    	#if NUM_CORES > 1
        buffer[pi_core_id()]=0;
        idx[pi_core_id()]=0;


		for(int j = start_NC; j < (start_NC + blockSize_NC); j++){

			buffer[pi_core_id()] = buffer[pi_core_id()] + v[idx[pi_core_id()]+start_NC]*v[idx[pi_core_id()]+start_NC];
			idx[pi_core_id()] = idx[pi_core_id()] + 1;
		}

		pi_cl_team_barrier();
        if(pi_core_id()==0){
            for(j=0; j<NUM_CORES; j++){
                temp += buffer[j];
            }
		}

		pi_cl_team_barrier();

		#else

		for (int j = 0; j < N_CHANNELS; j++){
			temp = temp + v[j]*v[j];
		}
		#endif

    	#if NUM_CORES > 1

		for(int j = start_NC; j < (start_NC + blockSize_NC); j++){
			if(j >= i){
				for(int k = i; k < N_CHANNELS; k++){
					H[j][k] = H[j][k] - two*v[k]*v[j]/temp;
				}
			}
		}
		pi_cl_team_barrier();

		#else

		for (j = i; j < N_CHANNELS; j++){
			for(int k = i; k < N_CHANNELS; k++){
				H[j][k] = H[j][k] - two*v[k]*v[j]/temp;
			}
		}

		#endif

		matMul(&H[0][0],&R[0][0],&R_temp[0][0],N_CHANNELS,N_CHANNELS,EV_WINDOWS_SIZE);
		matMul(&Q[0][0],&H[0][0],&Q_temp[0][0],N_CHANNELS,N_CHANNELS,N_CHANNELS);
		
		pi_cl_team_barrier();

		#if NUM_CORES > 1

		for(int j = start_NC; j < (start_NC + blockSize_NC); j++){
			for(int k = 0; k < EV_WINDOWS_SIZE; k++){
				R[j][k] = R_temp[j][k];
			}
		}
		pi_cl_team_barrier();

		#else

		for(j = 0;j < N_CHANNELS; j++){//??????
			for(int k = 0; k < EV_WINDOWS_SIZE; k++){
				R[j][k] = R_temp[j][k];
			}
		}

		#endif
		
		#if NUM_CORES > 1

		for(int j = start_NC; j < (start_NC + blockSize_NC); j++){
			for(int k = 0; k < N_CHANNELS; k++){
				Q[j][k] = Q_temp[j][k];
			}
		}
		pi_cl_team_barrier();

		#else

		for(j = 0;j < N_CHANNELS; j++){//??????
			for(int k = 0; k < N_CHANNELS; k++){
				Q[j][k] = Q_temp[j][k];
			}
		}

		#endif

	}
	return;
}

