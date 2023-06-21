#include "pmsis.h"
#include "stats.h" 
#include "inputData/qr40X40.h"
#include "inputData/data40X40.inc"
#include <math.h>

#define STATS 
#define STACK_SIZE 1024
#define MOUNT 1
#define UNMOUNT 0
#define CID 0
#define DEBUG

//Q is transposed
PI_L1 float Q[N_ROW][N_ROW]; 
PI_L1 float R[N_ROW][N_COL];

//PI_L1 float n;
//PI_L1 int i2;
PI_L1 static float buffer[NUM_CORES];
//PI_L1 static int idx[NUM_CORES];
PI_L1 float rk=0.0;
PI_L1 static float temp[NUM_CORES];
PI_L1 static float one = 1.0;

int cyclesArr[NUM_CORES];
int instrArr[NUM_CORES];
int activeArr[NUM_CORES];
int ldEXTArr[NUM_CORES];
int TCDMArr[NUM_CORES];
int ldstallArr[NUM_CORES];
int imissArr[NUM_CORES];



void cluster_main();
void qr_gramSmidt(float Q[][N_ROW], float R[][N_COL], float input[][N_ROW]);

	void pe_entry(void *arg){
		cluster_main();}

	void cluster_entry(void *arg){
		pi_cl_team_fork((NUM_CORES), pe_entry, 0);}

static int test_entry(){

	struct pi_device cluster_dev;
	struct pi_cluster_conf cl_conf;
	struct pi_cluster_task cl_task;

	pi_cluster_conf_init(&cl_conf);
	pi_open_from_conf(&cluster_dev, &cl_conf);
	if (pi_cluster_open(&cluster_dev)){
		return -1;
	}

	//pi_cluster_send_task_to_cl(&cluster_dev, pi_cluster_task(&cl_task, cluster_entry, NULL));

	pi_cluster_task(&cl_task, cluster_entry, NULL);
	cl_task.stack_size = STACK_SIZE;
	cl_task.slave_stack_size = STACK_SIZE;

	pi_cluster_send_task_to_cl(&cluster_dev, &cl_task);

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
				Q[j][i]=0;            
			}
		}


		for(int i=0; i<N_ROW; i++){
			for(int j=0; j<N_COL; j++){
				R[i][j]=0;            
			}
		}
	}
    
	pi_cl_team_barrier();
    
	pi_perf_reset();
	pi_perf_start();
    
    qr_gramSmidt(Q, R, input);

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
	
	#ifdef DEBUG
 	// if(pi_core_id()==0){
	// 	printf("\n\nQ = \n");
	// 	for(int i=0; i<N_ROW; i++){
	// 		for(int j=0; j<N_ROW; j++)
	// 			printf("%.20f ", Q[j][i]);
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
	// 				res[i][j] += Q[k][i] * R[k][j];
	// 			}
	// 			printf("%.20f ", res[i][j]);
	// 		}
	// 		printf("\n");
	// 	}

	// 	/*printf("\n\nResult = \n");
	// 	for(int i=0; i<N_ROW; i++){
	// 		for(int j=0; j<N_COL; j++)
	// 			printf("%f ", res[i][j]);
	// 		printf("\n");
	// 	}*/
	// }

	pi_cl_team_barrier();
	#endif
}

inline float Sqrt(float x) {
	float res;
	asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
	return res;
}

#pragma GCC push_options
#pragma GCC optimize ("-O3")
__attribute__ ((noinline))
void  norm(float *v, int row, int column, int k, int blockSize_column, int start_column) { 
	int i2 = 0;
	float n = 0.0f;
	int j, idx;

	#if NUM_CORES > 1
	int core_id = pi_core_id();

	int start_matrice = start_column;

	buffer[core_id]=0;
	idx=0;
	
	for(j = start_column; j<start_column + (blockSize_column & 0xfffffffe); j+=2){
		float t0 = v[idx+start_matrice];
		float t1 = v[idx+1+start_matrice];
		asm volatile ("":::"memory");
		buffer[core_id] = buffer[core_id] + t0*t0;
		buffer[core_id] = buffer[core_id] + t1*t1;
		idx+=2;
		asm volatile ("":::"memory");

	}
	if(blockSize_column & 0x1){
		buffer[core_id] = buffer[core_id] + v[idx-1+start_matrice]*v[idx-1+start_matrice];
	}
		
	pi_cl_team_barrier();
			
	if(core_id==0)
	{
		for(j=0; j<NUM_CORES; j++){
			n += buffer[j];
		}
		R[k][k] = Sqrt(n);
		rk = 1/R[k][k];
	}
	pi_cl_team_barrier();
	#else

	for(j=0; (j+1)<column; j+=2){
		float t0 = v[i2];float t1 = v[i2+1];
		float temp = t0*t0;
		temp = temp + t1*t1;
		n = n + temp;
		i2+=2*1;
	}
		
	if(j<column){
		float t0 = v[i2];
		float temp = t0*t0;
		n = n + temp;
		i2+=1*1;
	}

	R[k][k] = Sqrt(n);
        rk = 1/R[k][k];

	#endif
//	return sqrt(n);
	
}

void qr_gramSmidt(float Q[][N_ROW], float R[][N_COL], float input[][N_ROW]){
	int j;

	#if NUM_CORES > 1
	int blockSize_ROW = N_ROW/NUM_CORES;
	int start_ROW = pi_core_id()*blockSize_ROW;
	if(pi_core_id()==(NUM_CORES - 1)){
		blockSize_ROW = N_ROW - (NUM_CORES - 1)* blockSize_ROW;}
 
         int blockSize_column = N_COL/NUM_CORES;
         int start_column = pi_core_id()*blockSize_column;
         if(pi_core_id()==(NUM_CORES - 1)){
	          blockSize_column = N_COL - (NUM_CORES - 1)* blockSize_column;
         }
	
	#endif

	for(int k=0; k<N_COL; k++){
		#if NUM_CORES > 1
		for(j = start_ROW; ((j+1)<start_ROW + blockSize_ROW); j+=2){
			float in0 = input[k][j];float in1 = input[k][j+1];
			Q[k][j] = in0;
			Q[k][j+1] = in1;
		}
		if(j<start_ROW + blockSize_ROW){
			float in0 = input[k][j];
			Q[k][j] = in0;
		}
		#else
		for(j=0; (j+1)<N_ROW; j+=2){
			float in0 = input[k][j];float in1 = input[k][j+1];
			Q[k][j] = in0;
			Q[k][j+1] = in1;
		}
		if(j<N_ROW){
			float in0 = input[k][j];
			Q[k][j] = in0;
		}
		#endif		

		#if NUM_CORES > 1
		pi_cl_team_barrier();
		#endif

		for(int i=0; i<k; i++){
			#if NUM_CORES > 1
			temp[pi_core_id()]=0;
			for(j = start_ROW; ((j+1)<start_ROW + blockSize_ROW); j+=2){
				float ji0 = Q[i][j];float jk0 = Q[k][j];
				float ji1 = Q[i][j+1];float jk1 = Q[k][j+1];
				temp[pi_core_id()] = temp[pi_core_id()] + (ji0 * jk0) + (ji1 * jk1);				
			}
			if((j<start_ROW + blockSize_ROW)){
				float ji0 = Q[i][j];float jk0 = Q[k][j];
				temp[pi_core_id()] = temp[pi_core_id()] + (ji0 * jk0);
			}
			
			pi_cl_team_barrier();

			if(pi_core_id()==0)
				for(j=0; j<NUM_CORES; j++){
			    		R[i][k]+=temp[j];
			    	}
			pi_cl_team_barrier();

			#else

			for(j=0; (j+1)<N_ROW; j+=2){
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];float Qjk1 = Q[k][j+1];
				float Rik  = R[i][k];
				R[i][k] = Rik + (Qji0 * Qjk0) + (Qji1 * Qjk1);
            }
            if(j<N_ROW){				
				float Qji0 = Q[i][j]; float Qjk0 = Q[k][j];
				float Rik  = R[i][k];
				float temp = (Qji0 * Qjk0);
				R[i][k] = Rik + temp;
            }
			#endif
			
			#if NUM_CORES > 1
			for(j = start_ROW; (j+1<start_ROW + blockSize_ROW); j+=2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];float Qjk1 = Q[k][j+1];
				
				//hal_compiler_barrier();
                		
				Q[k][j] = Qjk0 - Rik*Qji0;
                Q[k][j+1] = Qjk1 - Rik*Qji1;
            }
			if(j<start_ROW + blockSize_ROW){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik*Qji0;
			}
			#else
			for(j=0; j+1<N_ROW; j+=2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];float Qjk1 = Q[k][j+1];

				Q[k][j] = Qjk0 - Rik*Qji0;
				Q[k][j+1] = Qjk1 - Rik*Qji1;
			}
			if(j<N_ROW){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik*Qji0;
			}
			#endif
			
			#if NUM_CORES > 1 
			pi_cl_team_barrier();
			#endif
		}
		#if NUM_CORES > 1 
		pi_cl_team_barrier();
		#endif
		
		// if(pi_core_id()==0){
		// 	printf("\n\nQ = \n");
		// 	for(int i=0; i<N_ROW; i++){
		// 		for(int j=0; j<N_ROW; j++)
		// 			printf("%.20f ", Q[j][i]);
		// 		printf("\n");
		// 	}

		// 	printf("\n\nR = \n");
		// 	for(int i=0; i<N_ROW; i++){
		// 		for(int j=0; j<N_COL; j++)
		// 			printf("%.20f ", R[i][j]);
		// 		printf("\n");
		// 	}
		// }

		#if NUM_CORES > 1
		norm(&Q[k][0],N_ROW,N_ROW,k, blockSize_column, start_column);
		#else
		norm(&Q[k][0],N_ROW,N_ROW,k, 0, 0);
		#endif
		//R[k][k] = norm(&Q[k][0],N_ROW,N_ROW);
		// printf("R[k][k] = %.20f \n", R[k][k]);
		//rk = one/R[k][k];

		
		#if NUM_CORES > 1
		for( j = start_ROW; (j+1)<start_ROW + (blockSize_ROW & 0xFFFFFFFE); j+=2){
			float Qjk0 = Q[k][j];float Qjk1 = Q[k][j+1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j+1] = Qjk1 * rk;
			
		} if (j<start_ROW + (blockSize_ROW & 0xFFFFFFFE)){
			float Qjk0 = Q[j][k];
			Q[k][j] = Qjk0 * rk;
		}
		#else
		for( j=0; (j+1)<N_ROW; j+=2){
			float Qjk0 = Q[k][j];float Qjk1 = Q[k][j+1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j+1] = Qjk1 * rk;
		} if (j<N_ROW){
			float Qjk0 = Q[k][j];
			Q[k][j] = Qjk0 * rk;
		}
		#endif
		
		#if NUM_CORES > 1 
		if(blockSize_ROW & 0x0001){
			Q[k][j] = Q[k][j]*rk;}
		pi_cl_team_barrier();
		#endif
		
	}
	return;
}

