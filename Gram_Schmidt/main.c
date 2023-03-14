#include "pmsis.h"
#include "stats.h" 
#include "qrGS.h"
#include "dataGS.inc"
#include <math.h>

#define STATS 
#define STACK_SIZE 2048
#define MOUNT 1
#define UNMOUNT 0
#define CID 0



PI_L1 float Q[N_CHANNELS][N_CHANNELS]; 
PI_L1 float R[N_CHANNELS][EV_WINDOWS_SIZE];

PI_L1 float n;
PI_L1 static float buffer[NUM_CORES];
PI_L1 float rk=0.0;
PI_L1 static float temp[NUM_CORES];
PI_L1 static float one = 1.0;

void cluster_main();
void qr_gramSmidt(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE], float input[][N_CHANNELS]);

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
		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<N_CHANNELS; j++){
				Q[j][i]=0;            
			}
		}


		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<EV_WINDOWS_SIZE; j++){
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
    
	pi_cl_team_barrier();
    
	if(pi_core_id()==0){
	printf("\n\nQ = \n");
		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<N_CHANNELS; j++)
				printf("%f ", Q[j][i]);
			printf("\n");
		}

	printf("\n\nR = \n");
		for(int i=0; i<N_CHANNELS; i++){
			for(int j=0; j<EV_WINDOWS_SIZE; j++)
				printf("%f ", R[i][j]);
			printf("\n");
		}
	}

	pi_cl_team_barrier();
}

inline float Sqrt(float x) {
	float res;
	asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
	return res;
}

float norm(float *v, int row, int column){ 
	n = 0.0f;
	int j;
	int idx;
	float res;

#if NUM_CORES > 1
	int core_id = pi_core_id();
	int blockSize_column = (column + NUM_CORES - 1) / NUM_CORES;
	int start_column = core_id * blockSize_column;
	int end_column = start_column + blockSize_column;
	if(end_column > column){
		end_column = column;
	}
	
	idx = 0;
	res = 0.0f;
	for(j = start_column; j < end_column; j++){
		res = res +  v[idx+start_column] * v[idx+start_column];
		idx++;
	}
	buffer[core_id] = res;
		
	pi_cl_team_barrier();
			
	if(core_id == 0){
		n = 0.0;
		float temp0 = buffer[0];
		float temp1 = buffer[1];
		for(int c=2; c<NUM_CORES; c+=2){
			float val0 = buffer[c];
			float val1 = buffer[c+1];
			hal_compiler_barrier();
			temp0 += val0;
			temp1 += val1;
		}
        	n = temp0 + temp1;
		n = Sqrt(n);
	}

	pi_cl_team_barrier();
#else
	idx = 0;
	for(j=0; j<column/2; j++){
		float t0 = v[idx];
		float t1 = v[idx+1];
		n += t0*t0;
		n += t1*t1;
		idx += 2;
	}
		
	if(column & 0x0001){
		float t0 = v[row-1];
		n = n + t0 * t0;
	}
	n = Sqrt(n);

#endif
	return n;
	
}


void qr_gramSmidt(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE], float input[][N_CHANNELS]){
	int j;
	float res;

	#if NUM_CORES > 1
	int core_id = pi_core_id();
	int blockSize_NC = (N_CHANNELS + NUM_CORES - 1) / NUM_CORES;
	int start_NC = core_id * blockSize_NC;
	if(core_id == (NUM_CORES-1)){
		blockSize_NC = N_CHANNELS - (NUM_CORES-1) * blockSize_NC;
	}
	int end_NC = start_NC + (blockSize_NC & 0xfffffffe);
	#endif

	for(int k=0; k<EV_WINDOWS_SIZE; k++){
		#if NUM_CORES > 1
		float *restrict  in = (float *) (&input[k][start_NC]);
		float *restrict out = (float *) (&Q[k][start_NC]);
		for(j = start_NC; j < end_NC; j+=2){
			float in0 = *(in++);
			float in1 = *(in++);
			*(out++) = in0;
			*(out++) = in1;
		}
		if(blockSize_NC & 0x1){
			float in0 = input[k][end_NC];
			Q[k][end_NC] = in0;
		}
		#else
		for(j=0; j<(N_CHANNELS & 0xfffffffe); j+=2){
			float in0 = input[k][j];
			float in1 = input[k][j+1];
			Q[k][j] = in0;
			Q[k][j+1] = in1;
		}
		if(N_CHANNELS & 0x1){
			float in0 = input[k][N_CHANNELS - 1];
			Q[k][N_CHANNELS - 1] = in0;
		}
		#endif		

		#if NUM_CORES > 1
		pi_cl_team_barrier();
		#endif

		for(int i=0; i<k; i++){
			#if NUM_CORES > 1
			temp[pi_core_id()]=0;
			for(j = start_NC; j < end_NC; j+=2){
				float ji0 = Q[i][j];
				float jk0 = Q[k][j];
				float ji1 = Q[i][j+1];
				float jk1 = Q[k][j+1];
				temp[core_id] = temp[core_id] + (ji0 * jk0) + (ji1 * jk1);				
			}
			if(blockSize_NC & 0x1){
				float ji0 = Q[i][end_NC];
				float jk0 = Q[k][end_NC];
				temp[core_id] = temp[core_id] + (ji0 * jk0);
			}
			
			pi_cl_team_barrier();

			if(core_id == 0) {
				float temp0 = temp[0];
				float temp1 = temp[1];
        			for(int c=2; c<NUM_CORES; c+=2){
					float val0 = temp[c];
					float val1 = temp[c+1];
					hal_compiler_barrier();
					temp0 += val0;
					temp1 += val1;
				}
				R[i][k] = temp0 + temp1;
			}			
			pi_cl_team_barrier();

			#else

			for(j=0; j<(N_CHANNELS & 0xfffffffe); j+=2){
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];
				float Qjk1 = Q[k][j+1];
				float Rik  = R[i][k];
				R[i][k] = Rik + (Qji0 * Qjk0) + (Qji1 * Qjk1);
			}
			if(N_CHANNELS & 0x1){				
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Rik  = R[i][k];
				float temp = (Qji0 * Qjk0);
				R[i][k] = Rik + temp;                	
			}
			#endif
			
			#if NUM_CORES > 1
			for(j = start_NC; j<end_NC; j+=2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];
				float Qjk1 = Q[k][j+1];
				
				//hal_compiler_barrier();
                		
				Q[k][j] = Qjk0 - Rik*Qji0;
				Q[k][j+1] = Qjk1 - Rik*Qji1;
			}
			if(blockSize_NC & 0x1){
				float Rik = R[i][k];
				float Qji0 = Q[i][end_NC];
				float Qjk0 = Q[k][end_NC];
				Q[k][end_NC] = Qjk0 - Rik*Qji0;
			}
			#else
			for(j=0; j<(N_CHANNELS & 0xfffffffe); j+=2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];
				float Qjk1 = Q[k][j+1];

				Q[k][j] = Qjk0 - Rik*Qji0;
				Q[k][j+1] = Qjk1 - Rik*Qji1;
			}
			if(N_CHANNELS & 0x1){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik*Qji0;
			}
			#endif
			
			#if NUM_CORES > 1 
			pi_cl_team_barrier();
			#endif
		}

		res = norm(&Q[k][0],N_CHANNELS,N_CHANNELS);

		#if NUM_CORES > 1
		if(core_id == 0) {
			R[k][k] = res;
	       	        rk = one / res;
		}
		pi_cl_team_barrier();
		#else
		R[k][k] = res;
		rk = one / res;
		#endif
		
		#if NUM_CORES > 1
		for( j = start_NC; j < end_NC; j+=2){
			float Qjk0 = Q[k][j];
			float Qjk1 = Q[k][j+1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j+1] = Qjk1 * rk;
			
		}
		
		if (blockSize_NC & 0x1){
			float Qjk0 = Q[k][end_NC];
			Q[k][end_NC] = Qjk0 * rk;
		}
		#else
		for( j=0; j<(N_CHANNELS & 0xfffffffe); j+=2){
			float Qjk0 = Q[k][j];
			float Qjk1 = Q[k][j+1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j+1] = Qjk1 * rk;
		}
		if (N_CHANNELS & 0x1){
			float Qjk0 = Q[k][j];
			Q[k][j] = Qjk0 * rk;
		}
		#endif
		
		#if NUM_CORES > 1 
		pi_cl_team_barrier();
		#endif
		
	}
	return;
}

