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
PI_L1 int i2;
PI_L1 static float buffer[NUM_CORES];
PI_L1 static int idx[NUM_CORES];
PI_L1 float rk=0.0;
PI_L1 static float temp[NUM_CORES];
PI_L1 static float one = 1.0;

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
    
	/*if(pi_core_id()==0){
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
	}*/

	pi_cl_team_barrier();
}

inline float Sqrt(float x) {
	float res;
	asm("fsqrt.s %0, %1":"=f"(res):"f"(x));
	return res;
}

float norm(float *v, int row, int column){ 
	i2=0;
	n = 0.0f;
	int j;
	#if NUM_CORES > 1

	int blockSize_column = column/NUM_CORES;
	int start_column = pi_core_id()*blockSize_column;
	int start_matrice = pi_core_id()*blockSize_column;

	if(pi_core_id()==(NUM_CORES - 1)){
		blockSize_column = column - (NUM_CORES - 1)* blockSize_column;
	}
	
	buffer[pi_core_id()]=0;
	idx[pi_core_id()]=0;
	
	for(j = start_column; (j<column) && (j<start_column + blockSize_column); j++){
		buffer[pi_core_id()] = buffer[pi_core_id()] + v[idx[pi_core_id()]+start_matrice]*v[idx[pi_core_id()]+start_matrice];
		idx[pi_core_id()]+=1;
	}
		
	pi_cl_team_barrier();
			
	if(pi_core_id()==0)
		for(j=0; j<NUM_CORES; j++){
			n += buffer[j];
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

	#endif
	return sqrt(n);
	
}


void qr_gramSmidt(float Q[][N_CHANNELS], float R[][EV_WINDOWS_SIZE], float input[][N_CHANNELS]){
	int j;

	#if NUM_CORES > 1
	int blockSize_NC = N_CHANNELS/NUM_CORES;
	int start_NC = pi_core_id()*blockSize_NC;
	if(pi_core_id()==(NUM_CORES - 1)){
		blockSize_NC = N_CHANNELS - (NUM_CORES - 1)* blockSize_NC;}
	#endif

	for(int k=0; k<EV_WINDOWS_SIZE; k++){
		#if NUM_CORES > 1
		for(j = start_NC; ((j+1)<start_NC + blockSize_NC); j+=2){
			float in0 = input[k][j];float in1 = input[k][j+1];
			Q[k][j] = in0;
			Q[k][j+1] = in1;
		}
		if(j<start_NC + blockSize_NC){
			float in0 = input[k][j];
			Q[k][j] = in0;
		}
		#else
		for(j=0; (j+1)<N_CHANNELS; j+=2){
			float in0 = input[k][j];float in1 = input[k][j+1];
			Q[k][j] = in0;
			Q[k][j+1] = in1;
		}
		if(j<N_CHANNELS){
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
			for(j = start_NC; ((j+1)<start_NC + blockSize_NC); j+=2){
				float ji0 = Q[i][j];float jk0 = Q[k][j];
				float ji1 = Q[i][j+1];float jk1 = Q[k][j+1];
				temp[pi_core_id()] = temp[pi_core_id()] + (ji0 * jk0) + (ji1 * jk1);				
			}
			if((j<start_NC + blockSize_NC)){
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

			for(j=0; (j+1)<N_CHANNELS; j+=2){
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];float Qjk1 = Q[k][j+1];
				float Rik  = R[i][k];
				R[i][k] = Rik + (Qji0 * Qjk0) + (Qji1 * Qjk1);
            }
            if(j<N_CHANNELS){				
				float Qji0 = Q[i][j]; float Qjk0 = Q[k][j];
				float Rik  = R[i][k];
				float temp = (Qji0 * Qjk0);
				R[i][k] = Rik + temp;                	
            }
			#endif
			
			#if NUM_CORES > 1
			for(j = start_NC; (j+1<start_NC + blockSize_NC); j+=2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];float Qjk1 = Q[k][j+1];
				
				//hal_compiler_barrier();
                		
				Q[k][j] = Qjk0 - Rik*Qji0;
                Q[k][j+1] = Qjk1 - Rik*Qji1;
            }
			if(j<start_NC + blockSize_NC){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik*Qji0;
			}
			#else
			for(j=0; j+1<N_CHANNELS; j+=2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j+1];float Qjk1 = Q[k][j+1];

				Q[k][j] = Qjk0 - Rik*Qji0;
				Q[k][j+1] = Qjk1 - Rik*Qji1;
			}
			if(j<N_CHANNELS){
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
		
		R[k][k] = norm(&Q[k][0],N_CHANNELS,N_CHANNELS);
		rk = one/R[k][k];
		
		#if NUM_CORES > 1
		for( j = start_NC; (j+1)<start_NC + (blockSize_NC & 0xFFFE); j+=2){
			float Qjk0 = Q[k][j];float Qjk1 = Q[k][j+1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j+1] = Qjk1 * rk;
			
		} if (j<start_NC + (blockSize_NC & 0xFFFE)){
			float Qjk0 = Q[j][k];
			Q[k][j] = Qjk0 * rk;
		}
		#else
		for( j=0; (j+1)<N_CHANNELS; j+=2){
			float Qjk0 = Q[k][j];float Qjk1 = Q[k][j+1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j+1] = Qjk1 * rk;
		} if (j<N_CHANNELS){
			float Qjk0 = Q[k][j];
			Q[k][j] = Qjk0 * rk;
		}
		#endif
		
		#if NUM_CORES > 1 
		if(blockSize_NC & 0x0001){
			Q[k][j] = Q[k][j]*rk;}
		pi_cl_team_barrier();
		#endif
		
	}
	return;
}

