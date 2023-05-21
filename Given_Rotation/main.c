#include "pmsis.h"
#include "qr.h"
#include "data.inc"
#include <math.h>
#include "stats.h" 


#define STATS 

#define MAX(a,b) ((a) > (b))? (a): (b)

#define STACK_SIZE 1024
#define MOUNT 1
#define UNMOUNT 0
#define CID 0

PI_L1 float a[MAX(EV_WINDOWS_SIZE,N_CHANNELS)];
PI_L1 float b[MAX(EV_WINDOWS_SIZE,N_CHANNELS)];

PI_L1 float Qt[N_CHANNELS][N_CHANNELS];
PI_L1 float Rt[N_CHANNELS][EV_WINDOWS_SIZE];

PI_L1 float one = 1.0f;
PI_L1 float zero = 0.0f;
PI_L1 float two = 2.0f;



PI_L1 unsigned long cycles = 0;
PI_L1 unsigned long instr = 0;
PI_L1 unsigned long active = 0;
PI_L1 unsigned long ldext = 0;
PI_L1 unsigned long tcdmcont = 0;
PI_L1 unsigned long ldstall = 0;
PI_L1 unsigned long imiss = 0;
PI_L1 unsigned long apu_cont = 0;
PI_L1 unsigned long apu_dep = 0;
PI_L1 unsigned long apu_type = 0;
PI_L1 unsigned long apu_wb = 0;

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
        //Q
	for(int i=0; i<N_CHANNELS; i++){
	    for(int j=0; j<N_CHANNELS; j++){
                if(i==j) {
                    Qt[i][j] = one;
		 } else Qt[i][j] = zero;
            }
        }
        
        //Rt
        for(int i=0; i<N_CHANNELS; i++){
            for(int j=0; j<EV_WINDOWS_SIZE; j++){
               Rt[i][j] = input[i][j];}

        }        
    }
	
    pi_cl_team_barrier();
    pi_perf_reset();
    pi_perf_start();

    factorization(Qt, Rt);

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
    printf("[%d] apu_cont = %lu\n", id, apu_cont/REPEAT);
    printf("[%d] apu_dep = %lu\n", id, apu_dep/REPEAT);
    printf("[%d] apu_type = %lu\n", id, apu_type/REPEAT);
    printf("[%d] apu_wb = %lu\n", id, apu_wb/REPEAT);

 
   
    pi_cl_team_barrier();
    
    /*if (pi_core_id()==0){
        printf("\n\nQ = \n");
        for(int i=0; i<N_CHANNELS; i++){
            for(int j=0; j<N_CHANNELS; j++){
        		printf("%f ", Qt[i][j]);
			}
        printf("\n");
		}
        
       
        printf("\n\nR = \n");
        for(int i=0; i<N_CHANNELS; i++){
            for(int j=0; j<EV_WINDOWS_SIZE; j++){
        		printf("%f ", Rt[i][j]);
			}
            printf("\n");
		}
    }*/

    pi_cl_team_barrier();
    
	return;
}

void inline givens_rotation(float x, float y, float *c, float *s){
    float r;
    
    if ((fabs(y)) >= (fabs(x))){ 
        r = x / y;
        *s = one/sqrtf(one + r*r);
        *c = (*s)*r;
	}
    else{
        r = y / x;
        *c = one / sqrtf(one + r*r);
        *s = (*c)*r;
	}
    return;
}

void factorization(float Qt[][N_CHANNELS], float Rt[][EV_WINDOWS_SIZE]){

	float c=0; 
	float s=0;
	float xi;
	float xj;
	
	int temp;

	
	#if NUM_CORES > 1
	int blockSize_EV = (EV_WINDOWS_SIZE + NUM_CORES - 1)/NUM_CORES;
	int start_EV = pi_core_id()*blockSize_EV;

	int blockSize_NC = (N_CHANNELS + NUM_CORES - 1)/ NUM_CORES;
	int start_NC = pi_core_id()*blockSize_NC; 

	if(pi_core_id()==(NUM_CORES - 1)){
		blockSize_EV = EV_WINDOWS_SIZE - (NUM_CORES - 1)* blockSize_EV;
		blockSize_NC = N_CHANNELS * (NUM_CORES - 1)*blockSize_NC;
	}
	
	#endif

	for(int k=0; k<EV_WINDOWS_SIZE; k++){
		for(int j=k+1; j<N_CHANNELS; j++){ 

	    	#if NUM_CORES > 1
			int blockSize_NC = ((N_CHANNELS - k) + (NUM_CORES) - 1)/ (NUM_CORES);
			int start_NC = (pi_core_id() % (NUM_CORES))*blockSize_NC;
			#endif

			if(Rt[j][k]!=0){
				xi=Rt[k][k];
				xj=Rt[j][k];
		       
				givens_rotation(xi, xj, &c, &s);
				
				pi_cl_team_barrier();	
							
				#if NUM_CORES > 1
				if(start_EV + blockSize_EV >= k){
					if(start_EV < k){
						temp = k;
						for(int i = k; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV + blockSize_EV); i+=2){
							float a0 = Rt[k][i];
							float a1 = Rt[k][i+1];						
							a[i] = a0;
							a[i+1] = a1;
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float a0 = Rt[k][temp];						
							a[temp] = a0;
						}
					} else {
						temp = start_EV;
						for(int i = start_EV; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV + blockSize_EV); i+=2){
							float a0 = Rt[k][i];
							float a1 = Rt[k][i+1];						
							a[i] = a0;
							a[i+1] = a1;
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float a0 = Rt[k][temp];						
							a[temp] = a0;
						}
					}
				}
				#else
			    	for(int i=k; i+1<EV_WINDOWS_SIZE; i+=2){
					float a0 = Rt[k][i];
					float a1 = Rt[k][i+1];
					a[i] = a0;
					a[i+1] = a1;	
					temp = i+two;
				}
				if(temp < EV_WINDOWS_SIZE){
					float a0 = Rt[k][temp];
					a[temp] = a0;
				}
			    #endif
				
				pi_cl_team_barrier();

				#if NUM_CORES > 1
				if(start_EV + blockSize_EV >= k){
					if(start_EV < k){
						temp = k;
						for(int i = k; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV+blockSize_EV); i+=2){
							float b0 = Rt[j][i];
							float b1 = Rt[j][i+1];
							b[i] = b0;
							b[i+1] = b1;
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float b0 = Rt[j][temp];
							b[temp] = b0;
						}								
					} else {
						temp = start_EV;
						for(int i = start_EV; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV+blockSize_EV); i+=2){
							float b0 = Rt[j][i];
							float b1 = Rt[j][i+1];
							b[i] = b0;
							b[i+1] = b1;
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float b0 = Rt[j][temp];
							b[temp] = b0;
						}
					}
				}
				#else
				for(int i=k; i+1<EV_WINDOWS_SIZE; i+=2){
					float b0 = Rt[j][i];
					float b1 = Rt[j][i+1];					
					b[i] = b0;
					b[i+1] = b1;
					temp = i+two;
				}
				if(temp < EV_WINDOWS_SIZE){
					float b0 = Rt[j][temp];
					b[temp] = b0;
				}
				#endif

				pi_cl_team_barrier();

				#if NUM_CORES > 1
				if(start_EV + blockSize_EV >= k){
					if(start_EV < k){
						temp = k;
						for(int i = k; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV + blockSize_EV); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[k][i]     = (c * a0) + (s * b0);
							Rt[k][i+1]     = (c * a1) + (s * b1);
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float a0 = a[temp];float b0 = b[temp];
							Rt[k][temp]     = (c * a0) + (s * b0);
						}									
					} else {
						temp = start_EV;
						for(int i = start_EV; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV + blockSize_EV); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[k][i]     = (c * a0) + (s * b0);
							Rt[k][i+1]     = (c * a1) + (s * b1);
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float a0 = a[temp];float b0 = b[temp];
							Rt[k][temp]     = (c * a0) + (s * b0);
						}
					}
				}
				#else
				for(int i=k; i+1<EV_WINDOWS_SIZE; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];
					Rt[k][i]     = (c * a0) + (s * b0);
					Rt[k][i+1]     = (c * a1) + (s * b1);
					temp = i+two;
				}
				if(temp < EV_WINDOWS_SIZE){
					float a0 = a[temp];float b0 = b[temp];
					Rt[k][temp]     = (c * a0) + (s * b0);
				}
				#endif
				
				#if NUM_CORES > 1
				if(start_EV + blockSize_EV >= k){
					if(start_EV < k){
						temp = k;
						for(int i = k; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV+blockSize_EV); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[j][i]     = (c * b0) - (s * a0);
							Rt[j][i+1]     = (c * b1) - (s * a1);
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float a0 = a[temp];float b0 = b[temp];
							Rt[j][temp]     = (c * b0) - (s * a0);
						}										
					} else {
						temp = start_EV;
						for(int i = start_EV; (i+1<EV_WINDOWS_SIZE)&&(i+1<start_EV+blockSize_EV); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[j][i] = (c * b0) - (s * a0);
							Rt[j][i+1] = (c * b1) - (s * a1);
							temp = i+two;
						}
						if((temp < EV_WINDOWS_SIZE)&&(temp<start_EV + blockSize_EV)){
							float a0 = a[temp];float b0 = b[temp];
							Rt[j][temp] = (c * b0) - (s * a0);
						}
					}
				}
				#else
				for(int i=k; i+1<EV_WINDOWS_SIZE; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];
					Rt[j][i]     = (c * b0) - (s * a0);
					Rt[j][i+1]     = (c * b1) - (s * a1);
					temp = i+two;
				}
				if(temp < EV_WINDOWS_SIZE){
					float a0 = a[temp];float b0 = b[temp];
					Rt[j][temp]     = (c * b0) - (s * a0);
				}
				#endif

				pi_cl_team_barrier();

				#if NUM_CORES > 1
				if(start_NC<N_CHANNELS){
					temp = start_NC;
					for(int i = start_NC; (i+1<N_CHANNELS)&&(i+1<start_NC + blockSize_NC); i+=2){
						float a0 = Qt[i][k];
						float a1 = Qt[i+1][k];
						a[i] = a0;
						a[i+1] = a1;
						temp = i+two;
					}
					if((temp < N_CHANNELS)&&(temp<start_NC + blockSize_NC)){
						float a0 = Qt[temp][k];
						a[temp] = a0;
					}
				}
				#else
				for(int i=0; i+1<N_CHANNELS; i+=2){
					float a0 = Qt[i][k];
					float a1 = Qt[i+1][k];
					a[i] = a0;
					a[i+1] = a1;
					temp = i+two;
					}
				if(temp < N_CHANNELS){
					float a0 = Qt[temp][k];
					a[temp] = a0;
				}
				#endif

				#if NUM_CORES > 1
				if(start_NC<N_CHANNELS){
					temp = start_NC;
					for(int i = start_NC; (i+1<N_CHANNELS)&&(i+1<start_NC + blockSize_NC); i+=2){
						float b0 = Qt[i][j];
						float b1 = Qt[i+1][j];
						b[i] = b0;
						b[i+1] = b1;
						temp = i+two;
					}
					if((temp < N_CHANNELS)&&(temp<start_NC + blockSize_NC)){
						float b0 = Qt[temp][j];
						b[temp] = b0;
					}
				}
				#else	
				for(int i=0; i+1<N_CHANNELS; i+=2){
					float b0 = Qt[i][j];
					float b1 = Qt[i+1][j];
					b[i] = b0;
					b[i+1] = b1;
				
					temp = i+two;
					}
				if(temp < N_CHANNELS){
					float b0 = Qt[temp][k];
					b[temp] = b0;
				}
				#endif

				pi_cl_team_barrier();				

				#if NUM_CORES > 1
				if(start_NC<N_CHANNELS){
					temp = start_NC;
					for(int i = start_NC; (i+1<N_CHANNELS)&&(i+1<start_NC + blockSize_NC); i+=2){
						float a0 = a[i];float b0 = b[i];
						float a1 = a[i+1];float b1 = b[i+1];						
						Qt[i][k] = (c * a0)+(s * b0);
						Qt[i+1][k] = (c * a1)+(s * b1);
						temp = i+two;					
					}
					if((temp < N_CHANNELS)&&(temp<start_NC + blockSize_NC)){
						float a0 = a[temp];float b0 = b[temp];
						Qt[temp][k] = (c * a0)+(s * b0);
					}
				}
				#else
				for(int i=0; i+1<N_CHANNELS; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];						
					Qt[i][k] = (c * a0)+(s * b0);
					Qt[i+1][k] = (c * a1)+(s * b1);
					temp = i+two;
				}
				if(temp < N_CHANNELS){
					float a0 = a[temp];float b0 = b[temp];
					Qt[temp][k] = (c * a0)+(s * b0);
				}
				#endif

				#if NUM_CORES > 1
				if(start_NC<N_CHANNELS){
					temp = start_NC;
					for(int i = start_NC; (i+1<N_CHANNELS)&&(i+1<start_NC + blockSize_NC); i+=2){
						float a0 = a[i];float b0 = b[i];
						float a1 = a[i+1];float b1 = b[i+1];						
						Qt[i][j] = (c * b0)-(s * a0);
						Qt[i+1][j] = (c * b1)-(s * a1);					
						temp = i+two;
					}
					if((temp < N_CHANNELS)&&(temp<start_NC + blockSize_NC)){
						float a0 = a[temp];float b0 = b[temp];
						Qt[temp][j] = (c * b0)-(s * a0);
					}
				}
				#else
				for(int i=0; i+1<N_CHANNELS; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];						
					Qt[i][j] = (c * b0)-(s * a0);
					Qt[i+1][j] = (c * b1)-(s * a1);					
					temp = i+two;
					}
				if(temp < N_CHANNELS){
					float a0 = a[temp];float b0 = b[temp];
					Qt[temp][j] = (c * b0)-(s * a0);
				}
				#endif
			}
        }
    }
    return;
}

