#include "pmsis.h"
#include "stats.h"
#include "inputData/qr40X40.h"
#include "inputData/data40X40.inc"
#include "math.h"



#define STATS 

#define MAX(a,b) ((a) > (b))? (a): (b)

#define STACK_SIZE 1024
#define MOUNT 1
#define UNMOUNT 0
#define CID 0

PI_L1 float a[MAX(N_COL,N_ROW)];
PI_L1 float b[MAX(N_COL,N_ROW)];

PI_L1 float Qt[N_ROW][N_ROW];
PI_L1 float Rt[N_ROW][N_COL];

PI_L1 float one = 1.0f;
PI_L1 float zero = 0.0f;
PI_L1 float two = 2.0f;

int cyclesArr[NUM_CORES];
int instrArr[NUM_CORES];
int activeArr[NUM_CORES];
int ldEXTArr[NUM_CORES];
int TCDMArr[NUM_CORES];
int ldstallArr[NUM_CORES];
int imissArr[NUM_CORES];



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
   
    pi_cl_team_barrier();

    if (pi_core_id()==0){
        //Q
	for(int i=0; i<N_ROW; i++){
	    for(int j=0; j<N_ROW; j++){
                if(i==j) {
                    Qt[i][j] = one;
		 } else Qt[i][j] = zero;
            }
        }
        
        //Rt
        for(int i=0; i<N_ROW; i++){
            for(int j=0; j<N_COL; j++){
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
    
    // if (pi_core_id()==0){
    //     printf("\n\nQ = \n");
    //     for(int i=0; i<N_ROW; i++){
    //         for(int j=0; j<N_ROW; j++){
    //     		printf("%.20f ", Qt[i][j]);
	// 		}
    //     printf("\n");
	// 	}
        
       
    //     printf("\n\nR = \n");
    //     for(int i=0; i<N_ROW; i++){
    //         for(int j=0; j<N_COL; j++){
    //     		printf("%.20f ", Rt[i][j]);
	// 		}
    //         printf("\n");
	// 	}
	// 	printf("\n\nRes = \n");
	// 	float res[N_ROW][N_COL];
	// 	for (int i = 0; i < N_ROW; i++) {
	// 		for (int j = 0; j < N_COL; j++) {
	// 			res[i][j] = 0;
	// 			for (int k = 0; k < N_ROW; k++) {
	// 				res[i][j] += Qt[i][k] * Rt[k][j];
	// 			}
	// 			printf("%.20f ", res[i][j]);
	// 		}
	// 		printf("\n");
	// 	}
    // }

    pi_cl_team_barrier();
    
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

void factorization(float Qt[][N_ROW], float Rt[][N_COL]){

	float c=0; 
	float s=0;
	float xi;
	float xj;
	
	int i;

	
	#if NUM_CORES > 1
	int blockSize_column = N_COL/NUM_CORES;
	int start_column = pi_core_id()*blockSize_column;

	int blockSize_row = N_ROW/ NUM_CORES;
	int start_row = pi_core_id()*blockSize_row; 

	if(pi_core_id()==(NUM_CORES - 1)){
		blockSize_column = N_COL - (NUM_CORES - 1) * blockSize_column;
		blockSize_row = N_ROW * (NUM_CORES - 1) * blockSize_row;
	}
	
	#endif

	for(int k=0; k<N_COL; k++){
		for(int j=k+1; j<N_ROW; j++){ 

	    	// #if NUM_CORES > 1
			// int blockSize_row = ((N_ROW - k) + (NUM_CORES) - 1)/ (NUM_CORES);
			// int start_row = pi_core_id() * blockSize_row;
			// #endif

			if(Rt[j][k]!=0){
				xi=Rt[k][k];
				xj=Rt[j][k];
		       
				givens_rotation(xi, xj, &c, &s);
				
				pi_cl_team_barrier();	
							
				#if NUM_CORES > 1
				if(start_column + blockSize_column >= k){
					if(start_column < k){
						for(i = k; (i+1<N_COL)&&(i+1<start_column + blockSize_column); i+=2){
							float a0 = Rt[k][i];
							float a1 = Rt[k][i+1];						
							a[i] = a0;
							a[i+1] = a1;
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float a0 = Rt[k][i];
							a[i] = a0;
						}
					} else {
						for(i = start_column; (i+1<N_COL)&&(i+1<start_column + blockSize_column); i+=2){
							float a0 = Rt[k][i];
							float a1 = Rt[k][i+1];						
							a[i] = a0;
							a[i+1] = a1;
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float a0 = Rt[k][i];
							a[i] = a0;
						}
					}
				}
				#else
				for(i=k; i+1<N_COL; i+=2){
					float a0 = Rt[k][i];
					float a1 = Rt[k][i+1];
					a[i] = a0;
					a[i+1] = a1;	
				}
				if(i < N_COL){
					float a0 = Rt[k][i];
					a[i] = a0;
				}
				
				// printf("a = \n");
				// for(int i = 0 ; i < 10; i++){
				// 	printf("%.20f\t",a[i]);
				// }
				// printf("\n");

				#endif
				
				pi_cl_team_barrier();

				#if NUM_CORES > 1
				if(start_column + blockSize_column >= k){
					if(start_column < k){
						for(i = k; (i+1<N_COL)&&(i+1<start_column+blockSize_column); i+=2){
							float b0 = Rt[j][i];
							float b1 = Rt[j][i+1];
							b[i] = b0;
							b[i+1] = b1;
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float b0 = Rt[j][i];
							b[i] = b0;
						}								
					} else {
						for(i = start_column; (i+1<N_COL)&&(i+1<start_column+blockSize_column); i+=2){
							float b0 = Rt[j][i];
							float b1 = Rt[j][i+1];
							b[i] = b0;
							b[i+1] = b1;
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float b0 = Rt[j][i];
							b[i] = b0;
						}
					}
				}
				#else
				for(i=k; i+1<N_COL; i+=2){
					float b0 = Rt[j][i];
					float b1 = Rt[j][i+1];					
					b[i] = b0;
					b[i+1] = b1;
				}
				if(i < N_COL){
					float b0 = Rt[j][i];
					b[i] = b0;
				}

				// printf("b = \n");
				// for(int i = 0 ; i < 10; i++){
				// 	printf("%.20f\t",b[i]);
				// }
				// printf("\n");

				#endif

				pi_cl_team_barrier();

				#if NUM_CORES > 1
				if(start_column + blockSize_column >= k){
					if(start_column < k){
						for(i = k; (i+1<N_COL)&&(i+1<start_column + blockSize_column); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[k][i]     = (c * a0) + (s * b0);
							Rt[k][i+1]     = (c * a1) + (s * b1);
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float a0 = a[i];float b0 = b[i];
							Rt[k][i]     = (c * a0) + (s * b0);
						}									
					} else {
						for(i = start_column; (i+1<N_COL)&&(i+1<start_column + blockSize_column); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[k][i]     = (c * a0) + (s * b0);
							Rt[k][i+1]     = (c * a1) + (s * b1);
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float a0 = a[i];float b0 = b[i];
							Rt[k][i]     = (c * a0) + (s * b0);
						}
					}
				}
				#else
				for(i=k; i+1<N_COL; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];
					Rt[k][i]     = (c * a0) + (s * b0);
					Rt[k][i+1]     = (c * a1) + (s * b1);
				}
				if(i < N_COL){
					float a0 = a[i];float b0 = b[i];
					Rt[k][i]     = (c * a0) + (s * b0);
				}
				#endif
				
				#if NUM_CORES > 1
				if(start_column + blockSize_column >= k){
					if(start_column < k){
						for(i = k; (i+1<N_COL)&&(i+1<start_column+blockSize_column); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[j][i]     = (c * b0) - (s * a0);
							Rt[j][i+1]     = (c * b1) - (s * a1);
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float a0 = a[i];float b0 = b[i];
							Rt[j][i]     = (c * b0) - (s * a0);
						}										
					} else {
						for(i = start_column; (i+1<N_COL)&&(i+1<start_column+blockSize_column); i+=2){
							float a0 = a[i];float b0 = b[i];
							float a1 = a[i+1];float b1 = b[i+1];
							Rt[j][i] = (c * b0) - (s * a0);
							Rt[j][i+1] = (c * b1) - (s * a1);
						}
						if((i < N_COL)&&(i < start_column + blockSize_column)){
							float a0 = a[i];float b0 = b[i];
							Rt[j][i] = (c * b0) - (s * a0);
						}
					}
				}
				#else
				for(i=k; i+1<N_COL; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];
					Rt[j][i]     = (c * b0) - (s * a0);
					Rt[j][i+1]     = (c * b1) - (s * a1);
				}
				if(i < N_COL){
					float a0 = a[i];float b0 = b[i];
					Rt[j][i]     = (c * b0) - (s * a0);
				}
				#endif

				pi_cl_team_barrier();

				#if NUM_CORES > 1
				if(start_row<N_ROW){
					for(i = start_row; (i+1<N_ROW)&&(i+1<start_row + blockSize_row); i+=2){
						float a0 = Qt[i][k];
						float a1 = Qt[i+1][k];
						a[i] = a0;
						a[i+1] = a1;
					}
					if((i < N_ROW)&&(i < start_row + blockSize_row)){
						float a0 = Qt[i][k];
						a[i] = a0;
					}
				}
				#else
				for(i=0; i+1<N_ROW; i+=2){
					float a0 = Qt[i][k];
					float a1 = Qt[i+1][k];
					a[i] = a0;
					a[i+1] = a1;
					}
				if(i < N_ROW){
					float a0 = Qt[i][k];
					a[i] = a0;
				}
				#endif

				#if NUM_CORES > 1
				if(start_row<N_ROW){
					for(i = start_row; (i+1<N_ROW)&&(i+1<start_row + blockSize_row); i+=2){
						float b0 = Qt[i][j];
						float b1 = Qt[i+1][j];
						b[i] = b0;
						b[i+1] = b1;
					}
					if((i < N_ROW)&&(i < start_row + blockSize_row)){
						float b0 = Qt[i][j];
						b[i] = b0;
					}
				}
				#else	
				for(i=0; i+1<N_ROW; i+=2){
					float b0 = Qt[i][j];
					float b1 = Qt[i+1][j];
					b[i] = b0;
					b[i+1] = b1;
				}
				if(i < N_ROW){
					float b0 = Qt[i][j];
					b[i] = b0;
				}
				#endif

				pi_cl_team_barrier();				

				#if NUM_CORES > 1
				if(start_row<N_ROW){
					for(i = start_row; (i+1<N_ROW)&&(i+1<start_row + blockSize_row); i+=2){
						float a0 = a[i];float b0 = b[i];
						float a1 = a[i+1];float b1 = b[i+1];						
						Qt[i][k] = (c * a0)+(s * b0);
						Qt[i+1][k] = (c * a1)+(s * b1);
					}
					if((i < N_ROW)&&(i < start_row + blockSize_row)){
						float a0 = a[i];float b0 = b[i];
						Qt[i][k] = (c * a0)+(s * b0);
					}
				}
				#else
				for(i=0; i+1<N_ROW; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];						
					Qt[i][k] = (c * a0)+(s * b0);
					Qt[i+1][k] = (c * a1)+(s * b1);
				}
				if(i < N_ROW){
					float a0 = a[i];float b0 = b[i];
					Qt[i][k] = (c * a0)+(s * b0);
				}
				#endif

				#if NUM_CORES > 1
				if(start_row<N_ROW){
					for(i = start_row; (i+1<N_ROW)&&(i+1<start_row + blockSize_row); i+=2){
						float a0 = a[i];float b0 = b[i];
						float a1 = a[i+1];float b1 = b[i+1];						
						Qt[i][j] = (c * b0)-(s * a0);
						Qt[i+1][j] = (c * b1)-(s * a1);					
					}
					if((i < N_ROW)&&(i < start_row + blockSize_row)){
						float a0 = a[i];float b0 = b[i];
						Qt[i][j] = (c * b0)-(s * a0);
					}
				}
				#else
				for(i=0; i+1<N_ROW; i+=2){
					float a0 = a[i];float b0 = b[i];
					float a1 = a[i+1];float b1 = b[i+1];						
					Qt[i][j] = (c * b0)-(s * a0);
					Qt[i+1][j] = (c * b1)-(s * a1);					
					}
				if(i < N_ROW){
					float a0 = a[i];float b0 = b[i];
					Qt[i][j] = (c * b0)-(s * a0);
				}
				#endif
			}
        }
    }
    return;
}
