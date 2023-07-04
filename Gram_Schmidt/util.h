#include "pmsis.h"

#define STATS 
#define STACK_SIZE 1024
#define MOUNT 1
#define UNMOUNT 0
#define CID 0
#define DEBUG


#define BarrierCounter //if(pi_core_id() == 0){numBarr++;};
#define definePrefArr 							\
int cyclesArr [NUM_CORES];						\
int instrArr  [NUM_CORES];						\
int activeArr [NUM_CORES];						\
int ldEXTArr  [NUM_CORES];						\
int TCDMArr   [NUM_CORES];						\
int ldstallArr[NUM_CORES];						\
int imissArr  [NUM_CORES];						\
long int AVGCycles = 0;   						\
long int AVGInstr  = 0;   						\
long int AVGActive = 0;   						\
long int AVGldEXT  = 0;   						\
long int AVGTCDM   = 0;   						\
long int AVGLdstall= 0;   						\
long int AVGImiss  = 0;   						\




#define printPerfCounters										\
						  										\
																\
unsigned long cycles   = pi_perf_read (PI_PERF_CYCLES);		 	\
unsigned long instr    = pi_perf_read (PI_PERF_INSTR);		 	\
unsigned long active   = pi_perf_read (PI_PERF_ACTIVE_CYCLES);	\
unsigned long ldext    = pi_perf_read (PI_PERF_LD_EXT);		 	\
unsigned long tcdmcont = pi_perf_read (PI_PERF_TCDM_CONT);		\
unsigned long ldstall  = pi_perf_read (PI_PERF_LD_STALL);		\
unsigned long imiss    = pi_perf_read (PI_PERF_IMISS);		 	\
																\
cyclesArr [pi_core_id()] = cycles;								\
instrArr  [pi_core_id()] = instr;								\
activeArr [pi_core_id()] = active;								\
ldEXTArr  [pi_core_id()] = ldext;								\
TCDMArr   [pi_core_id()] = tcdmcont;							\
ldstallArr[pi_core_id()] = ldstall;								\
imissArr  [pi_core_id()] = imiss;								\
																\
if(pi_core_id() == 0){											\
	for (int i = 0; i < NUM_CORES; i++){						\
		AVGCycles  += cyclesArr [i];							\
		AVGInstr   += instrArr  [i];							\
		AVGActive  += activeArr [i];							\
		AVGldEXT   += ldEXTArr  [i];							\
		AVGTCDM    += TCDMArr   [i];							\
		AVGLdstall += ldstallArr[i];							\
		AVGImiss   += imissArr  [i];							\
	}															\
	printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",				\
	AVGCycles/NUM_CORES,AVGInstr/NUM_CORES,		 				\
	AVGActive/NUM_CORES,AVGldEXT/NUM_CORES,		 				\
	AVGTCDM/NUM_CORES,AVGLdstall/NUM_CORES,		 				\
	AVGImiss/NUM_CORES);				 		 				\
}		 										 				\


void initialMatrixMatrix(float *inp,int rowSize,int colSize,float *initMatrix){
	for (int i = 0; i < rowSize; i++){
		for (int j = 0; j < colSize; j++){
			inp[i*rowSize+j] = initMatrix[i*rowSize+j];
		}
	}
}
void initialMatrixConst (float *inp,int rowSize,int colSize,float initValue){
	for (int i = 0; i < rowSize; i++){
		for (int j = 0; j < colSize; j++){
			inp[i*colSize+j] = initValue;
		}
	}
}

void initialMatrixEYE (float *inp,int rowSize,int colSize){
	for (int i = 0; i < rowSize; i++){
		for (int j = 0; j < colSize; j++){
			if(i == j){
				inp[i*colSize+j] = 1;
			}else{
				inp[i*colSize+j] = 0;
			}
		}
	}
}


void printMatrix(float *inp,int rowSize,int colSize,const char name[], int transposed){
	printf("\n\n%s = \n",name);
	for (int i = 0; i < rowSize; i++){
		for (int j = 0; j < colSize; j++){
			if(transposed){
				printf("%.20f ",inp[j*rowSize+i]);
			}else{
				printf("%.20f ",inp[i*colSize+j]);
			}
		}    
		printf("\n");
	}
}

