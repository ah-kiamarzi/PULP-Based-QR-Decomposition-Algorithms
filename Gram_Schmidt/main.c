#include "pmsis.h"
#include "util.h"
#include "inputData/T/qr40X40.h"
#include "inputData/T/data40X40.inc"
#include <math.h>

// Q is transposed
PI_L1 float Q[N_ROW][N_ROW];
PI_L1 float R[N_ROW][N_COL];

PI_L1 static float buffer[NUM_CORES];
PI_L1 float rk = 0.0;
PI_L1 static float temp[NUM_CORES];
PI_L1 static float one = 1.0;

float res[N_ROW][N_COL];
int numBarr = 0;
definePrefArr

void cluster_main();

void qr_gramSmidt(float Q[][N_ROW], float R[][N_COL], float input[][N_ROW]);

void pe_entry(void *arg)
{
	cluster_main();
}

void cluster_entry(void *arg)
{
	pi_cl_team_fork((NUM_CORES), pe_entry, 0);
}

static int test_entry()
{

	struct pi_device cluster_dev;
	struct pi_cluster_conf cl_conf;
	struct pi_cluster_task cl_task;

	pi_cluster_conf_init(&cl_conf);
	pi_open_from_conf(&cluster_dev, &cl_conf);
	if (pi_cluster_open(&cluster_dev))
	{
		return -1;
	}

	pi_cluster_task(&cl_task, cluster_entry, NULL);
	cl_task.stack_size = STACK_SIZE;
	cl_task.slave_stack_size = STACK_SIZE;

	pi_cluster_send_task_to_cl(&cluster_dev, &cl_task);

	pi_cluster_close(&cluster_dev);

	return 0;
}

static void test_kickoff(void *arg)
{
	int ret = test_entry();
	pmsis_exit(ret);
}

int main()
{
	return pmsis_kickoff((void *)test_kickoff);
}

void cluster_main()
{

	pi_perf_conf(
		(1 << PI_PERF_CYCLES) |
		(1 << PI_PERF_INSTR) |
		(1 << PI_PERF_ACTIVE_CYCLES) |
		(1 << PI_PERF_LD_EXT) |
		(1 << PI_PERF_TCDM_CONT) |
		(1 << PI_PERF_LD_STALL) |
		(1 << PI_PERF_IMISS) |
		(1 << 0x11) |
		(1 << 0x12) |
		(1 << 0x13) |
		(1 << 0x14));

	if (pi_core_id() == 0)
	{
		initialMatrixConst(&Q[0][0], N_ROW, N_ROW, 0);
		initialMatrixConst(&R[0][0], N_ROW, N_COL, 0);
	}

	pi_cl_team_barrier();

	pi_perf_reset();
	pi_perf_start();

	qr_gramSmidt(Q, R, input);

	pi_perf_stop();

	printPerfCounters

	if(pi_core_id() == 0){
		printf("NUMBARR = %d\n",numBarr);
	}

	pi_cl_team_barrier();

#ifdef DEBUG
	if (pi_core_id() == 0)
	{

		printMatrix(&Q[0][0],N_ROW,N_ROW,"Q",1);
		printMatrix(&R[0][0],N_ROW,N_COL,"R",0);

		for (int i = 0; i < N_ROW; i++)
		{
			for (int j = 0; j < N_COL; j++)
			{
				res[i][j] = 0;
				for (int k = 0; k < N_ROW; k++)
				{
					res[i][j] += Q[k][i] * R[k][j];
				}
			}
		}
		printMatrix(&res[0][0],N_ROW,N_COL,"Res",0);
	}

	pi_cl_team_barrier();
#endif
}

inline float Sqrt(float x)
{
	float res;
	asm("fsqrt.s %0, %1"
		: "=f"(res)
		: "f"(x));
	return res;
}

#pragma GCC push_options
#pragma GCC optimize("-O3")
__attribute__((noinline)) void norm(float *v, int row, int column, int k, int blockSize_ROW, int start_ROW)
{
	int i2 = 0;
	float n = 0.0f;
	int j, idx;

	#if NUM_CORES > 2
	int core_id = pi_core_id();

	int start_matrice = start_ROW;

	buffer[core_id] = 0;
	idx = 0;

	for (j = start_ROW; (j + 1) < start_ROW + (blockSize_ROW); j += 2){
		float t0 = v[idx + start_matrice];
		float t1 = v[idx + 1 + start_matrice];
		buffer[core_id] = buffer[core_id] + t0 * t0;
		buffer[core_id] = buffer[core_id] + t1 * t1;
		idx += 2;
	}

	if (j < start_ROW + blockSize_ROW){
		buffer[core_id] = buffer[core_id] + v[idx + start_matrice] * v[idx + start_matrice];
	}

	pi_cl_team_barrier();
	BarrierCounter

	if (core_id == 0){
		for (j = 0; j < NUM_CORES; j++){
			n += buffer[j];
		}
		R[k][k] = Sqrt(n);
		rk = 1 / R[k][k];
	}
	pi_cl_team_barrier();
	BarrierCounter

	#else

	for (j = 0; (j + 1) < column; j += 2){
		float t0 = v[i2];
		float t1 = v[i2 + 1];
		float temp = t0 * t0;
		temp = temp + t1 * t1;
		n = n + temp;
		i2 += 2 * 1;
	}

	if (j < column){
		float t0 = v[i2];
		float temp = t0 * t0;
		n = n + temp;
		i2 += 1 * 1;
	}

	R[k][k] = Sqrt(n);
	rk = one / R[k][k];
	#if NUM_CORES > 2
	pi_cl_team_barrier();
	BarrierCounter
	#endif

	#endif
}

void qr_gramSmidt(float Q[][N_ROW], float R[][N_COL], float input[][N_ROW]){
	int j;
	float RTEMP = 0;

#if NUM_CORES > 1
	int core_id = pi_core_id();

	int blockSize_ROW = N_ROW / NUM_CORES;
	int start_ROW = core_id * blockSize_ROW;
	if (core_id == (NUM_CORES - 1))
	{
		blockSize_ROW = N_ROW - (NUM_CORES - 1) * blockSize_ROW;
	}
	int end_ROW = start_ROW + blockSize_ROW;

#endif

	for (int k = 0; k < N_COL; k++){
		#if NUM_CORES > 1
		for (j = start_ROW; ((j + 1) < end_ROW); j += 2){
			float in0 = input[k][j];
			float in1 = input[k][j + 1]; 
			Q[k][j] = in0;
			Q[k][j + 1] = in1;
		}
		if (j < end_ROW){
			float in0 = input[k][j];
			Q[k][j] = in0;
		}
		#else
		for (j = 0; (j + 1) < N_ROW; j += 2){
			float in0 = input[k][j];
			float in1 = input[k][j + 1];
			Q[k][j] = in0;
			Q[k][j + 1] = in1;
		}
		if (j < N_ROW){
			float in0 = input[k][j];
			Q[k][j] = in0;
		}
		#endif

		#if NUM_CORES > 1
				pi_cl_team_barrier();
				BarrierCounter
		#endif

		for (int i = 0; i < k; i++){
			#if NUM_CORES > 1
			temp[core_id] = 0;
			for (j = start_ROW; ((j + 1) < end_ROW); j += 2){
				float ji0 = Q[i][j];
				float jk0 = Q[k][j];
				float ji1 = Q[i][j + 1];
				float jk1 = Q[k][j + 1];
				temp[core_id] = temp[core_id] + (ji0 * jk0) + (ji1 * jk1);
			}
			if ((j < end_ROW)){
				float ji0 = Q[i][j];
				float jk0 = Q[k][j];
				temp[core_id] = temp[core_id] + (ji0 * jk0);
			}

			pi_cl_team_barrier();
			BarrierCounter

			RTEMP = temp[0] + temp[1];
			for (j = 2; j < NUM_CORES; j++){
				RTEMP += temp[j];
			}

			R[i][k] = RTEMP;

			#else
			for (j = 0; (j + 1) < N_ROW; j += 2){
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j + 1];
				float Qjk1 = Q[k][j + 1];
				float Rik = R[i][k];
				R[i][k] = Rik + (Qji0 * Qjk0) + (Qji1 * Qjk1);
			}
			if (j < N_ROW){
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Rik = R[i][k];
				float temp = (Qji0 * Qjk0);
				R[i][k] = Rik + temp;
			}
			#endif

			#if NUM_CORES > 1
			for (j = start_ROW; (j + 1 < end_ROW); j += 2){
				float Rik = RTEMP;
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j + 1];
				float Qjk1 = Q[k][j + 1];

				Q[k][j] = Qjk0 - Rik * Qji0;
				Q[k][j + 1] = Qjk1 - Rik * Qji1;
			}
			if (j < end_ROW){
				float Rik = RTEMP;
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik * Qji0;
			}
			#else
			for (j = 0; j + 1 < N_ROW; j += 2){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j + 1];
				float Qjk1 = Q[k][j + 1];

				Q[k][j] = Qjk0 - Rik * Qji0;
				Q[k][j + 1] = Qjk1 - Rik * Qji1;
			}
			if (j < N_ROW){
				float Rik = R[i][k];
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik * Qji0;
			}
			#endif
		}

		#if NUM_CORES > 2
		norm(&Q[k][0], N_ROW, N_ROW, k, blockSize_ROW, start_ROW);
		#else
		norm(&Q[k][0], N_ROW, N_ROW, k, 0, 0);
		#endif

		#if NUM_CORES > 1
		for (j = start_ROW; (j + 1) < start_ROW + (blockSize_ROW & 0xFFFFFFFE); j += 2){
			float Qjk0 = Q[k][j];
			float Qjk1 = Q[k][j + 1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j + 1] = Qjk1 * rk;
		}
		if (j < start_ROW + (blockSize_ROW & 0xFFFFFFFE)){
			float Qjk0 = Q[k][j];
			Q[k][j] = Qjk0 * rk;
		}
		#else
		for (j = 0; (j + 1) < N_ROW; j += 2){
			float Qjk0 = Q[k][j];
			float Qjk1 = Q[k][j + 1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j + 1] = Qjk1 * rk;
		}
		if (j < N_ROW){
			float Qjk0 = Q[k][j];
			Q[k][j] = Qjk0 * rk;
		}
		#endif

		#if NUM_CORES > 1
		if (blockSize_ROW & 0x0001){
			Q[k][j] = Q[k][j] * rk;
		}
		#endif
	}
	return;
}
