/*
 * cuda_helper.cu
 *
 *  Created on: 13. Jan 2020
 *      Author: dixiw
 */

#include <stdio.h>
#include "value.h"
//#include "cuda_definition.h" // BUG why this include is not required when WORKLOAD is used here
#include "cuda_variables.h"

// ===== HELPER related functions =====
// === PRIVATE FOR THIS FILE - NOT IN HEADER ===

/* Function: set_zero
 * simple zeroing of given array in device memory (value-t variant)
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * set_array - pointer to array that should be cleared
 */
__global__ void set_zero(value_t *set_array)
{
	 set_array[blockIdx.x*blockDim.x+ threadIdx.x] = 0.0;
	 return;
}

/* Function: set_zero
 * simple zeroing of given array in device memory (int variant)
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * set_array - pointer to array that should be cleared
 */
__global__ void set_zero(int *set_array)
{
	 set_array[blockIdx.x*blockDim.x+ threadIdx.x] = 0;
	 return;
}

/* Function: set_zero
 * simple zeroing of given array in device memory unsigned int variant)
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * set_array - pointer to array that should be cleared
 */
__global__ void set_zero(unsigned int *set_array)
{
	 set_array[blockIdx.x*blockDim.x+ threadIdx.x] =(unsigned int) 0;
	 return;
}

/* Function: set_zero
 * simple zeroing of given array in device memory ushort2 variant)
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * set_array - pointer to array that should be cleared
 */
__global__ void set_zero(ushort2 *set_array)
{
	 set_array[blockIdx.x*blockDim.x+ threadIdx.x].x = (ushort) 0;
	 set_array[blockIdx.x*blockDim.x+ threadIdx.x].y = (ushort) 0;
	 return;
}

/* Function: set_false
 * simple setting given array to false in device memory
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * set_array - pointer to array that should be set to false
 */
__global__ void set_false(bool *set_array)
{
	 set_array[blockIdx.x*blockDim.x+ threadIdx.x] = false;
	 return;
}

/* Function: ceil_to_wokload
 * ceiling variables to WORKLOAD multiple if not already
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * set_array - pointer to array that should be cleared
 */
__global__ void ceil_to(unsigned int *set_array, unsigned int ceil_target)
{
	unsigned int tid = blockIdx.x*blockDim.x+ threadIdx.x;

	set_array[tid] = ((set_array[tid] + ceil_target - 1)/ceil_target)*ceil_target;
	return;
}

// === PUBLIC - IN HEADER ===

/* Function: writeout_1d
 * utility function for debug value writing on GPU (version with unsigned int)
 *
 * strided to reduce the amount of output data
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * write_variable - pointer to array that should be set written
 * mod_count	  - the size of the stride (for full use 1)
 */
__global__ void writeout_1d(unsigned int *write_variable, int mod_count)
{
//  simple output to console for unsigned int
//  mod count is modulo index which is written out - for complete set mod_count =1
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);


	if (idx%mod_count>0 ) return;

	printf("B:%02i T:%03i ID:%05i: %i \n", blockIdx.x, threadIdx.x, idx, write_variable[idx]);
}

__global__ void writeout_1d(ushort2 *write_variable, int mod_count)
{
//  simple output to console for unsigned int
//  mod count is modulo index which is written out - for complete set mod_count =1
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);


	if ( idx%mod_count>0 ) return;

	printf("B:%02i T:%03i ID:%05i: (%i, %i) \n", blockIdx.x, threadIdx.x, idx, write_variable[idx].x, write_variable[idx].y);

}

/* Function: writeout_1d
 * utility function for debug value writing on GPU (version with value_t)
 *
 * strided to reduce the amount of output data
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * write_variable - pointer to array that should be set written
 * mod_count	  - the size of the stride (for full use 1)
 */
__global__ void writeout_1d(value_t *write_variable, int mod_count)
{
/*
 * simple output to console
 * mod count is modulo index which is written out - for complete set mod_count =1
 */
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);


	if (idx%mod_count>0 ) return;

	//printf("B:%i T:%i ID:%i: %f \n", blockIdx.x, threadIdx.x, idx, write_variable[idx]);
	printf("B,%i,T,%i,ID,%i,%f\n", blockIdx.x, threadIdx.x, idx, write_variable[idx]);
}

/* Function: writeout_clock
 * utility function for debug write of execution time on GPU
 *
 * the time is calculated as (float)write_clock[idx]/clockRate/WORKLOAD)
 *
 * strided to reduce the amount of output data
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * write_clock - pointer to array that should be set written
 * clockRate   - the clock rate of GPU for correct time recalculation
 * mod_count   - the size of the stride (for full use 1)
 */
__global__ void writeout_clock(int *write_clock,int clockRate, int mod_count)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);

	if (idx%mod_count>0) return;

	printf("%i: %f ms \n", idx, (float)write_clock[idx]/clockRate/WORKLOAD);

//	for (int i = 0; i < sysp.atom_count_aligned; i +=mod_count)
//	{
//		printf("%i: %f micros\n",i, (float)(1000000*write_clock[i])/clockRate);
//	}
}

/* Function: writeout_3in1
 * utility function for debug value writing on GPU (version with value_t)
 *
 * strided to reduce the amount of output data
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 *
 * Parameter:
 * write_variable1 - pointer to array that should be set written at first position
 * write_variable2 - pointer to array that should be set written at second position
 * write_variable3 - pointer to array that should be set written at third position
 * mod_count	   - the size of the stride (for full use 1)
 */
__global__ void writeout_3in1(value_t *write_variable1, value_t *write_variable2, value_t *write_variable3, int mod_count)
 {
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);

	if (idx%mod_count>0) return; // modulo write-out not very effective
	printf("%i,%i: %f %f %f\n", blockIdx.x, threadIdx.x, write_variable1[idx], write_variable2[idx], write_variable3[idx]);
}
__device__ void writeout_1dM(value_t *write_variable, int mod_count)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);


	if (idx>=d_atom_count || idx%mod_count>0 ) return;

	//printf("B:%i T:%i ID:%i: %f \n", blockIdx.x, threadIdx.x, idx, write_variable[idx]);
	printf("B:%i T:%i ID:%i: %f\n", blockIdx.x, threadIdx.x, idx, write_variable[idx]);
}

