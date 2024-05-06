/*
 * cuda_helper.h
 *
 *  Created on: 13. Jan 2020
 *      Author: dixiw
 */

#ifndef CUDA_HELPER_H_
#define CUDA_HELPER_H_


#include "value.h"
/* ===== ===== CUDA CODE SECTION ===== ===== */
__global__ void set_zero(value_t *set_array);
__global__ void set_zero(int *set_array);
__global__ void set_zero(unsigned int *set_array);
__global__ void set_zero(ushort2 *set_array);

__global__ void set_false(bool *set_array);
__global__ void ceil_to(unsigned int *set_array, unsigned int ceil_target);

__global__ void writeout_1d(unsigned int *write_variable, int mod_count);
__global__ void writeout_1d(value_t *write_variable, int mod_count);
__global__ void writeout_1d(ushort2 *write_variable, int mod_count);
__global__ void writeout_clock(int *write_clock,int clockRate, int mod_count);
__global__ void writeout_3in1(value_t *write_variable1, value_t *write_variable2, value_t *write_variable3, int mod_count);
__device__ void writeout_1dM(value_t *write_variable, int mod_count);

#endif /* CUDA_HELPER_H_ */


