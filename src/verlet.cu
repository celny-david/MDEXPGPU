/*
 * verlet.cu
 *
 *  Created on: Dec 19, 2019
 *      Author: dixiw
 */
#include <math.h>

#include "err_check.h"
#include "value.h"
#include "cuda_helper.h"
#include "cuda_variables.h"
#include "cuda_definition.h"
#include "system_property.h"

/* ===== ===== =================== ===== ===== */
/* ===== ===== DEVICE CODE SECTION ===== ===== */
/* ===== ===== =================== ===== ===== */

/* Function: verlet_reconstruct
 * procedure for re/constructon of the Verlet list
 *
 * 2D block implementation,
 * exposes more parallelism at the cost of synchronisation
 *
 * Grid_requirements:
 * 2D grid in dimension of threadblock in x-dimension and n_mol/thread_per_block in y-dim
 *
 * Parameters:
 * d_position_x 	  - pointer to x position array
 * d_position_y 	  - pointer to y position array
 * d_position_z 	  - pointer to z position array
 * d_verlet			  - pointer to the Verlet list array
 * d_verlet_occupancy - pointer to the Verlet list occupancy array
 */
__global__ void verlet_continual_reconstruct(value_t* d_position_x, value_t* d_position_y,value_t* d_position_z,
											 ushort2* d_verlet, unsigned int* d_verlet_index)
{
	extern __shared__ value3_t j_positions[];
	value3_t i_position;
	value_t dx,dy,dz;
	uint i;
	uint row_id = blockIdx.x * blockDim.x + threadIdx.x;
	uint coll_id = blockIdx.y * blockDim.x + threadIdx.x;// each thread load one y element into shared memory

	i_position.x = d_position_x[row_id];
	i_position.y = d_position_y[row_id];
	i_position.z = d_position_z[row_id];

	// load the data into shared memory
	j_positions[threadIdx.x].x = d_position_x[coll_id];
	j_positions[threadIdx.x].y = d_position_y[coll_id];
	j_positions[threadIdx.x].z = d_position_z[coll_id];
	__syncthreads(); // can be a memory thread fence
	// __threadfence();

	for (i = 0; i < blockDim.x; i++)
	{// test if add idx to the verlet list

		dx = i_position.x - j_positions[i].x;
		dx -= d_lx*round(dx/d_lx);
		dy = i_position.y - j_positions[i].y;
		dy -= d_ly*round(dy/d_ly);
		dz = i_position.z - j_positions[i].z;
		dz -= d_lz*round(dz/d_lz);
		if( (dx*dx+dy*dy+dz*dz < d_cutoff_verlet*d_cutoff_verlet) && (dx*dx+dy*dy+dz*dz > 0.0) )
		{
			// on the actual position save the index
			// !idx does not correspond to correct column as it is valid only for loading into shared memory
			// ! because tile*blockDim.x <= idx <= (tile+1)*blockDim.x
			// this problem exhibit data-race situation -> data locking
			// because atomicInc perform controll ((old >= VERLET_WIDTH+1) ? 0 : (old+1))
			// old BUG should fist increment occupancy and then save the index on the old position which sould be empty
			// create the pair and save it to memory
			ushort2 pair;
			pair.x = row_id;
			pair.y = blockIdx.y * blockDim.x + i;
			// required to obtain the correct position in the list:
			// the proper position in verlet list is calculated as
			// skip for i-th row in d_verlet_index[row_id] 
			// increment of already evaluated pairs -> if index is atomic ncremented then if fulfill that
			// unsigned int ver_pos = atomicInc(&d_verlet_index[row_id],UINT_MAX); // BEWARE atomicInc return old value
			unsigned int ver_pos = atomicAdd(&d_verlet_index[row_id],1); // BEWARE atomicInc return old value
			d_verlet[ver_pos] = pair;
		}
	}
}
//TODO reimplement for the dynamic verlet
__global__ void verlet_reconstruct_exp(value_t* d_position_x, value_t* d_position_y, value_t* d_position_z,
		                           		ushort2* d_verlet, unsigned int* d_verlet_index, 
		                           		value_t*d_lx_buffer, unsigned long ii)
{
	extern __shared__ value3_t j_positions[];
	value3_t i_position;
	value_t dx,dy,dz;
	int i;
	int row_id = blockIdx.x * blockDim.x + threadIdx.x;
	int coll_id = blockIdx.y * blockDim.x + threadIdx.x;// each thread load one y element into shared memory

	i_position.x = d_position_x[row_id];
	i_position.y = d_position_y[row_id];
	i_position.z = d_position_z[row_id];

	// load the data into shared memory
	j_positions[threadIdx.x].x = d_position_x[coll_id];
	j_positions[threadIdx.x].y = d_position_y[coll_id];
	j_positions[threadIdx.x].z = d_position_z[coll_id];
	__syncthreads();

	for (i = 0; i < blockDim.x; i++)
	{// test if add idx to the verlet list

		dx = i_position.x - j_positions[i].x;
		dx -= d_lx_buffer[2*ii]*round(dx/ d_lx_buffer[2*ii]);
		dy = i_position.y - j_positions[i].y;
		dy -=  d_lx_buffer[2*ii]*round(dy/ d_lx_buffer[2*ii]);
		dz = i_position.z - j_positions[i].z;
		dz -=  d_lx_buffer[2*ii]*round(dz/ d_lx_buffer[2*ii]);
		//__threadfence();
		if( (dx*dx+dy*dy+dz*dz < d_cutoff_verlet*d_cutoff_verlet) && (dx*dx+dy*dy+dz*dz > 0.0) )
		{
			// on the actual position save the index
			// !idx does not correspond to correct column as it is valid only for loading into shared memory
			// ! because tile*blockDim.x <= idx <= (tile+1)*blockDim.x
			// this problem exhibit data-race situation -> data locking
			// because atomicInc perform controll ((old >= VERLET_WIDTH+1) ? 0 : (old+1))
			// old BUG should fist increment occupancy and then save the index on the old position which sould be empty
			// create the pair and save it to memory
			ushort2 pair;
			pair.x = row_id;
			pair.y = blockIdx.y * blockDim.x + i;
			// required to obtain the correct position in the list:
			// the proper position in verlet list is calculated as
			// skip for i-th row in d_verlet_index[row_id] 
			// increment of already evaluated pairs -> if index is atomic ncremented then if fulfill that
			// unsigned int ver_pos = atomicInc(&d_verlet_index[row_id],UINT_MAX); // BEWARE atomicInc return old value
			unsigned int ver_pos = atomicAdd(&d_verlet_index[row_id],1); // BEWARE atomicInc return old value
			d_verlet[ver_pos] = pair;
		}
	}
}

/*
 * the idea of the two pass verlet reconstruct algorithm is to first count the neighbors
 * and then allocate appropriate pitched space
 * the d_verlet_occupancy is then used for access in cumulative sum fashion
 * NOTE verlet_occupancy needs to be +1 larger to start from zero for ease of access
 * TODO consider forward calculation of cumulative sum here
 * BUG the increment support only int therefore no ushort for d_verlet_occupancy
 */
__global__ void verlet_dynamic_count(value_t* d_position_x, value_t* d_position_y, value_t* d_position_z,
	   	   	   	   	   	   	   	     unsigned int* d_verlet_occupancy)
{
	extern __shared__ value3_t j_positions[];
	value3_t i_position;
	value_t dx,dy,dz;
	int i;
	int row_id = blockIdx.x * blockDim.x + threadIdx.x;
	int coll_id = blockIdx.y * blockDim.x + threadIdx.x;// each thread load one y element into shared memory

	i_position.x = d_position_x[row_id];
	i_position.y = d_position_y[row_id];
	i_position.z = d_position_z[row_id];
	// DEBUG - need bigger printf buffer
	// printf("row_id: %i coll_id: %i i_pos: ( %f, %f, %f)\n",row_id, coll_id, i_position.x, i_position.y, i_position.z );

	// load the data into shared memory
	j_positions[threadIdx.x].x = d_position_x[coll_id];
	j_positions[threadIdx.x].y = d_position_y[coll_id];
	j_positions[threadIdx.x].z = d_position_z[coll_id];
	// DEBUG - need bigger printf buffer
	// printf("row_id: %i coll_id: %i i_pos: ( %f, %f, %f)\n",row_id, coll_id, j_positions[threadIdx.x].x, j_positions[threadIdx.x].y, j_positions[threadIdx.x].z );
	// __threadfence();
	__syncthreads();
	for (i = 0; i < blockDim.x; i++)
	{// test if add idx to the verlet list

		dx = i_position.x - j_positions[i].x;
		dx -= d_lx*round(dx/d_lx);
		dy = i_position.y - j_positions[i].y;
		dy -= d_ly*round(dy/d_ly);
		dz = i_position.z - j_positions[i].z;
		dz -= d_lz*round(dz/d_lz);

		if( (dx*dx+dy*dy+dz*dz < d_cutoff_verlet*d_cutoff_verlet) && (dx*dx+dy*dy+dz*dz > 0.0) )
		{
			// atomicInc(&d_verlet_occupancy[row_id],UINT_MAX);
			atomicAdd(&d_verlet_occupancy[row_id],1);
		}
	}
}

/* ===== ===== ================= ===== ===== */
/* ===== ===== HOST CODE SECTION ===== ===== */
/* ===== ===== ================= ===== ===== */

/* Function: verlet_const_memory_copy
 * procedure to copy verlet data to the const memory device
 *
 * Notes:
 * > wrapper intended to be callable from copy to symbol overlay method
 * >! emits and catch cudaError bound with improper alloc
 *
 */
void verlet_const_memory_copy()
{
	cudaSafeCall( cudaMemcpyToSymbol(d_cutoff_verlet,&sysp.cutoff_verlet,sizeof(value_t)) );
}

/* Function: verlet_device_alloc
 * procedure to allocate verlet data on the device
 *
 * Notes:
 * >! emits and catch cudaError bound with improper alloc
 * > wrapper intended to be callable from interface method
 *
 */
void verlet_device_alloc(ushort verlet_size)
{
	cudaSafeCall( cudaMalloc( (void**)&d_verlet, verlet_size * sizeof(ushort2)) );
	cudaSafeCall( cudaMalloc( (void**)&d_verlet_occupancy, sysp.molecule_count_aligned * sizeof(unsigned int)) );
	cudaSafeCall( cudaMalloc( (void**)&d_verlet_index, sysp.molecule_count_aligned * sizeof(unsigned int)) );
}


/* Function: verlet_device_free
 * procedure to deallocate verlet data on the device
 *
 * Notes:
 * >! emits and catch cudaError bound with improper free
 * > wrapper intended to be callable from interface method
 *
 */
void verlet_device_free()
{
	cudaSafeCall( cudaFree(d_verlet) );
	cudaSafeCall( cudaFree(d_verlet_occupancy) );
	cudaSafeCall( cudaFree(d_verlet_index) );
}

/* Function: verlet_continual_accumulate
 * accumulate the prefix and total sum of input array into the memory
 *
 * Notes:
 * >! this is host code - consider writing at gpu in case of bottleneck
 *
 */
void verlet_continual_accumulate(unsigned int *array_in, unsigned int length, unsigned int *prefix_sum, unsigned int *total_sum)
{
 	unsigned int curr_count;

 	curr_count = 0;
	for (int i = 0; i < length; ++i)
	{
		prefix_sum[i] = curr_count; // first index is zero
		curr_count += array_in[i];		
	}
	total_sum[0] = curr_count;			
}

int once = 10; // DEBUG- how many lines of verlet refresh sizes should be printed

/* Function: verlet_refresh
 * manage the refresh of the verlet list
 *
 * Notes:
 * >! emits and catch cudaError bound with improper kernell run
 * > wrapper intended to be callable from iterative schema
 * > TODO reallocates the size only when needed
 * > TODO determine optimal expansion factor so it is not reallocated every reconstruct
 *
 */
unsigned int verlet_refresh(value_t *d_position_x, value_t *d_position_y, value_t *d_position_z,
							value_t *d_lx_buffer, unsigned long ii)
{
	// const float expansion_factor = 1.0;
	unsigned int verlet_size;
	unsigned int verlet_size_new;
	unsigned int *verlet_occupancy = (unsigned int*) malloc(sysp.molecule_count_aligned * sizeof(unsigned int));
	unsigned int *verlet_index = (unsigned int*) malloc(sysp.molecule_count_aligned * sizeof(unsigned int));
	unsigned const int n_block_construct = sysp.molecule_count_aligned/THREAD_P_BLOCK_RECONSTRUCT;	
	dim3 block_dim_construct(n_block_construct,n_block_construct);
	
	// cudaSafeCall( cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 10*1024*1024));
	// writeout_3d<<<sysp.molecule_count_aligned/64,64>>>(d_position_x,d_position_y,d_position_z,1);							
	// cudaDeviceSynchronize();
	// printf("n_block_construct: %i, THREAD_P_BLOCK_RECONSTRUCT: %i\n", n_block_construct, THREAD_P_BLOCK_RECONSTRUCT);

	// cudaDeviceSynchronize();
	set_zero <<< n_block_construct,
				 THREAD_P_BLOCK_RECONSTRUCT >>> (d_verlet_occupancy);
	// cudaSafeKernell();

	verlet_dynamic_count <<<block_dim_construct,
							THREAD_P_BLOCK_RECONSTRUCT,
							THREAD_P_BLOCK_RECONSTRUCT * sizeof(value3_t)>>>
							(d_position_x, d_position_y, d_position_z, d_verlet_occupancy);	
	cudaSafeKernell();
	
	// TODO the continual count does not fill occupancy - index cant be constructed and verlet size cant be constructed from ocupancy[end] + index[end] as index is exclusive
	cudaSafeCall( cudaMemcpy(verlet_occupancy, d_verlet_occupancy,sizeof(unsigned int)*sysp.molecule_count_aligned,cudaMemcpyDeviceToHost)); // move the size device
	verlet_continual_accumulate(verlet_occupancy, sysp.molecule_count_aligned, verlet_index, &verlet_size_new);	
	
	// if (verlet_size <= verlet_size_new) { // if verlet size is not enough make it larger
	if (once)
	{
		printf("DEBUG ver_size: %i, ver_size_new: %i, ", verlet_size,
													     verlet_size_new);
	}

	// padding of the verlet 											
	verlet_size = ((verlet_size_new + THREAD_P_BLOCK_RECONSTRUCT - 1)/ THREAD_P_BLOCK_RECONSTRUCT)*THREAD_P_BLOCK_RECONSTRUCT;

	if (once)
	{
		printf("ver_size_new_align: %i, Int/mol: %f\n",  verlet_size,
									    			(1.0*verlet_size)/sysp.molecule_count_aligned);
		once -= 1;
	}
	cudaSafeCall( cudaFree(d_verlet)); // clear the verlet list TODO not efficient
	
	if (verlet_size == 0) // Handle the empty neighbour list
	{
		return verlet_size;
	}

	// allocate the verlet again
	// BUG allocate does not fill the array with zeros
	cudaSafeCall( cudaMalloc( (void**)&d_verlet, verlet_size*sizeof(ushort2) )); 
	// TODO rework reconstruct so that it fills the leftover size with zeros
	// TODO consider padding of individual molecule interaction for comapcted work over same data block -> workloading
	dim3 block_dim_set(verlet_size/THREAD_P_BLOCK_RECONSTRUCT);	
	
	set_zero <<<block_dim_set,
				THREAD_P_BLOCK_RECONSTRUCT>>> (d_verlet); // DEBUG zero out the verlet list
	cudaSafeKernell();
	// }

	cudaSafeCall( cudaMemcpy(d_verlet_index, verlet_index,sizeof(unsigned int)*sysp.molecule_count_aligned,cudaMemcpyHostToDevice)); // move the size device
	
#ifdef EXPANSION
	verlet_reconstruct_exp<<<block_dim,THREAD_P_BLOCK_RECONSTRUCT,THREAD_P_BLOCK_RECONSTRUCT*sizeof(value3_t)>>>(d_position_x, d_position_y, d_position_z,
																							       	   	   	   	     d_verlet, d_verlet_occupancy,d_lx_buffer,ii);
#else
	verlet_continual_reconstruct<<< block_dim_construct,
									THREAD_P_BLOCK_RECONSTRUCT,
									THREAD_P_BLOCK_RECONSTRUCT*sizeof(value3_t)>>>
									(d_position_x, d_position_y, d_position_z,
									 d_verlet, d_verlet_index);
#endif
	cudaSafeKernell();
	return verlet_size;
}
