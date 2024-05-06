/*
 * verlet.cu
 *
 *  Created on: Dec 19, 2019
 *      Author: dixiw
 */
#include <stdio.h>
#include <math.h>

#include "value.h"
#include "cuda_definition.h"
#include "system_property.h"


/* ===== ===== ==++=============== ===== ===== */
/* ===== ===== DEVICE CODE SECTION ===== ===== */
/* ===== ===== ++================= ===== ===== */

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
void verlet_continual_reconstruct(value_t* d_position_x, value_t* d_position_y,value_t* d_position_z,
											 ushort2** d_verlet, unsigned int* d_verlet_index,
											 unsigned int row_id, unsigned int coll_id)
{
	value3_t j_positions;
	value3_t i_position;
	value_t dx,dy,dz;

	if (row_id == coll_id)	return;

	i_position.x = d_position_x[row_id];
	i_position.y = d_position_y[row_id];
	i_position.z = d_position_z[row_id];

	// load the data into shared memory
	j_positions.x = d_position_x[coll_id];
	j_positions.y = d_position_y[coll_id];
	j_positions.z = d_position_z[coll_id];


	dx = i_position.x - j_positions.x;
	dx -= sysp.lx*round(dx/sysp.lx);
	dy = i_position.y - j_positions.y;
	dy -= sysp.ly*round(dy/sysp.ly);
	dz = i_position.z - j_positions.z;
	dz -= sysp.lz*round(dz/sysp.lz);

	if( (dx*dx+dy*dy+dz*dz < sysp.cutoff_verlet*sysp.cutoff_verlet) && (dx*dx+dy*dy+dz*dz > 0.0) )
	// if( (dx*dx+dy*dy+dz*dz > 0.0) ) // full list test
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
		// pair.y = blockIdx.y * blockDim.x + i;
		pair.y = coll_id;
		// required to obtain the correct position in the list:
		// the proper position in verlet list is calculated as
		// skip for i-th row in d_verlet_index[row_id] 
		// increment of already evaluated pairs -> if index is atomic ncremented then if fulfill that
		// unsigned int ver_pos = atomicInc(&d_verlet_index[row_id],UINT_MAX); // BEWARE atomicInc return old value
		unsigned int ver_pos = d_verlet_index[row_id]; // BEWARE atomicInc return old value
		d_verlet_index[row_id] += 1; // DC NOTE this is how the atomic increment look like
		(*d_verlet)[ver_pos] = pair;
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
void verlet_dynamic_count(value_t* d_position_x, value_t* d_position_y, value_t* d_position_z,
   	   	   	   	   	     unsigned int* d_verlet_occupancy,
   	   	   	   	   	     unsigned int row_id, unsigned int coll_id)
{
	value3_t j_positions;
	value3_t i_position;
	value_t dx,dy,dz;

	if (row_id == coll_id) return;

	i_position.x = d_position_x[row_id];
	i_position.y = d_position_y[row_id];
	i_position.z = d_position_z[row_id];

	// DEBUG - need bigger printf buffer
	// printf("row_id: %i coll_id: %i i_pos: ( %f, %f, %f)\n",row_id, coll_id, i_position.x, i_position.y, i_position.z );

	// load the data into shared memory
	j_positions.x = d_position_x[coll_id];
	j_positions.y = d_position_y[coll_id];
	j_positions.z = d_position_z[coll_id];

	// DEBUG - need bigger printf buffer
	// printf("row_id: %i coll_id: %i i_pos: ( %f, %f, %f)\n",row_id, coll_id, j_positions[threadIdx.x].x, j_positions[threadIdx.x].y, j_positions[threadIdx.x].z );

	dx = i_position.x - j_positions.x;
	dx -= sysp.lx*round(dx/sysp.lx);
	dy = i_position.y - j_positions.y;
	dy -= sysp.ly*round(dy/sysp.ly);
	dz = i_position.z - j_positions.z;
	dz -= sysp.lz*round(dz/sysp.lz);

	if( (dx*dx+dy*dy+dz*dz < sysp.cutoff_verlet*sysp.cutoff_verlet) && (dx*dx+dy*dy+dz*dz > 0.0) )
	// if( (dx*dx+dy*dy+dz*dz > 0.0) ) // full list test
	{
		d_verlet_occupancy[row_id] += 1;
	}
}

/* Function: verlet_continual_accumulate_seq
 * accumulate the prefix and total sum of input array into the memory
 *
 * Notes:
 *
 */
void verlet_continual_accumulate_seq(unsigned int *array_in, unsigned int length, unsigned int *prefix_sum, unsigned int *total_sum)
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


int once_seq = 10; // DEBUG- how many lines of verlet refresh sizes should be printed

/* Function: verlet_continual_refresh
 * manage the refresh of the dynamic verlet list
 *
 * Notes:
 * >! emits and catch cudaError bound with improper kernell run
 * > wrapper intended to be callable from iterative schema
 * > TODO reallocates the size only when needed
 * > TODO determine optimal expansion factor so it is not reallocated every reconstruct
 *
 */
unsigned int verlet_continual_refresh(value_t *position_x, value_t *position_y, value_t *position_z, unsigned int verlet_size, ushort2 **verlet)
{
	const float expansion_factor = 1.0;
	unsigned int verlet_size_new;
	unsigned int *verlet_occupancy = (unsigned int*) malloc(sysp.molecule_count_aligned * sizeof(unsigned int));
	unsigned int *verlet_index = (unsigned int*) malloc(sysp.molecule_count_aligned * sizeof(unsigned int));
	
	for (int i = 0; i < sysp.molecule_count_aligned; ++i)
	{
		verlet_occupancy[i] = 0.0;

		for (int j = 0; j < sysp.molecule_count_aligned; ++j)
		{
			verlet_dynamic_count (position_x, position_y, position_z,
								 verlet_occupancy,								 
								 i, j);
		}
	}
	
	verlet_continual_accumulate_seq(verlet_occupancy, sysp.molecule_count_aligned, verlet_index, &verlet_size_new);	

	if (once_seq)
	{
		printf("DEBUG ver_size: %i, ver_size_new: %i, ", verlet_size,
													     verlet_size_new);
	}

	// padding of the verlet 											
	verlet_size = ((verlet_size_new + THREAD_P_BLOCK_RECONSTRUCT - 1)/ THREAD_P_BLOCK_RECONSTRUCT)*THREAD_P_BLOCK_RECONSTRUCT;

	if (once_seq)
	{
		printf("ver_size_new_align: %i, Int/mol: %f\n",  verlet_size,
									    			(1.0*verlet_size)/sysp.molecule_count_aligned);
		once_seq -= 1;
	}
	
	if (verlet_size == 0) // Handle the empty neighbour list
	{
		return verlet_size;
	}
	
	free(*verlet);
	*verlet = (ushort2*) malloc( verlet_size*sizeof(ushort2) ); // allocate the verlet again expect it is filled with zeros
	for (int i = 0; i < verlet_size; ++i)
	{
		(*verlet)[i].x = 0.0;
		(*verlet)[i].y = 0.0;
	}

	for (int i = 0; i < sysp.molecule_count_aligned; ++i)
	{
		for (int j = 0; j < sysp.molecule_count_aligned; ++j)
		{
			verlet_continual_reconstruct(position_x, position_y, position_z,
										 verlet, verlet_index,
										 i, j);
		}
	}
	
	free(verlet_occupancy);
	free(verlet_index);
	return verlet_size;
}
