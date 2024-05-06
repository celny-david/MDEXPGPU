/*
 * forcefield.cu
 * provides the standard force field implementations for the MD
 * {12-6 LJ, }
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#include "atomic_add_double.h" // for cuda capable devices <6 provides double atomicAddOld
#include "err_check.h"
#include "lj_spline.h"
#include "lj_elspline.h"
#include "value.h"
#include "cuda_helper.h"
#include "cuda_variables.h"
#include "cuda_definition.h"
#include "system_property.h"

// === =========================== ===
// === SPLINE evaluation functions ===
// === =========================== ===

/* Function: lj_spline2mmk
 * device function for calculation of the standard potential interaction
 *
 * between two atoms given with the distance differences
 *
 * Notes:
 * >the implementation has checks for the exceed of cutoff distance
 * >optimised for minimal local variable use
 *
 * Parameter:
 * dx_in  - value_t position difference between atoms
 * dy_in  - value_t position difference between atoms
 * dz_in  - value_t position difference between atoms
 * i_in   - integer index of corresponding output force location
 * output - pointer to the output array [epot, force_x, force_y, force_z]
 */
__device__ void lj_spline2mmk (value_t dx_in, value_t dy_in, value_t dz_in,
							   int i_in, value_t *output)
{
	value_t dr_fr;

#if FOR_TYPE == 0
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot requires r
	if (dr_fr<d_cutoff)
	{
		output[0] += lj_pot(dr_fr);
		dr_fr = lj_for(dr_fr);

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif FOR_TYPE == 1
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot requires r
	if (dr_fr < d_cutoff)
	{
		output[0] += lj_pot(dr_fr);
		dr_fr = lj_for(dr_fr)/(dr_fr);

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif FOR_TYPE == 2
	dr_fr = dx_in*dx_in +dy_in*dy_in +dz_in*dz_in; // present version lj_force, lj_pot requires r
	if (dr_fr < (d_cutoff*d_cutoff))
	{
#if SUBSTANCE == ARGON
		// 16.07.2020 edit rescale for argon case for usage of LJF potential 
		// sigma_ar = 3.4[Angstrom]
		// epsilon_ar = 120[K] 
		dr_fr /= 3.4*3.4; // not possible to put to width - x eval in polynomial would be wrong
		output[0] += 120*lj_pot(dr_fr);
		dr_fr = 120*lj_for(dr_fr)/(3.4*3.4);
		// dr_fr = 120*lj_for(dr_fr);
#else
		output[0] += lj_pot(dr_fr);
		dr_fr = lj_for(dr_fr);
#endif
		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#endif
	return;
}

#if SUBSTANCE ==SPCE
/* Function: lj_spline2mmk
 * device function for calculation of the electrostatic potential interaction
 *
 * between two atoms given with the distance differences
 *
 * Notes:
 * >the implementation has checks for the exceed of cutoff distance
 * >optimised for minimal local variable use
 * > factor works such that instead of [epot, force_x, force_y, force_z]
 * >                     we get factor*[epot, force_x, force_y, force_z]
 *
 * Parameter:
 * dx_in  - value_t position difference between atoms
 * dy_in  - value_t position difference between atoms
 * dz_in  - value_t position difference between atoms
 * i_in   - integer index of corresponding output force location
 * factor - signed char number representing to modification of the force and epot calculation
 * output - pointer to the output array [epot, force_x, force_y, force_z]
 */
__device__ void lj_elspline2mmk (value_t dx_in, value_t dy_in, value_t dz_in,
							     int i_in, signed char factor, value_t *output)
{
	value_t dr_fr;

#if ELFOR_TYPE == 0
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot requires r
	if (dr_fr<d_cutoff_elstat )
	{
		output[0] += K_E_INV*factor*lj_elpot(dr_fr);
		dr_fr = K_E_INV*factor*lj_elfor(dr_fr);

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif ELFOR_TYPE == 1
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot requires r
	if (dr_fr<d_cutoff_elstat )
	{		
		output[0] += K_E_INV*factor*lj_elpot(dr_fr);
		dr_fr = K_E_INV*factor*lj_elfor(dr_fr)/(dr_fr);

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif ELFOR_TYPE == 2
	dr_fr = dx_in*dx_in +dy_in*dy_in +dz_in*dz_in; // present version lj_force, lj_pot requires r
	if (dr_fr<d_cutoff_elstat*d_cutoff_elstat )
	{
		output[0] += K_E_INV*factor*lj_elpot(dr_fr);
		dr_fr = K_E_INV*factor*lj_elfor(dr_fr); 

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#endif
	return;
}
#endif


// === ============================== ===
// === MOLECULE interaction functions ===
// === ============================== ===
 
#if SUBSTANCE == SPCE
/* Function: mol_interaction
 * device function for calculation of the molecule-molecule interactions
 * expected to be inlined 
 * the calls of the potential calculations are done from here
 *
 * Notes:
 * >the implementation expect only valid mol-mol pair
 * >performs double loop over the interaction pairs
 * > has a preprocessor controlled variant for the SPCE and OTHER
 *
 * Parameter:
 * id_i   			   - unsigned integer id of i-th moelcule
 * id_j   			   - unsigned integer id of j-th moelcule
 * dx_in 			   - value_t position difference between molecule centres
 * dy_in 			   - value_t position difference between molecule centres
 * dz_in  			   - value_t position difference between molecule centres
 * d_position_ref_x_in - pointer to x atomal position array
 * d_position_ref_y_in - pointer to y atomal position array
 * d_position_ref_z_in - pointer to z atomal position array
 * output - pointer to the output array [epot, force_x, force_y, force_z]
 */
__device__ void mol_interaction(unsigned int id_i,unsigned int id_j,
								value_t dx_in, value_t dy_in, value_t dz_in,
		                        value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
		                        value_t *output)
{
/*
 * the special case for SPCE
 */
	// normal pot for O-O interaction
	lj_spline2mmk(dx_in+ (d_position_ref_x_in[SUBSTANCE*id_i] - d_position_ref_x_in[SUBSTANCE*id_j]),
				  dy_in+ (d_position_ref_y_in[SUBSTANCE*id_i] - d_position_ref_y_in[SUBSTANCE*id_j]),
				  dz_in+ (d_position_ref_z_in[SUBSTANCE*id_i] - d_position_ref_z_in[SUBSTANCE*id_j]),
				  0, output);
	// electrostatic interaction
#pragma unroll
	for (int ii = 0; ii < SUBSTANCE; ii++)
	{
#pragma unroll
		for (int jj = 0; jj < SUBSTANCE; jj++)
		{
			lj_elspline2mmk(dx_in+ (d_position_ref_x_in[SUBSTANCE*id_i+ii] - d_position_ref_x_in[SUBSTANCE*id_j+jj]),
						    dy_in+ (d_position_ref_y_in[SUBSTANCE*id_i+ii] - d_position_ref_y_in[SUBSTANCE*id_j+jj]),
						    dz_in+ (d_position_ref_z_in[SUBSTANCE*id_i+ii] - d_position_ref_z_in[SUBSTANCE*id_j+jj]),
						    ii,d_mol_charge_template[ii]*d_mol_charge_template[jj] , output);
		}
	}
}
#else
__device__ void mol_interaction(unsigned int id_i,unsigned int id_j,
								value_t dx_in, value_t dy_in, value_t dz_in,
		                        value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
		                        value_t *output)
{

#pragma unroll
	for (int ii = 0; ii < SUBSTANCE; ii++)
	{
#pragma unroll
		for (int jj = 0; jj < SUBSTANCE; jj++)
		{
			lj_spline2mmk(dx_in+ (d_position_ref_x_in[SUBSTANCE*id_i+ii] - d_position_ref_x_in[SUBSTANCE*id_j+jj]),
						  dy_in+ (d_position_ref_y_in[SUBSTANCE*id_i+ii] - d_position_ref_y_in[SUBSTANCE*id_j+jj]),
						  dz_in+ (d_position_ref_z_in[SUBSTANCE*id_i+ii] - d_position_ref_z_in[SUBSTANCE*id_j+jj]),
						  ii, output);
		}
	}
}
#endif

// === ================================ ===
// === DATASTRUCTURE evaluation methods ===
// === ================================ ===
 
/* Function: n_body_dline
 * kernel function for molecule pair interaction evaluation
 *
 * naive "sequential" approach for dynamic Verlet list
 *
 * Notes:
 * >expect the 1D call with the threads are distributed into:
 * >operates with single thread per neighbour list interaction pair
 *
 * Parameter:
 * d_position_x_in		- pointer to x molecular position array
 * d_position_y_in		- pointer to y molecular position array
 * d_position_z_in		- pointer to z molecular position array
 * d_position_ref_x_in	- pointer to x atomal position array
 * d_position_ref_y_in	- pointer to y atomal position array
 * d_position_ref_z_in	- pointer to z atomal position array
 * d_force_ref_x_in		- pointer to x atomal force array
 * d_force_ref_y_in		- pointer to y atomal force array
 * d_force_ref_z_in		- pointer to z atomal force array
 * d_epot_in			- pointer to array of potential energy per atom
 * d_verlet_cont	    - pointer to dynamic Verlet list with pair interactions
 */
__global__ void n_body_dline (value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					  		  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							  value_t *d_epot_in,
							  ushort2 *d_verlet_cont)
{
/*
 * calculates the forces from dynamic verlet pair datastructure that is not padded
 * it is not workloaded with single dimension expected call thread per pair
 */
	uint tid = (blockDim.x*blockIdx.x + threadIdx.x); // the ith atom index in block // BEWARE it has to be larger than the verlet list
	value_t output[3*SUBSTANCE+1] = {0.0};//default init all atomal forces +energy
	value_t dx,dy,dz;
	ushort id_i =d_verlet_cont[tid].x;
	ushort id_j =d_verlet_cont[tid].y;

	// DEBUG section
	// printf("%d, %d", id_i, id_j);
	// required for factoring out the empty 0-0 alignment interactions in dynamic verlet list
	if( id_i == id_j) 
	{
		// NOTE catch the aligment 0-0 invalid interaction and i==j interaction
		// printf("Error id_i == id_j = %d\n",id_i );
		return;
	}

	dx = d_position_x_in[id_i] - d_position_x_in[id_j];
	dx -= d_lx*round(dx/d_lx);
	dy = d_position_y_in[id_i] - d_position_y_in[id_j];
	dy -= d_ly*round(dy/d_ly);
	dz = d_position_z_in[id_i] - d_position_z_in[id_j];
	dz -= d_lz*round(dz/d_lz);

	mol_interaction(id_i,id_j,
					dx, dy, dz,
					d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
					output);

#if PREC == DOUBLE && __CUDA_ARCH__ <600
#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		atomicAddOld(&d_force_ref_x_in[SUBSTANCE*id_i+i], output[1+3*i]);
		atomicAddOld(&d_force_ref_y_in[SUBSTANCE*id_i+i], output[2+3*i]);
		atomicAddOld(&d_force_ref_z_in[SUBSTANCE*id_i+i], output[3+3*i]);
	}
	atomicAddOld(&d_epot_in[id_i], output[0]);
#else
#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		atomicAdd(&d_force_ref_x_in[SUBSTANCE*id_i+i], output[1+3*i]);
		atomicAdd(&d_force_ref_y_in[SUBSTANCE*id_i+i], output[2+3*i]);
		atomicAdd(&d_force_ref_z_in[SUBSTANCE*id_i+i], output[3+3*i]);
	}
	atomicAdd(&d_epot_in[id_i], output[0]);
#endif

	return;
}

/* Function: n_body_wdline
 * kernel function for molecule pair interaction evaluation with workloading
 * 
 * naive "sequential" approach for dynamic Verlet list
 *
 * BEWARE: UNFINISHED IMPLEMENTATION
 * Notes:
 * NOT READY YET 18-07-2020
 * >expect the 1D call with the threads are distributed into:
 * >operates with single thread per neighbour list interaction pair
 *
 * Parameter:
 * d_position_x_in		- pointer to x molecular position array
 * d_position_y_in		- pointer to y molecular position array
 * d_position_z_in		- pointer to z molecular position array
 * d_position_ref_x_in	- pointer to x atomal position array
 * d_position_ref_y_in	- pointer to y atomal position array
 * d_position_ref_z_in	- pointer to z atomal position array
 * d_force_ref_x_in		- pointer to x atomal force array
 * d_force_ref_y_in		- pointer to y atomal force array
 * d_force_ref_z_in		- pointer to z atomal force array
 * d_epot_in			- pointer to array of potential energy per atom
 * d_verlet_cont	    - pointer to dynamic Verlet list with pair interactions
 */
__global__ void n_body_wdline (value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					  		  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							  value_t *d_epot_in,
							  ushort2 *d_verlet_cont)
{
/*
 * calculates the forces from dynamic verlet pair datastructure that is not padded
 * it is workloaded
 */
	ushort tid = WORKLOAD*(blockDim.x*blockIdx.x + threadIdx.x); // the ith atom index in block
	value_t output[3*SUBSTANCE+1] = {0.0};//default init all atomal forces +energy
	value_t dx,dy,dz;
	ushort id_i =d_verlet_cont[tid].x;
	ushort id_j =d_verlet_cont[tid].y;

// TODO finish the implementation
#pragma unroll
	for (int ii= 0; ii < WORKLOAD; ii++)
	{
		id_i =d_verlet_cont[tid].x;
		id_j =d_verlet_cont[tid].y;
		dx = d_position_x_in[id_i] - d_position_x_in[id_j];
		dx -= d_lx*round(dx/d_lx);
		dy = d_position_y_in[id_i] - d_position_y_in[id_j];
		dy -= d_ly*round(dy/d_ly);
		dz = d_position_z_in[id_i] - d_position_z_in[id_j];
		dz -= d_lz*round(dz/d_lz);

		mol_interaction(id_i,id_j,
						dx, dy, dz,
						d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
						output);
		tid++;
	}

#if PREC == DOUBLE && __CUDA_ARCH__ <600
#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		atomicAddOld(&d_force_ref_x_in[SUBSTANCE*id_i+i], output[1+3*i]);
		atomicAddOld(&d_force_ref_y_in[SUBSTANCE*id_i+i], output[2+3*i]);
		atomicAddOld(&d_force_ref_z_in[SUBSTANCE*id_i+i], output[3+3*i]);
	}
	atomicAddOld(&d_epot_in[id_i], output[0]);
#else
#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		atomicAdd(&d_force_ref_x_in[SUBSTANCE*id_i+i], output[1+3*i]);
		atomicAdd(&d_force_ref_y_in[SUBSTANCE*id_i+i], output[2+3*i]);
		atomicAdd(&d_force_ref_z_in[SUBSTANCE*id_i+i], output[3+3*i]);
	}
	atomicAdd(&d_epot_in[id_i], output[0]);
#endif

	return;
}

// === ========================= ===
// === FORCE CALCULATION methods ===
// === ========================= ===

/* Function: n_body_dline
 * kernel function for molecule pair interaction evaluation
 * incorporate the expansion which requires the variable box size input
 *
 * naive "sequential" approach for dynamic Verlet list
 *
 * Notes:
 * >expect the 1D call with the threads are distributed into:
 * >operates with single thread per neighbour list interaction pair
 * > 
 *
 * Parameter:
 * d_position_x_in		- pointer to x molecular position array
 * d_position_y_in		- pointer to y molecular position array
 * d_position_z_in		- pointer to z molecular position array
 * d_position_ref_x_in	- pointer to x atomal position array
 * d_position_ref_y_in	- pointer to y atomal position array
 * d_position_ref_z_in	- pointer to z atomal position array
 * d_force_ref_x_in		- pointer to x atomal force array
 * d_force_ref_y_in		- pointer to y atomal force array
 * d_force_ref_z_in		- pointer to z atomal force array
 * d_epot_in			- pointer to array of potential energy per atom
 * d_verlet_cont	    - pointer to dynamic Verlet list with pair interactions
 * d_lx_buffer		    - pointer to the calculation volume sizes (assume the sizes are same)
 * d_lx_possition	    - the index of the current size in use
 */
__global__ void n_body_dline_exp (value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							  	  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
						  		  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
								  value_t *d_epot_in,
								  ushort2 *d_verlet_cont,
								  value_t *d_lx_buffer, unsigned long lx_position)
{
/*
 * calculates the forces from dynamic verlet pair datastructure that is not padded
 * it is not workloaded with single dimension expected call thread per pair
 * employs the expansion -> require the sizes buffer and the index into the buffer
 */
	uint tid = (blockDim.x*blockIdx.x + threadIdx.x); // the ith atom index in block // BEWARE it has to be larger than the verlet list
	value_t output[3*SUBSTANCE+1] = {0.0};//default init all atomal forces +energy
	value_t dx,dy,dz;
	ushort id_i =d_verlet_cont[tid].x;
	ushort id_j =d_verlet_cont[tid].y;

	// DEBUG section
	// printf("%d, %d", id_i, id_j);
	// required for factoring out the empty 0-0 alignment interactions in dynamic verlet list
	if( id_i == id_j) 
	{
		// NOTE catch the aligment 0-0 invalid interaction and i==j interaction
		// printf("Error id_i == id_j = %d\n",id_i );
		return;
	}

	dx = d_position_x_in[id_i] - d_position_x_in[id_j];
	dx -= d_lx_buffer[2*lx_position]*round(dx/d_lx_buffer[2*lx_position]);
	dy = d_position_y_in[id_i] - d_position_y_in[id_j];
	dy -= d_lx_buffer[2*lx_position]*round(dy/d_lx_buffer[2*lx_position]);
	dz = d_position_z_in[id_i] - d_position_z_in[id_j];
	dz -= d_lx_buffer[2*lx_position]*round(dz/d_lx_buffer[2*lx_position]);


	mol_interaction(id_i,id_j,
					dx, dy, dz,
					d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
					output);

#if PREC == DOUBLE && __CUDA_ARCH__ <600
#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		atomicAddOld(&d_force_ref_x_in[SUBSTANCE*id_i+i], output[1+3*i]);
		atomicAddOld(&d_force_ref_y_in[SUBSTANCE*id_i+i], output[2+3*i]);
		atomicAddOld(&d_force_ref_z_in[SUBSTANCE*id_i+i], output[3+3*i]);
	}
	atomicAddOld(&d_epot_in[id_i], output[0]);
#else
#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		atomicAdd(&d_force_ref_x_in[SUBSTANCE*id_i+i], output[1+3*i]);
		atomicAdd(&d_force_ref_y_in[SUBSTANCE*id_i+i], output[2+3*i]);
		atomicAdd(&d_force_ref_z_in[SUBSTANCE*id_i+i], output[3+3*i]);
	}
	atomicAdd(&d_epot_in[id_i], output[0]);
#endif

	return;
}

/* Function: force_calculation
 * procedure for calculation of the forces
 *
 * naive "sequential" workloaded approach for Verlet list
 *
 * Notes:
 * >has check for Verlet list width overflow -> then invalid calculation
 * >should implement the sparse matrix format
 * >operates with single thread per neighbour list Workload block
 *
 * Parameter:
 * d_position_x_in		- pointer to x molecular position array
 * d_position_y_in		- pointer to y molecular position array
 * d_position_z_in		- pointer to z molecular position array
 * d_position_ref_x_in	- pointer to x atomal position array
 * d_position_ref_y_in	- pointer to y atomal position array
 * d_position_ref_z_in	- pointer to z atomal position array
 * d_force_ref_x_in		- pointer to x atomal force array
 * d_force_ref_y_in		- pointer to y atomal force array
 * d_force_ref_z_in		- pointer to z atomal force array
 * d_ekin_in			- pointer to array of kinetic energy per atom
 * d_verlet		   		- pointer to Verlet list
 * d_d_verlet_occupancy - pointer to Verlet occupancy list
 */
 void force_calculation(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
 					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
 					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
 					   value_t *d_epot_in,
 					   unsigned int verlet_size_in,
 					   ushort2 *d_verlet,
 					   value_t*d_lx_buffer, unsigned long lx_postition)
{
	// TODO properly set the kernel launch parameter to cover the vector dimension
	unsigned int set_tpb = THREAD_P_BLOCK_RECONSTRUCT; //the mol count has to be divisible with this
	unsigned int set_blk = sysp.molecule_count_aligned/set_tpb;

	// printf("DEBUG blk, tpb: %d %d ", set_blk, set_tpb );	

	set_zero<<<set_blk*SUBSTANCE, set_tpb>>>(d_force_ref_x_in); // atom count lenght
	set_zero<<<set_blk*SUBSTANCE, set_tpb>>>(d_force_ref_y_in); // atom count lenght
	set_zero<<<set_blk*SUBSTANCE, set_tpb>>>(d_force_ref_z_in); // atom count lenght

	set_zero<<<set_blk, set_tpb>>>(d_epot_in); // molecule count length
	
	cudaSafeKernell();

	// DEBUG write the verlet list 
	// writeout_1d<<<ceil(verlet_size_in/64),64>>>(d_verlet,1);
	// printf("DEBUG verlet_size %d \n", verlet_size_in/THREAD_P_BLOCK_RECONSTRUCT);
	// cudaSafeKernell();
#ifdef EXPANSION
	// TODO reimplement call and parameters
	n_body_dline_exp<<<verlet_size_in/THREAD_P_BLOCK_RECONSTRUCT,
					THREAD_P_BLOCK_RECONSTRUCT>>>
					(d_position_x_in, d_position_y_in, d_position_z_in,
					 d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
					 d_force_ref_x_in, d_force_ref_y_in, d_force_ref_z_in,
					 d_epot_in,
					 d_verlet,
					 d_lx_buffer,lx_postition);

#else	
	n_body_dline<<< verlet_size_in/THREAD_P_BLOCK_RECONSTRUCT,
					THREAD_P_BLOCK_RECONSTRUCT>>>
					(d_position_x_in, d_position_y_in, d_position_z_in,
					 d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
					 d_force_ref_x_in, d_force_ref_y_in, d_force_ref_z_in,
					 d_epot_in,
					 d_verlet);
#endif

 	cudaSafeKernell();
 	return;
}
