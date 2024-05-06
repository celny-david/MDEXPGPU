/*
 * forcefield_seq.cu
 * provides the standard force field implementations for the MD in sequential manner
 * {12-6 LJ, }
 *  Created on: 18 07 2020
 *      Author: dixiw
 */
#include <stdio.h>
#include "lj_spline.h"
#include "lj_elspline.h"
#include "value.h"
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
void lj_spline2mmk_seq (value_t dx_in, value_t dy_in, value_t dz_in,
				  		int i_in, value_t *output)
{
	value_t dr_fr;

#if FOR_TYPE == 0
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot_seq requires r
	if (dr_fr<sysp.cutoff)
	{
		output[0] += lj_pot_seq(dr_fr);
		dr_fr = lj_for_seq(dr_fr);

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif FOR_TYPE == 1
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot_seq requires r
	if (dr_fr<sysp.cutoff)
	{
		output[0] += lj_pot_seq(dr_fr);
		dr_fr = lj_for_seq(dr_fr)/(dr_fr);

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif FOR_TYPE == 2
	dr_fr = dx_in*dx_in +dy_in*dy_in +dz_in*dz_in; // present version lj_force, lj_pot_seq requires r
	if (dr_fr < (sysp.cutoff*sysp.cutoff))
	{
#if SUBSTANCE == ARGON 
		// 16.07.2020 edit rescale for argon case for usage of LJF potential case -1,6
		// sigma_ar = 3.4[Angstrom]
		// epsilon_ar = 120[K] 
		dr_fr /= 3.4*3.4; 
		output[0] += 120*lj_pot_seq(dr_fr);
		dr_fr = 120*lj_for_seq(dr_fr)/(3.4*3.4);
#else
		output[0] += lj_pot_seq(dr_fr);
		dr_fr = lj_for_seq(dr_fr);
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
void lj_elspline2mmk_seq (value_t dx_in, value_t dy_in, value_t dz_in,
  				  		  int i_in, signed char factor, value_t *output)
{
	value_t dr_fr;

#if ELFOR_TYPE == 0
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot_seq requires r
	if (dr_fr<sysp.cutoff_elstat )
	{
		output[0] += lj_elpot_seq(dr_fr)*factor;
		dr_fr = lj_elfor_seq(dr_fr)*factor; // BUG DC edit if the elstat forces are flipped as + + charges should

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif ELFOR_TYPE == 1
	dr_fr = sqrt(dx_in*dx_in +dy_in*dy_in +dz_in*dz_in); // present version lj_force, lj_pot_seq requires r
	if (dr_fr<sysp.cutoff_elstat )
	{
		output[0] += K_E_INV*factor*lj_elpot_seq(dr_fr);
		dr_fr = K_E_INV*factor*lj_elfor_seq(dr_fr)/(dr_fr); // BUG DC edit if the elstat forces are not flipped as + + should do force opposite

		output[1+3*i_in] += dr_fr*dx_in;
		output[2+3*i_in] += dr_fr*dy_in;
		output[3+3*i_in] += dr_fr*dz_in;
	}
#elif ELFOR_TYPE == 2
	dr_fr = dx_in*dx_in +dy_in*dy_in +dz_in*dz_in; // present version lj_force, lj_pot_seq requires r
	if (dr_fr<sysp.cutoff_elstat*sysp.cutoff_elstat )
	{
		output[0] += lj_elpot_seq(dr_fr)*factor;
		dr_fr = lj_elfor_seq(dr_fr)*factor;

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
 *
 * the calls of the potential calculations are done from here
 *
 * Notes:
 * >the implementation expect only valid mol-mol pair
 * >performs double loop over the interaction pairs
 * > has a preprocessor controlled variant for the SPCE
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
void mol_interaction_seq(unsigned int id_i,unsigned int id_j,
					value_t dx_in, value_t dy_in, value_t dz_in,
                    value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
                    value_t *output)
{
/*
 * the special case for SPCE
 */
	// normal pot for O-O interaction
	lj_spline2mmk_seq(dx_in + (d_position_ref_x_in[SUBSTANCE*id_i] - d_position_ref_x_in[SUBSTANCE*id_j]),
					  dy_in + (d_position_ref_y_in[SUBSTANCE*id_i] - d_position_ref_y_in[SUBSTANCE*id_j]),
					  dz_in + (d_position_ref_z_in[SUBSTANCE*id_i] - d_position_ref_z_in[SUBSTANCE*id_j]),
					  0, output);

	// electrostatic interaction
#pragma unroll
	for (int ii = 0; ii < SUBSTANCE; ii++)
	{
#pragma unroll
		for (int jj = 0; jj < SUBSTANCE; jj++)
		{
			lj_elspline2mmk_seq(dx_in+ (d_position_ref_x_in[SUBSTANCE*id_i+ii] - d_position_ref_x_in[SUBSTANCE*id_j+jj]),
							    dy_in+ (d_position_ref_y_in[SUBSTANCE*id_i+ii] - d_position_ref_y_in[SUBSTANCE*id_j+jj]),
							    dz_in+ (d_position_ref_z_in[SUBSTANCE*id_i+ii] - d_position_ref_z_in[SUBSTANCE*id_j+jj]),
							    ii, sysp.mol_mass_template[ii]*sysp.mol_mass_template[jj], output);
		}
	}	
}
#else
void mol_interaction_seq(unsigned int id_i,unsigned int id_j,
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
			lj_spline2mmk_seq(dx_in+ (d_position_ref_x_in[SUBSTANCE*id_i+ii] - d_position_ref_x_in[SUBSTANCE*id_j+jj]),
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
 * position_x_in		- pointer to x molecular position array
 * position_y_in		- pointer to y molecular position array
 * position_z_in		- pointer to z molecular position array
 * position_ref_x_in	- pointer to x atomal position array
 * position_ref_y_in	- pointer to y atomal position array
 * position_ref_z_in	- pointer to z atomal position array
 * force_ref_x_in		- pointer to x atomal force array
 * force_ref_y_in		- pointer to y atomal force array
 * force_ref_z_in		- pointer to z atomal force array
 * epot_in				- pointer to array of potential energy per atom
 * verlet_cont	 	    - pointer to dynamic Verlet list with pair interactions
 * tid 					- the substitute for thread index used in iteration
 */

int once_seq2 = -1;
int fcount = 0;

void n_body_dline (value_t *position_x_in, value_t *position_y_in, value_t *position_z_in,
				  value_t *position_ref_x_in, value_t *position_ref_y_in, value_t *position_ref_z_in,
		  		  value_t *force_ref_x_in, value_t *force_ref_y_in, value_t *force_ref_z_in,
				  value_t *epot_in,
				  ushort2 *verlet_cont,
				  unsigned int tid)
{
/*
 * calculates the forces from dynamic verlet pair datastructure that is not padded
 * it is not workloaded with single dimension expected call thread per pair
 */
	value_t output[3*SUBSTANCE+1] = {0.0};//default init all atomal forces +energy
	value_t dx,dy,dz;
	ushort id_i =verlet_cont[tid].x;
	ushort id_j =verlet_cont[tid].y;

	if( id_i == id_j) 
	{
		// NOTE catch the aligment 0-0 invalid interaction and i==j interaction
		// printf("Error id_i == id_j = %d\n",id_i );
		return;
	}

	dx = position_x_in[id_i] - position_x_in[id_j];
	dx -= sysp.lx*round(dx/sysp.lx);
	dy = position_y_in[id_i] - position_y_in[id_j];
	dy -= sysp.ly*round(dy/sysp.ly);
	dz = position_z_in[id_i] - position_z_in[id_j];
	dz -= sysp.lz*round(dz/sysp.lz);
		
	mol_interaction_seq(id_i,id_j,
						dx, dy, dz,
						position_ref_x_in, position_ref_y_in, position_ref_z_in,
						output);

	// DEBUG print out the dx,dy,dz and partial forces for interaction pair - expect full verle list for correct delimiting
	// if (once_seq2  == fcount)
	// {	
	// 	printf("[ %+f, %+f, %+f] [ %+f, %+f, %+f]", dx, dy, dz, output[1], output[2], output[3]);
	// 	if((tid+1)%(sysp.molecule_count_aligned-1)==0) printf("\n");
	// }

#pragma unroll
	for(int i=0; i<SUBSTANCE; i++)
	{
		force_ref_x_in[SUBSTANCE*id_i+i] += output[1+3*i];
		force_ref_y_in[SUBSTANCE*id_i+i] += output[2+3*i];
		force_ref_z_in[SUBSTANCE*id_i+i] += output[3+3*i];
	}
	epot_in[id_i] +=  output[0];

	return;
}

// === ========================= ===
// === FORCE CALCULATION methods ===
// === ========================= ===

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
 * position_x_in		- pointer to x molecular position array
 * position_y_in		- pointer to y molecular position array
 * position_z_in		- pointer to z molecular position array
 * position_ref_x_in	- pointer to x atomal position array
 * position_ref_y_in	- pointer to y atomal position array
 * position_ref_z_in	- pointer to z atomal position array
 * force_ref_x_in		- pointer to x atomal force array
 * force_ref_y_in		- pointer to y atomal force array
 * force_ref_z_in		- pointer to z atomal force array
 * ekin_in				- pointer to array of kinetic energy per atom
 * verlet		   		- pointer to Verlet list
 * verlet_occupancy 	- pointer to Verlet occupancy list
 */
 void force_calculation_seq(value_t *position_x_in, value_t *position_y_in, value_t *position_z_in,
	 					   value_t *position_ref_x_in, value_t *position_ref_y_in, value_t *position_ref_z_in,
	 					   value_t *force_ref_x_in, value_t *force_ref_y_in, value_t *force_ref_z_in,
	 					   value_t *epot_in,
	 					   unsigned int verlet_size_in,
	 					   ushort2 *verlet)
{
	// DEBUG - print the verlet list 
	// if (once_seq2  == fcount)
	// {
	// 	for (int j = 0; j < verlet_size_in; ++j)
	// 	{
	// 		if(j%(sysp.molecule_count_aligned-1) == 0) printf("\n");
	// 		printf("[%2d,%2d], ", verlet[j].x, verlet[j].y);
	// 	}
	// 	printf("\n");
	// // 	// once_seq2 = false;
	// }

	for (int i = 0; i < sysp.atom_count_aligned; ++i)
	{
		force_ref_x_in[i] = 0.0;
		force_ref_y_in[i] = 0.0;
		force_ref_z_in[i] = 0.0;
	}
	for (int i = 0; i < sysp.molecule_count_aligned; ++i)
	{
		epot_in[i] = 0.0;
	}

	for (int i = 0; i < verlet_size_in; ++i)
	{
		n_body_dline(position_x_in, position_y_in, position_z_in,
				 position_ref_x_in, position_ref_y_in, position_ref_z_in,
				 force_ref_x_in, force_ref_y_in, force_ref_z_in,
				 epot_in,
				 verlet,
				 i);
		// DEBUG print the final force contribution for molecule - it has to be sum uf partial printed from n_body_line
		// if (once_seq2  == fcount)
		// {
		// 	if((i+1)%(sysp.molecule_count_aligned-1)==0)
		// 	{
		// 		int tmp = i/(sysp.molecule_count_aligned-1);
		// 		printf("                                   = %+f, %+f, %+f= \n", force_ref_x_in[tmp], force_ref_y_in[tmp], force_ref_z_in[tmp]);
		// 	}
		// 	printf("\n");
		// }
	}
	// if (once_seq2  == fcount)
	// {			
	// 	once_seq2 = false;
	// }
	// fcount += 1;
 	return;
}
