/*
 * iterators.cu
 *
 *  Created on: 19 Dec 2019
 *      Author: dixiw
 */
#include "err_check.h"
#include "value.h"
#include "cuda_variables.h"
#include "cuda_definition.h"

// === ========= ===
// === SHAKE ALG ===
// === ========= ===

/* Function: shake
 * device function for preserving the bonds of bonded molecules via SHAKE iterative algorithm
 *
 * operate on the molecular level of CoM implementation
 *
 * Notes:
 * >the implementation is in form of triple loop
 * >the most outer loop is for shake iteration and inner loops are for atom pairs in the considered molecule
 *
 *
 * Parameter:
 * d_position_ref_x_in - pointer to x atomal position array
 * d_position_ref_y_in - pointer to y atomal position array
 * d_position_ref_z_in - pointer to z atomal position array
 * position_old		   - pointer to array with old positions of the moleucle
 * rij				   - pointer to the distance array for reuse purposes
 * mol_ind_ref		   - stride index pointing to the selected molecule in atomal arrays
 */
__device__ void shake(value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					  value_t *position_old,
					  value_t *rij,
					  int mol_ind_ref)
{
	// TODO do the position manipulation in local array to save repeated memory access into d_position
	value_t alpha_m; // alpha variable with mass addition
	value_t rij_old[3];
#pragma unroll
	for (int it = 0; it < d_n_shake; it++)
	{
		for (int ii = 0; ii < SUBSTANCE; ii++)
		{
			for (int jj = ii+1; jj < SUBSTANCE; jj++)
			{
				// r_ij calculation as vector subtraction
				rij[0] = d_position_ref_x_in[mol_ind_ref+ii] - d_position_ref_x_in[mol_ind_ref+jj];
				rij[1] = d_position_ref_y_in[mol_ind_ref+ii] - d_position_ref_y_in[mol_ind_ref+jj];
				rij[2] = d_position_ref_z_in[mol_ind_ref+ii] - d_position_ref_z_in[mol_ind_ref+jj];
				// r_ij calculation as vector subtraction
				// NOTE this is independent on iteration and can be precalculated in advance 
				rij_old[0] = position_old[  3*ii] - position_old[  3*jj];
				rij_old[1] = position_old[1+3*ii] - position_old[1+3*jj];
				rij_old[2] = position_old[2+3*ii] - position_old[2+3*jj];

				// alpha initial computation
				alpha_m = 1.28* // the bulgarian constant that enhances efficiency of shake
						 ( rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]
#if SUBSTANCE==NITROGEN
						   -1.08892776*1.08892776)/
#elif SUBSTANCE==SPCE
				           -d_bond_template[ii]*d_bond_template[ii])/ //indexing is taken from the first atom that determines type of bond O-H, H-H
#else // beware this can lead to drift in bonds - espetially in case of lower shake iteration count
				           -rij_old[0]*rij_old[0] - rij_old[1]*rij_old[1] - rij_old[2]*rij_old[2])/
#endif
				          (rij[0]*rij_old[0] + rij[1]*rij_old[1] + rij[2]*rij_old[2]);
				// the mass addition with leftover 1/2.0
				alpha_m /= 2.0*(d_mol_mass_template[ii] + d_mol_mass_template[jj]);
				// TODO chceck if alpha is already small enough -> return to prevent overiteration -> numerical precision induced erros
				// TODO OPTIMIZATION perform in local variable and update at the end - can be checked for increments as well
				// update of position of i atom
				d_position_ref_x_in[mol_ind_ref+ii] -= alpha_m*d_mol_mass_template[jj]*rij_old[0];
				d_position_ref_y_in[mol_ind_ref+ii] -= alpha_m*d_mol_mass_template[jj]*rij_old[1];
				d_position_ref_z_in[mol_ind_ref+ii] -= alpha_m*d_mol_mass_template[jj]*rij_old[2];
				// update of position of j atom
				d_position_ref_x_in[mol_ind_ref+jj] += alpha_m*d_mol_mass_template[ii]*rij_old[0];
				d_position_ref_y_in[mol_ind_ref+jj] += alpha_m*d_mol_mass_template[ii]*rij_old[1];
				d_position_ref_z_in[mol_ind_ref+jj] += alpha_m*d_mol_mass_template[ii]*rij_old[2];
			}
		}
	}
}

// === ========== ===
// === KOLAFA ALG ===
// === ========== ===

/* Function: kolafa_adjust
 * kernel function for iteration step of the leap-frog scheme
 *
 * updates positions, velocities and kinetic energies
 *
 * direct call to the <shake> device function
 *
 * Notes:
 * >operates with single thread per molecule
 * >the implementation is targeting molecules for easy shake operation
 * >kolafa centre of mass algorithm from 15.06.2018
 *
 * beware kernel does not have overflow control for molecule count
 * has to be controlled externally with launch configuration
 * (it is the case with number of molecules as power of 2)
 *
 * Parameter:
 * d_position_x_in     - pointer to x molecular position array
 * d_position_y_in     - pointer to y molecular position array
 * d_position_z_in 	   - pointer to z molecular position array
 * d_position_ref_x_in - pointer to x atomal position array
 * d_position_ref_y_in - pointer to y atomal position array
 * d_position_ref_z_in - pointer to z atomal position array
 * d_velocity_x_in     - pointer to x molecular velocity array
 * d_velocity_y_in     - pointer to y molecular velocity array
 * d_velocity_z_in 	   - pointer to z molecular velocity array
 * d_velocity_ref_x_in - pointer to x atomal velocity array
 * d_velocity_ref_y_in - pointer to y atomal velocity array
 * d_velocity_ref_z_in - pointer to z atomal velocity array
 * d_force_ref_x_in	   - pointer to x atomal force array
 * d_force_ref_y_in    - pointer to y atomal force array
 * d_force_ref_z_in    - pointer to z atomal force array
 * d_ekin_in		   - pointer to array of kinetic energy per atom
 */

__global__ void kolafa_adjust(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
							  value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
							  value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
							  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							  value_t *d_ekin_in)
{
/*
 * compute value_t the ekin
 * !!! requires the edit of initial velocity *dt
 * !!! correct energy fluctuation is obtained only for mean kinetic energy before and after velocity update

 * this is version with thread per molecule
 *  -> consider thread per atom approach that uses shared memory to hold the
 *   intermediate variables (should count the warp%substances unused threads in warp)
 */
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);
	int i;
	value_t tmp_a[3];// = {0.0,0.0,0.0};
	value_t tmp_ekin = 0.0;
	value_t tmp_pos[3*SUBSTANCE];

	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;
#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// acceleration weighted sum computation
		tmp_a[0] += d_force_ref_x_in[idx*SUBSTANCE+i];
		tmp_a[1] += d_force_ref_y_in[idx*SUBSTANCE+i];
		tmp_a[2] += d_force_ref_z_in[idx*SUBSTANCE+i];
		//kinetic energy per molecule part 1 calculation
		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
				    +(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
				    +(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * d_mol_mass_template[i];
	}
	// the weight of force sum
	tmp_a[0] *= d_mol_mass_template[SUBSTANCE]; // d_mol_mass_template[substance] is already inverted
	tmp_a[1] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[2] *= d_mol_mass_template[SUBSTANCE];
	// molecular velocity update
	d_velocity_x_in[idx] += d_dt2*tmp_a[0];
	d_velocity_y_in[idx] += d_dt2*tmp_a[1];
	d_velocity_z_in[idx] += d_dt2*tmp_a[2];
	//molecular position update
	d_position_x_in[idx] += d_velocity_x_in[idx];
	d_position_y_in[idx] += d_velocity_y_in[idx];
	d_position_z_in[idx] += d_velocity_z_in[idx];
	// PCB crop per molecule
	d_position_x_in[idx] -= d_lx*floor(d_position_x_in[idx]/d_lx);
	d_position_y_in[idx] -= d_ly*floor(d_position_y_in[idx]/d_ly);
	d_position_z_in[idx] -= d_lz*floor(d_position_z_in[idx]/d_lz);

#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// update atom velocity
		d_velocity_ref_x_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_x_in[idx*SUBSTANCE+i]/d_mol_mass_template[i] - tmp_a[0]);
		d_velocity_ref_y_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_y_in[idx*SUBSTANCE+i]/d_mol_mass_template[i] - tmp_a[1]);
		d_velocity_ref_z_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_z_in[idx*SUBSTANCE+i]/d_mol_mass_template[i] - tmp_a[2]);
		// save the previous iteration atom position
		tmp_pos[  3*i] = d_position_ref_x_in[idx*SUBSTANCE+i];
		tmp_pos[1+3*i] = d_position_ref_y_in[idx*SUBSTANCE+i];
		tmp_pos[2+3*i] = d_position_ref_z_in[idx*SUBSTANCE+i];
		// update the atom position
		d_position_ref_x_in[idx*SUBSTANCE+i] += d_velocity_ref_x_in[idx*SUBSTANCE+i];
		d_position_ref_y_in[idx*SUBSTANCE+i] += d_velocity_ref_y_in[idx*SUBSTANCE+i];
		d_position_ref_z_in[idx*SUBSTANCE+i] += d_velocity_ref_z_in[idx*SUBSTANCE+i];
	}

	// shake algorithm
	shake(d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
		  tmp_pos, tmp_a, idx*SUBSTANCE);

	// recycle this for the new center of molecule calculation
	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;
	
	#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// update atom velocity
		tmp_a[0] += d_position_ref_x_in[idx*SUBSTANCE+i]*d_mol_mass_template[i];
		tmp_a[1] += d_position_ref_y_in[idx*SUBSTANCE+i]*d_mol_mass_template[i];
		tmp_a[2] += d_position_ref_z_in[idx*SUBSTANCE+i]*d_mol_mass_template[i];
	}

	// update of the molecular centre of mas with the shift in atomal CoM
	tmp_a[0] *= d_mol_mass_template[SUBSTANCE]; // d_mol_mass_template[substance] is already inverted
	tmp_a[1] *= d_mol_mass_template[SUBSTANCE]; 
	tmp_a[2] *= d_mol_mass_template[SUBSTANCE];

#pragma unroll // valid only when shaking - only kinetic energy
	for (i = 0; i < SUBSTANCE; i++)
	{
		// modification of the position with regards to new CoM
		d_position_ref_x_in[idx*SUBSTANCE+i] -= tmp_a[0];
		d_position_ref_y_in[idx*SUBSTANCE+i] -= tmp_a[1];
		d_position_ref_z_in[idx*SUBSTANCE+i] -= tmp_a[2];
		// velocity recalculation after shake
		d_velocity_ref_x_in[idx*SUBSTANCE+i] = d_position_ref_x_in[idx*SUBSTANCE+i] - tmp_pos[  3*i];
		d_velocity_ref_y_in[idx*SUBSTANCE+i] = d_position_ref_y_in[idx*SUBSTANCE+i] - tmp_pos[1+3*i];
		d_velocity_ref_z_in[idx*SUBSTANCE+i] = d_position_ref_z_in[idx*SUBSTANCE+i] - tmp_pos[2+3*i];
		//kinetic energy per molecule part 2 calculation
		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
					+(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
					+(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * d_mol_mass_template[i];
	}
	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy
}

__global__ void kolafa_adjust_correction(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
										 value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
										 value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
										 value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
										 value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
										 value_t *d_ekin_in)
{
/*
 * OPTIM if the registers are issue rework it with tmp_a[3] and two loops first over position as in adjust and then for velocities
 */
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);
	int i;
	value_t tmp_a[6];// = {0.0,0.0,0.0};
	value_t tmp_ekin = 0.0;
	value_t tmp_pos[3*SUBSTANCE];

	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;
#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// acceleration weighted sum computation
		tmp_a[0] += d_force_ref_x_in[idx*SUBSTANCE+i];
		tmp_a[1] += d_force_ref_y_in[idx*SUBSTANCE+i];
		tmp_a[2] += d_force_ref_z_in[idx*SUBSTANCE+i];
		//kinetic energy per molecule part 1 calculation
		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
				    +(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
				    +(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * d_mol_mass_template[i];
	}
	// the weight of force sum
	tmp_a[0] *= d_mol_mass_template[SUBSTANCE]; // d_mol_mass_template[substance] is already inverted
	tmp_a[1] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[2] *= d_mol_mass_template[SUBSTANCE];
	// molecular velocity update
	d_velocity_x_in[idx] += d_dt2*tmp_a[0];
	d_velocity_y_in[idx] += d_dt2*tmp_a[1];
	d_velocity_z_in[idx] += d_dt2*tmp_a[2];
	//molecular position update
	d_position_x_in[idx] += d_velocity_x_in[idx];
	d_position_y_in[idx] += d_velocity_y_in[idx];
	d_position_z_in[idx] += d_velocity_z_in[idx];
	// PCB crop per molecule
	d_position_x_in[idx] -= d_lx*floor(d_position_x_in[idx]/d_lx);
	d_position_y_in[idx] -= d_ly*floor(d_position_y_in[idx]/d_ly);
	d_position_z_in[idx] -= d_lz*floor(d_position_z_in[idx]/d_lz);

#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// update atom velocity
		d_velocity_ref_x_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_x_in[idx*SUBSTANCE+i]/d_mol_mass_template[i] - tmp_a[0]);
		d_velocity_ref_y_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_y_in[idx*SUBSTANCE+i]/d_mol_mass_template[i] - tmp_a[1]);
		d_velocity_ref_z_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_z_in[idx*SUBSTANCE+i]/d_mol_mass_template[i] - tmp_a[2]);
		// save the previous iteration atom position
		tmp_pos[  3*i] = d_position_ref_x_in[idx*SUBSTANCE+i];
		tmp_pos[1+3*i] = d_position_ref_y_in[idx*SUBSTANCE+i];
		tmp_pos[2+3*i] = d_position_ref_z_in[idx*SUBSTANCE+i];
		// update the atom position
		d_position_ref_x_in[idx*SUBSTANCE+i] += d_velocity_ref_x_in[idx*SUBSTANCE+i];
		d_position_ref_y_in[idx*SUBSTANCE+i] += d_velocity_ref_y_in[idx*SUBSTANCE+i];
		d_position_ref_z_in[idx*SUBSTANCE+i] += d_velocity_ref_z_in[idx*SUBSTANCE+i];
	}

	// shake algorithm
	shake(d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
		  tmp_pos, tmp_a, idx*SUBSTANCE);

	// recycle this for the new center of molecule calculation
	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;
	tmp_a[3] = 0.0;
	tmp_a[4] = 0.0;
	tmp_a[5] = 0.0;


	#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// velocity recalculation after shake
		d_velocity_ref_x_in[idx*SUBSTANCE+i] = d_position_ref_x_in[idx*SUBSTANCE+i] - tmp_pos[  3*i];
		d_velocity_ref_y_in[idx*SUBSTANCE+i] = d_position_ref_y_in[idx*SUBSTANCE+i] - tmp_pos[1+3*i];
		d_velocity_ref_z_in[idx*SUBSTANCE+i] = d_position_ref_z_in[idx*SUBSTANCE+i] - tmp_pos[2+3*i];
		// the rounding errors correction for velocity
		tmp_a[0] += d_mol_mass_template[i]*d_velocity_ref_x_in[idx*SUBSTANCE+i];
		tmp_a[1] += d_mol_mass_template[i]*d_velocity_ref_y_in[idx*SUBSTANCE+i];
		tmp_a[2] += d_mol_mass_template[i]*d_velocity_ref_z_in[idx*SUBSTANCE+i];
		// the rounding errors correction for position
		tmp_a[3] += d_mol_mass_template[i]*d_position_ref_x_in[idx*SUBSTANCE+i];
		tmp_a[4] += d_mol_mass_template[i]*d_position_ref_y_in[idx*SUBSTANCE+i];
		tmp_a[5] += d_mol_mass_template[i]*d_position_ref_z_in[idx*SUBSTANCE+i];
	}


		// printf("idx: %d  tmp_a[3]: %.16g\n",idx,tmp_a[3]); // DEBUG index and tmp_a printout


	// update of the molecular centre of mas with the shift in atomal CoM
	tmp_a[0] *= d_mol_mass_template[SUBSTANCE]; // d_mol_mass_template[substance] is already inverted
	tmp_a[1] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[2] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[3] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[4] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[5] *= d_mol_mass_template[SUBSTANCE];

#pragma unroll // valid only when shaking - only kinetic energy
	for (i = 0; i < SUBSTANCE; i++)
	{
		// velocity correction for rounding errors
		d_velocity_ref_x_in[idx*SUBSTANCE+i] -= tmp_a[0];
		d_velocity_ref_y_in[idx*SUBSTANCE+i] -= tmp_a[1];
		d_velocity_ref_z_in[idx*SUBSTANCE+i] -= tmp_a[2];
		// position correction for rounding errors
		d_position_ref_x_in[idx*SUBSTANCE+i] -= tmp_a[3];
		d_position_ref_y_in[idx*SUBSTANCE+i] -= tmp_a[4];
		d_position_ref_z_in[idx*SUBSTANCE+i] -= tmp_a[5];
		//kinetic energy per molecule part 2 calculation
		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
					+(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
					+(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * d_mol_mass_template[i];
	}
	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy
}
