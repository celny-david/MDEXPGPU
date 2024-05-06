/*
 * iterators.cu
 *
 *  Created on: 19 Dec 2019
 *      Author: dixiw
 */
#include "value.h"
#include "cuda_definition.h"
#include "system_property.h"
#include "iterative_property.h"
#include <stdio.h>
// ===== LEAP FROG scheme related functions =====

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
void shake_seq(value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
			  value_t *position_old,
			  value_t *rij,
			  int mol_ind_ref)
{
	value_t alpha_m; // alpha variable with mass addition
	value_t rij_old[3];
#if SUBSTANCE==SPCE
	static const value_t bond_template[2] = {1.0,2.0*sin(0.5*acos(-1.0/3.0))};
#endif
#pragma unroll
	for (int it = 0; it < itep.n_shake; it++)
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
				rij_old[0] = position_old[  3*ii] - position_old[  3*jj];
				rij_old[1] = position_old[1+3*ii] - position_old[1+3*jj];
				rij_old[2] = position_old[2+3*ii] - position_old[2+3*jj];

				// alpha initial computation
				alpha_m = 1.28* // the bulgarian constant that enhances efficiency of shake
						 ( rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]
#if SUBSTANCE==NITROGEN
						   -1.08892776*1.08892776)/
#elif SUBSTANCE==SPCE
				           -bond_template[ii]*bond_template[ii])/ //indexing is taken from the first atom that determines type of bond O-H, H-H
#else // beware this can lead to drift in bonds - espetially in case of lower shake iteration count
				           -rij_old[0]*rij_old[0] - rij_old[1]*rij_old[1] - rij_old[2]*rij_old[2])/
#endif
				          (rij[0]*rij_old[0] + rij[1]*rij_old[1] + rij[2]*rij_old[2]);
				// the mass addition with leftover 1/2.0
				alpha_m /= 2.0*(sysp.mol_mass_template[ii] + sysp.mol_mass_template[jj]);
				// TODO chceck if alpha is already small enough -> return to prevent overiteration -> numerical precision induced erros
				// TODO OPTIMIZATION perform in local variable and update at the end - can be checked for increments as well
				// update of position of i atom
				d_position_ref_x_in[mol_ind_ref+ii] -= alpha_m*sysp.mol_mass_template[jj]*rij_old[0];
				d_position_ref_y_in[mol_ind_ref+ii] -= alpha_m*sysp.mol_mass_template[jj]*rij_old[1];
				d_position_ref_z_in[mol_ind_ref+ii] -= alpha_m*sysp.mol_mass_template[jj]*rij_old[2];
				// update of position of j atom
				d_position_ref_x_in[mol_ind_ref+jj] += alpha_m*sysp.mol_mass_template[ii]*rij_old[0];
				d_position_ref_y_in[mol_ind_ref+jj] += alpha_m*sysp.mol_mass_template[ii]*rij_old[1];
				d_position_ref_z_in[mol_ind_ref+jj] += alpha_m*sysp.mol_mass_template[ii]*rij_old[2];
			}
		}
	}
}

void kolafa_adjust_seq(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					  value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					  value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					  value_t *d_ekin_in,
					  unsigned int idx)
{
/*
 * compute value_t the ekin
 * !!! requires the edit of initial velocity *dt
 * !!! correct energy fluctuation is obtained only for mean kinetic energy before and after velocity update

 * this is version with thread per molecule
 *  -> consider thread per atom approach that uses shared memory to hold the
 *   intermediate variables (should count the warp%substances unused threads in warp)
 */
	// int idx = (blockIdx.x*blockDim.x+ threadIdx.x);
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
				   * sysp.mol_mass_template[i];
	}
	// the weight of force sum
	tmp_a[0] *= sysp.mol_mass_template[SUBSTANCE]; // d_mol_mass_template[substance] is already inverted
	tmp_a[1] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[2] *= sysp.mol_mass_template[SUBSTANCE];
	// molecular velocity update
	d_velocity_x_in[idx] += itep.dt*itep.dt*tmp_a[0];
	d_velocity_y_in[idx] += itep.dt*itep.dt*tmp_a[1];
	d_velocity_z_in[idx] += itep.dt*itep.dt*tmp_a[2];
	//molecular position update
	d_position_x_in[idx] += d_velocity_x_in[idx];
	d_position_y_in[idx] += d_velocity_y_in[idx];
	d_position_z_in[idx] += d_velocity_z_in[idx];
	// PCB crop per molecule
	d_position_x_in[idx] -= sysp.lx*floor(d_position_x_in[idx]/sysp.lx);
	d_position_y_in[idx] -= sysp.ly*floor(d_position_y_in[idx]/sysp.ly);
	d_position_z_in[idx] -= sysp.lz*floor(d_position_z_in[idx]/sysp.lz);

#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// update atom velocity
		d_velocity_ref_x_in[idx*SUBSTANCE+i] += itep.dt*itep.dt*(d_force_ref_x_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i] - tmp_a[0]);
		d_velocity_ref_y_in[idx*SUBSTANCE+i] += itep.dt*itep.dt*(d_force_ref_y_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i] - tmp_a[1]);
		d_velocity_ref_z_in[idx*SUBSTANCE+i] += itep.dt*itep.dt*(d_force_ref_z_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i] - tmp_a[2]);
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
	 shake_seq(d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
	 	  	  tmp_pos, tmp_a, idx*SUBSTANCE);

	// recycle this for the new center of molecule calculation
	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;
	
 #pragma unroll
 	for (i = 0; i < SUBSTANCE; i++)
 	{
 		// update atom velocity
 		tmp_a[0] += d_position_ref_x_in[idx*SUBSTANCE+i]*sysp.mol_mass_template[i];
 		tmp_a[1] += d_position_ref_y_in[idx*SUBSTANCE+i]*sysp.mol_mass_template[i];
 		tmp_a[2] += d_position_ref_z_in[idx*SUBSTANCE+i]*sysp.mol_mass_template[i];
 	}

 	// update of the molecular centre of mas with the shift in atomal CoM
 	tmp_a[0] *= sysp.mol_mass_template[SUBSTANCE]; // d_mol_mass_template[substance] is already inverted
 	tmp_a[1] *= sysp.mol_mass_template[SUBSTANCE];
 	tmp_a[2] *= sysp.mol_mass_template[SUBSTANCE];

	d_position_x_in[idx] += tmp_a[0];
	d_position_y_in[idx] += tmp_a[1];
	d_position_z_in[idx] += tmp_a[2];


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
					   * sysp.mol_mass_template[i];
	}
	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy
}

void get_mh_positions_seq(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					  value_t *d_position_x_mh_in, value_t *d_position_y_mh_in, value_t *d_position_z_mh_in,
					  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					  value_t *d_position_ref_x_mh_in, value_t *d_position_ref_y_mh_in, value_t *d_position_ref_z_mh_in,
					  value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					  value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					  unsigned int idx)
{
//NOTE: works just for ARGON now, the SHAKE inclusion is required for other molecules
	int i;
	value_t tmp_a[3];
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
	}
	// the weight of force sum
	tmp_a[0] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[1] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[2] *= sysp.mol_mass_template[SUBSTANCE];


	d_position_x_mh_in[idx]=d_position_x_in[idx]-d_velocity_x_in[idx]+0.5*itep.dt*itep.dt*tmp_a[0];
	d_position_y_mh_in[idx]=d_position_y_in[idx]-d_velocity_y_in[idx]+0.5*itep.dt*itep.dt*tmp_a[1];
	d_position_z_mh_in[idx]=d_position_z_in[idx]-d_velocity_z_in[idx]+0.5*itep.dt*itep.dt*tmp_a[2];

	// PCB crop per molecule- just for certainity
	d_position_x_mh_in[idx] -= sysp.lx*floor(d_position_x_in[idx]/sysp.lx);
	d_position_y_mh_in[idx] -= sysp.ly*floor(d_position_y_in[idx]/sysp.ly);
	d_position_z_mh_in[idx] -= sysp.lz*floor(d_position_z_in[idx]/sysp.lz);


#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		d_position_ref_x_mh_in[idx*SUBSTANCE+i] = d_position_ref_x_in[idx*SUBSTANCE+i]-d_velocity_ref_x_in[idx*SUBSTANCE+i]+0.5*itep.dt*itep.dt*d_force_ref_x_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i];
		d_position_ref_y_mh_in[idx*SUBSTANCE+i] = d_position_ref_y_in[idx*SUBSTANCE+i]-d_velocity_ref_y_in[idx*SUBSTANCE+i]+0.5*itep.dt*itep.dt*d_force_ref_y_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i];
		d_position_ref_z_mh_in[idx*SUBSTANCE+i] = d_position_ref_z_in[idx*SUBSTANCE+i]-d_velocity_ref_z_in[idx*SUBSTANCE+i]+0.5*itep.dt*itep.dt*d_force_ref_z_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i];

	}

}

void verlet_step_seq_argon( value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							value_t *d_position_x_mh_in, value_t *d_position_y_mh_in, value_t *d_position_z_mh_in,
							value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
							value_t *d_position_ref_x_mh_in, value_t *d_position_ref_y_mh_in, value_t *d_position_ref_z_mh_in,
							value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
							value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
							value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							value_t *d_ekin_in,
							unsigned int idx)
{
	value_t tmp_a[3];// = {0.0,0.0,0.0};
	value_t tmp_ekin = 0.0;
	value_t old_pos[3];

	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;

	// acceleration weighted sum computation
	tmp_a[0] += d_force_ref_x_in[idx*SUBSTANCE];
	tmp_a[1] += d_force_ref_y_in[idx*SUBSTANCE];
	tmp_a[2] += d_force_ref_z_in[idx*SUBSTANCE];
	//kinetic energy per molecule part 1 calculation
	tmp_ekin +=( d_velocity_x_in[idx] *d_velocity_x_in[idx]
			    +d_velocity_y_in[idx] *d_velocity_y_in[idx]
			    +d_velocity_z_in[idx] *d_velocity_z_in[idx] )
				   * sysp.mol_mass_template[0];

	// the weight of force sum
	tmp_a[0] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[1] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[2] *= sysp.mol_mass_template[SUBSTANCE];

	//saving the positions in t
	old_pos[0]=d_position_x_in[idx];
	old_pos[1]=d_position_y_in[idx];
	old_pos[2]=d_position_z_in[idx];

	//calculation of position in t+h
	d_position_x_in[idx]= 2.0*d_position_x_in[idx] - d_position_x_mh_in[idx] + itep.dt*itep.dt*tmp_a[0];
	d_position_y_in[idx]= 2.0*d_position_y_in[idx] - d_position_y_mh_in[idx] + itep.dt*itep.dt*tmp_a[1];
	d_position_z_in[idx]= 2.0*d_position_z_in[idx] - d_position_z_mh_in[idx] + itep.dt*itep.dt*tmp_a[2];

	//calculation of velocity in t+h, v(t)=(r(t+h)-r(t))/h, h is assumed 1 in our units
	d_velocity_x_in[idx]= (d_position_x_in[idx] - old_pos[0]);
	d_velocity_y_in[idx]= (d_position_y_in[idx] - old_pos[1]);
	d_velocity_z_in[idx]= (d_position_z_in[idx] - old_pos[2]);

	// PCB for velocities
	d_velocity_x_in[idx] -= sysp.lx*round(d_velocity_x_in[idx]/sysp.lx);
	d_velocity_y_in[idx] -= sysp.ly*round(d_velocity_y_in[idx]/sysp.ly);
	d_velocity_z_in[idx] -= sysp.lz*round(d_velocity_z_in[idx]/sysp.lz);

	// PCB crop per molecule
	d_position_x_in[idx] -= sysp.lx*floor(d_position_x_in[idx]/sysp.lx);
	d_position_y_in[idx] -= sysp.ly*floor(d_position_y_in[idx]/sysp.ly);
	d_position_z_in[idx] -= sysp.lz*floor(d_position_z_in[idx]/sysp.lz);

	//saving the positions in t into positions t-h for future step
	d_position_x_mh_in[idx]=old_pos[0];
	d_position_y_mh_in[idx]=old_pos[1];
	d_position_z_mh_in[idx]=old_pos[2];

	tmp_ekin +=( d_velocity_x_in[idx] *d_velocity_x_in[idx]
			    +d_velocity_y_in[idx] *d_velocity_y_in[idx]
			    +d_velocity_z_in[idx] *d_velocity_z_in[idx] )
			   * sysp.mol_mass_template[0];

	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy
}


void verlet_step_seq_shaked( value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							 value_t *d_position_x_mh_in, value_t *d_position_y_mh_in, value_t *d_position_z_mh_in,
							 value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
							 value_t *d_position_ref_x_mh_in, value_t *d_position_ref_y_mh_in, value_t *d_position_ref_z_mh_in,
							 value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
							 value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
							 value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							 value_t *d_ekin_in,
							 unsigned int idx)
{
//NOTE: works for NITROGEN and SPCE, using the SHAKE
	int i;
	value_t tmp_a[3];// = {0.0,0.0,0.0};
	value_t tmp_ekin = 0.0;
	value_t old_pos[3];
	value_t old_pos_ref[3*SUBSTANCE];
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
				   * sysp.mol_mass_template[i];
	}
	// the weight of force sum
	tmp_a[0] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[1] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[2] *= sysp.mol_mass_template[SUBSTANCE];

	//saving molecular positions in t
	old_pos[0]=d_position_x_in[idx];
	old_pos[1]=d_position_y_in[idx];
	old_pos[2]=d_position_z_in[idx];

	//calculation of molecular position in t+h
	d_position_x_in[idx]=2*d_position_x_in[idx]-d_position_x_mh_in[idx]+itep.dt*itep.dt*tmp_a[0];
	d_position_y_in[idx]=2*d_position_y_in[idx]-d_position_y_mh_in[idx]+itep.dt*itep.dt*tmp_a[1];
	d_position_z_in[idx]=2*d_position_z_in[idx]-d_position_z_mh_in[idx]+itep.dt*itep.dt*tmp_a[2];

#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		//saving atom positions in t
		old_pos_ref[  3*i]=d_position_ref_x_in[idx*SUBSTANCE+i];
		old_pos_ref[1+3*i]=d_position_ref_y_in[idx*SUBSTANCE+i];
		old_pos_ref[2+3*i]=d_position_ref_z_in[idx*SUBSTANCE+i];

		//calculation of atom positions in t+h
		d_position_ref_x_in[idx*SUBSTANCE+i]=2*d_position_ref_x_in[idx*SUBSTANCE+i]-d_position_ref_x_mh_in[idx*SUBSTANCE+i]+itep.dt*itep.dt*d_force_ref_x_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i];
		d_position_ref_y_in[idx*SUBSTANCE+i]=2*d_position_ref_y_in[idx*SUBSTANCE+i]-d_position_ref_y_mh_in[idx*SUBSTANCE+i]+itep.dt*itep.dt*d_force_ref_y_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i];
		d_position_ref_z_in[idx*SUBSTANCE+i]=2*d_position_ref_z_in[idx*SUBSTANCE+i]-d_position_ref_z_mh_in[idx*SUBSTANCE+i]+itep.dt*itep.dt*d_force_ref_z_in[idx*SUBSTANCE+i]/sysp.mol_mass_template[i];

		// save the atom positions before shake
		tmp_pos[  3*i] = d_position_ref_x_in[idx*SUBSTANCE+i];
		tmp_pos[1+3*i] = d_position_ref_y_in[idx*SUBSTANCE+i];
		tmp_pos[2+3*i] = d_position_ref_z_in[idx*SUBSTANCE+i];
	}
	// shake algorithm
	shake_seq(d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
			  tmp_pos, tmp_a, idx*SUBSTANCE);

	// recycle this for the new center of molecule calculation
	tmp_a[0] = 0.0;
	tmp_a[1] = 0.0;
	tmp_a[2] = 0.0;

#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// center of mass calculaton part one
		tmp_a[0] += d_position_ref_x_in[idx*SUBSTANCE+i]*sysp.mol_mass_template[i];
		tmp_a[1] += d_position_ref_y_in[idx*SUBSTANCE+i]*sysp.mol_mass_template[i];
		tmp_a[2] += d_position_ref_z_in[idx*SUBSTANCE+i]*sysp.mol_mass_template[i];
	}

	// center of mass calculaton part two
	tmp_a[0] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[1] *= sysp.mol_mass_template[SUBSTANCE];
	tmp_a[2] *= sysp.mol_mass_template[SUBSTANCE];

		// adding center of mass to molecular position
	d_position_x_in[idx] += tmp_a[0];
	d_position_y_in[idx] += tmp_a[1];
	d_position_z_in[idx] += tmp_a[2];

	//calculation of molecular velocity in t+h, v(t)=(r(t+h)-r(t))/h, h is assumed 1 in our units
	d_velocity_x_in[idx]=d_position_x_in[idx]-old_pos[0];
	d_velocity_y_in[idx]=d_position_y_in[idx]-old_pos[1];
	d_velocity_z_in[idx]=d_position_z_in[idx]-old_pos[2];

	// PCB for velocities
	// BUG is the PCB for velocities really required ???
	d_velocity_x_in[idx] -= sysp.lx*round(d_velocity_x_in[idx]/sysp.lx);
	d_velocity_y_in[idx] -= sysp.ly*round(d_velocity_y_in[idx]/sysp.ly);
	d_velocity_z_in[idx] -= sysp.lz*round(d_velocity_z_in[idx]/sysp.lz);

	// PCB crop per molecule
	d_position_x_in[idx] -= sysp.lx*floor(d_position_x_in[idx]/sysp.lx);
	d_position_y_in[idx] -= sysp.ly*floor(d_position_y_in[idx]/sysp.ly);
	d_position_z_in[idx] -= sysp.lz*floor(d_position_z_in[idx]/sysp.lz);

	//saving the positions in t into positions t-h for future step
	d_position_x_mh_in[idx]=old_pos[0];
	d_position_y_mh_in[idx]=old_pos[1];
	d_position_z_mh_in[idx]=old_pos[2];

#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// modification of the atom position with regards to new CoM
		d_position_ref_x_in[idx*SUBSTANCE+i] -= tmp_a[0];
		d_position_ref_y_in[idx*SUBSTANCE+i] -= tmp_a[1];
		d_position_ref_z_in[idx*SUBSTANCE+i] -= tmp_a[2];

//		calculation of velocity in t+h, v(t)=(r(t+h)-r(t-h))/2h, h is assumed 1 in our units
//		d_velocity_ref_x_in[idx*SUBSTANCE+i]=0.5*(d_position_ref_x_in[idx*SUBSTANCE+i]-d_position_ref_x_mh_in[idx*SUBSTANCE+i]);
//		d_velocity_ref_y_in[idx*SUBSTANCE+i]=0.5*(d_position_ref_y_in[idx*SUBSTANCE+i]-d_position_ref_y_mh_in[idx*SUBSTANCE+i]);
//		d_velocity_ref_z_in[idx*SUBSTANCE+i]=0.5*(d_position_ref_z_in[idx*SUBSTANCE+i]-d_position_ref_z_mh_in[idx*SUBSTANCE+i]);

		//calculation of atom velocity in t+h, v(t)=(r(t+h)-r(t))/h, h is assumed 1 in our units
		d_velocity_ref_x_in[idx*SUBSTANCE+i]=d_position_ref_x_in[idx*SUBSTANCE+i]-old_pos_ref[  3*i];
		d_velocity_ref_y_in[idx*SUBSTANCE+i]=d_position_ref_y_in[idx*SUBSTANCE+i]-old_pos_ref[1+3*i];
		d_velocity_ref_z_in[idx*SUBSTANCE+i]=d_position_ref_z_in[idx*SUBSTANCE+i]-old_pos_ref[2+3*i];

		//saving the positions in t into positions t-h for future step
		d_position_ref_x_mh_in[idx*SUBSTANCE+i]=old_pos_ref[  3*i];
		d_position_ref_y_mh_in[idx*SUBSTANCE+i]=old_pos_ref[1+3*i];
		d_position_ref_z_mh_in[idx*SUBSTANCE+i]=old_pos_ref[2+3*i];

		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
					+(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
					+(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * sysp.mol_mass_template[i];
	}

	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy
}
