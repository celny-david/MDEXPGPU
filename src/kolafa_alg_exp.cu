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

// === ============ ===
// === TRVP HELPERS ===
// === ============ ===

/* Function: TRVP2_right_sight
 * time-reversible velocity predictor of velocity in time t
 * order of predictor is k=2, paper DOI: 10.1021/ct200108g
 * history of positions is used for prediction
 *
 * Notes:
 * >returns normal velocity[angstrom/ps], NOT h*velocity[angstrom]
 *
 * Parameter:
 * d_position_in 		- input position in time t
 * d_position_m1h_in 	- input position in time t-h
 * d_position_m2h_in 	- input position in time t-2*h
 * d_position_m3h_in	- input position in time t-3*h
 * mol_ind_ref	 		- identification number of particle
 */
__device__ value_t TRVP2_right_sight(value_t *d_position_in, value_t *d_position_m1h_in, value_t *d_position_m2h_in,
									 value_t *d_position_m3h_in, unsigned int mol_ind_ref)
{
	return ( ( 10*d_position_in[mol_ind_ref]
	         - 15*d_position_m1h_in[mol_ind_ref]
	         + 6*d_position_m2h_in[mol_ind_ref]
	         - d_position_m3h_in[mol_ind_ref]) / d_dt) / 0.6e1;
}

/* Function: TRVP2Vel_right_sight
 * time-reversible velocity predictor of velocity in time t
 * order of predictor is k=2, paper DOI: 10.1021/ct200108g
 * history of velocities is used for prediction
 *
 * Notes:
 * >returns normal velocity[angstrom/ps], NOT h*velocity[angstrom]
 *
 * Parameter:
 * d_velocity_m05h_in 	- input velocity in time t-0.5*h
 * d_velocity_m15h_in 	- input velocity in time t-1.5*h
 * d_velocity_m25h_in 	- input velocity in time t-2.5*h
 * mol_ind_ref	 		- identification number of particle
 */
__device__ value_t TRVP2Vel_right_sight(value_t *d_velocity_m05h_in, value_t *d_velocity_m15h_in, value_t *d_velocity_m25h_in,
										unsigned int mol_ind_ref)
{
	return (0.5e1 / 0.3e1 * d_velocity_m05h_in[mol_ind_ref]
	      - 0.5e1 / 0.6e1 * d_velocity_m15h_in[mol_ind_ref]
	      + d_velocity_m25h_in[mol_ind_ref] / 0.6e1) / d_dt;
}

__device__ value_t calc_lambda_buffer(value_t *d_lx_buffer, unsigned long i)
{
	unsigned long actual=2*i;
	if(actual<=1)
	{
		actual=1;
	}
	return (d_lx_buffer[actual+1]/d_lx_buffer[actual-1]);
}

__device__ value_t calc_vf_buffer(value_t *d_lx_buffer, unsigned long i)
{
	unsigned long actual=2*i;
	if(actual<=1)
	{
		actual=1;
	}
	return log(d_lx_buffer[actual+1]/d_lx_buffer[actual-1])/d_dt;
}

// === ========= ===
// === TRVP ALGS ===
// === ========= ===

/* Function: TRVP2
 * time-reversible velocity predictor of molecular velocity in time t
 * order of predictor is k=2, paper DOI: 10.1021/ct200108g
 * history of molecular positions is used for prediction
 *
 * Notes:
 * >calls device function TRVP2_right_sight
 *
 * Parameter:
 * d_velocityTRVP2_x_out - output array of predicted molecular x velocities in time t
 * d_velocityTRVP2_y_out - output array of predicted molecular y velocities in time t
 * d_velocityTRVP2_z_out - output array of predicted molecular z velocities in time t
 * d_position_x_in		 - input array of molecular x positions in time t
 * d_position_y_in	 	 - input array of molecular y positions in time t
 * d_position_z_in	     - input array of molecular z positions in time t
 * d_position_x_m1h_in	 - input array of molecular x positions in time t-h
 * d_position_y_m1h_in	 - input array of molecular y positions in time t-h
 * d_position_z_m1h_in	 - input array of molecular z positions in time t-h
 * d_position_x_m2h_in	 - input array of molecular x positions in time t-2h
 * d_position_y_m2h_in	 - input array of molecular y positions in time t-2h
 * d_position_z_m2h_in	 - input array of molecular z positions in time t-2h
 * d_position_x_m3h_in 	 - input array of molecular x positions in time t-3h
 * d_position_y_m3h_in	 - input array of molecular y positions in time t-3h
 * d_position_z_m3h_in	 - input array of molecular z positions in time t-3h
 */
__global__ void TRVP2( value_t *d_velocityTRVP2_x_out, value_t *d_velocityTRVP2_y_out, value_t *d_velocityTRVP2_z_out,
					   value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_x_m1h_in, value_t *d_position_y_m1h_in, value_t *d_position_z_m1h_in,
					   value_t *d_position_x_m2h_in, value_t *d_position_y_m2h_in, value_t *d_position_z_m2h_in,
					   value_t *d_position_x_m3h_in, value_t *d_position_y_m3h_in, value_t *d_position_z_m3h_in)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);

	d_velocityTRVP2_x_out[idx] =  TRVP2_right_sight(d_position_x_in, d_position_x_m1h_in ,d_position_x_m2h_in,d_position_x_m3h_in, idx);
	d_velocityTRVP2_y_out[idx] =  TRVP2_right_sight(d_position_y_in, d_position_y_m1h_in ,d_position_y_m2h_in,d_position_y_m3h_in, idx);
	d_velocityTRVP2_z_out[idx] =  TRVP2_right_sight(d_position_z_in, d_position_z_m1h_in ,d_position_z_m2h_in,d_position_z_m3h_in, idx);

}

/* Function: TRVP2_ref
 * time-reversible velocity predictor of atom relative velocity in time t
 * order of predictor is k=2, paper DOI: 10.1021/ct200108g
 * history of atom positions is used for prediction
 *
 * Notes:
 * >calls device function TRVP2_right_sight
 *
 * Parameter:
 * d_velocityTRVP2_ref_x_out - output array of predicted atom x velocities in time t
 * d_velocityTRVP2_ref_y_out - output array of predicted atom y velocities in time t
 * d_velocityTRVP2_ref_z_out - output array of predicted atom z velocities in time t
 * d_position_ref_x_in		 - input array of atom x positions in time t
 * d_position_ref_y_in	 	 - input array of atom y positions in time t
 * d_position_ref_z_in	     - input array of atom z positions in time t
 * d_position_ref_x_m1h_in	 - input array of atom x positions in time t-h
 * d_position_ref_y_m1h_in	 - input array of atom y positions in time t-h
 * d_position_ref_z_m1h_in	 - input array of atom z positions in time t-h
 * d_position_ref_x_m2h_in	 - input array of atom x positions in time t-2h
 * d_position_ref_y_m2h_in	 - input array of atom y positions in time t-2h
 * d_position_ref_z_m2h_in	 - input array of atom z positions in time t-2h
 * d_position_ref_x_m3h_in 	 - input array of atom x positions in time t-3h
 * d_position_ref_y_m3h_in	 - input array of atom y positions in time t-3h
 * d_position_ref_z_m3h_in	 - input array of atom z positions in time t-3h
 */
__global__ void TRVP2_ref( value_t *d_velocityTRVP2_ref_x_out, value_t *d_velocityTRVP2_ref_y_out, value_t *d_velocityTRVP2_ref_z_out,
						   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
						   value_t *d_position_ref_x_m1h_in, value_t *d_position_ref_y_m1h_in, value_t *d_position_ref_z_m1h_in,
						   value_t *d_position_ref_x_m2h_in, value_t *d_position_ref_y_m2h_in, value_t *d_position_ref_z_m2h_in,
						   value_t *d_position_ref_x_m3h_in, value_t *d_position_ref_y_m3h_in, value_t *d_position_ref_z_m3h_in)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);

#pragma unroll
	for (int i = 0; i < SUBSTANCE; i++)
	{
		d_velocityTRVP2_ref_x_out[idx*SUBSTANCE+i] =  TRVP2_right_sight(d_position_ref_x_in, d_position_ref_x_m1h_in ,d_position_ref_x_m2h_in,d_position_ref_x_m3h_in , idx*SUBSTANCE+i);
		d_velocityTRVP2_ref_y_out[idx*SUBSTANCE+i] =  TRVP2_right_sight(d_position_ref_y_in, d_position_ref_y_m1h_in ,d_position_ref_y_m2h_in,d_position_ref_y_m3h_in , idx*SUBSTANCE+i);
		d_velocityTRVP2_ref_z_out[idx*SUBSTANCE+i] =  TRVP2_right_sight(d_position_ref_z_in, d_position_ref_z_m1h_in ,d_position_ref_z_m2h_in,d_position_ref_z_m3h_in , idx*SUBSTANCE+i);
	}
}

/* Function: TRVP2_vel
 * time-reversible velocity predictor of molecular velocity in time t
 * order of predictor is k=2, paper DOI: 10.1021/ct200108g
 * history of molecular velocities is used for prediction
 *
 * Notes:
 * >calls device function TRVP2Vel_right_sight
 *
 * Parameter:
 * d_velocityTRVP2_x_out - output array of predicted molecular x velocities in time t
 * d_velocityTRVP2_y_out - output array of predicted molecular y velocities in time t
 * d_velocityTRVP2_z_out - output array of predicted molecular z velocities in time t
 * d_velocity_m05h_x_in	 - input array of molecular x velocities in time t-0.5*h
 * d_velocity_m05h_y_in	 - input array of molecular y velocities in time t-0.5*h
 * d_velocity_m05h_z_in	 - input array of molecular z velocities in time t-0.5*h
 * d_velocity_m15h_x_in	 - input array of molecular x velocities in time t-1.5*h
 * d_velocity_m15h_y_in	 - input array of molecular y velocities in time t-1.5*h
 * d_velocity_m15h_z_in	 - input array of molecular z velocities in time t-1.5*h
 * d_velocity_m25h_x_in	 - input array of molecular x velocities in time t-2.5*h
 * d_velocity_m25h_y_in	 - input array of molecular y velocities in time t-2.5*h
 * d_velocity_m25h_z_in	 - input array of molecular z velocities in time t-2.5*h
 */
__global__ void TRVP2_vel( value_t *d_velocityTRVP2_x_out, value_t *d_velocityTRVP2_y_out, value_t *d_velocityTRVP2_z_out,
					   value_t *d_velocity_m05h_x_in, value_t *d_velocity_m05h_y_in, value_t *d_velocity_m05h_z_in,
					   value_t *d_velocity_m15h_x_in, value_t *d_velocity_m15h_y_in, value_t *d_velocity_m15h_z_in,
					   value_t *d_velocity_m25h_x_in, value_t *d_velocity_m25h_y_in, value_t *d_velocity_m25h_z_in)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);

	d_velocityTRVP2_x_out[idx] =  TRVP2Vel_right_sight(d_velocity_m05h_x_in, d_velocity_m15h_x_in,d_velocity_m25h_x_in,idx);
	d_velocityTRVP2_y_out[idx] =  TRVP2Vel_right_sight(d_velocity_m05h_y_in, d_velocity_m15h_y_in,d_velocity_m25h_y_in,idx);
	d_velocityTRVP2_z_out[idx] =  TRVP2Vel_right_sight(d_velocity_m05h_z_in, d_velocity_m15h_z_in,d_velocity_m25h_z_in,idx);
}

/* Function: TRVP2_ref_vel
 * time-reversible velocity predictor of atom relative velocity in time t
 * order of predictor is k=2, paper DOI: 10.1021/ct200108g
 * history of atom velocities is used for prediction
 *
 * Notes:
 * >calls device function TRVP2Vel_right_sight
 *
 * Parameter:
 * d_velocityTRVP2_ref_x_out - output array of predicted atom x velocities in time t
 * d_velocityTRVP2_ref_y_out - output array of predicted atom y velocities in time t
 * d_velocityTRVP2_ref_z_out - output array of predicted atom z velocities in time t
 * d_velocity_m05h_ref_x_in	 - input array of atom x velocities in time t-0.5*h
 * d_velocity_m05h_ref_y_in	 - input array of atom y velocities in time t-0.5*h
 * d_velocity_m05h_ref_z_in	 - input array of atom z velocities in time t-0.5*h
 * d_velocity_m15h_ref_x_in	 - input array of atom x velocities in time t-1.5*h
 * d_velocity_m15h_ref_y_in	 - input array of atom y velocities in time t-1.5*h
 * d_velocity_m15h_ref_z_in	 - input array of atom z velocities in time t-1.5*h
 * d_velocity_m25h_ref_x_in	 - input array of atom x velocities in time t-2.5*h
 * d_velocity_m25h_ref_y_in	 - input array of atom y velocities in time t-2.5*h
 * d_velocity_m25h_ref_z_in	 - input array of atom z velocities in time t-2.5*h
 */
__global__ void TRVP2_ref_vel( value_t *d_velocityTRVP2_ref_x_out, value_t *d_velocityTRVP2_ref_y_out, value_t *d_velocityTRVP2_ref_z_out,
		 	 	 	 	 	   value_t *d_velocity_m05h_ref_x_in, value_t *d_velocity_m05h_ref_y_in, value_t *d_velocity_m05h_ref_z_in,
							   value_t *d_velocity_m15h_ref_x_in, value_t *d_velocity_m15h_ref_y_in, value_t *d_velocity_m15h_ref_z_in,
							   value_t *d_velocity_m25h_ref_x_in, value_t *d_velocity_m25h_ref_y_in, value_t *d_velocity_m25h_ref_z_in)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);

#pragma unroll
	for (int i = 0; i < SUBSTANCE; i++)
	{
		d_velocityTRVP2_ref_x_out[idx*SUBSTANCE+i] =  TRVP2Vel_right_sight(d_velocity_m05h_ref_x_in, d_velocity_m15h_ref_x_in,d_velocity_m25h_ref_x_in,idx*SUBSTANCE+i);
		d_velocityTRVP2_ref_y_out[idx*SUBSTANCE+i] =  TRVP2Vel_right_sight(d_velocity_m05h_ref_y_in, d_velocity_m15h_ref_y_in,d_velocity_m25h_ref_y_in,idx*SUBSTANCE+i);
		d_velocityTRVP2_ref_z_out[idx*SUBSTANCE+i] =  TRVP2Vel_right_sight(d_velocity_m05h_ref_z_in, d_velocity_m15h_ref_z_in,d_velocity_m25h_ref_z_in,idx*SUBSTANCE+i);
	}
}

__global__ void kinetic_energy_initial_calculation(value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
												   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
												   value_t *d_ekin_in)
{
	int idx = (blockIdx.x*blockDim.x+ threadIdx.x);
	value_t tmp_ekin = 0.0;

#pragma unroll
	for (int i = 0; i < SUBSTANCE; i++)
	{
		//kinetic energy per molecule part 1 calculation
		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
				    +(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
				    +(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * d_mol_mass_template[i]*2;// multiplycation by 2 is just for same "scale" as usual
	}

	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy
}

// === ========= ===
// === SHAKE ALG ===
// === ========= ===

/* Function: shake_expansion
 * device function for preserving the bonds of bonded molecules via SHAKE iterative algorithm
 * modification of shake for expansion
 * was found that is NOT required for expansion!!!!
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
__device__ void shake_expansion(value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					  value_t *position_old,
					  value_t *rij,
					  int mol_ind_ref,
					  value_t m_lambda_tphh)
{
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
				rij_old[0] = position_old[  3*ii] - position_old[  3*jj];
				rij_old[1] = position_old[1+3*ii] - position_old[1+3*jj];
				rij_old[2] = position_old[2+3*ii] - position_old[2+3*jj];

				// alpha initial computation
				alpha_m = 1.28* // the bulgarian constant that enhances efficiency of shake
						 ( rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]
#if SUBSTANCE==NITROGEN
						   -1.08892776*1.08892776/(m_lambda_tphh*m_lambda_tphh))/
#elif SUBSTANCE==SPCE
				           -d_bond_template[ii]*d_bond_template[ii]/(m_lambda_tphh*m_lambda_tphh))/ //indexing is taken from the first atom that determines type of bond O-H, H-H
#else // beware this can lead to drift in bonds - espetially in case of lower shake iteration count
				           -rij_old[0]*rij_old[0] - rij_old[1]*rij_old[1] - rij_old[2]*rij_old[2])/
#endif
				          (rij[0]*rij_old[0] + rij[1]*rij_old[1] + rij[2]*rij_old[2]);
				// the mass addition with leftover 1/2.0
				alpha_m /= 2.0*(d_mol_mass_template[ii] + d_mol_mass_template[jj]);

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

/* Function: kolafa_expansion
 * kernel function for iteration step of the leap-frog scheme with box rescaling using TRVP
 *
 * updates positions, velocities and kinetic energies
 *
 * direct call to the <shake> device function
 *
 *.box file required
 *
 * with correction for the numeric error
 *
 * Notes:
 * >the implementation is in form of single loop for number of atoms in molecule
 * >kolafa centre of mass algorithm from 15.06.2018
 * >operates with single thread per molecule
 *
 * beware kernel does not have overflow control
 * has to be controlled externally with launch configuration
 * (it is the case with number of molecules as power of 2)
 *
 * Parameter:
 * d_position_x_in     		- pointer to x molecular position array
 * d_position_y_in     		- pointer to y molecular position array
 * d_position_z_in 	   		- pointer to z molecular position array
 * d_position_ref_x_in 		- pointer to x atom position array
 * d_position_ref_y_in 		- pointer to y atom position array
 * d_position_ref_z_in 		- pointer to z atom position array
 * d_velocity_x_in     		- pointer to x molecular velocity array
 * d_velocity_y_in     		- pointer to y molecular velocity array
 * d_velocity_z_in 	   		- pointer to z molecular velocity array
 * d_velocity_ref_x_in 		- pointer to x atom velocity array
 * d_velocity_ref_y_in 		- pointer to y atom velocity array
 * d_velocity_ref_z_in 		- pointer to z atom velocity array
 * d_force_ref_x_in	   		- pointer to x atom force array
 * d_force_ref_y_in    		- pointer to y atom force array
 * d_force_ref_z_in    		- pointer to z atom force array
 * d_velocity_x_TRVP_in		- pointer to x molecular velocity array predicted by TRVP
 * d_velocity_y_TRVP_in		- pointer to y molecular velocity array predicted by TRVP
 * d_velocity_z_TRVP_in		- pointer to z molecular velocity array predicted by TRVP
 * d_velocity_ref_x_TRVP_in	- pointer to x atom velocity array predicted by TRVP
 * d_velocity_ref_y_TRVP_in	- pointer to y atom velocity array predicted by TRVP
 * d_velocity_ref_z_TRVP_in	- pointer to z atom velocity array predicted by TRVP
 * d_ekin_in		   		- pointer to array of atom energy per atom
 * m_lambda_tphh			- pointer to value of box volume ratio lambda
 * m_vft					- pointer to value volume change logarithm vf
 * m_lxyz					- pointer to array of box edge lengths
 */

__global__ void kolafa_expansion_NoMM(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					   value_t *d_velocity_x_TRVP_in, value_t *d_velocity_y_TRVP_in,
					   value_t *d_velocity_z_TRVP_in, value_t *d_velocity_ref_x_TRVP_in, value_t *d_velocity_ref_y_TRVP_in,
					   value_t *d_velocity_ref_z_TRVP_in,
					   value_t *d_ekin_in, 
					   value_t *d_lx_buffer, unsigned long lx_position)
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
	value_t lambda=calc_lambda_buffer(d_lx_buffer,lx_position);
	value_t vf=calc_vf_buffer(d_lx_buffer,lx_position);

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
	tmp_a[0] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[1] *= d_mol_mass_template[SUBSTANCE];
	tmp_a[2] *= d_mol_mass_template[SUBSTANCE];

	// molecular velocity update - modified for EXPANSION
	d_velocity_x_in[idx] += d_dt2*tmp_a[0]-d_dt2*(vf*d_velocity_x_TRVP_in[idx]);//-(m_vft*d_velocity_x_TRVP_in[idx]) added for expansion
	d_velocity_y_in[idx] += d_dt2*tmp_a[1]-d_dt2*(vf*d_velocity_y_TRVP_in[idx]);//-||-
	d_velocity_z_in[idx] += d_dt2*tmp_a[2]-d_dt2*(vf*d_velocity_z_TRVP_in[idx]);//-||-

	//molecular position update
	d_position_x_in[idx] += d_velocity_x_in[idx];
	d_position_y_in[idx] += d_velocity_y_in[idx];
	d_position_z_in[idx] += d_velocity_z_in[idx];

	//molecular position update
	d_position_x_in[idx] *= lambda;
	d_position_y_in[idx] *= lambda;
	d_position_z_in[idx] *= lambda;

	// PCB crop per molecule
	d_position_x_in[idx] -= d_lx_buffer[2*lx_position]*floor(d_position_x_in[idx]/d_lx_buffer[2*lx_position]);
	d_position_y_in[idx] -= d_lx_buffer[2*lx_position]*floor(d_position_y_in[idx]/d_lx_buffer[2*lx_position]);
	d_position_z_in[idx] -= d_lx_buffer[2*lx_position]*floor(d_position_z_in[idx]/d_lx_buffer[2*lx_position]);


#pragma unroll
	for (i = 0; i < SUBSTANCE; i++)
	{
		// update atom velocity - modified for EXPANSION
		d_velocity_ref_x_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_x_in[idx*SUBSTANCE+i]/d_mol_mass_template[i]- tmp_a[0])-d_dt2*(vf*d_velocity_ref_x_TRVP_in[idx*SUBSTANCE+i]);//-(m_vft*d_velocity_ref_x_TRVP_in[idx*SUBSTANCE+i]) added for expansion
		d_velocity_ref_y_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_y_in[idx*SUBSTANCE+i]/d_mol_mass_template[i]- tmp_a[1])-d_dt2*(vf*d_velocity_ref_y_TRVP_in[idx*SUBSTANCE+i]);//-||-
		d_velocity_ref_z_in[idx*SUBSTANCE+i] += d_dt2*(d_force_ref_z_in[idx*SUBSTANCE+i]/d_mol_mass_template[i]- tmp_a[2])-d_dt2*(vf*d_velocity_ref_z_TRVP_in[idx*SUBSTANCE+i]);//-||-

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
	shake_expansion(d_position_ref_x_in, d_position_ref_y_in, d_position_ref_z_in,
			  		tmp_pos, tmp_a,
			  		idx*SUBSTANCE, 
			  		lambda);

#pragma unroll // valid only when shaking - only kinetic energy
	for (i = 0; i < SUBSTANCE; i++)
	{

//		// velocity recalculation after shake
		d_velocity_ref_x_in[idx*SUBSTANCE+i] = d_position_ref_x_in[idx*SUBSTANCE+i] - tmp_pos[  3*i];
		d_velocity_ref_y_in[idx*SUBSTANCE+i] = d_position_ref_y_in[idx*SUBSTANCE+i] - tmp_pos[1+3*i];
		d_velocity_ref_z_in[idx*SUBSTANCE+i] = d_position_ref_z_in[idx*SUBSTANCE+i] - tmp_pos[2+3*i];
//		//kinetic energy per molecule part 2 calculation
		tmp_ekin +=( (d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])*(d_velocity_x_in[idx] + d_velocity_ref_x_in[idx*SUBSTANCE+i])
					+(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])*(d_velocity_y_in[idx] + d_velocity_ref_y_in[idx*SUBSTANCE+i])
					+(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i])*(d_velocity_z_in[idx] + d_velocity_ref_z_in[idx*SUBSTANCE+i]) )
				   * d_mol_mass_template[i];

	}
	d_ekin_in[idx] = tmp_ekin; // the 1/4.0 is performed on processor at final summation of kinetic energy

}
