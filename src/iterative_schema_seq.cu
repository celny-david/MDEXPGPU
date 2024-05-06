/*
 *  iterative_schema.cu
 *	iterative schemas collection
 *	{leap_frog, }
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#include "value.h"
#include "triplet.h"
#include "output_property.h"
#include "iterative_property.h"
#include "system_property.h"
#include "verlet_seq.h"
#include "forcefield_seq.h"
#include "kolafa_alg_seq.h"
#include "output.h"

/* Function: leap_frog
 * the procedure responsible for the time iteration
 *
 * done according to the leap-frog scheme
 *
 * contains the Verlet reconstruction force calculation and update and data output
 *
 * verlet reconstructed with <verlet_reconstruct>
 * force is calculate with <force_calculation>
 * the positions are updated with <kolafa> & <kolafa_c>
 *
 * Notes:
 * >the implementation is in form of single loop over timesteps
 * >! untested for high iteration count
 * > the block, thread dimensions have to correspond to the solved task size
 *
 * Parameter:
 * n_block     		- unsigned int of how many block to use
 * thread_per_block - unsigned int of how many threads per block to use
 */
void leap_frog_seq(triplet *mol_pos_tmp, triplet *atm_pos_tmp, triplet *mol_vel_tmp, triplet *atm_vel_tmp, triplet *atm_fcs_tmp)
{
/*
 * it is made not responsible for allocation/deallocation of device memory
 * this task initial move of data is left in main
 *
 * * leap frog schema: v(t+dt/2) = v(t) + a(t)*dt/2
 * * 				   r(t+dt)   = r(t) + v(t+dt/2)*dt
 * * 				   a(t+dt)   = F( r(t+dt) )/m
 * * 				   v(t+dt)   = v(t+dt/2) + a(t+dt)*dt/2
 */
	uint j; // for loop counter
	double t = 0.0; // [picosecond]
	unsigned long max_step = ceil(itep.total_time/itep.dt); // number of steps
	unsigned int verlet_size = 0; // to start things up 

	// tmp variables for energy logging
	value_t *ekin_tmp = (value_t*) malloc(sysp.atom_count_aligned*sizeof(value_t));
	value_t *epot_tmp = (value_t*) malloc(sysp.atom_count_aligned*sizeof(value_t));
	unsigned int *verlet_occupancy_tmp = (unsigned int*) malloc(sysp.molecule_count_aligned*sizeof(unsigned int));
	ushort2 *verlet = (ushort2 *) malloc(10*sizeof(value_t)); // just for initial startup

	// ===== ITERATION =====
	for (itep.step = 0; itep.step <= max_step; itep.step++)
	{
		if(itep.step % itep.verlet_refresh == 0)
		{
			verlet_size = verlet_continual_refresh(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z, verlet_size, &verlet);
		}

		// force calculation DYNAMIC verlet list
		if (verlet_size != 0)
		{
			force_calculation_seq(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
								  atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
								  atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
								  epot_tmp,
								  verlet_size,
								  verlet);
		}

		// leap_forg iteration		
		for (j = 0; j < sysp.molecule_count_aligned; ++j)
		{
			kolafa_adjust_seq(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
							atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
							mol_vel_tmp->x, mol_vel_tmp->y, mol_vel_tmp->z,
							atm_vel_tmp->x, atm_vel_tmp->y, atm_vel_tmp->z,
							atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
							ekin_tmp,
							j); 
		}
		
		if(itep.printout==true)
		{
			print_into_files(mol_pos_tmp,atm_pos_tmp,
							 mol_vel_tmp,atm_vel_tmp,
							 atm_fcs_tmp,
							 ekin_tmp,epot_tmp);
		}

		t += itep.dt;
	}

	print_into_files(mol_pos_tmp,atm_pos_tmp,
					 mol_vel_tmp,atm_vel_tmp,
					 atm_fcs_tmp,
					 ekin_tmp,epot_tmp);
}

void verlet_integration_seq(triplet *mol_pos_tmp, triplet *atm_pos_tmp, triplet *mol_vel_tmp, triplet *atm_vel_tmp, triplet *atm_fcs_tmp)
{
	uint j; // for loop counter
	double t = 0.0; // [picosecond]
	unsigned long max_step = ceil(itep.total_time/itep.dt); // number of steps

	unsigned int verlet_size = 0; // to start things up

	// tmp variables for energy logging
	value_t *ekin_tmp = (value_t*) malloc(sysp.atom_count_aligned*sizeof(value_t));
	value_t *epot_tmp = (value_t*) malloc(sysp.atom_count_aligned*sizeof(value_t));
	unsigned int *verlet_occupancy_tmp = (unsigned int*) malloc(sysp.molecule_count_aligned*sizeof(unsigned int));
	ushort2 *verlet = (ushort2 *) malloc(10*sizeof(value_t)); // just for initial startup

	//alocation of positions in t-h
	triplet *mol_pos_mh_tmp = (triplet*) malloc(sizeof(triplet));
	triplet *atm_pos_mh_tmp = (triplet*) malloc(sizeof(triplet));
	triplet_alloc(mol_pos_mh_tmp,sysp.molecule_count_aligned);
	triplet_alloc(atm_pos_mh_tmp,sysp.atom_count_aligned);


	//obtaining positions in t-h
	verlet_size = verlet_continual_refresh(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z, verlet_size, &verlet);

	force_calculation_seq(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
						  atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
						  atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
						  epot_tmp,
						  verlet_size,
						  verlet);
	for (j = 0; j < sysp.molecule_count_aligned; ++j)
	{
		get_mh_positions_seq( mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
							  mol_pos_mh_tmp->x, mol_pos_mh_tmp->y, mol_pos_mh_tmp->z,
							  atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
							  atm_pos_mh_tmp->x, atm_pos_mh_tmp->y, atm_pos_mh_tmp->z,
							  mol_vel_tmp->x, mol_vel_tmp->y, mol_vel_tmp->z,
							  atm_vel_tmp->x, atm_vel_tmp->y, atm_vel_tmp->z,
							  atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
							  j);
//		printf("%f  %f \n",mol_pos_tmp->z[k],mol_pos_mh_tmp->z[k]); // DEBUG
	}

	// ===== ITERATION =====
	for (itep.step = 0; itep.step <= max_step; itep.step++)
	{
		if(itep.step % itep.verlet_refresh == 0)
		{
			verlet_size = verlet_continual_refresh(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z, verlet_size, &verlet);
		}

		if (verlet_size != 0)
		{
			force_calculation_seq(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
								  atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
								  atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
								  epot_tmp,
								  verlet_size,
								  verlet);
		}

		//integration step
		for (j = 0; j < sysp.molecule_count_aligned; ++j)
		{
#if SUBSTANCE == ARGON
			verlet_step_seq_argon(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
								  mol_pos_mh_tmp->x, mol_pos_mh_tmp->y, mol_pos_mh_tmp->z,
								  atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
								  atm_pos_mh_tmp->x, atm_pos_mh_tmp->y, atm_pos_mh_tmp->z,
								  mol_vel_tmp->x, mol_vel_tmp->y, mol_vel_tmp->z,
								  atm_vel_tmp->x, atm_vel_tmp->y, atm_vel_tmp->z,
								  atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
								  ekin_tmp,j);
#else
			verlet_step_seq_shaked(mol_pos_tmp->x, mol_pos_tmp->y, mol_pos_tmp->z,
								   mol_pos_mh_tmp->x, mol_pos_mh_tmp->y, mol_pos_mh_tmp->z,
								   atm_pos_tmp->x, atm_pos_tmp->y, atm_pos_tmp->z,
								   atm_pos_mh_tmp->x, atm_pos_mh_tmp->y, atm_pos_mh_tmp->z,
								   mol_vel_tmp->x, mol_vel_tmp->y, mol_vel_tmp->z,
								   atm_vel_tmp->x, atm_vel_tmp->y, atm_vel_tmp->z,
								   atm_fcs_tmp->x, atm_fcs_tmp->y, atm_fcs_tmp->z,
								   ekin_tmp,j);
#endif
		}

		if(itep.printout == true)
		{
			print_into_files(mol_pos_tmp,atm_pos_tmp,
							 mol_vel_tmp,atm_vel_tmp,
							 atm_fcs_tmp,
							 ekin_tmp,epot_tmp);
		}

		t += itep.dt;
	}

	print_into_files(mol_pos_tmp,atm_pos_tmp,
					 mol_vel_tmp,atm_vel_tmp,
					 atm_fcs_tmp,
					 ekin_tmp,epot_tmp);
	
	triplet_dealloc(mol_pos_mh_tmp);
	triplet_dealloc(atm_pos_mh_tmp);
}
