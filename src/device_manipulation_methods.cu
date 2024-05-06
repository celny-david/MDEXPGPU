/*
 *  device_manipulation_methods.cu
 *	collects all the methods that operate on Host-Device layer
 *	{ }
 *  Created on: 19. Dec 2019
 *      Author: dixiw
 */

#include "err_check.h"
#include "value.h"
#include "triplet.h"
#include "verlet.h"
#include "cuda_variables.h"
#include "cuda_definition.h"
#include "system_property.h"
#include "iterative_property.h"

// == CPU<->GPU data handling ==

/* Function: prop2device
 * procedure to transfer the constant memory data on the device
 *
 * intended as an interface for the forcefield function which operates with device <forcefield_prop2device>
 *
 */
void prop2device()
{
/*
 * wrapper due to the requirement of file locality for __constant__ using functions
 */
	value_t dt2 = itep.dt*itep.dt;
#if SUBSTANCE== SPCE // this implies ELPOt_TYPE = 20
	signed char mol_charge_template[SUBSTANCE] = {-2,1,1};
	value_t bond_template[2] = {1.0,2.0*sin(0.5*acos(-1.0/3.0))}; // the bond length between O-H, H-H
	cudaSafeCall( cudaMemcpyToSymbol(d_mol_charge_template,&mol_charge_template,(SUBSTANCE)*sizeof(signed char)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_bond_template,&bond_template,2*sizeof(value_t)) );
#endif

	cudaSafeCall( cudaMemcpyToSymbol(d_dt,&itep.dt,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_dt2,&dt2,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_mol_mass_template,&sysp.mol_mass_template,(SUBSTANCE+1)*sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_n_shake,&itep.n_shake,sizeof(unsigned int)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_lx,&sysp.lx,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_ly,&sysp.ly,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_lz,&sysp.lz,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_cutoff,&sysp.cutoff,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_cutoff_elstat,&sysp.cutoff_elstat,sizeof(value_t)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_atom_count,&sysp.atom_count_aligned,sizeof(unsigned int)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_molecule_count,&sysp.molecule_count_aligned,sizeof(unsigned int)) );
	cudaSafeCall( cudaMemcpyToSymbol(d_atom_per_molecule,&sysp.atom_per_molecule,sizeof(unsigned int)) );

	verlet_const_memory_copy();
}

/* Function: data_device_alloc
 * procedure to alocate data on the device
 *
 * Notes:
 * >! emits and catch cudaError bound with improper alloc
 *
 */
void data_device_alloc()
{
/*
 * malloc the data on device
 * wrapper intended to be callable from main/scheme
 */
	cudaSafeCall( cudaMalloc( (void**)&d_position_x, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_position_y, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_position_z, sysp.molecule_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_position_ref_x, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_position_ref_y, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_position_ref_z, sysp.atom_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_x, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_y, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_z, sysp.molecule_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_x, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_y, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_z, sysp.atom_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_force_ref_x, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_force_ref_y, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_force_ref_z, sysp.atom_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_epot_tmp, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_ekin_tmp, sysp.molecule_count_aligned * sizeof(value_t)) );

	// TODO NOTE - this is ment to cause error in case there is not enough space on the GPU to hold maximal size of verlet list
	verlet_device_alloc(sysp.molecule_count_aligned*sysp.molecule_count_aligned); // TODO this is inapropriate set up correctly

#ifdef EXPANSION
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_x_TRVP_out, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_y_TRVP_out, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_z_TRVP_out, sysp.molecule_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_x_TRVP_out, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_y_TRVP_out, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_z_TRVP_out, sysp.atom_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_x_m3l2h, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_y_m3l2h, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_z_m3l2h, sysp.molecule_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_x_m5l2h, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_y_m5l2h, sysp.molecule_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_z_m5l2h, sysp.molecule_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_x_m3l2h, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_y_m3l2h, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_z_m3l2h, sysp.atom_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_x_m5l2h, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_y_m5l2h, sysp.atom_count_aligned * sizeof(value_t)) );
	cudaSafeCall( cudaMalloc( (void**)&d_velocity_ref_z_m5l2h, sysp.atom_count_aligned * sizeof(value_t)) );

	cudaSafeCall( cudaMalloc( (void**)&d_lx_buffer, 2*(ceil(itep.total_time/itep.dt)+1) * sizeof(value_t)) );

#endif
}

/* Function: data_device_free
 * procedure to dealocate data on the device
 *
 * Notes:
 * >! emits and catch cudaError bound with improper free
 *
 */
void data_device_free()
{
/*
 * free the data on device
 * wrapper intended to be callable from main/scheme
 */
	cudaSafeCall( cudaFree(d_position_x) );
	cudaSafeCall( cudaFree(d_position_y) );
	cudaSafeCall( cudaFree(d_position_z) );
	cudaSafeCall( cudaFree(d_position_ref_x) );
	cudaSafeCall( cudaFree(d_position_ref_y) );
	cudaSafeCall( cudaFree(d_position_ref_z) );
	cudaSafeCall( cudaFree(d_velocity_x) );
	cudaSafeCall( cudaFree(d_velocity_y) );
	cudaSafeCall( cudaFree(d_velocity_z) );
	cudaSafeCall( cudaFree(d_velocity_ref_x) );
	cudaSafeCall( cudaFree(d_velocity_ref_y) );
	cudaSafeCall( cudaFree(d_velocity_ref_z) );

	cudaSafeCall( cudaFree(d_force_ref_x) );
	cudaSafeCall( cudaFree(d_force_ref_y) );
	cudaSafeCall( cudaFree(d_force_ref_z) );

	cudaSafeCall( cudaFree(d_epot_tmp) );
	cudaSafeCall( cudaFree(d_ekin_tmp) );

	verlet_device_free();

#ifdef EXPANSION
	cudaFree(d_velocity_x_TRVP_out);
	cudaFree(d_velocity_y_TRVP_out);
	cudaFree(d_velocity_z_TRVP_out);

	cudaFree(d_velocity_ref_x_TRVP_out);
	cudaFree(d_velocity_ref_y_TRVP_out);
	cudaFree(d_velocity_ref_z_TRVP_out);

	cudaFree(d_velocity_x_m3l2h);
	cudaFree(d_velocity_y_m3l2h);
	cudaFree(d_velocity_z_m3l2h);

	cudaFree(d_velocity_x_m5l2h);
	cudaFree(d_velocity_y_m5l2h);
	cudaFree(d_velocity_z_m5l2h);

	cudaFree(d_velocity_ref_x_m3l2h);
	cudaFree(d_velocity_ref_y_m3l2h);
	cudaFree(d_velocity_ref_z_m3l2h);

	cudaFree(d_velocity_ref_x_m5l2h);
	cudaFree(d_velocity_ref_y_m5l2h);
	cudaFree(d_velocity_ref_z_m5l2h);

	cudaFree(d_lx_buffer);
#endif

}

/* Function: data_host2device
 * data transfer from host to device
 *
 * used to transfer the initial configuration data
 *
 * Notes:
 * >! emits and catch cudaError bound with improper transfer
 *
 * Parameter:
 * pos_in - pointer to the host molecular-size position triplet
 * pos_ref_in - pointer to the host atomal-size position triplet
 * vel_in - pointer to the host molecular-size velocity triplet
 * vel_ref_in - pointer to the host atomal-size velocity triplet
 * for_ref_in - pointer to the host atomal-size force triplet
 */
void data_host2device(triplet *pos_in, triplet *pos_ref_in, triplet *vel_in, triplet *vel_ref_in, triplet *for_ref_in)
{
/*
 * memcpy the initial data to device
 * wrapper intended to be callable from main function
 * move only the necessary - epot/ekin are not moved
 */
	cudaSafeCall( cudaMemcpy(d_position_x,pos_in->x,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_position_y,pos_in->y,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_position_z,pos_in->z,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );

	cudaSafeCall( cudaMemcpy(d_position_ref_x,pos_ref_in->x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_position_ref_y,pos_ref_in->y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_position_ref_z,pos_ref_in->z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );

	cudaSafeCall( cudaMemcpy(d_velocity_x,vel_in->x,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_y,vel_in->y,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_z,vel_in->z,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );

	cudaSafeCall( cudaMemcpy(d_velocity_ref_x,vel_ref_in->x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_ref_y,vel_ref_in->y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_ref_z,vel_ref_in->z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );

	cudaSafeCall( cudaMemcpy(d_force_ref_x,for_ref_in->x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_force_ref_y,for_ref_in->y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
	cudaSafeCall( cudaMemcpy(d_force_ref_z,for_ref_in->z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyHostToDevice) );
}

/* Function: data_device2host_position
 * specialised version of data transfer from device to host
 *
 * this is only used for transferring position x,y,z arrays for atoms and molecule centres
 *
 * Notes:
 * >! emits and catch cudaError bound with improper transfer
 *
 * Parameter:
 * mol_pos_in - pointer to the host molecular-size position triplet
 * atm_pos_in - pointer to the host atomal-size position triplet
 */
void data_device2host_position(triplet *mol_pos_in,triplet *atm_pos_in)
{
/*
 * memcpy all the data from device to host
 * wrapper intended to be callable from main function
 */
	cudaSafeCall( cudaMemcpy(mol_pos_in->x,d_position_x,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(mol_pos_in->y,d_position_y,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(mol_pos_in->z,d_position_z,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );

	cudaSafeCall( cudaMemcpy(atm_pos_in->x,d_position_ref_x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(atm_pos_in->y,d_position_ref_y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(atm_pos_in->z,d_position_ref_z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
}

/* Function: data_device2host_velocity
 * specialised version of data transfer from device to host
 *
 * this is only used for transferring velocity x,y,z arrays for atoms and molecule centres
 *
 * Notes:
 * >! emits and catch cudaError bound with improper transfer
 *
 * Parameter:
 * mol_vel_in - pointer to the host molecular-size velocity triplet
 * atm_vel_in - pointer to the host atomal-size velocity triplet
 */
void data_device2host_velocity(triplet *mol_vel_in,triplet *atm_vel_in)
{
/*
 * memcpy all the data from device to host
 * wrapper intended to be callable from main function
 */
	cudaSafeCall( cudaMemcpy(mol_vel_in->x,d_velocity_x,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(mol_vel_in->y,d_velocity_y,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(mol_vel_in->z,d_velocity_z,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );

	cudaSafeCall( cudaMemcpy(atm_vel_in->x,d_velocity_ref_x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(atm_vel_in->y,d_velocity_ref_y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(atm_vel_in->z,d_velocity_ref_z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
}

/* Function: data_device2host_force
 * specialised version of data transfer from device to host
 *
 * this is only used for transferring force x,y,z array
 *
 * Notes:
 * >! emits and catch cudaError bound with improper transfer
 *
 * Parameter:
 * for_in - pointer to the host force triplet
 */
void data_device2host_force(triplet *for_in)
{
/*
 * memcpy all the data from device to host
 * wrapper intended to be callable from main function
 */
	cudaSafeCall( cudaMemcpy(for_in->x,d_force_ref_x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(for_in->y,d_force_ref_y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(for_in->z,d_force_ref_z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
}

/* Function: data_device2host_energy
 * specialised version of data transfer from device to host
 *
 * this is only used for transferring partial energies vectors
 *
 * Notes:
 * >! emits and catch cudaError bound with improper transfer
 *
 * Parameter:
 * ekin_in - pointer to the host kinetic energy array of appropriate size
 * epot_in - pointer to the host potential energy array of appropriate size
 */
void data_device2host_energy(value_t *ekin_in, value_t *epot_in)
{
/*
 * memcpy the ekin/epot vecotor from device to host
 * wrapper intended to be callable from main function
 */
	cudaSafeCall( cudaMemcpy(ekin_in,d_ekin_tmp,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
	cudaSafeCall( cudaMemcpy(epot_in,d_epot_tmp,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToHost) );
}
