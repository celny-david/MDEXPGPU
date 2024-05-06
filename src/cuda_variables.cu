/*
 * cuda_variables.cu
 *
 *  Created on: 13. Jan 2020
 *      Author: dixiw
 */


#include "value.h"
#include "cuda_definition.h"

// == GPU residing variables ==
// for external manage-ability (data_2device) declared as global

/* Variables: GPU global memory core data

	d_position_x 		- molecular centre of mass x positions
	d_position_y 		- molecular centre of mass y positions
	d_position_z 		- molecular centre of mass z positions
	d_position_ref_x	- atomal x positions
	d_position_ref_y	- atomal y positions
	d_position_ref_z	- atomal z positions
	d_velocity_x 		- molecular centre of mass x velocities
	d_velocity_y 		- molecular centre of mass y velocities
	d_velocity_z 		- molecular centre of mass z velocities
	d_velocity_ref_x 	- atomal x velocities
	d_velocity_ref_y 	- atomal y velocities
	d_velocity_ref_z 	- atomal z velocities
	d_force_ref_x 		- atomal x forces (molecular forces are not needed)
	d_force_ref_y 		- atomal y forces (molecular forces are not needed)
	d_force_ref_z 		- atomal z forces (molecular forces are not needed)
	d_ekin_tmp 			- kinetic energy per individual atom (for total perform sum)
	d_epot_tmp 			- potential energy per individual atom (for total perform sum)
*/

value_t *d_position_x;
value_t *d_position_y;
value_t *d_position_z;
value_t *d_position_ref_x;
value_t *d_position_ref_y;
value_t *d_position_ref_z;
value_t *d_velocity_x;
value_t *d_velocity_y;
value_t *d_velocity_z;
value_t *d_velocity_ref_x;
value_t *d_velocity_ref_y;
value_t *d_velocity_ref_z;

value_t *d_force_ref_x;
value_t *d_force_ref_y;
value_t *d_force_ref_z;

value_t *d_ekin_tmp;
value_t *d_epot_tmp;

/* Variables: GPU global memory expansion data

	d_velocity_x_TRVP_out 		- molecular centre of mass x positions
	d_velocity_y_TRVP_out 		- molecular centre of mass y positions
	d_velocity_z_TRVP_out 		- molecular centre of mass z positions
	d_velocity_ref_x_TRVP_out	- atomal x positions
	d_velocity_ref_y_TRVP_out	- atomal y positions
	d_velocity_ref_z_TRVP_out	- atomal z positions
	d_velocity_x_m3l2h 			- molecular centre of mass x velocities
	d_velocity_y_m3l2h 			- molecular centre of mass y velocities
	d_velocity_z_m3l2h 			- molecular centre of mass z velocities
	d_velocity_x_m5l2h 			- atomal x velocities
	d_velocity_y_m5l2h 			- atomal y velocities
	d_velocity_z_m5l2h 			- atomal z velocities
	d_velocity_ref_x_m3l2h 		- atomal x forces (molecular forces are not needed)
	d_velocity_ref_y_m3l2h 		- atomal y forces (molecular forces are not needed)
	d_velocity_ref_z_m3l2h 		- atomal z forces (molecular forces are not needed)
	d_velocity_ref_x_m5l2h 		- kinetic energy per individual atom (for total perform sum)
	d_velocity_ref_y_m5l2h 		- kinetic energy per individual atom (for total perform sum)
	d_velocity_ref_z_m5l2h 		- potential energy per individual atom (for total perform sum)
	d_lx_buffer 				- potential energy per individual atom (for total perform sum)
*/
#ifdef EXPANSION

value_t *d_velocity_x_TRVP_out;
value_t *d_velocity_y_TRVP_out;
value_t *d_velocity_z_TRVP_out;
value_t *d_velocity_ref_x_TRVP_out;
value_t *d_velocity_ref_y_TRVP_out;
value_t *d_velocity_ref_z_TRVP_out;

value_t *d_velocity_x_m3l2h;
value_t *d_velocity_y_m3l2h;
value_t *d_velocity_z_m3l2h;

value_t *d_velocity_x_m5l2h;
value_t *d_velocity_y_m5l2h;
value_t *d_velocity_z_m5l2h;

value_t *d_velocity_ref_x_m3l2h;
value_t *d_velocity_ref_y_m3l2h;
value_t *d_velocity_ref_z_m3l2h;

value_t *d_velocity_ref_x_m5l2h;
value_t *d_velocity_ref_y_m5l2h;
value_t *d_velocity_ref_z_m5l2h;

value_t *d_lx_buffer;
#endif

/* Variables: GPU global memory verlet data

	d_verlet 			- ushort2 pointer - the dynamic verlet list (stores the indice pairs)
	d_verlet_occupancy 	- uint pointer - number of neighbours of each atom i.e. occupied space in verlet lists
	d_verlet_index 		- uint pointer - access indices containing staritng points into neighbor list
*/
ushort2 *d_verlet; // neighbouring atoms indices
unsigned int *d_verlet_occupancy; // number of neighbours of each atom i.e. occupied space in verlet lists
unsigned int *d_verlet_index; // access indices containing staritng points into neighbor list

// == GPU global/constant parameters ==
/* Constants: GPU constant memory data

	d_atom_count 			- computational total atom count [#]
	d_molecule_count 		- computational total molecule count [#]
	d_atom_per_molecule 	- atoms per molecule count [#]
	d_cutoff 				- standard cutoff [Angstrom]
	d_cutoff_elstat 		- electrostatic cutoff [Angstrom]
	d_cutoff_verlet 		- verlet list construction cutoff [Angstrom]
	d_lx 					- x dimension of computational volume [Angstrom]
	d_ly 					- y dimension of computational volume [Angstrom]
	d_lz 					- z dimension of computational volume [Angstrom]
	d_dt 					- time step [picosecond]
	d_dt2 					- time step squared [picosecond]
	d_mol_mass_template 	- array (size SUBSTANCE+1) representing molar masses of atoms in molecule - !!! VALID FOR SPCE ONLY [macsimus mass units]
	d_bond_template 		- array representing bonds within the molecule
	d_n_shake 				- number of shake iteration [#]
	d_mol_charge_template 	- array representing molar charges of atoms in molecule [charge number]
*/

__constant__ unsigned int d_atom_count;
__constant__ unsigned int d_molecule_count;
__constant__ unsigned int d_atom_per_molecule;
__constant__ value_t d_cutoff;
__constant__ value_t d_cutoff_elstat;
__constant__ value_t d_cutoff_verlet;
__constant__ value_t d_lx;
__constant__ value_t d_ly;
__constant__ value_t d_lz;
__constant__ value_t d_dt;
__constant__ value_t d_dt2;
__constant__ value_t d_mol_mass_template[SUBSTANCE+1];

#if SUBSTANCE==SPCE
__constant__ value_t d_bond_template[2];
#endif
__constant__ unsigned int d_n_shake;

#if ELPOT_TYPE >0
__constant__ signed char d_mol_charge_template[SUBSTANCE];
#endif
