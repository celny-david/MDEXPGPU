/*
 * cuda_variables.h
 *
 *  Created on: 13. Jan 2020
 *      Author: dixiw
 */

#ifndef CUDA_VARIABLES_H_
#define CUDA_VARIABLES_H_


#include "value.h"
#include "cuda_definition.h"

extern value_t *d_position_x;
extern value_t *d_position_y;
extern value_t *d_position_z;
extern value_t *d_position_ref_x;
extern value_t *d_position_ref_y;
extern value_t *d_position_ref_z;
extern value_t *d_velocity_x;
extern value_t *d_velocity_y;
extern value_t *d_velocity_z;
extern value_t *d_velocity_ref_x;
extern value_t *d_velocity_ref_y;
extern value_t *d_velocity_ref_z;
extern value_t *d_force_ref_x;
extern value_t *d_force_ref_y;
extern value_t *d_force_ref_z;
extern value_t *d_ekin_tmp;
extern value_t *d_epot_tmp;

#ifdef EXPANSION
extern value_t *d_velocity_x_TRVP_out;
extern value_t *d_velocity_y_TRVP_out;
extern value_t *d_velocity_z_TRVP_out;
extern value_t *d_velocity_ref_x_TRVP_out;
extern value_t *d_velocity_ref_y_TRVP_out;
extern value_t *d_velocity_ref_z_TRVP_out;

extern value_t *d_velocity_x_m3l2h;
extern value_t *d_velocity_y_m3l2h;
extern value_t *d_velocity_z_m3l2h;

extern value_t *d_velocity_x_m5l2h;
extern value_t *d_velocity_y_m5l2h;
extern value_t *d_velocity_z_m5l2h;

extern value_t *d_velocity_ref_x_m3l2h;
extern value_t *d_velocity_ref_y_m3l2h;
extern value_t *d_velocity_ref_z_m3l2h;

extern value_t *d_velocity_ref_x_m5l2h;
extern value_t *d_velocity_ref_y_m5l2h;
extern value_t *d_velocity_ref_z_m5l2h;

extern value_t *d_lx_buffer;
#endif

extern ushort2 *d_verlet; // new implementation of dynamic verlet list
extern unsigned int *d_verlet_occupancy; //
extern unsigned int *d_verlet_index; // access indices containing staritng points into neighbor list

extern __constant__ unsigned int d_atom_count;
extern __constant__ unsigned int d_molecule_count;
extern __constant__ unsigned int d_atom_per_molecule;
extern __constant__ value_t d_cutoff;
extern __constant__ value_t d_cutoff_elstat;
extern __constant__ value_t d_cutoff_verlet;
extern __constant__ value_t d_lx;
extern __constant__ value_t d_ly;
extern __constant__ value_t d_lz;
extern __constant__ value_t d_dt;
extern __constant__ value_t d_dt2;
extern __constant__ value_t d_mol_mass_template[SUBSTANCE+1];

#if SUBSTANCE==SPCE
extern __constant__ value_t d_bond_template[2];
#endif
extern __constant__ unsigned int d_n_shake;

#if ELPOT_TYPE >0
extern __constant__ signed char d_mol_charge_template[SUBSTANCE];
#endif

#endif /* CUDA_VARIABLES_H_ */
