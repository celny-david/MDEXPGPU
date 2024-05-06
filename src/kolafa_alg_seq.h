/*
 * forcefield.h
 *
 *  Created on: 19 Dec 2019
 *      Author: dixiw
 */

#ifndef KOLAFA_ALG_SEQ_H_
#define KOLAFA_ALG_SEQ_H_

#include "value.h"

void kolafa_adjust_seq(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_x_in, value_t *d_force_y_in, value_t *d_force_z_in,
					   value_t *d_ekin_in,
					   unsigned int idx);

void get_mh_positions_seq(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
						  value_t *d_position_x_mh_in, value_t *d_position_y_mh_in, value_t *d_position_z_mh_in,
						  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
						  value_t *d_position_ref_x_mh_in, value_t *d_position_ref_y_mh_in, value_t *d_position_ref_z_mh_in,
						  value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
						  value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
						  value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
						  unsigned int idx);

void verlet_step_seq_argon( value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							value_t *d_position_x_mh_in, value_t *d_position_y_mh_in, value_t *d_position_z_mh_in,
							value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
							value_t *d_position_ref_x_mh_in, value_t *d_position_ref_y_mh_in, value_t *d_position_ref_z_mh_in,
							value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
							value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
							value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							value_t *d_ekin_in,
							unsigned int idx);

void verlet_step_seq_shaked( value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							 value_t *d_position_x_mh_in, value_t *d_position_y_mh_in, value_t *d_position_z_mh_in,
							 value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
							 value_t *d_position_ref_x_mh_in, value_t *d_position_ref_y_mh_in, value_t *d_position_ref_z_mh_in,
							 value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
							 value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
							 value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
							 value_t *d_ekin_in,
							 unsigned int idx);

#endif /* KOLAFA_ALG_SEQ_H_ */
