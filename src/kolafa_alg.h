/*
 * forcefield.h
 *
 *  Created on: 19 Dec 2019
 *      Author: dixiw
 */

#ifndef KOLAFA_ALG_H_
#define KOLAFA_ALG_H_

#include "value.h"

__global__ void kolafa_adjust(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
							  value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
							  value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
							  value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
							  value_t *d_force_x_in, value_t *d_force_y_in, value_t *d_force_z_in,
							  value_t *d_ekin_in);

__global__ void kolafa_adjust_correction(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
										 value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
										 value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
										 value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
										 value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
										 value_t *d_ekin_in);

#endif /* KOLAFA_ALG_H_ */
