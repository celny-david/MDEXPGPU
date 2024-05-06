/*
 * iterative_schema.h
 *
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#ifndef ITERATIVE_SCHEMA_H_
#define ITERATIVE_SCHEMA_H_

#include "value.h"

//calculation methods
void leap_frog(unsigned int n_block, unsigned int thread_per_block);
value_t L_in_time(value_t t);
void calc_vf(value_t t);//for expansion algorithm, see kolafa_expansion() documentation
void TRVP_execute_and_copy_velocities( unsigned int n_block, unsigned int thread_per_block,unsigned long i,
									   value_t *d_velocityTRVP2_x_out, value_t *d_velocityTRVP2_y_out, value_t *d_velocityTRVP2_z_out,
									   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
									   value_t *d_velocity_x_m3l2h_in, value_t *d_velocity_y_m3l2h_in, value_t *d_velocity_z_m3l2h_in,
									   value_t *d_velocity_x_m5l2h_in, value_t *d_velocity_y_m5l2h_in, value_t *d_velocity_z_m5l2h_in,
									   value_t *d_velocityTRVP2_ref_x_out, value_t *d_velocityTRVP2_ref_y_out, value_t *d_velocityTRVP2_ref_z_out,
									   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
									   value_t *d_velocity_ref_x_m3l2h_in, value_t *d_velocity_ref_y_m3l2h_in, value_t *d_velocity_ref_z_m3l2h_in,
									   value_t *d_velocity_ref_x_m5l2h_in, value_t *d_velocity_ref_y_m5l2h_in, value_t *d_velocity_ref_z_m5l2h_in);

#endif /* ITERATIVE_SCHEMA_H_ */
