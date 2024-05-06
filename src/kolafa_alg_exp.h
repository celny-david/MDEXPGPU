/*
 * forcefield.h
 *
 *  Created on: 19 Dec 2019
 *      Author: dixiw
 */

#ifndef KOLAFA_ALG_EXP_H_
#define KOLAFA_ALG_EXP_H_

#include "value.h"

__global__ void kolafa(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_x_in, value_t *d_force_y_in, value_t *d_force_z_in,
					   value_t *d_ekin_in);

__global__ void kolafa_c(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					     value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					     value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					     value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					     value_t *d_force_x_in, value_t *d_force_y_in, value_t *d_force_z_in,
					     value_t *d_ekin_in);

__global__ void kolafa_expansion(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					   value_t * d_velocity_x_TRVP_in,value_t * d_velocity_y_TRVP_in,
					   value_t * d_velocity_z_TRVP_in,value_t * d_velocity_ref_x_TRVP_in,value_t * d_velocity_ref_y_TRVP_in,
					   value_t * d_velocity_ref_z_TRVP_in,value_t *d_ekin_in,value_t m_lambda_tphh,
					   value_t m_vft,value_t *m_lxyz);

__global__ void kolafa_expansion_NoMM(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					   value_t * d_velocity_x_TRVP_in,value_t * d_velocity_y_TRVP_in,
					   value_t * d_velocity_z_TRVP_in,value_t * d_velocity_ref_x_TRVP_in,value_t * d_velocity_ref_y_TRVP_in,
					   value_t * d_velocity_ref_z_TRVP_in,value_t *d_ekin_in,value_t *d_lx_buffer,
					   unsigned long ii);

__global__ void kolafa_graviti_center_corr(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					   value_t *d_ekin_in);

__global__ void kolafa_expansion_NoMM_graviti_center_corr(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
					   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
					   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
					   value_t * d_velocity_x_TRVP_in,value_t * d_velocity_y_TRVP_in,
					   value_t * d_velocity_z_TRVP_in,value_t * d_velocity_ref_x_TRVP_in,value_t * d_velocity_ref_y_TRVP_in,
					   value_t * d_velocity_ref_z_TRVP_in,value_t *d_ekin_in,value_t *d_lx_buffer,
					   unsigned long ii);

__global__ void TRVP2( value_t *d_velocityTRVP2_x_out, value_t *d_velocityTRVP2_y_out, value_t *d_velocityTRVP2_z_out,
		 	 	 	   value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
					   value_t *d_position_x_m1h_in, value_t *d_position_y_m1h_in, value_t *d_position_z_m1h_in,
					   value_t *d_position_x_m2h_in, value_t *d_position_y_m2h_in, value_t *d_position_z_m2h_in,
					   value_t *d_position_x_m3h_in, value_t *d_position_y_m3h_in, value_t *d_position_z_m3h_in);

__global__ void TRVP2_ref( value_t *d_velocityTRVP2_ref_x_out, value_t *d_velocityTRVP2_ref_y_out, value_t *d_velocityTRVP2_ref_z_out,
						   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
						   value_t *d_position_ref_x_m1h_in, value_t *d_position_ref_y_m1h_in, value_t *d_position_ref_z_m1h_in,
						   value_t *d_position_ref_x_m2h_in, value_t *d_position_ref_y_m2h_in, value_t *d_position_ref_z_m2h_in,
						   value_t *d_position_ref_x_m3h_in, value_t *d_position_ref_y_m3h_in, value_t *d_position_ref_z_m3h_in);

__global__ void TRVP2_vel( value_t *d_velocityTRVP2_x_out, value_t *d_velocityTRVP2_y_out, value_t *d_velocityTRVP2_z_out,
					       value_t *d_velocity_m05h_x_in, value_t *d_velocity_m05h_y_in, value_t *d_velocity_m05h_z_in,
					       value_t *d_velocity_m15h_x_in, value_t *d_velocity_m15h_y_in, value_t *d_velocity_m15h_z_in,
					       value_t *d_velocity_m25h_x_in, value_t *d_velocity_m25h_y_in, value_t *d_velocity_m25h_z_in);

__global__ void TRVP2_ref_vel( value_t *d_velocityTRVP2_ref_x_out, value_t *d_velocityTRVP2_ref_y_out, value_t *d_velocityTRVP2_ref_z_out,
		 	 	 	 	 	   value_t *d_velocity_m05h_ref_x_in, value_t *d_velocity_m05h_ref_y_in, value_t *d_velocity_m05h_ref_z_in,
							   value_t *d_velocity_m15h_ref_x_in, value_t *d_velocity_m15h_ref_y_in, value_t *d_velocity_m15h_ref_z_in,
							   value_t *d_velocity_m25h_ref_x_in, value_t *d_velocity_m25h_ref_y_in, value_t *d_velocity_m25h_ref_z_in);

__global__ void kinetic_energy_initial_calculation(value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
												   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
												   value_t *d_ekin_in);



#endif /* KOLAFA_ALG_EXP_H_ */
