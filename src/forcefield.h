/*
 * forcefield.h
 *
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#ifndef FORCEFIELD_H_
#define FORCEFIELD_H_

#include "value.h"

// variant with dynamic verlet list
void force_calculation(value_t *d_position_x_in, value_t *d_position_y_in, value_t *d_position_z_in,
 					   value_t *d_position_ref_x_in, value_t *d_position_ref_y_in, value_t *d_position_ref_z_in,
 					   value_t *d_force_ref_x_in, value_t *d_force_ref_y_in, value_t *d_force_ref_z_in,
 					   value_t *d_epot_in,
 					   unsigned int verlet_size_in,
 					   ushort2 *d_verlet,
 					   value_t *d_lx_buffer, unsigned long lx_postition);

#endif /* FORCEFIELD_H_ */
