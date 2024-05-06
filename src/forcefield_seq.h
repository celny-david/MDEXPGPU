/*
 * forcefield.h
 *
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#ifndef FORCEFIELD_SEQ_H_
#define FORCEFIELD_SEQ_H_

#include "value.h"

// variant with dynamic verlet list

void force_calculation_seq(value_t *position_x_in, value_t *position_y_in, value_t *position_z_in,
	 					   value_t *position_ref_x_in, value_t *position_ref_y_in, value_t *position_ref_z_in,
	 					   value_t *force_ref_x_in, value_t *force_ref_y_in, value_t *force_ref_z_in,
	 					   value_t *epot_in,
	 					   unsigned int verlet_size_in,
	 					   ushort2 *verlet);
#endif /* FORCEFIELD_SEQ_H_ */
