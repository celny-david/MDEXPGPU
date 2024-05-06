/*
 * verlet.h
 *
 *  Created on: 19. Dec 2019
 *      Author: dixiw
 */

#ifndef VERLET_SEQ_H_
#define VERLET_SEQ_H_


#include "value.h"

unsigned int verlet_continual_refresh(value_t *d_position_x, value_t *d_position_y, value_t *d_position_z,
									 unsigned int verlet_size,
									 ushort2 **verlet);

#endif /* VERLET_SEQ_H_ */
