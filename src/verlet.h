/*
 * verlet.h
 *
 *  Created on: 19. Dec 2019
 *      Author: dixiw
 */

#ifndef VERLET_H_
#define VERLET_H_


#include "value.h"

void verlet_const_memory_copy();
void verlet_device_alloc(ushort verlet_size);
void verlet_device_free();
unsigned int verlet_refresh(value_t *d_position_x, value_t *d_position_y, value_t *d_position_z,
							value_t *d_lx_buffer, unsigned long ii);

#endif /* VERLET_H_ */
