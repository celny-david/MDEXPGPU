/*
 * iterative_schema.h
 *
*  Created on: 19. Dec 2019
 *      Author: dixiw
 */

#ifndef DEVICE_MANIPULATION_METHODS_H_
#define DEVICE_MANIPULATION_METHODS_H_

#include "value.h"
#include "triplet.h"
#include "system_property.h"

// manipulation methods to/on/from device
void prop2device();
void data_device_alloc();
void data_device_free();
void data_host2device(triplet *pos_in,
		              triplet *pos_ref_in,
		              triplet *vel_in,
		              triplet *vel_ref_in,
		              triplet *for_ref_in);
void data_device2host_position(triplet *mol_pos_in,triplet *atm_pos_in);
void data_device2host_velocity(triplet *mol_vel_in,triplet *atm_vel_in);
void data_device2host_force(triplet *for_in);
void data_device2host_energy(value_t *ekin_in, value_t *epot_in);

#endif /* DEVICE_MANIPULATION_METHODS_H_ */
