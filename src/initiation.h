/*
 * initiation.h
 *
 *  Created on: 9 Aug 2017
 *      Author: dixiw
 */

#ifndef INITIATION_H_
#define INITIATION_H_

#include "triplet.h"

// === pull/obtain functions ===
void fill_starting_properties(triplet *molecule_position, triplet *atom_position,
							  triplet *molecule_velocity, triplet *atom_velocity,
							  triplet *atom_force);

#endif /* INITIATION_H_ */
