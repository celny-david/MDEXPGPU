/*
 * triplet.h
 *	define structure for position/velocity/force type of data
 *	data contain x,y,z triplet per element and in theory form array of float3
 *	due to the requirement of memeory coalescent the AoS is changed to SoA named triplet
 *  Created on: 26 Jan 2018
 *      Author: dixiw
 */

#ifndef TRIPLET_H_
#define TRIPLET_H_

#include "value.h"

/* Structures: triplet
 * the enveloping structure of three <value_t> vectors
 *
 * x - the value_t pointer to x data array
 * y - the value_t pointer to y data array
 * z - the value_t pointer to z data array
 *
 */
typedef struct
{
	value_t *x;
	value_t *y;
	value_t *z;
}triplet;


void triplet_alloc(triplet *what_to_alloc,int alloc_count);
void triplet_dealloc(triplet *what_to_dealloc);
#endif /* TRIPLET_H */
