/*
 * triplet.cu
 *
 *  Created on: 26 Jan 2018
 *      Author: dixiw
 */

#include "triplet.h"


/* Function: triplet_alloc
 * function for triplet structure allocation
 *
 * <triplet> is the <value_t> pointers to x,y,z data
 *
 * Notes:
 * > initial declaration of the triplet is left to the user
 *
 * Parameter:
 * what_to_alloc - triplet pointer what is to be allocated
 * alloc_count 	 - the number of ellocated elements
 */
void triplet_alloc(triplet *what_to_alloc,int alloc_count)
{
	what_to_alloc->x = (value_t*) malloc(alloc_count*sizeof(value_t));
	what_to_alloc->y = (value_t*) malloc(alloc_count*sizeof(value_t));
	what_to_alloc->z = (value_t*) malloc(alloc_count*sizeof(value_t));
	return;
}

/* Function: triplet_dealloc
 * function for triplet structure deallocation
 *
 * Parameter:
 * what_to_dealloc - triplet pointer what is to be deallocated
 */
void triplet_dealloc(triplet *what_to_dealloc)
{
	if(what_to_dealloc->x!= NULL) free(what_to_dealloc->x);
	if(what_to_dealloc->y!= NULL) free(what_to_dealloc->y);
	if(what_to_dealloc->z!= NULL) free(what_to_dealloc->z);
	if(what_to_dealloc != NULL) free(what_to_dealloc);
}

