/*
 * lj_spline.h
 *
 *  Created on: 2 Oct 2017
 *      Author: dixiw
 */

#ifndef LJ_SPLINE_H_
#define LJ_SPLINE_H_

#include "value.h"

void set_splines();

__device__ value_t lj_for(value_t x);
__device__ value_t lj_pot(value_t x);

/* !!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!*/
/* !!!!! CODE FOR SEQUENTIAL RUN !!!! */
/* !!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!*/
value_t lj_for_seq(value_t x);
value_t lj_pot_seq(value_t x);

#endif /* LJ_SPLINE_H_ */
