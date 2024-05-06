/*
 * lj_elspline.h
 *
 *  Created on: 31 Aug 2018
 *      Author: dixiw
 */

#ifndef LJ_ELSPLINE_H_
#define LJ_ELSPLINE_H_

#include "value.h"
#include "cuda_definition.h"

// POT_TYPE= 20 -> qq4r-12.6875-32 // FOR_POT = 1


void set_elsplines();

#if ELPOT_TYPE >0

__device__ value_t lj_elfor(value_t x);
__device__ value_t lj_elpot(value_t x);

/* !!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!*/
/* !!!!! CODE FOR SEQUENTIAL RUN !!!! */
/* !!!!! !!!!!!!!!!!!!!!!!!!!!!! !!!!!*/
value_t lj_elfor_seq(value_t x);
value_t lj_elpot_seq(value_t x);

#endif

#endif /* LJ_ELSPLINE_H_ */
