/*
 * value.h
 *
 *  Created on: Mar 13, 2018
 *      Author: celnyd
 */

#ifndef VALUE_H_
#define VALUE_H_

#include "cuda_definition.h"

/* Defines: value_t
 * typedef of the used precision
 *
 * automatically handled by <PREC> flag and can be either float or double
 */
/* Defines: value3_t
 * typedef of the used precision for derived type
 *
 * automatically handled by <PREC> flag and can be either float3 or double3
 */
// == typedef handling ==
#if PREC == FLOAT
	typedef float value_t;
	typedef float3 value3_t;
//	typedef float4 value4_t;
#elif PREC == DOUBLE
	typedef double value_t;
	typedef double3 value3_t;
//	typedef double4 value4_t;
#endif


#endif /* VALUE_H_ */
