/*
 * iterative_schema.h
 *
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#ifndef ITERATIVE_SCHEMA_SEQ_H_
#define ITERATIVE_SCHEMA_SEQ_H_

#include "triplet.h"

//calculation methods
// void leap_frog_seq(unsigned int n_block, unsigned int thread_per_block);//MK_modified
void leap_frog_seq(triplet *mol_pos_tmp, triplet *atm_pos_tmp, triplet *mol_vel_tmp, triplet *atm_vel_tmp, triplet *atm_fcs_tmp);
void verlet_integration_seq(triplet *mol_pos_tmp, triplet *atm_pos_tmp, triplet *mol_vel_tmp, triplet *atm_vel_tmp, triplet *atm_fcs_tmp);
#endif /* ITERATIVE_SCHEMA_SEQ_H_ */
