/*
 * output.h
 *
 *  Created on: Dec 16, 2019
 *      Author: klimam
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <stdio.h>
#include "value.h"
#include "triplet.h"

// helper method
void position_combiner(triplet *position_in,triplet *molecule_postition_in,triplet *atom_position_in);

//********* Used just in output.cu
// void save_to_energy(value_t *ekin_in, value_t *epot_in);
// void save_to_history(triplet *molecule_postition_in, triplet *atom_position_in);
// void save_to_history_separate(triplet *molecule_postition_in, triplet *atom_position_in);
// void save_to_fcs(triplet *atom_forces_in);

// void save_to_vel(triplet *molecule_velocity_in, triplet *atom_velocity_in);
// void save_to_oob(triplet *molecule_postition_in, triplet *atom_position_in);
// void save_to_test4(triplet *molecule_postition_in, triplet *atom_position_in);

// void save_to_plb(triplet *molecule_postition_in, triplet *atom_position_in);
// void save_to_cp(triplet *molecule_postition_in, triplet *atom_position_in,triplet *atom_forces_in,value_t* ekin_in, value_t* epot_in);
// void save_to_cfa(triplet *molecule_postition_in, triplet *atom_position_in,triplet *molecule_velocity_in, triplet *atom_velocity_in, value_t* ekin_in, value_t* epot_in);
// void save_to_cfg(triplet *molecule_postition_in, triplet *atom_position_in,triplet *molecule_velocity_in, triplet *atom_velocity_in, value_t* ekin_in, value_t* epot_in);

// save methods (only top level methods)
void print_into_files(triplet *mol_pos_tmp,triplet *atm_pos_tmp,triplet *mol_vel_tmp,triplet *atm_vel_tmp,triplet *atm_fcs_tmp,value_t *ekin_tmp,value_t *epot_tmp);


#endif /* OUTPUT_H_ */
