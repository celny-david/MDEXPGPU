/*
 * property.h
 *
 *  Created on: 22 Nov 2017
 *      Author: dixiw
 */

#ifndef PROPERTY_H_
#define PROPERTY_H_

#include <stdio.h>
#include "value.h"
#include "triplet.h"
#include "cuda_definition.h"

// == global variable parameter settings ==
const unsigned int p_mol_count = 1024; // <- set this
const float p_boxsize_x = 100; // <- set this
const float p_boxsize_y = 100; // <- set this
const float p_boxsize_z = 100; // <- set this
const float p_density=100;     // <- set this in kg*m^(-3)
const unsigned int p_temperature = 300; // <- set this
const unsigned int p_verlet_refresh = 10; // <- set this for the frequency of Verlet refresh
									// TODO develop the dynamic way to set this
// = default variable parameters =
const double p_d_dt = 0.001;         // <- set if you want to change default behaviour[ps]
const double p_d_total_time =1;     // <- set if you want to change default behaviour[ps]
const unsigned long int p_d_cp_timemark_steps = 1000;// <- how often[in printout_steps] the time mark(containing real time and date) will be placed into .cp file
#ifndef TESTING
const int p_d_printout_steps = ceil(0.001/p_d_dt);   // <- set if you want to change default behaviour of energy, cp ... outputting[ps]
const int p_d_plb_hist_printout_steps = ceil(1/p_d_dt);   // <- set if you want to change default behaviour of history and .plb outputting[ps]
const char p_d_output_type[32] = "th"; // <- set if you want to change default behaviour
								// e = energy file
								// h = history file (positions only)
								// t = test file (velocities only)
								// 2 = vel file (forces only)
								// p = plb file (positions only, Macsimus standard)
								// 3 = outlying molecules file
								// c = convergence profile .cp file
								// f = binary file with positions and velocities .cfg file
								// m = MACSIMUS standart output- .cp .plb+.mol and .cfg
								// a = text file with positions and velocities .cfa file
#else
const int p_d_printout_steps = 1;   // <- set if you want to change default behaviour of history/energy/... outputting[ps]
const int p_d_printout_steps =1;    // <- set if you want to change default behaviour of history and .plb outputting[ps]
const char p_d_output_type[32] = "eht2p3cfa";
#endif
const unsigned int p_d_nshake = 10; // <- set if you want to change default behaviour


// === functions returning some property ===
value_t Return_Etot(value_t* ekin_in, value_t* epot_in);
value_t Return_Epot(value_t* epot_in);
value_t Return_Ekin(value_t* ekin_in);
value_t Return_Pcfg(triplet *molecule_postition_in, triplet *atom_position_in,triplet *atm_fcs_tmp, value_t* ekin_in);
value_t Return_Temperature(triplet *atom_velocity_in);


#endif /* PROPERTY_H_ */
