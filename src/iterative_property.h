/*
 * iterative_property.h
 *
 *  Created on: Jan 18, 2020
 *      Author: klimam
 */

#ifndef ITERATIVE_PROPERTY_H_
#define ITERATIVE_PROPERTY_H_

#include "value.h"

// == ==================================== ==
// == iterative scheme relevant parameters ==

/* Structures: Iterative_prop
 * enveloping of the iteration properties
 *
 * useful mainly for iteration schema
 *
 * step						- unsigned long actual step of calculation [#]
 * dt 						- delta t time value [picosecond]
 * total_time 				- total time value [picosecond]
 * printout 				- if any data should be printed [true/false]
 * printout_steps 			- number of steps between position printout [# steps]
 * verlet_refresh 			- number of steps between verlet reconstruct [# steps]
 * n_shake 					- number of shake iterations [# steps]
 * plb_hist_printout_steps 	- number of steps between .plb and history printout [# steps]
 * vft 						- for expansion algorithm [MK TODO]
 * lambda_tphh 				- for expansion algorithm, means lambda(t+0.5*h) t p(lus)hh(alf) [MK TODO]
 */
struct Iterative_prop
{
	unsigned long step; // current step of
	unsigned long max_step; // maximal reachable step
	value_t dt; // delta t time value_t [picosecond]
	value_t total_time; // total time value_t [picosecond]
	bool printout; // if any data should be printed [true/false]
	unsigned int printout_steps; // number of steps between printout
	unsigned int verlet_refresh; // number of steps between verlet reconstruct
	unsigned int n_shake; //number of shake iteration in the shake algorithm
	unsigned int plb_hist_printout_steps; // number of steps between .plb and history printout

	value_t vft;//for expansion algorithm, see kolafa_expansion() documentation
	value_t lambda_tphh;//for expansion algorithm, see kolafa_expansion() documentation, means lambda(t+0.5*h) t p(lus)hh(alf)

};
extern Iterative_prop itep;

// === set functions iterative scheme ===
void set_dt(value_t dt_in);
void set_total_time (value_t total_time_in);
void set_printout_steps (int printout_steps_in);
void set_plb_hist_printout_steps (int plb_hist_printout_steps_in);//MK_function
void set_verlet_refresh (int printout_steps_in);
void set_n_shake(unsigned int n_shake_in);
void set_max_step();


#endif /* ITERATIVE_PROPERTY_H_ */
