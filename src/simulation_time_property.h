/*
 * Simulation_time_property.h
 *
 *  Created on: Jan 23, 2020
 *      Author: dixiw
 */

#ifndef SIMULATION_TIME_PROPERTY_H_
#define SIMULATION_TIME_PROPERTY_H_

/* Structures: Simulation_time_prop
 * enveloping of the timing functionality of program
 *
 * used in main for holding the timing values
 *
 * start 	- clock_t time corresponding to the start of program
 * setup 	- clock_t time after the setup section of program (input parse, setup start configuration )
 * init 	- clock_t time after the initialization section of the program (param setting, gpu param init)
 * iter 	- clock_t time containing the whole run of iteration schema (depends on the # of steps)
 * end 		- clock_t time corresponding to the end of program (after the dealloc of variables)
 *
 */
struct Simulation_time_prop
{
	clock_t start;
	clock_t setup;
	clock_t init;
	clock_t iter;
	clock_t end;
};
extern Simulation_time_prop simtp;

// === functions for structure Simulation_time_prop type  ===
void print_info(bool also_log_to_file);
void print_time_main(bool also_log_to_file);

#endif /* SIMULATION_TIME_PROPERTY_H_ */
