/*
 * input.h
 *
 *  Created on: Dec 16, 2019
 *      Author: klimam
 */

#ifndef INPUT_PROPERTY_H_
#define INPUT_PROPERTY_H_

#include "value.h"
#include "triplet.h"

void set_parameters_from_property(void);
char* removeLastN(char* str, int n );
void set_parameters_from_input_file(void);
void print_input_error(void);
void set_outputs_and_periods(void);
void L_from_density_file(value_t t);

/* Structures: Input_prop
 * enveloping of the basic input capacity of the program
 *
 * useful mainly for setup input side of IO
 *
 * input_source		- int flag if where the source of parameters comes from [0 =property.h ,1= *.cdef]
 *
 * input_name		- char array[256] name buffer for input file (*.cdef)
 * box_file_name	- char array[256] name buffer for box file (*.box)
 * input_cfa_name	- char array[256] name buffer for cfa file (*.cfa)
 * input_cfg_name	- char array[256] name buffer for cfg file (*.cfg)
 *
 */

struct Input_prop
{
	int input_source=0; //0= source is property.h
						//1= source is external .cdef file
	char input_name[256]; //with .cdef suffix
	char box_file_name[256];//with .box suffix
						    //1= source is external .cdef file
	char input_cfa_name[256]; //with .cfa suffix
	char input_cfg_name[256]; //with .cfg suffix

};
extern Input_prop inpp;

void set_input_from_file(char *input_file_whole);
void set_box_file_name();
void set_input_cfa_name(char input_cfa_name_buffer[]);
void set_input_cfg_name(char input_cfa_name_buffer[]);
void set_parameters_from_cdef_cfa_file(triplet *molecule_position, triplet *atom_position,
		   	   	   	   	   	   	   	   triplet *molecule_velocity, triplet *atom_velocity,
		   	   	   	   	   	   	   	   triplet *atom_force);
void set_parameters_from_cdef_cfg_file(triplet *molecule_position, triplet *atom_position,
		   	   	   	   	   	   	   	   triplet *molecule_velocity, triplet *atom_velocity,
		   	   	   	   	   	   	   	   triplet *atom_force);
void position_decombiner(triplet *position_in,triplet *molecule_postition_in,triplet *atom_position_in);
value_t get_lx_from_density_file(value_t t, FILE* rhof);
void fill_lx_buffer(value_t *lx_buffer,unsigned int buffer_size);

#endif /* INPUT_PROPERTY_H_ */
