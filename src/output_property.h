/*
 * output_property.h
 *
 *  Created on: Jan 23, 2020
 *      Author: dixiw
 */

#ifndef OUTPUT_PROPERTY_H_
#define OUTPUT_PROPERTY_H_

/* Structures: Output_prop
 * enveloping of the basic output capacity of the program
 *
 * useful mainly for setup output side of IO
 *
 * energy_out			- bool flag if energy output is produced
 * history_out			- bool flag if history output is produced
 * fcs_out				- bool flag if fcs output is produced
 * vel_out				- bool flag if vel output is produced
 * oob_out				- bool flag if oob output is produced
 * test4_out			- bool flag if test4 output is produced
 * plb_out				- bool flag if plb output is produced
 * cp_out				- bool flag if cp output is produced
 * cfg_out				- bool flag if cfg output is produced
 * cfa_out				- bool flag if cfa output is produced
 * log_out				- bool flag if log output is produced
 *
 * output_files_name	- char array[256] name buffer for output file
 * energy_name			- char array[256] name buffer for energy file
 * history_name			- char array[256] name buffer for history file
 * fcs_name				- char array[256] name buffer for fcs file
 * vel_name				- char array[256] name buffer for vel file
 * oob_name				- char array[256] name buffer for oob file
 * test4_name			- char array[256] name buffer for test4 file
 * plb_name				- char array[256] name buffer for plb file
 * mol_name				- char array[256] name buffer for mol file
 * cp_name				- char array[256] name buffer for cp file
 * cfg_name				- char array[256] name buffer for cfg file
 * cfa_name				- char array[256] name buffer for cfa file
 * log_name				- char array[256] name buffer for log file
 * stp_name				- char array[256] name buffer for stp file
 *
 * energy_print_period	- uint frequency of energy printout if requested
 * history_print_period	- uint frequency of history printout if requested
 * vel_print_period		- uint frequency of vel printout if requested
 * oob_print_period		- uint frequency of oob printout if requested
 * test4_print_period	- uint frequency of test4 printout if requested
 * fcs_print_period		- uint frequency of fcs printout if requested
 * plb_print_period		- uint frequency of plb printout if requested
 * cp_print_period		- uint frequency of cp printout if requested
 * cfg_print_period		- uint frequency of cfg printout if requested
 * cfa_print_period		- uint frequency of cfa printout if requested
 * cp_print_timemark	- uint frequency of cp_p printout if requested (DEFAULT value=1000)
 *
 */

struct Output_prop
{
	bool energy_out = false;
	bool history_out = false;
	bool fcs_out = false;
	bool vel_out = false;
	bool oob_out = false;
	bool test4_out = false;
	bool plb_out = false;
	bool cp_out = false;
	bool cfg_out = false;
	bool cfa_out = false;
	bool log_out = true;

	char output_files_name[256];
	char energy_name[256];
	char history_name[256];
	char fcs_name[256];
	char vel_name[256];
	char oob_name[256];
	char test4_name[256];
	char plb_name[256];
	char mol_name[256];
	char cp_name[256];
	char cfg_name[256];
	char cfa_name[256];
	char log_name[256];
	char stp_name[256];

	unsigned int energy_print_period=UINT_MAX; //given in steps
	unsigned int history_print_period=UINT_MAX;
	unsigned int vel_print_period=UINT_MAX;
	unsigned int oob_print_period=UINT_MAX;
	unsigned int test4_print_period=UINT_MAX;
	unsigned int fcs_print_period=UINT_MAX;
	unsigned int plb_print_period=UINT_MAX;
	unsigned int cp_print_period=UINT_MAX;
	unsigned int cfg_print_period=UINT_MAX;
	unsigned int cfa_print_period=UINT_MAX;
	unsigned int cp_print_timemark=1000;//UINT_MAX;  //time mark is written every 1000 data printouts
};
extern Output_prop outp;

// === set functions for structure type Output_properties ===
void set_output_name(char output_name_buffer[]);
void set_output_name(void);
void create_files(void);
void create_log_file(void);
void set_all_printouts (unsigned int printout_steps_in);//given in steps

#endif /* OUTPUT_PROPERTY_H_ */
