/* About: Main program
 *	 Name        - Mackimus_GPU_module
 *	 Author      - David Celny (Dixiw)
 *	 Version     - TBA
 *	 Copyright   - TBA
 *	 Description - Molecular simulation software macsimus GPU compatible module
 *				   for the heterogeneous molecular dynamics simulation
 */

#include <stdio.h>
#include <unistd.h>

#include "err_check.h"
#include "triplet.h"
#include "initiation.h"
#include "input_property.h"
#include "system_property.h"
#include "iterative_property.h"
#include "simulation_time_property.h"
#include "output_property.h"
#include "iterative_schema.h"
#include "iterative_schema_seq.h"
#include "cuda_definition.h"
#include "device_manipulation_methods.h"

System_prop sysp;
Iterative_prop itep;
Output_prop outp;
Simulation_time_prop simtp;
Input_prop inpp;

/* Function: Main
 * the main function
 *
 *  + initialisation
 *  + setting up the parameters
 *  + execution
 *  + timing
 *  + output
 *  + cleanup
 *
 * Parameters:
 *
 * 	  argc - the command line argument count.
 *    argv[] - the array of command line arguments.
 */

long int seed = 123456789;

int main( int argc, char *argv[] )
{

	// === initialisation section ===
	triplet *molecule_position = (triplet*) malloc(sizeof(triplet));
	triplet *molecule_velocity = (triplet*) malloc(sizeof(triplet));
	triplet *atom_position = (triplet*) malloc(sizeof(triplet));
	triplet *atom_velocity = (triplet*) malloc(sizeof(triplet));
	triplet *atom_force = (triplet*) malloc(sizeof(triplet));
	simtp.start = clock();

	// == system specific parameters ==
	if(argc==2)
	{
		set_input_from_file(argv[1]);
		create_log_file();
		if( access( inpp.input_cfg_name , R_OK ) != -1 )
		{
			set_parameters_from_cdef_cfg_file(molecule_position, atom_position,
					   	   	   	   	   	   	  molecule_velocity, atom_velocity,
					   	   	   	   	   	   	  atom_force);
			create_files();
		}
		else if( access( inpp.input_cfa_name , R_OK ) != -1 )
		{
			set_parameters_from_cdef_cfa_file(molecule_position, atom_position,
											  molecule_velocity, atom_velocity,
											  atom_force);
			create_files();			
		}
		else
		{
			set_parameters_from_input_file();
			create_files();
			
			// === create atom positions types & velocities ===
			srand(seed); // seed set up
			fill_starting_properties(molecule_position, atom_position,
									 molecule_velocity, atom_velocity,
									 atom_force);
		}

	}
	else
	{
		print_input_error();
	}

	// === info display ===
	print_info(outp.log_out);
	simtp.setup = clock();

	// === copying system properties to GPU device ===
	// the properties have to be fully available at this time
	prop2device();

	// == device alloc + transfer ==
	data_device_alloc();
	data_host2device(molecule_position, atom_position, molecule_velocity, atom_velocity, atom_force);
	check_settings(); //TODO move to apropriate location durin input setting
	simtp.init = clock();

	// === run the computation ===
	// cudaSafeCall( cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
	// GPU PARALLEL
	leap_frog(sysp.molecule_count_aligned/THREAD_P_BLOCK,THREAD_P_BLOCK);//falls for Argon and Nitrogen here in Float and Double
	// CPU SEQUENTIAL
	// leap_frog_seq(molecule_position, atom_position, molecule_velocity, atom_velocity, atom_force);//falls for Argon and Nitrogen here in Float and Double
	// verlet_integration_seq(molecule_position, atom_position, molecule_velocity, atom_velocity, atom_force);
	// cudaSafeCall(cudaDeviceSynchronize());
	simtp.iter = clock();

	// === clear the garbage ===
	triplet_dealloc(molecule_position);
	triplet_dealloc(molecule_velocity);
	triplet_dealloc(atom_position);
	triplet_dealloc(atom_velocity);
	triplet_dealloc(atom_force);
	data_device_free();
	simtp.end = clock();

	// === info display ===
	print_time_main(outp.log_out);

	// === finish, reset of GPU ===
	// cudaSafeCall( cudaDeviceReset());
	return 0;
}

