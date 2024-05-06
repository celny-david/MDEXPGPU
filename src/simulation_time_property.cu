/*
 * simulation_time_property.cu
 *
 *  Created on: Jan 23, 2020
 *      Author: dixiw
 */

#include <stdio.h>
#include <math.h>
#include <unistd.h> // hostname

#include "err_check.h"
#include "cuda_definition.h"
#include "simulation_time_property.h"
#include "system_property.h"
#include "iterative_property.h"
#include "output_property.h"

/* Function: print_info
 * print to console of basic info about system settings, host and GPU
 *
 * if also_log_to_file is true than same is printed to runlog_file
 *
 * Notes:
 * > number of molecules aligned to be divisible by number of threads
 * > number of Blocks, number of Thread per block, number of work assigned per threa
 * > volume settings Lx * Ly * Lz [Angstrom]
 * > cutoff settings cutoff [Angstrom]
 * > velet width, verlet reset and # of shake iteration info [#]
 * > total time of simulation, timestep [ms]
 * > temperature of system []
 * > time string [hh:mm:ss] - date string [dd.mm.yyyy]
 * > name of the host computer
 * > name of the GPU the first GPU
 * > number of SM available in GPU, Global memory, Shared memory per block, Constant memory
 *
 * Parameter:
 * also_log_to_file - bool value if the info should be written into runlog_file
 */
void print_info(bool also_log_to_file)
{
	time_t current_time;
	char time_string[80];
	char hostname[1024];
	FILE * File_bufer;

	// = get host specific info =
	current_time = time(NULL);
	gethostname(hostname, sizeof(hostname)-1);

	// = get GPU specific info =
	int current_device;
	cudaDeviceProp device_prop;
	cudaSafeCall( cudaGetDevice(&current_device));
	cudaSafeCall( cudaGetDeviceProperties(&device_prop,current_device));

	// = print info to console =
	strftime(time_string,80,"%X - %d.%m.%Y",localtime(&current_time));
	switch(SUBSTANCE)
	{
		case ARGON:
			printf("*** %i molecules of argon***\n",sysp.molecule_count_aligned);
			break;
		case NITROGEN:
			printf("*** %i molecules of nitrogen***\n",sysp.molecule_count_aligned);
			break;
		case SPCE:
			printf("*** %i molecules of SPCE water***\n",sysp.molecule_count_aligned);
			break;
		case TIP4P:
			printf("*** %i molecules of TIP4P water***\n",sysp.molecule_count_aligned);
			break;
	}
#if PREC == FLOAT
	printf("*** potential: %i ,in FLOAT ***\n",POT_TYPE);
#elif PREC == DOUBLE
	printf("*** potential: %i ,in DOUBLE ***\n",POT_TYPE);
#endif
	printf("*** BLOCK: %i, THREAD: %i, WORKLOAD: %i ***\n",sysp.atom_count_aligned/THREAD_P_BLOCK,THREAD_P_BLOCK, WORKLOAD);
	//printf("*** %4.2f*%4.2f*%4.2f ***\n",sysp.lx,sysp.ly,sysp.lz);
	printf("*** %f*%f*%f ***\n",sysp.lx,sysp.ly,sysp.lz);
	printf("*** cutoffs: norm= %10.6f",sysp.cutoff);
	if(sysp.cutoff_verlet>0) printf(", verlet= %10.6f",sysp.cutoff_verlet);
	if(sysp.cutoff_elstat>0) printf(", elstat= %10.6f",sysp.cutoff_elstat);
	printf(" ***\n");
	printf("*** V_width= %i, V_reset= %i, n_shake= %i ***\n",0,itep.verlet_refresh ,itep.n_shake);
	printf("*** t= %4.2e, dt= %8.4e ***\n",itep.total_time,itep.dt);
	printf("*** T = %6.2f ***\n",sysp.temperature);
#ifdef EXPANSION
	printf("*** EXPANSION ***\n");
#else
	printf("*** Rho = %6.2f ***\n",sysp.density);
#endif
	printf("*** %s ***\n",time_string);
	printf("*** output: E=%d, H=%d, T=%d, T2=%d ***\n",outp.energy_out,outp.history_out,outp.fcs_out, outp.vel_out);
	printf("*** %s ***\n", hostname);
	printf("*** %s ***\n", device_prop.name);
	printf("*** SM: %i, Gmem: %d MB, Smem/B: %d kB, Cmem: %d kB ***\n",(int)device_prop.multiProcessorCount
																	  ,(int)device_prop.totalGlobalMem/1024/1024
																	  ,(int)device_prop.sharedMemPerBlock/1024
																	  ,(int)device_prop.totalConstMem/1024);

	if(also_log_to_file)
	{
		File_bufer= fopen(outp.log_name,"a");
		// = print info to console =
		fprintf(File_bufer,"\n ==================== \n");
		switch(SUBSTANCE)
		{
			case ARGON:
				fprintf(File_bufer,"*** %i molecules of argon***\n",sysp.molecule_count_aligned);
				break;
			case NITROGEN:
				fprintf(File_bufer,"*** %i molecules of nitrogen***\n",sysp.molecule_count_aligned);
				break;
			case SPCE:
				fprintf(File_bufer,"*** %i molecules of SPCE water***\n",sysp.molecule_count_aligned);
				break;
			case TIP4P:
				fprintf(File_bufer,"*** %i molecules of TIP4P water***\n",sysp.molecule_count_aligned);
				break;
		}
#if PREC == FLOAT
		fprintf(File_bufer,"*** potential: %i ,in FLOAT ***\n",POT_TYPE);
#elif PREC == DOUBLE
		fprintf(File_bufer,"*** potential: %i ,in DOUBLE ***\n",POT_TYPE);
#endif

		fprintf(File_bufer,"*** BLOCK: %i, THREAD: %i, WORKLOAD: %i ***\n",sysp.molecule_count_aligned/THREAD_P_BLOCK,THREAD_P_BLOCK, WORKLOAD);
		fprintf(File_bufer,"*** %4.2f*%4.2f*%4.2f ***\n",sysp.lx,sysp.ly,sysp.lz);
		fprintf(File_bufer,"*** cutoffs: norm= %10.6f",sysp.cutoff);
		if(sysp.cutoff_verlet>0) fprintf(File_bufer,", verlet= %10.6f",sysp.cutoff_verlet);
		if(sysp.cutoff_elstat>0) fprintf(File_bufer,", elstat= %10.6f",sysp.cutoff_elstat);
		fprintf(File_bufer," ***\n");
		fprintf(File_bufer,"*** V_width= %i, V_reset= %i, n_shake= %i ***\n",0,itep.verlet_refresh ,itep.n_shake);
		fprintf(File_bufer,"*** t= %4.2f, dt= %8.4f ***\n",itep.total_time,itep.dt);
		fprintf(File_bufer,"*** T = %6.2f ***\n",sysp.temperature);
#ifdef EXPANSION
	fprintf(File_bufer,"*** EXPANSION ***\n");
#else
	fprintf(File_bufer,"*** Rho = %6.2f ***\n",sysp.density);
#endif
		fprintf(File_bufer,"*** %s ***\n",time_string);
		fprintf(File_bufer,"*** output: E=%d, H=%d, T=%d, T2=%d ***\n",outp.energy_out,outp.history_out,outp.fcs_out, outp.vel_out);
		fprintf(File_bufer,"*** %s ***\n", hostname);
		fprintf(File_bufer,"*** %s ***\n", device_prop.name);
		fprintf(File_bufer,"*** SM: %i, Gmem: %d MB, Smem/B: %d kB, Cmem: %d kB ***\n",(int)device_prop.multiProcessorCount
																		  ,(int)device_prop.totalGlobalMem/1024/1024
																		  ,(int)device_prop.sharedMemPerBlock/1024
																		  ,(int)device_prop.totalConstMem/1024);
		fclose(File_bufer);
	}

	return;
}

/* Function: print_info
 * print the run time values for the specified time profile of the main function
 *
 * if also_log_to_file is true than same is printed to runlog_file
 *
 * Notes:
 * > total time [s]
 * > setup time [s], percentage of total
 * > iteration time [s], percentage of total
 * >    single iteration time [ms], percentage of total
 * > cleanup time [s], percentage of total
 *
 * Parameter:
 * t_start - clock_t time of startup beginning
 * t_setup - clock_t time of setup beginning
 * t_init - clock_t time of initialisation beginning
 * t_iter - clock_t time of iteration beginning
 * t_end - clock_t time of end
 * also_log_to_file - bool value if the info should be written into runlog_file
 */
void print_time_main(bool also_log_to_file)
{
	FILE * File_buffer;


	// === print statistics section ===
	printf("=== TIMING PROFILE ===\n");
	printf("t_total = %f s\n",(float)(simtp.end - simtp.start)/CLOCKS_PER_SEC );
	printf("t_setup = %f s, %09.6f %%\n",(float)(simtp.setup - simtp.start)/CLOCKS_PER_SEC, (float) 100*(simtp.setup - simtp.start)/(simtp.end - simtp.start) );
	printf("t_init  = %f s, %09.6f %%\n",(float)(simtp.init - simtp.setup)/CLOCKS_PER_SEC, (float) 100*(simtp.init - simtp.setup)/(simtp.end - simtp.start) );
	printf("t_iter  = %f s, %09.6f %%\n",(float)(simtp.iter - simtp.init)/CLOCKS_PER_SEC, (float) 100*(simtp.iter - simtp.init)/(simtp.end - simtp.start) );
	printf(" 1_iter ~ %f ms \n",(float)(simtp.iter - simtp.init)/CLOCKS_PER_SEC/(itep.total_time/itep.dt)*1000 );
	printf("t_clean = %f s, %09.6f %%\n",(float)(simtp.end - simtp.iter)/CLOCKS_PER_SEC, (float) 100*(simtp.end - simtp.iter)/(simtp.end - simtp.start) );

	if(also_log_to_file)
	{
		File_buffer=fopen(outp.log_name,"a");
		fprintf(File_buffer,"=== TIMING PROFILE ===\n");
		fprintf(File_buffer,"t_total = %f s\n",(float)(simtp.end - simtp.start)/CLOCKS_PER_SEC );
		fprintf(File_buffer,"t_setup = %f s, %09.6f %%\n",(float)(simtp.setup - simtp.start)/CLOCKS_PER_SEC, (float) 100*(simtp.setup - simtp.start)/(simtp.end - simtp.start) );
		fprintf(File_buffer,"t_init  = %f s, %09.6f %%\n",(float)(simtp.init - simtp.setup)/CLOCKS_PER_SEC, (float) 100*(simtp.init - simtp.setup)/(simtp.end - simtp.start) );
		fprintf(File_buffer,"t_iter  = %f s, %09.6f %%\n",(float)(simtp.iter - simtp.init)/CLOCKS_PER_SEC, (float) 100*(simtp.iter - simtp.init)/(simtp.end - simtp.start) );
		fprintf(File_buffer," 1_iter ~ %f ms \n",(float)(simtp.iter - simtp.init)/CLOCKS_PER_SEC/(itep.total_time/itep.dt)*1000 );
		fprintf(File_buffer,"t_clean = %f s, %09.6f %%\n",(float)(simtp.end - simtp.iter)/CLOCKS_PER_SEC, (float) 100*(simtp.end - simtp.iter)/(simtp.end - simtp.start) );
		fclose(File_buffer);
	}

	return;
}

