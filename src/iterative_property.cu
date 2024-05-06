/*
 * iterative_property.cu
 *
 *  Created on: Jan 18, 2020
 *      Author: klimam
 */

#include "iterative_property.h"

// ===== ===== ========= ===== =====
// ===== ITERATIVE SCHEME PROP =====
// ===== ===== ========= ===== =====

//struct Iterative_prop
//{
//	value_t dt; // delta t time value [picosecond]
//	value_t total_time; // total time value [picosecond]
//  bool printout; // if any data should be printed [true/false]
//	unsigned short printout_steps; // number of steps between position printout
//  unsigned int verlet_refresh; // number of steps between verlet reconstruct
//}itep;

// === set functions iterative scheme ===

/* Function: set_dt
 * the method for setting the time step
 *
 * Notes:
* >!no input correctness check is performed
 *
 * Parameter:
 * dt_in - the timestep value [ps]
 */
void set_dt(value_t dt_in)
{
	itep.dt = dt_in;
	return;
}

/* Function: set_total_time
 * the method for setting the total time of simulation
 *
 * Notes:
 * >no input correctness check is performed
 *
 * Parameter:
 * total_time_in - the total time [ps]
 */
void set_total_time (value_t total_time_in)
{
	itep.total_time = total_time_in;
	set_max_step();
	return;
}

/* Function: set_printout_steps
 * the method for setting the frequency of printout
 *
 * this is given by printout steps count (i.e. each 100 step)
 *
 * Notes:
* >if the printout step is <=0 then no printout is performed
 *
 * Parameter:
 * printout_steps_in - number of steps per printout
 */
void set_printout_steps (int printout_steps_in)
{
	if( printout_steps_in <= 0 )
	{
		itep.printout = false;
		itep.printout_steps = 0;
	}else
	{
		itep.printout = true;
		itep.printout_steps = printout_steps_in;
	}
	return;
}


/* Function: set_plb_hist_printout_steps
 * the method for setting the frequency of printout of history and .plb file
 *
 * this is given by printout steps count (i.e. each 100 step)
 *
 * Notes:
* >if the printout step is <=0 then no printout is performed
 *
 * Parameter:
 * printout_steps_in - number of steps per printout
 */
void set_plb_hist_printout_steps (int plb_hist_printout_steps_in)
{
	if( plb_hist_printout_steps_in <= 0 )
		{
			itep.printout = false;
			itep.plb_hist_printout_steps = 0;
		}
		else
		{
			itep.printout = true;
			itep.plb_hist_printout_steps = plb_hist_printout_steps_in;
		}
		return;
}

/* Function: set_verlet_refresh
 * the method for setting the frequency of verlet list reconstruct
 *
 * Notes:
 * > depends on the system and timestep - usually (~20)
 * >!no input correctness check is performed
 *
 * Parameter:
 * refresh_steps_in - number of steps per refresh
 */
void set_verlet_refresh (int refresh_steps_in)
{
	itep.verlet_refresh = refresh_steps_in;
	return;
}

/* Function: set_n_shake
 * set the number of iterations in shake algorithm
 *
 * Notes:
 * > the number should be positive integer larger than zero
 *
 * Parameter:
 * n_shake_in - unsigned int with number of shake iterations
 */
void set_n_shake(unsigned int n_shake_in)
{
	itep.n_shake=n_shake_in;
	return;
}

void set_max_step()
{
	itep.max_step = ceil(itep.total_time/itep.dt);
}
