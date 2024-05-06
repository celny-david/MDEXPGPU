/*
 *  iterative_schema.cu
 *	iterative schemas collection
 *	{leap_frog, }
 *  Created on: 31 Aug 2017
 *      Author: dixiw
 */

#include <unistd.h>

#include "err_check.h"
#include "value.h"
#include "triplet.h"
#include "device_manipulation_methods.h"
#include "verlet.h"
#include "kolafa_alg.h"
#include "kolafa_alg_exp.h"
#include "cuda_variables.h"
#include "output.h"
#include "output_property.h"
#include "forcefield.h"
#include "iterative_property.h"
#include "input_property.h"
#include "cuda_definition.h"
#include "iterative_schema.h"

/* Function: leap_frog
 * the procedure responsible for the time iteration
 *
 * done according to the leap-frog scheme
 *
 * contains the Verlet reconstruction force calculation and update and data output
 *
 * verlet reconstructed with <verlet_reconstruct>
 * force is calculate with <force_calculation>
 * the positions are updated with <kolafa> & <kolafa_c>
 *
 * Notes:
 * >the implementation is in form of single loop over timesteps
 * >! untested for high iteration count
 * > the block, thread dimensions have to correspond to the solved task size
 *
 * Parameter:
 * n_block     		- unsigned int of how many block to use
 * thread_per_block - unsigned int of how many threads per block to use
 */
void leap_frog(unsigned int n_block, unsigned int thread_per_block)
{
/*
 * it is made not responsible for allocation/deallocation of device memory
 * this task initial move of data is left in main
 *
 * * leap frog schema: v(t+dt/2) = v(t) + a(t)*dt/2
 * * 				   r(t+dt)   = r(t) + v(t+dt/2)*dt
 * * 				   a(t+dt)   = F( r(t+dt) )/m
 * * 				   v(t+dt)   = v(t+dt/2) + a(t+dt)*dt/2
 */
	double t = 0.0; // [picosecond]
//	unsigned long max_step = ceil(itep.total_time/itep.dt); // number of steps
	FILE * File_buffer; // for premature termination purposes - detect the existence of the file
	unsigned int verlet_size = 0; // NOTE initialization that verlet list is reconstructed at the start

	// tmp variables for energy logging
	value_t *ekin_tmp = (value_t*) malloc(sysp.atom_count_aligned*sizeof(value_t));
	value_t *epot_tmp = (value_t*) malloc(sysp.atom_count_aligned*sizeof(value_t));
	uint *verlet_occupancy_tmp = (uint*) malloc(sysp.molecule_count_aligned*sizeof(uint));
	triplet *mol_pos_tmp = (triplet*) malloc(sizeof(triplet));
	triplet *atm_pos_tmp = (triplet*) malloc(sizeof(triplet));
	triplet *mol_vel_tmp = (triplet*) malloc(sizeof(triplet));
	triplet *atm_vel_tmp = (triplet*) malloc(sizeof(triplet));
	triplet *atm_fcs_tmp = (triplet*) malloc(sizeof(triplet));

	triplet_alloc(mol_pos_tmp,sysp.molecule_count_aligned);
	triplet_alloc(atm_pos_tmp,sysp.atom_count_aligned);
	triplet_alloc(mol_vel_tmp,sysp.molecule_count_aligned);
	triplet_alloc(atm_vel_tmp,sysp.atom_count_aligned);
	triplet_alloc(atm_fcs_tmp,sysp.atom_count_aligned);

#ifdef EXPANSION


	unsigned int boxLx_buffer_size=2*(itep.max_step+1);
	value_t * boxLx_buffer=(value_t*) malloc(boxLx_buffer_size*sizeof(value_t));
	fill_lx_buffer(boxLx_buffer,boxLx_buffer_size);
	cudaSafeCall( cudaMemcpy(d_lx_buffer,boxLx_buffer,boxLx_buffer_size*sizeof(value_t),cudaMemcpyHostToDevice) );

	cudaSafeCall( cudaMemcpy(d_velocity_x_m3l2h,d_velocity_x,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_y_m3l2h,d_velocity_y,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_z_m3l2h,d_velocity_z,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

	cudaSafeCall( cudaMemcpy(d_velocity_x_m5l2h,d_velocity_x,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_y_m5l2h,d_velocity_y,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_z_m5l2h,d_velocity_z,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

	cudaSafeCall( cudaMemcpy(d_velocity_ref_x_m3l2h,d_velocity_ref_x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_ref_y_m3l2h,d_velocity_ref_y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_ref_z_m3l2h,d_velocity_ref_z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

	cudaSafeCall( cudaMemcpy(d_velocity_ref_x_m5l2h,d_velocity_ref_x,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_ref_y_m5l2h,d_velocity_ref_y,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
	cudaSafeCall( cudaMemcpy(d_velocity_ref_z_m5l2h,d_velocity_ref_z,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

#else
	value_t *d_lx_buffer;
#endif

	// ===== SAVE INITIAL CONFIGURATION =====
	// TODO resolve the pre iteration set up and calls to verlet refresh
	verlet_size = verlet_refresh(d_position_x, d_position_y, d_position_z, d_lx_buffer, itep.step);

	force_calculation(d_position_x, d_position_y, d_position_z,
					  d_position_ref_x, d_position_ref_y, d_position_ref_z,
					  d_force_ref_x, d_force_ref_y, d_force_ref_z,
					  d_epot_tmp,
					  verlet_size,
					  d_verlet,
					  d_lx_buffer, 0);

	kinetic_energy_initial_calculation<<<n_block,thread_per_block>>>(d_velocity_x, d_velocity_y, d_velocity_z,
																	 d_velocity_ref_x, d_velocity_ref_y, d_velocity_ref_z,
																	 d_ekin_tmp);

	print_into_files(mol_pos_tmp, atm_pos_tmp,
					 mol_vel_tmp, atm_vel_tmp,
					 atm_fcs_tmp,
					 ekin_tmp, epot_tmp);

	// ===== ITERATION =====
	for (itep.step = 0; itep.step <= itep.max_step; itep.step++)
	{
		if(itep.step % itep.verlet_refresh == 0)
		{
			// TODO change the call to expanding verlet refresh
			verlet_size = verlet_refresh(d_position_x, d_position_y, d_position_z, d_lx_buffer, itep.step);
			//check if the .stp file exists
			if( access( outp.stp_name , R_OK ) != -1 )
			{
				File_buffer=fopen(outp.log_name,"a");
				printf("*** Simulation interupted by user! ***\n");
				fprintf(File_buffer,"*** Simulation interupted by user! ***\n");
				fclose(File_buffer);
				break;
			}
		}
// NOTE velocity predictor for expansion evaluation
#ifdef EXPANSION
		TRVP_execute_and_copy_velocities(n_block, thread_per_block,i,
										 d_velocity_x_TRVP_out, d_velocity_y_TRVP_out, d_velocity_z_TRVP_out,
										 d_velocity_x, d_velocity_y, d_velocity_z,
										 d_velocity_x_m3l2h, d_velocity_y_m3l2h, d_velocity_z_m3l2h,
										 d_velocity_x_m5l2h, d_velocity_y_m5l2h, d_velocity_z_m5l2h,
										 d_velocity_ref_x_TRVP_out, d_velocity_ref_y_TRVP_out, d_velocity_ref_z_TRVP_out,
										 d_velocity_ref_x, d_velocity_ref_y, d_velocity_ref_z,
										 d_velocity_ref_x_m3l2h, d_velocity_ref_y_m3l2h, d_velocity_ref_z_m3l2h,
										 d_velocity_ref_x_m5l2h, d_velocity_ref_y_m5l2h, d_velocity_ref_z_m5l2h);
		//cudaDeviceSynchronize();
		//L_from_density_file(t);
		//printf("Buffer:%lf, direct:%lf, diff:%.16g \n",boxLx_buffer[2*i],sysp.lx,boxLx_buffer[2*i]-sysp.lx);
		//calc_vf((value_t)t);
		//value_t pred=itep.vft;
		//calc_vf_buffer_host((value_t)t, boxLx_buffer);
		//printf("Buffer:%.16g, direct:%.16g, diff:%.16g \n",itep.vft,pred,pred-itep.vft);
		//calc_lambda(t+0.5*itep.dt);//lambda in t+0.5h
		//value_t pred=itep.lambda_tphh;
		//calc_lambda_buffer_host(t+0.5*itep.dt,boxLx_buffer);
		//printf("Buffer:%.16g, direct:%.16g, diff:%.16g \n",itep.lambda_tphh,pred,pred-itep.lambda_tphh);
		
		set_l(boxLx_buffer[2*i],boxLx_buffer[2*i],boxLx_buffer[2*i]);
		set_density_from_l(boxLx_buffer[2*i]);
#endif

		// force calculation DYNAMIC verlet list
		if (verlet_size != 0)
		{
			force_calculation(d_position_x, d_position_y, d_position_z,
							  d_position_ref_x, d_position_ref_y, d_position_ref_z,
							  d_force_ref_x, d_force_ref_y, d_force_ref_z,
							  d_epot_tmp,
							  verlet_size,
							  d_verlet,
							  d_lx_buffer, itep.step);
		}
		//cudaDeviceSynchronize(); // TODO test the necessity of the synchronization

// Kolafa alg section

#if defined EXPANSION
		kolafa_expansion_NoMM_graviti_center_corr<<<n_block,thread_per_block>>>(d_position_x, d_position_y, d_position_z,
													d_position_ref_x, d_position_ref_y, d_position_ref_z,
													d_velocity_x, d_velocity_y,  d_velocity_z,
													d_velocity_ref_x, d_velocity_ref_y, d_velocity_ref_z,
													d_force_ref_x, d_force_ref_y, d_force_ref_z,
													d_velocity_x_TRVP_out, d_velocity_y_TRVP_out,
													d_velocity_z_TRVP_out,d_velocity_ref_x_TRVP_out, d_velocity_ref_y_TRVP_out,
													d_velocity_ref_z_TRVP_out, d_ekin_tmp,d_lx_buffer,i);

#else
		// BUG may be incorrect for different GPU or larger number of atoms
//		kolafa_adjust<<<n_block,thread_per_block>>>(d_position_x, d_position_y, d_position_z,
//												 	d_position_ref_x, d_position_ref_y, d_position_ref_z,
//												 	d_velocity_x, d_velocity_y, d_velocity_z,
//												 	d_velocity_ref_x, d_velocity_ref_y, d_velocity_ref_z,
//												 	d_force_ref_x, d_force_ref_y, d_force_ref_z,
//												 	d_ekin_tmp); 
		
		kolafa_adjust_correction<<<n_block,thread_per_block>>>(d_position_x, d_position_y, d_position_z,
															   d_position_ref_x, d_position_ref_y, d_position_ref_z,
															   d_velocity_x, d_velocity_y, d_velocity_z,
															   d_velocity_ref_x, d_velocity_ref_y, d_velocity_ref_z,
															   d_force_ref_x, d_force_ref_y, d_force_ref_z,
															   d_ekin_tmp); // BUG not correct for different GPU or larger number of atoms
#endif

		if(itep.printout == true)
		{
			print_into_files(mol_pos_tmp, atm_pos_tmp,
							mol_vel_tmp, atm_vel_tmp,
							atm_fcs_tmp,
							ekin_tmp, epot_tmp);
		}

		t += itep.dt;
	}

	print_into_files(mol_pos_tmp, atm_pos_tmp,
					 mol_vel_tmp, atm_vel_tmp,
					 atm_fcs_tmp,
					 ekin_tmp, epot_tmp);


#ifdef EXPANSION
	free(boxLx_buffer);
#endif
	triplet_dealloc(mol_pos_tmp);
	triplet_dealloc(atm_pos_tmp);
	triplet_dealloc(mol_vel_tmp);
	triplet_dealloc(atm_vel_tmp);
	triplet_dealloc(atm_fcs_tmp);

}

/* Function: L_in_time
 * reads density from .box file, interpolates between points, calculates and returns box edhe length
 *
 * Notes:
 * >returns cubic box edge length [angstrom]
 *
 * Parameter:
 * t - value_t - input actual time [ps]
 */
value_t L_in_time(value_t t)
{
  typedef double vector[3];

  static struct history_s {
    double t;
    vector L;
  } one,hist[3];
  int k=0;
  static int nl=0;
  static FILE *rhof=NULL;
  char line[128];
  double rho;
  int i;


  rhof=fopen(inpp.box_file_name,"rt");

while (t>=hist[2].t || nl<3) {
    if (!fgets(line,128,rhof)) break;
    if (-1==-1) {
      if (sscanf(line,"%lf%lf",&one.t,&rho)<2) continue;
      one.L[0]=one.L[1]=one.L[2]=cbrt(sysp.system_real_mass /rho); }
    else
      if (sscanf(line,"%lf%lf%lf%lf",&one.t,one.L,one.L+1,one.L+2)<4) continue;
    hist[0]=hist[1]; hist[1]=hist[2]; hist[2]=one;
    nl++; }

  i=t>hist[1].t;
  fclose(rhof);
  return 0.1e11*(hist[i].L[k]+(t-hist[i].t)/(hist[i+1].t-hist[i].t)*(hist[i+1].L[k]-hist[i].L[k]));
}

/* Function: calc_vf
 * calculates v_f[ps^-1], for more info see paper, DOI: 10.1021/acs.jctc.8b00066
 *
 * Notes:
 * > calls L_in_time two times
 * > sets itep.vft
 *
 * Parameter:
 * t - value_t - input actual time [ps]
 */
void calc_vf(value_t t)
{
	itep.vft=log(L_in_time(t+0.500*itep.dt)/L_in_time(t-0.500*itep.dt))/itep.dt;
}

void calc_vf_buffer_host(value_t time, value_t* lx_buffer)
{
	unsigned int actual;
	actual=round(time/(0.5*itep.dt));
	if(actual<=1)
	{
		actual=1;
	}
	itep.vft=log(lx_buffer[actual+1]/lx_buffer[actual-1])/itep.dt;
}

/* Function: calc_lambda
 * calculates lambda[dimensionless], for more info see paper, DOI: 10.1021/acs.jctc.8b00066
 *
 * Notes:
 * > calls L_in_time two times
 * > sets itep.lambda_tphh
 *
 * Parameter:
 * t - value_t - input actual time [ps]
 */
void calc_lambda(value_t t)
{
	itep.lambda_tphh=L_in_time(t+0.500*itep.dt)/L_in_time(t-0.500*itep.dt);
	return;
}

void calc_lambda_buffer_host(value_t time, value_t* lx_buffer)
{
	unsigned int actual;
	actual=round(time/(0.5*itep.dt));
	if(actual<=1)
	{
		actual=1;
	}
	itep.lambda_tphh=lx_buffer[actual+1]/lx_buffer[actual-1];
}
/* Function: TRVP_execute_and_copy_velocities
 * performs TRVP and reduces number of copying history vectors.
 *
 * Notes:
 * > calls L_in_time two times
 *
 * Parameter:
 * n_block 					- unsigned int 	- number of blocks
 * thread_per_block 		- unsigned int 	- threads per block
 * i	 					- unsigned int	- actual step
 * d_velocityTRVP2_x_out	- * value_t		- array for storing molecular x velocity predicted by TRVP [angstrom]
 * d_velocityTRVP2_y_out	- * value_t		- array for storing molecular y velocity predicted by TRVP [angstrom]
 * d_velocityTRVP2_z_out	- * value_t		- array for storing molecular z velocity predicted by TRVP [angstrom]
 * d_velocity_x_in			- * value_t		- array of molecular x velocity in time t-0.5h [angstrom]
 * d_velocity_y_in			- * value_t		- array of molecular y velocity in time t-0.5h [angstrom]
 * d_velocity_z_in			- * value_t		- array of molecular z velocity in time t-0.5h [angstrom]
 * d_velocity_x_m3l2h_in	- * value_t		- array of molecular x velocity in time t-1.5h [angstrom]
 * d_velocity_y_m3l2h_in	- * value_t		- array of molecular y velocity in time t-1.5h [angstrom]
 * d_velocity_z_m3l2h_in	- * value_t		- array of molecular z velocity in time t-1.5h [angstrom]
 * d_velocity_x_m5l2h_in	- * value_t		- array of molecular x velocity in time t-2.5h [angstrom]
 * d_velocity_y_m5l2h_in	- * value_t		- array of molecular y velocity in time t-2.5h [angstrom]
 * d_velocity_z_m5l2h_in	- * value_t		- array of molecular z velocity in time t-2.5h [angstrom]
 * d_velocityTRVP2_ref_x_out- * value_t		- array for storing atom x velocity predicted by TRVP [angstrom]
 * d_velocityTRVP2_ref_y_out- * value_t		- array for storing atom y velocity predicted by TRVP [angstrom]
 * d_velocityTRVP2_ref_z_out- * value_t		- array for storing atom z velocity predicted by TRVP [angstrom]
 * d_velocity_ref_x_in		- * value_t		- array of atom x velocity in time t-0.5h [angstrom]
 * d_velocity_ref_y_in		- * value_t		- array of atom y velocity in time t-0.5h [angstrom]
 * d_velocity_ref_z_in		- * value_t		- array of atom z velocity in time t-0.5h [angstrom]
 * d_velocity_x_m3l2h_in	- * value_t		- array of atom x velocity in time t-1.5h [angstrom]
 * d_velocity_y_m3l2h_in	- * value_t		- array of atom y velocity in time t-1.5h [angstrom]
 * d_velocity_z_m3l2h_in	- * value_t		- array of atom z velocity in time t-1.5h [angstrom]
 * d_velocity_x_m5l2h_in	- * value_t		- array of atom x velocity in time t-2.5h [angstrom]
 * d_velocity_y_m5l2h_in	- * value_t		- array of atom y velocity in time t-2.5h [angstrom]
 * d_velocity_z_m5l2h_in	- * value_t		- array of atom z velocity in time t-2.5h [angstrom]
 */
void TRVP_execute_and_copy_velocities(unsigned int n_block, unsigned int thread_per_block,unsigned long i,
									   value_t *d_velocityTRVP2_x_out, value_t *d_velocityTRVP2_y_out, value_t *d_velocityTRVP2_z_out,
									   value_t *d_velocity_x_in, value_t *d_velocity_y_in, value_t *d_velocity_z_in,
									   value_t *d_velocity_x_m3l2h_in, value_t *d_velocity_y_m3l2h_in, value_t *d_velocity_z_m3l2h_in,
									   value_t *d_velocity_x_m5l2h_in, value_t *d_velocity_y_m5l2h_in, value_t *d_velocity_z_m5l2h_in,
									   value_t *d_velocityTRVP2_ref_x_out, value_t *d_velocityTRVP2_ref_y_out, value_t *d_velocityTRVP2_ref_z_out,
									   value_t *d_velocity_ref_x_in, value_t *d_velocity_ref_y_in, value_t *d_velocity_ref_z_in,
									   value_t *d_velocity_ref_x_m3l2h_in, value_t *d_velocity_ref_y_m3l2h_in, value_t *d_velocity_ref_z_m3l2h_in,
									   value_t *d_velocity_ref_x_m5l2h_in, value_t *d_velocity_ref_y_m5l2h_in, value_t *d_velocity_ref_z_m5l2h_in)
{
	if(i%2==0)
	{
		TRVP2_vel<<<n_block,thread_per_block>>>(d_velocityTRVP2_x_out, d_velocityTRVP2_y_out, d_velocityTRVP2_z_out,
												d_velocity_x_in, d_velocity_y_in, d_velocity_z_in,
											    d_velocity_x_m3l2h_in, d_velocity_y_m3l2h_in,d_velocity_z_m3l2h_in,
											    d_velocity_x_m5l2h_in, d_velocity_y_m5l2h_in, d_velocity_z_m5l2h_in);

		TRVP2_ref_vel<<<n_block,thread_per_block>>>(d_velocityTRVP2_ref_x_out, d_velocityTRVP2_ref_y_out, d_velocityTRVP2_ref_z_out,
													d_velocity_ref_x_in, d_velocity_ref_y_in, d_velocity_ref_z_in,
													d_velocity_ref_x_m3l2h_in, d_velocity_ref_y_m3l2h_in, d_velocity_ref_z_m3l2h_in,
													d_velocity_ref_x_m5l2h_in, d_velocity_ref_y_m5l2h_in, d_velocity_ref_z_m5l2h_in);

		cudaSafeCall( cudaMemcpy(d_velocity_x_m5l2h_in,d_velocity_x_in,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_y_m5l2h_in,d_velocity_y_in,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_z_m5l2h_in,d_velocity_z_in,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

		cudaSafeCall( cudaMemcpy(d_velocity_ref_x_m5l2h_in,d_velocity_ref_x_in,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_ref_y_m5l2h_in,d_velocity_ref_y_in,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_ref_z_m5l2h_in,d_velocity_ref_z_in,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

	}
	else
	{
		TRVP2_vel<<<n_block,thread_per_block>>>(d_velocityTRVP2_x_out, d_velocityTRVP2_y_out, d_velocityTRVP2_z_out,
												d_velocity_x_in, d_velocity_y_in, d_velocity_z_in,
											    d_velocity_x_m5l2h_in, d_velocity_y_m5l2h_in, d_velocity_z_m5l2h_in,
											    d_velocity_x_m3l2h_in, d_velocity_y_m3l2h_in,d_velocity_z_m3l2h_in);

		TRVP2_ref_vel<<<n_block,thread_per_block>>>(d_velocityTRVP2_ref_x_out, d_velocityTRVP2_ref_y_out, d_velocityTRVP2_ref_z_out,
													d_velocity_ref_x_in, d_velocity_ref_y_in, d_velocity_ref_z_in,
													d_velocity_ref_x_m5l2h_in, d_velocity_ref_y_m5l2h_in, d_velocity_ref_z_m5l2h_in,
													d_velocity_ref_x_m3l2h_in, d_velocity_ref_y_m3l2h_in, d_velocity_ref_z_m3l2h_in);

		cudaSafeCall( cudaMemcpy(d_velocity_x_m3l2h_in,d_velocity_x_in,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_y_m3l2h_in,d_velocity_y_in,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_z_m3l2h_in,d_velocity_z_in,sysp.molecule_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

		cudaSafeCall( cudaMemcpy(d_velocity_ref_x_m3l2h_in,d_velocity_ref_x_in,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_ref_y_m3l2h_in,d_velocity_ref_y_in,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );
		cudaSafeCall( cudaMemcpy(d_velocity_ref_z_m3l2h_in,d_velocity_ref_z_in,sysp.atom_count_aligned*sizeof(value_t),cudaMemcpyDeviceToDevice) );

	}

}
