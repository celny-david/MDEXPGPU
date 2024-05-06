/*
 ============================================================================
 Name        : output.cu
 Author      : klimam
 Version     :
 Description : Collects all parameters and functions relevant for creating and writing outputs

 NOTES : *
 ============================================================================
 */

#include <stdio.h>
#include <string.h> // strcat
#include <math.h>

#include "value.h"
#include "property.h"
#include "constants.h"
#include "system_property.h"
#include "iterative_property.h"
#include "output_property.h"
#include "device_manipulation_methods.h"

// required for mymalloc in save_to_cfg_directly 
#define TYPEOF(X)

typedef double vector[3];

// === LOCAL STRUCTURES FOR MACSIMUS OUTPUT ===
struct rec_s
{
    int key;
    int intval;
    vector vecval;
} rec;

struct a_s
{
    int  size;    // whole struct in bytes
    int  padd;    // padded to 8 bytes
    double logs;  // log of the Nose variable s
    vector lambda;// log(box)
    vector shape; // (reserved for box shape)
    vector rp[1]; // contiguous array of all r[ns] (and p[ns])
} *a[3];

struct spec_s
{
  int N;           // number of molecules
  int ns;          // number of sites
} *spec;           // [nspec]


// == ============== ==
// == SAVE TO ENERGY ==
// == ============== ==

/* Function: save_to_energy
 * the formated output of energy file (overloaded)
 * TOP level call method
 * 
 * write the kinetic energy, potential energy, total energy for one timestep
 *
 * Notes:
 * >overloaded with the number of molecules for size
 *
 * Parameter:
 * ekin_in - kinetic energy vector for one iteration input
 * epot_in - potential energy vector for one iteration input
 */
void save_to_energy(value_t* ekin_in, value_t* epot_in)
{
	FILE* File_bufer;
	File_bufer = fopen(outp.energy_name,"a"); // for clearing content of existing file

	value_t Ekin=Return_Ekin(ekin_in);
	value_t Epot=Return_Epot(epot_in);

#if PREC == FLOAT
	fprintf(File_bufer," %.8e %.8e %.8e %.8e\n",Ekin, Epot, Ekin + Epot,Ekin/sysp.molecule_count_aligned);
#elif PREC == DOUBLE
	fprintf(File_bufer," %.16e %.16e %.16e %.16e\n",Ekin, Epot, Ekin + Epot,Ekin/sysp.molecule_count_aligned);
#endif
	fclose(File_bufer);
	return;
}

// == =============== ==
// == SAVE TO HISTORY ==
// == =============== ==

/* Function: position_combiner
 * the helper tool for pre-formating output for history write-out
 *
 * combines data form molecule and atom position into position data
 *
 * Notes:
 * > assume the input arrays are of length sysp.atom_count
 * > performs simple addition of CoM to all corresponding atoms in molecule:
 *   this is done by integer division of index by atoms per molecule
 *
 * Parameter:
 * position_in 			 - triplet array of positions used to store output final position of atoms
 * molecule_postition_in - triplet array of molecular CoM
 * atom_position_in 	 -  triplet array of atomal positions
 */
void position_combiner(triplet *position_in,triplet *molecule_postition_in,triplet *atom_position_in)
{
	unsigned int i; //counters
	for(i=0; i<sysp.atom_count_aligned; i++ ) // NOTE assume the input arrays are of length sysp.atom_count
	{// expectation that atom position holds 0.0 at position of reference atom
	 // or that the molecule position is CoM and atom positions are differences of CoM
		position_in->x[i] = molecule_postition_in->x[i/sysp.atom_per_molecule] + atom_position_in->x[i];
		position_in->y[i] = molecule_postition_in->y[i/sysp.atom_per_molecule] + atom_position_in->y[i];
		position_in->z[i] = molecule_postition_in->z[i/sysp.atom_per_molecule] + atom_position_in->z[i];
	}
}

/* Function: save_to_history
 * the formated output of history file (overloaded for CPU styled input)
 * TOP LEVEL CALL method
 * write the positions of the atoms
 *
 * Notes:
 * > has check for existence of history file
 *
 * Parameter:
 * molecule_postition_in - triplet array of molecular CoM
 * atom_position_in 	 - triplet array of atomal positions
 * timestep		 		 - the iteration time of the printed snapshot
 */
void save_to_history(triplet *molecule_postition_in, triplet *atom_position_in)
{
	triplet* position_in = (triplet*) malloc(sizeof(triplet));
	FILE* File_buffer;
	File_buffer = fopen(outp.history_name,"a");
	triplet_alloc(position_in,sysp.atom_count_aligned);
	position_combiner(position_in,molecule_postition_in,atom_position_in);
	
#if SUBSTANCE == ARGON
	const char *def_name[1] = {"Ar"};
#elif SUBSTANCE == NITROGEN
	const char *def_name[2] = {"N ", "N "};
#elif SUBSTANCE == SPCE
	const char *def_name[3] = {"O ", "H ", "H "};
#endif
	fprintf(File_buffer,"%i\n",sysp.atom_count_aligned);
	fprintf(File_buffer,"#=== time: %f [ps] ===\n",itep.step*itep.dt);
	for(unsigned int i=0; i<sysp.atom_count_aligned; i++ )
	{
		fprintf(File_buffer,"%s %f %f %f\n",def_name[i%SUBSTANCE],position_in->x[i],position_in->y[i],position_in->z[i]);
	}
	fclose(File_buffer);	
	triplet_dealloc(position_in);
}

// == ======================== ==
// == SAVE TO HISTORY SEPARATE ==
// == ======================== ==

/* Function: save_to_history_separate
 * the formated output of history file (overloaded for CPU styled input)
 * TOP LEVEL CALL method
 * write the positions of the atoms
 *
 * Notes:
 * > has check for existence of history file
 *
 * Parameter:
 * molecule_postition_in - triplet array of molecular CoM
 * atom_position_in 	 -  triplet array of atomal positions
 * timestep		 		 - the iteration time of the printed snapshot
 */
void save_to_history_separate(triplet *molecule_postition_in, triplet *atom_position_in)
{

	FILE* File_buffer;
	File_buffer = fopen(outp.history_name,"a");

#if SUBSTANCE == ARGON
	const char *def_name[1] = {"Ar"};
#elif SUBSTANCE == NITROGEN
	const char *def_name[2] = {"N ", "N "};
#elif SUBSTANCE == SPCE
	const char *def_name[3] = {"O ", "H ", "H "};
#endif
	fprintf(File_buffer,"%i\n",sysp.atom_count_aligned);
	fprintf(File_buffer,"#=== time: %f [ps] ===\n",itep.step*itep.dt);
	for(unsigned int i=0; i<sysp.atom_count_aligned; i++ )
	{
		fprintf(File_buffer,"%s %f %f %f %f %f %f\n", def_name[i%SUBSTANCE],
													  molecule_postition_in->x[i/sysp.atom_per_molecule],
													  molecule_postition_in->y[i/sysp.atom_per_molecule],
													  molecule_postition_in->z[i/sysp.atom_per_molecule],
													  atom_position_in->x[i],
													  atom_position_in->y[i],
													  atom_position_in->z[i]);

	}
	fclose(File_buffer);		
}

// == ============ ==
// == SAVE TO TEST ==
// == ============ ==

/* Function: save_to_test
 * the formated output of test file
 *
 * write the values of the atoms i.e. speeds/forces
 *
 * Notes:
 * > macro for the molecule names
 *
 * Parameter:
 * where	  - the file pointer to where it should be save
 * value_in_x - the x value i.e. speed of the atom written to the file
 * value_in_y - the y value i.e. speed of the atom written to the file
 * value_in_z - the z value i.e. speed of the atom written to the file
 * timestep	  - the iteration time of the printed snapshot
 */
void general_file_output_internal(FILE* where, value_t* value_in_x, value_t* value_in_y, value_t* value_in_z)
{
#if SUBSTANCE == ARGON
	const char *def_name[1] = {"Ar"};
#elif SUBSTANCE == NITROGEN
	const char *def_name[2] = {"N ", "N "};
#elif SUBSTANCE == SPCE
	const char *def_name[3] = {"O ", "H ", "H "};
#endif
	value_t mean_norm=0.0;

	fprintf(where,"%i\n",sysp.atom_count_aligned);
	fprintf(where,"#=== time: %f [ps] ===\n",itep.step*itep.dt);
	for(unsigned int i=0; i<sysp.atom_count_aligned; i++ ) // NOTE assume the input arrays are of length sysp.atom_count
	{
		fprintf(where,"%s %d %f %f %f\n",def_name[i%SUBSTANCE],i,value_in_x[i],value_in_y[i],value_in_z[i]);

		mean_norm += sqrt(value_in_x[i]*value_in_x[i]
					 	 +value_in_y[i]*value_in_y[i]
						 +value_in_z[i]*value_in_z[i]);
	}
	fprintf(where,"#=== mean norm of previous time-step: %f ===\n",mean_norm/sysp.atom_count_aligned);
}

/* Function: save_to_test
 * the formated output of test file (overloaded for CPU styled input)
 * TOP LEVEL CALL
 * write the values of the atoms i.e. speeds/forces
 *
 * this version does not account for CoM approach (does not implicitly call the combiner)
 *
 * Notes:
 * > macro for the molecule names
 *
 * Parameter:
 * where	  - the file pointer to where it should be save
 * value_in_x - the x value i.e. speed of the atom written to the file
 * value_in_y - the y value i.e. speed of the atom written to the file
 * value_in_z - the z value i.e. speed of the atom written to the file
 * timestep	  - the iteration time of the printed snapshot
 */
void general_file_output(FILE* where, triplet *atom_property_in)
{
	if (where != NULL)
	{
		general_file_output_internal(where,atom_property_in->x, atom_property_in->y, atom_property_in->z);
	}
}

// == ============ ==
// == SAVE TO REST ==
// == ============ ==

/* Function: save_to_vel
 * todo write the description
 */
void save_to_vel(triplet *molecule_velocity_in, triplet *atom_velocity_in)
{
	FILE * File_buffer;
	strcpy(outp.vel_name,outp.output_files_name);
	strcat(outp.vel_name,".vel");
	File_buffer = fopen(outp.vel_name,"a");

	triplet* velocity_in = (triplet*) malloc(sizeof(triplet));
	triplet_alloc(velocity_in,sysp.atom_count_aligned);
	position_combiner(velocity_in,molecule_velocity_in,atom_velocity_in);
	general_file_output_internal(File_buffer,velocity_in->x, velocity_in->y, velocity_in->z);
	triplet_dealloc(velocity_in);

	fclose(File_buffer);
}


/* Function: save_to_plb_directly
 * appends one frame of atom positions to .plb file
 *
 * Notes:
 * > The .plb file is MACSIMUS standard binary file containing atom positions.
 * > It is important to have .mol file to show it using
 *
 * Parameter:
 * boxsize_x - float size of simulation box in x axis [Angstrom]
 * boxsize_y - float size of simulation box in y axis [Angstrom]
 * boxsize_z - float size of simulation box in z axis [Angstrom]
 * position_in_x - value_t* pointer to vector of atom positions in x axis [Angstrom]
 * position_in_y - value_t* pointer to vector of atom positions in y axis [Angstrom]
 * position_in_z - value_t* pointer to vector of atom positions in z axis [Angstrom]
 */

void save_to_plb_directly(float boxsize_x , float boxsize_y, float boxsize_z, value_t* position_in_x, value_t* position_in_y, value_t* position_in_z)
{

	FILE * File_buffer;
	File_buffer = fopen(outp.plb_name,"ab");


#if	PREC == DOUBLE
	float buffer_position_in_x,buffer_position_in_y,buffer_position_in_z;
#endif

	fwrite(&boxsize_x,1,sizeof(boxsize_x),File_buffer);
	fwrite(&boxsize_y,1,sizeof(boxsize_y),File_buffer);
	fwrite(&boxsize_z,1,sizeof(boxsize_z),File_buffer);

	for(int i=0;i<sysp.atom_count_aligned;i++)
	{
#if	PREC == DOUBLE
		buffer_position_in_x=(float)position_in_x[i];
		buffer_position_in_y=(float)position_in_y[i];
		buffer_position_in_z=(float)position_in_z[i];

		fwrite(&buffer_position_in_x,1,sizeof(float),File_buffer);
		fwrite(&buffer_position_in_y,1,sizeof(float),File_buffer);
		fwrite(&buffer_position_in_z,1,sizeof(float),File_buffer);
#else
		fwrite(&position_in_x[i],1,sizeof(float),File_buffer);
		fwrite(&position_in_y[i],1,sizeof(float),File_buffer);
		fwrite(&position_in_z[i],1,sizeof(float),File_buffer);
#endif
	}

	fclose(File_buffer);
	return;
}

/* Function: save_to_plb
 * obtains atom positions calling position_combiner
 * calls save_to_plb_directly to append frame to .plb file
 *
 * Notes:
 * > The .plb file is MACSIMUS standard binary file containing atom positions.
 * > It is important to have .mol file to show it using
 * > Uses pointers to structure triplet(contains three value_t* pointers to vectors of positions in x,y,z).
 *
 * Parameter:
 * molecule_postition_in - triplet* molecular center of gravity coordinates [Angstrom]
 * atom_position_in - triplet* atom coordinates in molecule [Angstrom]
 */
void save_to_plb(triplet *molecule_postition_in, triplet *atom_position_in)
{
	triplet* position_in = (triplet*) malloc(sizeof(triplet));

	triplet_alloc(position_in, sysp.atom_count_aligned);
	position_combiner(position_in, molecule_postition_in, atom_position_in);
	save_to_plb_directly((float) sysp.lx, (float) sysp.ly, (float)sysp.lz, position_in->x, position_in->y, position_in->z);
	triplet_dealloc(position_in);
	return;
}

/* Function: find_outlying_molecules
 * Test if atoms are in the simulation box.
 * If the atom is out, following data are written:
 * 		serial number of atom
 * 		atom coordinates
 * 		serial number of its molecule
 * 		coordinates of molecular center of gravity
 * 		test if the molecular center of gravity is in the box
 *
 * Notes:
 * > Uses pointers to structure triplet(contains three value_t* pointers to vectors of positions in x,y,z).
 * > Calls position_combiner
 *
 * Parameter:
 * molecule_postition_in - triplet* molecular center of gravity coordinates [Angstrom]
 * atom_position_in - triplet* atom coordinates in molecule [Angstrom]
 * frame_number - int serial number of actual frame
 */

/* Function: save_to_oob
 * todo write/edit the description
 */
void save_to_oob(triplet *molecule_postition_in, triplet *atom_position_in)
{	
	FILE * File_buffer;
	File_buffer=fopen(outp.oob_name,"a");

	int molecule_number=0, first_atom_number=0;
	float gravity_center_x, gravity_center_y, gravity_center_z;
	triplet* position_in = (triplet*) malloc(sizeof(triplet));
	triplet_alloc(position_in,sysp.atom_count_aligned);
	position_combiner(position_in,molecule_postition_in,atom_position_in);

	fprintf(File_buffer,"***** FRAME NUMBER: %lu *****\n",itep.step);

	for(int i=0; i<sysp.atom_count_aligned; i++ )
	{
		if(position_in->x[i] > sysp.lx || position_in->y[i] > sysp.ly || position_in->z[i] > sysp.lz)
		{
			molecule_number=(i-(i%sysp.atom_per_molecule))/sysp.atom_per_molecule;

#if SUBSTANCE == SPCE

			first_atom_number=molecule_number*3;
			gravity_center_x=((position_in->x[first_atom_number]*MASS_O)+(position_in->x[first_atom_number+1]*MASS_H)+(position_in->x[first_atom_number+2]*MASS_H))/(MASS_O+2*MASS_H);
			gravity_center_y=((position_in->y[first_atom_number]*MASS_O)+(position_in->y[first_atom_number+1]*MASS_H)+(position_in->y[first_atom_number+2]*MASS_H))/(MASS_O+2*MASS_H);
			gravity_center_z=((position_in->z[first_atom_number]*MASS_O)+(position_in->z[first_atom_number+1]*MASS_H)+(position_in->z[first_atom_number+2]*MASS_H))/(MASS_O+2*MASS_H);

			if(i%sysp.atom_per_molecule==0)
			{
				fprintf(File_buffer," O  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",
										i,
										molecule_number,
										atom_position_in->x[i],
										atom_position_in->y[i],
										atom_position_in->z[i],
										molecule_postition_in->x[molecule_number],
										molecule_postition_in->y[molecule_number],
										molecule_postition_in->z[molecule_number],
										gravity_center_x,
										gravity_center_y,
										gravity_center_z);
			}
			else
			{
				fprintf(File_buffer," H  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",
										i,
										molecule_number,
										atom_position_in->x[i],
										atom_position_in->y[i],
										atom_position_in->z[i],
										molecule_postition_in->x[molecule_number],
										molecule_postition_in->y[molecule_number],
										molecule_postition_in->z[molecule_number],
										gravity_center_x,
										gravity_center_y,
										gravity_center_z);
			}
			for(int printed_atom=first_atom_number; printed_atom<(first_atom_number+sysp.atom_per_molecule); printed_atom++)
			{
				if(printed_atom%sysp.atom_per_molecule==0)
				{
					fprintf(File_buffer,"*O  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",
											printed_atom,
											molecule_number,
											atom_position_in->x[i],
											atom_position_in->y[i],
											atom_position_in->z[i],
											molecule_postition_in->x[molecule_number],
											molecule_postition_in->y[molecule_number],
											molecule_postition_in->z[molecule_number],
											gravity_center_x,
											gravity_center_y,
											gravity_center_z);
				}
				else
				{
					fprintf(File_buffer,"*H  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",
											printed_atom,
											molecule_number,
											atom_position_in->x[i],
											atom_position_in->y[i],
											atom_position_in->z[i],
											molecule_postition_in->x[molecule_number],
											molecule_postition_in->y[molecule_number],
											molecule_postition_in->z[molecule_number],
											gravity_center_x,
											gravity_center_y,
											gravity_center_z);
				}
			}

#elif	SUBSTANCE == ARGON

			fprintf(File_buffer,"Ar  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f\n",
									i,
									molecule_number,
									position_in->x[i],
									position_in->y[i],
									position_in->z[i],
									molecule_postition_in->x[molecule_number],
									molecule_postition_in->y[molecule_number],
									molecule_postition_in->z[molecule_number]);

#elif	SUBSTANCE == NITROGEN

			first_atom_number = molecule_number*2;
			gravity_center_x = ((position_in->x[first_atom_number]*MASS_N) + (position_in->x[first_atom_number+1]*MASS_N))/(2*MASS_N);
			gravity_center_y = ((position_in->y[first_atom_number]*MASS_N) + (position_in->y[first_atom_number+1]*MASS_N))/(2*MASS_N);
			gravity_center_z = ((position_in->z[first_atom_number]*MASS_N) + (position_in->z[first_atom_number+1]*MASS_N))/(2*MASS_N);

			fprintf(File_buffer," N  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",
									i,
									molecule_number,
									position_in->x[i],
									position_in->y[i],
									position_in->z[i],
									molecule_postition_in->x[molecule_number],
									molecule_postition_in->y[molecule_number],
									molecule_postition_in->z[molecule_number],
									gravity_center_x,
									gravity_center_y,
									gravity_center_z);

			for(int printed_atom=first_atom_number; printed_atom<(first_atom_number+sysp.atom_per_molecule); printed_atom++)
			{
				fprintf(File_buffer,"*N  %4d  %5d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",
										printed_atom,
										molecule_number,
										position_in->x[printed_atom],
										position_in->y[printed_atom],
										position_in->z[printed_atom],
										molecule_postition_in->x[molecule_number],
										molecule_postition_in->y[molecule_number],
										molecule_postition_in->z[molecule_number],
										gravity_center_x,
										gravity_center_y,
										gravity_center_z);
			}

#endif
		} // END if(position_in->x[i] > sysp.lx || position_in->y[i] > sysp.ly || position_in->z[i] > sysp.lz)
	} // END for(int i=0; i<sysp.atom_count_aligned; i++ )

	triplet_dealloc(position_in);
	fclose(File_buffer);
	return;
}

void save_to_test4(triplet *molecule_postition_in, triplet *atom_position_in)
{
	FILE * File_buffer;
	File_buffer=fopen(outp.test4_name,"a");

	char character[5];
	int molecule_number=0, first_atom_number=0;
	float gravity_center_x, gravity_center_y, gravity_center_z;
	triplet* position_in = (triplet*) malloc(sizeof(triplet));
	triplet_alloc(position_in,sysp.atom_count_aligned);
	position_combiner(position_in,molecule_postition_in,atom_position_in);
	float owerlap=1.5f;
	float enabledError=1e-4;

#ifdef EXPANSION
	fprintf(File_buffer,"***** TIME: %f, STEP: %lu, OVERLAP: %f, TOLERANCE: %f, BOX: L_x=%f  L_y=%f  L_z=%f, DENSITY: %f, LAMBDA: %f, v_ft: %f *****\n",itep.step*itep.dt, itep.step, owerlap, enabledError, sysp.lx, sysp.ly, sysp.lz, sysp.density, itep.lambda_tphh, itep.vft);
#else
	fprintf(File_buffer,"***** TIME: %f, STEP: %lu, OVERLAP: %f, TOLERANCE: %f, BOX: L_x=%f  L_y=%f  L_z=%f, DENSITY: %f *****\n",itep.step*itep.dt, itep.step, owerlap, enabledError, sysp.lx, sysp.ly, sysp.lz, sysp.density);
#endif


	for(int i=0; i<sysp.atom_count_aligned; i++ )
	{
			molecule_number=(i-(i%sysp.atom_per_molecule))/sysp.atom_per_molecule;
			strcpy(character,"    ");
			if(abs(position_in->x[i]) > (sysp.lx+owerlap) || abs(position_in->y[i]) > (sysp.ly+owerlap) || abs(position_in->z[i]) > (sysp.lz+owerlap))
			{
				character[0]='#';
			}
			if(abs(position_in->x[i]) > sysp.lx || abs(position_in->y[i]) > sysp.ly || abs(position_in->z[i]) > sysp.lz)
			{
				character[1]='*';
			}
#if SUBSTANCE == SPCE
			if((MASS_O*atom_position_in->x[SUBSTANCE*molecule_number]+MASS_H*atom_position_in->x[SUBSTANCE*molecule_number+1]+MASS_H*atom_position_in->x[SUBSTANCE*molecule_number+2]) > enabledError || (MASS_O*atom_position_in->y[SUBSTANCE*molecule_number]+MASS_H*atom_position_in->y[SUBSTANCE*molecule_number+1]+MASS_H*atom_position_in->y[SUBSTANCE*molecule_number+2]) > enabledError || (MASS_O*atom_position_in->z[SUBSTANCE*molecule_number]+MASS_H*atom_position_in->z[SUBSTANCE*molecule_number+1]+MASS_H*atom_position_in->z[SUBSTANCE*molecule_number+2]) > enabledError)
			{
				character[2]='x';
			}
#elif SUBSTANCE == NITROGEN
			if((MASS_N*atom_position_in->x[SUBSTANCE*molecule_number]+MASS_N*atom_position_in->x[SUBSTANCE*molecule_number+1]) > enabledError || (MASS_N*atom_position_in->y[SUBSTANCE*molecule_number]+MASS_N*atom_position_in->y[SUBSTANCE*molecule_number+1]) > enabledError || (MASS_N*atom_position_in->z[SUBSTANCE*molecule_number]+MASS_N*atom_position_in->z[SUBSTANCE*molecule_number+1]) > enabledError)
			{
				character[2]='x';
			}

#endif
			else
			{
				strcpy(character,"    ");
			}

#if SUBSTANCE == SPCE

			first_atom_number=molecule_number*3;
			gravity_center_x=((position_in->x[first_atom_number]*MASS_O)+(position_in->x[first_atom_number+1]*MASS_H)+(position_in->x[first_atom_number+2]*MASS_H))/(MASS_O+2*MASS_H);
			gravity_center_y=((position_in->y[first_atom_number]*MASS_O)+(position_in->y[first_atom_number+1]*MASS_H)+(position_in->y[first_atom_number+2]*MASS_H))/(MASS_O+2*MASS_H);
			gravity_center_z=((position_in->z[first_atom_number]*MASS_O)+(position_in->z[first_atom_number+1]*MASS_H)+(position_in->z[first_atom_number+2]*MASS_H))/(MASS_O+2*MASS_H);

			if(i%sysp.atom_per_molecule==0)
			{
				fprintf(File_buffer,"%sO  %4d  %5d  %10f  %10f  %10f  %10f %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n",character,i,molecule_number,
															atom_position_in->x[i],
															atom_position_in->y[i],
															atom_position_in->z[i],
															position_in->x[i],
															position_in->y[i],
															position_in->z[i],
															molecule_postition_in->x[molecule_number],
															molecule_postition_in->y[molecule_number],
															molecule_postition_in->z[molecule_number],
															gravity_center_x,gravity_center_y,gravity_center_z);
			}
			else
			{
				fprintf(File_buffer,"%sH  %4d  %5d  %10f  %10f  %10f  %10f  %10f %10f  %10f  %10f %10f  %10f  %10f  %10f\n",character,i,molecule_number,
																			atom_position_in->x[i],
																			atom_position_in->y[i],
																			atom_position_in->z[i],
																			position_in->x[i],
																			position_in->y[i],
																			position_in->z[i],
																			molecule_postition_in->x[molecule_number],
																			molecule_postition_in->y[molecule_number],
																			molecule_postition_in->z[molecule_number],
																			gravity_center_x,gravity_center_y,gravity_center_z);
			}



#elif	SUBSTANCE == ARGON

			fprintf(File_buffer,"%sAr  %4d  %5d  %10f  %10f  %10f  %10f  %10f %10f  %10f  %10f %10f  %10f  %10f  %10f\n",character,i,molecule_number,
																		atom_position_in->x[i],
																		atom_position_in->y[i],
																		atom_position_in->z[i],
																		position_in->x[i],
																		position_in->y[i],
																		position_in->z[i],
																		molecule_postition_in->x[molecule_number],
																		molecule_postition_in->y[molecule_number],
																		molecule_postition_in->z[molecule_number],
																		0.0f,0.0f,0.0f);


#elif	SUBSTANCE == NITROGEN

			first_atom_number=molecule_number*2;
			gravity_center_x=((position_in->x[first_atom_number]*MASS_N)+(position_in->x[first_atom_number+1]*MASS_N))/(2*MASS_N);
			gravity_center_y=((position_in->y[first_atom_number]*MASS_N)+(position_in->y[first_atom_number+1]*MASS_N))/(2*MASS_N);
			gravity_center_z=((position_in->z[first_atom_number]*MASS_N)+(position_in->z[first_atom_number+1]*MASS_N))/(2*MASS_N);


			fprintf(File_buffer,"%sN  %4d  %5d  %10f  %10f  %10f  %10f  %10f %10f  %10f  %10f %10f  %10f  %10f  %10f\n",character,i,molecule_number,
																		atom_position_in->x[i],
																		atom_position_in->y[i],
																		atom_position_in->z[i],
																		position_in->x[i],
																		position_in->y[i],
																		position_in->z[i],
																		molecule_postition_in->x[molecule_number],
																		molecule_postition_in->y[molecule_number],
																		molecule_postition_in->z[molecule_number],
																		gravity_center_x,gravity_center_y,gravity_center_z);

#endif

	}

	triplet_dealloc(position_in);
	fclose(File_buffer);
	return;
}

/* Function: VarPut
 * writes wariable into .cfg file
 *
 * Notes:
 * > Function was taken from MACSIMUS utilities
 *
 * Parameters:
 * file_name - char * name of .cfg file including sufix
 * v - void * data to write in file
 * size - int size of v
 */
void VarPut(char * file_name, void *v, int size)
{
	int i,s;
	unsigned char a;
	FILE* File_buffer;
	File_buffer= fopen(file_name,"ab");

	if (!size) printf("0 record to write");

	fwrite(&size, sizeof(size), 1, File_buffer);
	s = 0x4A4B;
	
	for(i=0;i<size;i++) //loop(I,FROM,TO) for ((I)=(FROM); (I)<(TO); (I)++)
	{
		a = *((unsigned char*)v+i);
		s += a;
		putc(a,File_buffer);
	};
	if (fwrite(&s, sizeof(s), 1, File_buffer)!=1) printf("write");
	
	fclose(File_buffer);
}

/* Function: VarClose
 * closes output .cfg file
 *
 * Notes:
 * > Function was taken from MACSIMUS utilities
 *
 * Parameters:
 * file_name - char* name of .cfg file
 */
void VarClose(char * file_name)
{
	int z=0;
	FILE* File_buffer;
	File_buffer= fopen(file_name,"ab");
	if (fwrite(&z, sizeof(z), 1, File_buffer)!=1) printf("write EOF");
	if (fclose(File_buffer)) printf("close");
}


/* Function: mymalloc
 * alocates structure a
 *
 * Notes:
 * > Function was taken from MACSIMUS utilities without deeper understanding !!!
 *
 * Parameter:
 * size - int size of allocated structure a in bytes
 * zero - int unknown
 */
a_s *mymalloc(int size,int zero)
{
	a_s *a;

	if (size>0x40000000)
	{
		printf("AllocSizeLim exceeded");
		return NULL;
	}
	if ((a=(a_s*)malloc(size)) == NULL)
	{
		printf("malloc failed");
		return NULL;
	}
	if (zero) memset(a,0,size);

	return a;
}

/* Function: save_to_cfg_directly
 * writes configuration into .cfg file.
 * File contains:
 * 		Number of molecules
 * 		Number of sights in molecule
 * 		Number of atoms
 * 		Total energy [K]
 * 		Running simulation time [ps]
 * 		Time step [ps]
 * 		Box sizes in x,y,z (x=y=z default) [Angstrom]
 * 		Thermostat type (0 default)
 * 		Atom positions x,y,z [Angstrom]
 * 		Atom velocities multiplied by time step x,y,z [UNKNOWN]
 *		Atom accelerations multiplied by time step squared x,y,z [UNKNOWN] (0 default)
 *
 * Notes:
 * > The .cfg file is MACSIMUS standard binary file.
 * > Uses pointers to structure triplet(contains three value_t* pointers to vectors of positions in x,y,z).
 *
 * Parameters:
 * position_in - triplet* pointer to atom positions [Angstrom]
 * velocity_in - triplet* pointer to atom velocities [UNKNOWN]
 * Etot - double total energy [K]
 * h - double time step [ps]
 * sim_time - double running simulation time [ps]
 * box_edge_length - double length of simulation box edge [Angstrom] (cubic box expected)
 */
void save_to_cfg_directly(char *file_name,triplet* position_in, triplet* velocity_in, double Etot, double h, double sim_time, double box_edge_length)
{
    //h= time step [ps]
    //sim_time= running simulation time [ps]
    //box_edge_length= box sizes [Angstrom]
    //Etot= total energy [K]

    int ns=0;
    int size=0;
    int np=0;
    int i;

////First part of header
    int cfgkey=1; 						//1 for nonpolar version, 2 for polar version
    int options[32]={-9,100,1,9,
    				 0,-1,0,-32333,
    				 0,2,0,-9,
    				 0,2,0,1,
    				 229,100,1,0,
    				 0,0,3,1,
    				 6,2147483647,0,0,
    				 -32333,0,0,-32333}; //unknown meaning
    int nspec=1; 						 //number of molecule types
    np = cfgkey&2 ? 2 : 1;

    int specie_specification[2]; 							 //for each molecule type
        specie_specification[0]=sysp.molecule_count_aligned; //number of molecules
        specie_specification[1]=sysp.atom_per_molecule; //number of sites in molecule

    VarPut(file_name,&cfgkey,sizeof(cfgkey));
    VarPut(file_name,&options,sizeof(options));
    VarPut(file_name,&nspec,sizeof(nspec));
    VarPut(file_name,&specie_specification,sizeof(specie_specification));

////First type record specifications
	rec.key=1;      						//key for record type specification
	rec.intval=sysp.molecule_count_aligned; //Number of molecules
	rec.vecval[0]=sim_time; 				//running simulation time [ps]
	rec.vecval[1]=h; 						//time step [ps]
	rec.vecval[2]=Etot; 					//total energy [K]

	VarPut(file_name,&rec,sizeof(rec));

////Second type record specifications
	rec.key=2; 							//key for record type specification
	rec.intval=sysp.atom_count_aligned; //Number of atoms
	rec.vecval[0]=box_edge_length;		//box sizes [Angstrom]
	rec.vecval[1]=box_edge_length;		//box sizes [Angstrom]
	rec.vecval[2]=box_edge_length;		//box sizes [Angstrom]
	ns=rec.intval;
	VarPut(file_name,&rec,sizeof(rec));

////Thermostat specifications
	rec.key=4;		 //key for record type specification
	rec.intval=0; 	 //Thermostat type
	rec.vecval[0]=0; //RvdW (for van der Waals radius setup)
	rec.vecval[1]=0; //tau.T
	rec.vecval[2]=0; //tau.P

    VarPut(file_name,&rec,sizeof(rec));

////alocation of structures a[i]
	size=sizeof(*a[0])+sizeof(vector)*(ns*np-1);
	for (i=0;i<3;i++)
	{
		(a[i])=TYPEOF((a[i])) mymalloc(size,1);
		*(int*)(a[i])=size;
	}
////logos
    a[0]->logs=0;
    a[1]->logs=0;
    a[2]->logs=0;
////ln(L)=lambda
    a[0]->lambda[0]=log(box_edge_length);
    a[0]->lambda[1]=log(box_edge_length);
    a[0]->lambda[2]=log(box_edge_length);
    a[1]->lambda[0]=0;
    a[1]->lambda[1]=0;
    a[1]->lambda[2]=0;
    a[2]->lambda[0]=0;
    a[2]->lambda[1]=0;
    a[2]->lambda[2]=0;
////setting of positions, velocities and accelerations
for (i=0;i<ns*np;i++) {
//// positions x,y,z
    a[0]->rp[i][0]=(double) position_in->x[i];
    a[0]->rp[i][1]=(double) position_in->y[i];
    a[0]->rp[i][2]=(double) position_in->z[i];
//// h*velocities x,y,z
    a[1]->rp[i][0]=(double) velocity_in->x[i];
    a[1]->rp[i][1]=(double) velocity_in->y[i];
    a[1]->rp[i][2]=(double) velocity_in->z[i];
//// h*h*accelerations x,y,z
    a[2]->rp[i][0]=(double) 0.f;
    a[2]->rp[i][1]=(double) 0.f;
    a[2]->rp[i][2]=(double) 0.f;
  }

//// write into cfg file
	VarPut(file_name,a[0],size);
	VarPut(file_name,a[1],size);
	VarPut(file_name,a[2],size);
	VarClose(file_name);
	free(a[2]);
	free(a[1]);
	free(a[0]);
	return;
}

/* Function: save_to_cfa_directly
 * writes configuration into .cfg file.
 * File contains:
 * 		Number of molecules
 * 		Number of sights in molecule
 * 		Number of atoms
 * 		Total energy [K]
 * 		Running simulation time [ps]
 * 		Time step [ps]
 * 		Box sizes in x,y,z (x=y=z default) [Angstrom]
 * 		Thermostat type (0 default)
 * 		Atom positions x,y,z [Angstrom]
 * 		Atom velocities multiplied by time step x,y,z [UNKNOWN]
 *		Atom accelerations multiplied by time step squared x,y,z [UNKNOWN] (0 default)
 *
 * Notes:
 * > The .cfa file is MACSIMUS standard text file equivalent to .cfg binary file
 * > Uses pointers to structure triplet(contains three value_t* pointers to vectors of positions in x,y,z).
 *
 *
 * position_in - triplet* pointer to atom positions [Angstrom]
 * velocity_in - triplet* pointer to atom velocities [UNKNOWN]
 * Etot - double total energy [K]
 * h - double time step [ps]
 * sim_time - double running simulation time [ps]
 * box_edge_length - double length of simulation box edge [Angstrom] (cubic box expected)
 */
 void save_to_cfa_directly(FILE* File_buffer,triplet* position_in, triplet* velocity_in, double Etot, double h, double sim_time, double box_edge_length)
 {
	int optionlist[32]={-9,100,1,9,
						0,-1,0,-32333,
						0,2,0,-9,
						0,2,0,1,
						229,100,1,0,
						0,0,3,1,
						6,2147483647,0,0,
						-32333,0,0,-32333};
	char a='a';

	fprintf(File_buffer,"1 cfgkey (nonpolar)\n");
	fprintf(File_buffer,"-@%d\n",optionlist[0]);
	for(int i=0;i<26;i++)
	{
		fprintf(File_buffer,"-%c%d\n",a+i,optionlist[i+1]);
	}
	a='[';
	for(int i=0;i<5;i++)
	{
		fprintf(File_buffer,"-%c%d\n",a+i,optionlist[i+27]);
	}


	fprintf(File_buffer,"1 nspec\n");
	fprintf(File_buffer,"%d %d N ns (of species 0)\n",sysp.molecule_count_aligned,sysp.atom_per_molecule);//molecules, sights in molecule
	fprintf(File_buffer,"%d %d %.16g %.16g %.16g key N t h En.tot\n",1, sysp.molecule_count_aligned, sim_time, h, Etot);
	fprintf(File_buffer,"%d %d %.16g %.16g %.16g key ns L[3]\n",2, sysp.atom_count_aligned, box_edge_length, box_edge_length, box_edge_length);
	fprintf(File_buffer,"%d %d %.16g %.16g %.16g key intval vecval[3]\n",4, 0, 0.f, 0.f, 0.f);
	fprintf(File_buffer,". end of header: cfg, vel, acc follow, VarFile.size=%d\n",3664);
	fprintf(File_buffer,"! value h*velocity h^2*acceleration[3]; for leap-frog, h*velocity=value(t)-value(t-h)\n");
	fprintf(File_buffer,"%20.16f %.16g %.16g logs\n",0.0, 0.0, 0.0);
	fprintf(File_buffer,"%20.16f %20.16f %20.16f  %.16g %.16g %.16g  %.16g %.16g %.16g  lambda=ln(L)\n",log(box_edge_length),log(box_edge_length),log(box_edge_length), 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
	fprintf(File_buffer,"! r[3] h*velocity[3] h^2*acceleration[3] site\n");

	for(int i=0;i<sysp.atom_count_aligned;i++)
	{
		fprintf(File_buffer,"%20.16f %20.16f %20.16f  ",(double)position_in->x[i],(double)position_in->y[i],(double)position_in->z[i]); //positions
		fprintf(File_buffer,"%.16g %.16g %.16g  ",(double)velocity_in->x[i],(double)velocity_in->y[i],(double)velocity_in->z[i]); //velocities
		fprintf(File_buffer,"%.16g %.16g %.16g",0.f,0.f,0.f);												  							  //accelerations ignored
		fprintf(File_buffer," site[%d]\n",i); //end of line
	}

	return;
}

/* Function: save_to_cp
 * calculates configuration pressure, kinetic and potential energy of the system
 *
 *
 * Notes:
 * > The .cp file is MACSIMUS standard binary file containing measured data from system.
 *
 * Parameter:
 * molecule_postition_in - triplet* molecular center of gravity coordinates [Angstrom]
 * atom_position_in - triplet* atom coordinates in molecule [Angstrom]
 * atom_forces_in - triplet* force acting on each atom [UNKNOWN]
 * ekin_in - value_t* pointer to vector of kinetic energies of molecules
 * epot_in - value_t* pointer to vector of potential energies of molecules
 * printout_step - unsigned long int Frequency of real computer time writing into .cp file
 */
void save_to_cp(triplet *molecule_postition_in, triplet *atom_position_in,triplet *atom_forces_in,value_t* ekin_in, value_t* epot_in)
{

	//.cp file structure is described in MACSIMUS manual paragraph 16.2 and units of are in paragraph 14.2
	FILE* File_buffer;
	float into_cp[collumns_in_cp];
	char bufr[256];
	time_t curtime;
	struct tm *loc_time;

	value_t Pcfg=Return_Pcfg(molecule_postition_in,atom_position_in,atom_forces_in, ekin_in);
	value_t Ekin=Return_Ekin(ekin_in);
	value_t Epot=Return_Epot(epot_in);

	File_buffer= fopen(outp.cp_name,"ab");

	if (File_buffer != NULL)
	{																				// units were chosen using MACSIMUS standard
		into_cp[0]=(float) (Ekin + Epot);										// Total energy[K]
		into_cp[1]=(float) Ekin/sysp.molecule_count_aligned;						// Temperature [K]

		into_cp[2]=(float) Epot*kB_real*N_Avogadro/sysp.molecule_count_aligned; 	// Potential energy[J/mol]
		into_cp[3]=(float) Ekin*kB_real*N_Avogadro/sysp.molecule_count_aligned;  // Kinetic energy[J/mol]
		into_cp[4]=(float) Pcfg;										            //Pressure [Pa]
		into_cp[5]=(float) sysp.density;  											// System density[kg*m^(-3)]
		into_cp[6]=(float) 5.254;  							    //Now unknown     Translation temperature [K]
		into_cp[7]=(float) (Epot+Ekin)*kB_real*N_Avogadro/sysp.molecule_count_aligned+(Pcfg*sysp.lx*sysp.ly*sysp.lz)*1e-30; //Enthalpy [J/mol]														// Enthalpy without cutoff correction

		if((unsigned long int)ceil(itep.step/outp.cp_print_period) % outp.cp_print_timemark == 0 || itep.step*itep.dt>= itep.total_time)
		{
			curtime = time (NULL);
			loc_time = localtime (&curtime);
			strcpy(bufr,asctime (loc_time));

			fwrite(&CPmark,1,sizeof(CPmark),File_buffer);
			fwrite(&bufr,sizeof(char),(4*collumns_in_cp-4),File_buffer);

			double CPmarkW_buffer[2];
			CPmarkW_buffer[0]= (double) itep.step*itep.dt;// actual time [ps]
			CPmarkW_buffer[1]=(double) outp.cp_print_period*itep.dt;
			fwrite(&CPmarkW, sizeof(float), 1, File_buffer);
			fwrite(&CPmarkW_buffer, sizeof(float), (collumns_in_cp-1), File_buffer);

		}

		fwrite(&into_cp, 1, sizeof(into_cp), File_buffer);
		fclose(File_buffer);
	}
	else
	{
		printf(".cp file cant be opened.");
	}
	return;
}

/* Function: save_to_cfa
 * Prepares data for writing .cfa file and calls save_to_cfa_directly
 *
 * Notes:
 * > The .cfa file is MACSIMUS standard text file equivalent to .cfg binary file
 * > Uses pointers to structure triplet(contains three value_t* pointers to vectors of positions in x,y,z).
 * > Calls position_combiner
 *
 * Parameter:
 * t - double simulation time [ps]
 * molecule_postition_in - triplet* molecular center of gravity coordinates [Angstrom]
 * atom_position_in - triplet* atom coordinates in molecule [Angstrom]
 * molecule_velocity_in - triplet* molecular center of gravity velocities [UNKNOWN]
 * atom_velocity_in - triplet* atom velocities in molecule [UNKNOWN]
 * ekin_in - value_t* pointer to vector of molecule kinetic energies
 * epot_in - value_t* pointer to vector of molecule potential energies
 */
 void save_to_cfa(triplet *molecule_postition_in, triplet *atom_position_in,triplet *molecule_velocity_in, triplet *atom_velocity_in, value_t* ekin_in, value_t* epot_in)
 {
	 triplet* position_in = (triplet*) malloc(sizeof(triplet));
	 triplet* velocity_in = (triplet*) malloc(sizeof(triplet));
 
	 FILE* File_buffer;
	 char name_bufer[256];
	 char time_bufer[256];

	 strcpy(name_bufer,outp.output_files_name);

	 if(itep.step*itep.dt==0.f)
	 {
		 sprintf(time_bufer,"_t0");
	 }
	 else
	 {
		 sprintf(time_bufer,"_t%.3f",itep.step*itep.dt);
	 }

	 strcat(name_bufer,time_bufer);
	 strcat(name_bufer,".cfa");
 
	 File_buffer= fopen(name_bufer,"w");
	 fclose(File_buffer);
	 File_buffer= fopen(name_bufer,"a");
 
	 triplet_alloc(position_in,sysp.atom_count_aligned);
	 triplet_alloc(velocity_in,sysp.atom_count_aligned);
 
	 position_combiner(position_in,molecule_postition_in,atom_position_in);
	 position_combiner(velocity_in,molecule_velocity_in,atom_velocity_in);
 
	 save_to_cfa_directly(File_buffer,position_in,velocity_in,(double) Return_Etot(ekin_in, epot_in),(double) p_d_dt,(double) p_d_total_time, (double)sysp.lx);
 
	 fclose(File_buffer);
 
	 triplet_dealloc(position_in);
	 triplet_dealloc(velocity_in);
 
	 return;
 }
 
 /* Function: save_to_cfg
  * Prepares data for writing .cfg file and calls save_to_cfg_directly
  *
  * Notes:
  * > The .cfg file is MACSIMUS standard binary file
  * > Uses pointers to structure triplet(contains three value_t* pointers to vectors of positions in x,y,z).
  * > Calls position_combiner
  *
  * Parameter:
  * t - double simulation time [ps]
  * molecule_postition_in - triplet* molecular center of gravity coordinates [Angstrom]
  * atom_position_in - triplet* atom coordinates in molecule [Angstrom]
  * molecule_velocity_in - triplet* molecular center of gravity velocities [UNKNOWN]
  * atom_velocity_in - triplet* atom velocities in molecule [UNKNOWN]
  * ekin_in - value_t* pointer to vector of molecule kinetic energies
  * epot_in - value_t* pointer to vector of molecule potential energies
  */
 void save_to_cfg(triplet *molecule_postition_in, triplet *atom_position_in,triplet *molecule_velocity_in, triplet *atom_velocity_in, value_t* ekin_in, value_t* epot_in)
 {
	 triplet* position_in = (triplet*) malloc(sizeof(triplet));
	 triplet* velocity_in = (triplet*) malloc(sizeof(triplet));
 
	 FILE* File_buffer;
	 char name_bufer[256];
	 char time_bufer[256];
 
	 strcpy(name_bufer,outp.output_files_name);
 
	 if(itep.step*itep.dt==0.f)
	 {
		 sprintf(time_bufer,"_t0");
	 }
	 else
	 {
		 sprintf(time_bufer,"_t%.3f",itep.step*itep.dt);
	 }
	 strcat(name_bufer,time_bufer);
	 strcat(name_bufer,".cfg");
 
	 File_buffer= fopen(name_bufer,"w");
	 fclose(File_buffer);
 
	 triplet_alloc(position_in,sysp.atom_count_aligned);
	 triplet_alloc(velocity_in,sysp.atom_count_aligned);
 
	 position_combiner(position_in,molecule_postition_in,atom_position_in);
	 position_combiner(velocity_in,molecule_velocity_in,atom_velocity_in);
 
	 save_to_cfg_directly(name_bufer, position_in,velocity_in,(double) Return_Etot(ekin_in, epot_in),(double) itep.dt,(double) itep.total_time, (double)sysp.lx);
 
	 triplet_dealloc(position_in);
	 triplet_dealloc(velocity_in);
 
	 return;
 }
 
 void save_to_fcs(triplet* atom_forces_in)
 {
	 FILE * File_buffer;
	 File_buffer = fopen(outp.fcs_name,"a");

	 fprintf(File_buffer,"***** TIME: %f, STEP: %lu, BOX: L_x=%f  L_y=%f  L_z=%f, DENSITY: %f *****\n",itep.step*itep.dt,itep.step,sysp.lx,sysp.ly,sysp.lz,sysp.density);

	 for(unsigned int i=0;i<sysp.atom_count_aligned;i++)
	 {
#if SUBSTANCE == SPCE
		if(i%SUBSTANCE==0)
		{
			fprintf(File_buffer,"O %u \t %f \t %f \t %f\n",i,atom_forces_in->x[i],atom_forces_in->y[i],atom_forces_in->z[i]);

		}
		else
		{
			fprintf(File_buffer,"H %u \t %f \t %f \t %f\n",i,atom_forces_in->x[i],atom_forces_in->y[i],atom_forces_in->z[i]);
		}
#elif SUBSTANCE == NITROGEN
		fprintf(File_buffer,"N %u %f \t %f \t %f\n",i,atom_forces_in->x[i],atom_forces_in->y[i],atom_forces_in->z[i]);
#elif SUBSTANCE == ARGON
		fprintf(File_buffer,"Ar %u %f \t %f \t %f\n",i,atom_forces_in->x[i],atom_forces_in->y[i],atom_forces_in->z[i]);
#else
		fprintf(File_buffer,"NON %u %f \t %f \t %f\n",i,atom_forces_in->x[i],atom_forces_in->y[i],atom_forces_in->z[i]);
#endif
	 }
	 fclose(File_buffer);
 }

 /* Function: copy_data_from_GPU
  * copies data from GPU memory into appropriate triplet* structures
  *
  * Notes:
  * > Calls functions data_device2host_*
  *
  * Parameter:
  * mol_pos_tmp - triplet* molecular center of gravity coordinates [Angstrom]
  * atm_pos_tmp - triplet* atom coordinates in molecule [Angstrom]
  * mol_vel_tmp - triplet* molecular center of gravity velocities [UNKNOWN]
  * atm_vel_tmp - triplet* atom velocities in molecule [UNKNOWN]
  * atm_fcs_tmp - triplet* force acting on each atom [UNKNOWN]
  * ekin_tmp - value_t* pointer to vector of molecule kinetic energies
  * epot_tmp - value_t* pointer to vector of molecule potential energies
  * energy - bool sets if energy file should be written
  * history - bool sets if history file should be written
  * test - bool sets if test file should be written
  * test2 - bool sets if test2 file should be written
  * test3 - bool sets if test3 file should be written
  * plb - bool sets if plb file should be written
  * cp - bool sets if cp file should be written
  * cfg - bool sets if cfg file should be written
  * cfa - bool sets if cfa file should be written
  */
 void copy_data_from_GPU(triplet *mol_pos_tmp,triplet *atm_pos_tmp,triplet *mol_vel_tmp,triplet *atm_vel_tmp,triplet *atm_fcs_tmp,value_t *ekin_tmp,value_t *epot_tmp,
						 bool energy ,bool history,bool vel,bool oob,bool test4,bool fcs,bool plb,bool cp,bool cfg,bool cfa)
 {
	 if(energy || cp || cfg || cfa)
	 {
		 data_device2host_energy(ekin_tmp, epot_tmp);
	 }
	 if(history || oob || test4 || plb || cp || cfg || cfa)
	 {
		 data_device2host_position(mol_pos_tmp,atm_pos_tmp);
	 }
	 if(cfg || cfa || vel)
	 {
		 data_device2host_velocity(mol_vel_tmp,atm_vel_tmp);
	 }
	 if(cp || fcs)
	 {
		 data_device2host_force(atm_fcs_tmp);
	 }
 
	 return;
 }
 
 void print_vectors(FILE * output_file, value_t multiplycation_parameter, value_t * vector_x,value_t * vector_y,value_t * vector_z, unsigned int vector_size, unsigned int printout_period)
 {
	 if(itep.step % printout_period == 0)
	 {
		 fprintf(output_file,"#=== time: %f  step: %lu ===\n",itep.step*itep.dt,itep.step);
		 for(unsigned int i=0;i<vector_size;i++)
		 {
			 fprintf(output_file,"Ar %f %f %f\n",multiplycation_parameter*vector_x[i],multiplycation_parameter*vector_y[i],multiplycation_parameter*vector_z[i]);
		 }

	 }
	 return;
 }

  void print_into_files(triplet *mol_pos_tmp,triplet *atm_pos_tmp,triplet *mol_vel_tmp,triplet *atm_vel_tmp,triplet *atm_fcs_tmp,value_t *ekin_tmp,value_t *epot_tmp)
  {
	 bool energy_now = false;
	 bool history_now = false;
	 bool vel_now = false;
	 bool oob_now = false;
	 bool test4_now = false;
	 bool fcs_now = false;
	 bool plb_now = false;
	 bool cp_now = false;
	 bool cfg_now = false;
	 bool cfa_now = false;


	 if(((itep.step % outp.energy_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.energy_out)
	 {
		energy_now=true;
	 }
	 if(((itep.step % outp.history_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.history_out)
	 {
		history_now=true;
	 }
	 if(((itep.step % outp.vel_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.vel_out)
	 {
		vel_now=true;
     }
	 if(((itep.step % outp.oob_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.oob_out)
	 {
		oob_now=true;
	 }
	 if(((itep.step % outp.test4_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.test4_out)
	 {
		test4_now=true;
	 }
	 if(((itep.step % outp.fcs_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.fcs_out)
	 {
		fcs_now=true;
	 }
	 if(((itep.step % outp.plb_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.plb_out)
	 {
		plb_now=true;
	 }
	 if(((itep.step % outp.cp_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.cp_out)
	 {
		cp_now=true;
	 }
	 if(((itep.step % outp.cfg_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.cfg_out)
	 {
		cfg_now=true;
	 }
	 if(((itep.step % outp.cfa_print_period == 0) || itep.step == 0 || itep.step == itep.max_step) && outp.cfa_out)
	 {
		cfa_now=true;
	 }

	 copy_data_from_GPU(mol_pos_tmp, atm_pos_tmp, mol_vel_tmp, atm_vel_tmp, atm_fcs_tmp, ekin_tmp, epot_tmp, energy_now, history_now, vel_now, oob_now, test4_now, fcs_now, plb_now,cp_now, cfg_now, cfa_now);

   	 if(energy_now)
 	 {
 		 save_to_energy(ekin_tmp, epot_tmp);
 	 }
 	 if(history_now)
 	 {
 		 save_to_history(mol_pos_tmp,atm_pos_tmp);
 	 }
 	 if(vel_now)
 	 {
 		 save_to_vel(mol_vel_tmp,atm_vel_tmp);
 	 }
 	 if(oob_now)
 	 {
 		 save_to_oob(mol_pos_tmp, atm_pos_tmp);
 	 }
	 if(test4_now)
 	 {
 		 save_to_test4(mol_pos_tmp,atm_pos_tmp);
 	 }
	 if(fcs_now)
	 {
		 save_to_fcs(atm_fcs_tmp);
	 }
 	 if(plb_now)
 	 {
 		 save_to_plb(mol_pos_tmp,atm_pos_tmp);
 	 }
 	 if(cp_now)
 	 {
 		 save_to_cp(mol_pos_tmp, atm_pos_tmp,atm_fcs_tmp,ekin_tmp,epot_tmp);
 	 }
 	 if(cfg_now)
 	 {
 		 save_to_cfg(mol_pos_tmp,atm_pos_tmp,mol_vel_tmp,atm_vel_tmp,ekin_tmp, epot_tmp );
 	 }
 	 if(cfa_now)
 	 {
 		 save_to_cfa(mol_pos_tmp,atm_pos_tmp,mol_vel_tmp,atm_vel_tmp,ekin_tmp, epot_tmp );
 	 }

  }
