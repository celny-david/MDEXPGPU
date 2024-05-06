/*
 ============================================================================
 Name        : property.cu
 Author      : Dixiw
 Version     :
 Description : Collects all parameters relevant for simulation

 NOTES : * contain methods for setting global parameters
 	 	 * for additional safety the structs can become classes with
 	 	 only constructor setting the private parameters and public getters
 ============================================================================
 */

#include <string.h> // strcat

#include "property.h"
#include "output.h"
#include "constants.h"
#include "cuda_definition.h"
#include "system_property.h"
#include "iterative_property.h"


void create_custom_file(FILE* the_file,const char* name,const char* suffix)
{
	char name_buffer[256];

	if(strlen(name)>=(256-4) )
	{
		fprintf (stderr,"Error name length exceeded allowed limit\n");
		return;
	}
	if(strlen(suffix)>=4 )
	{
		fprintf (stderr,"Error suffix length exceeded allowed limit\n");
		return;
	}

	strcpy (name_buffer,name);

	strcat (name_buffer,".");
	strcat (name_buffer,suffix);
	the_file = fopen(name_buffer,"w"); // for clearing content of existing file
	fclose(the_file);
	the_file = fopen(name_buffer,"a");
}
/* Function: Return_Etot
 * returns value_t total energy of system [K*Molecules]
 *
 * Notes:
 * > Run serial on CPU
 *
 * Parameter:
 * ekin_in - value_t*  pointer to vector of kinetic energies of molecules
 * epot_in - value_t*  pointer to vector of potential energies of molecules
 */
value_t Return_Etot(value_t* ekin_in, value_t* epot_in)
{
	value_t ekin_tmp =0.0;
	value_t epot_tmp =0.0;

	// variant with addition on CPU
	for(unsigned int i=0; i<sysp.molecule_count_aligned; i++ ) // NOTE assume the input arrays are of length sysp.atom_count
	{
		ekin_tmp += ekin_in[i];
		epot_tmp += epot_in[i];
	}
	return ((ekin_tmp/(4.0*itep.dt*itep.dt)) + (epot_tmp/2.0));
	// from 14.08.2020 the dt^2 division is performed here
}

value_t Return_Epot(value_t* epot_in)
{
	value_t epot_tmp =0.0;

	// variant with addition on CPU
	for(unsigned int i=0; i<sysp.molecule_count_aligned; i++ ) // NOTE assume the input arrays are of length sysp.atom_count
	{
		epot_tmp += epot_in[i];
	}
	return (epot_tmp/2.0);
}

value_t Return_Ekin(value_t* ekin_in)
{
	value_t ekin_tmp =0.0;

	// variant with addition on CPU
	for(unsigned int i=0; i<sysp.molecule_count_aligned; i++ ) // NOTE assume the input arrays are of length sysp.atom_count
	{
		ekin_tmp += ekin_in[i];
	}
	return (ekin_tmp/(4.0*itep.dt*itep.dt));

}

/* Function: Return_Pcfg
 * returns the configuration pressure in Pascals
 * Pressure is calculated using following formula:
 * 		Pcfg=1/(3*V)*(2*Ekin+sum(r_i*f_i))
 * 		V=simulation box volume
 * 		Ekin=total kinetic energy of system
 * 		r_i=position vector of atom i
 * 		f_i=vector of forces acting to atom i
 *
 * Notes:
 * > internal using of sysp.lx, sysp.ly, and sysp.lz all in [Angstrom]
 * > internal using of number of molecules and atoms in simulation
 * > internal calling of position_combiner()
 *
 * Parameter:
 * molecule_postition_in - triplet* molecular center of gravity coordinates [Angstrom]
 * atom_position_in - triplet* atom coordinates in molecule [Angstrom]
 * atom_forces_in - triplet* force acting on each atom [UNKNOWN]
 * ekin_in - triplet*  kinetic energy of each atom
 */
value_t Return_Pcfg(triplet *molecule_postition_in, triplet *atom_position_in,triplet *atom_forces_in, value_t* ekin_in)
{
	value_t Pressure=0.0, E_kin=0.0;

	// obtaining atom positions in position_in
	triplet* position_in = (triplet*) malloc(sizeof(triplet));
	triplet_alloc(position_in,sysp.atom_count_aligned);
	position_combiner(position_in,molecule_postition_in,atom_position_in);

	// obtaining total kinetic energy
	for(unsigned int i=0; i<sysp.molecule_count_aligned; i++ )
	{
		E_kin += ekin_in[i];
	}
	E_kin=E_kin/(4.0*itep.dt*itep.dt);
	// from 14.08.2020 the dt^2 division is performed here

	// calculation of sum
	for(unsigned int i=0; i<sysp.atom_count_aligned; i++ )
	{
		Pressure +=position_in->x[i]*atom_forces_in->x[i];
		Pressure +=position_in->y[i]*atom_forces_in->y[i];
		Pressure +=position_in->z[i]*atom_forces_in->z[i];
	}
	triplet_dealloc(position_in);

	Pressure+=2*E_kin;
	Pressure=Pressure/(3*sysp.lx*sysp.ly*sysp.lz);
	Pressure=Pressure*1e30*kB_real/sysp.molecule_count_aligned;//recalculation to Pascals
	return Pressure;
}

value_t Return_Temperature(triplet *atom_velocity_in)
{
	value_t ekin_tmp=0.0f;
	value_t result;
	unsigned int u;

	for(unsigned int i=0;i<sysp.atom_count_aligned;i++)
	{
		u=i%SUBSTANCE;
		ekin_tmp += (atom_velocity_in->x[i]*atom_velocity_in->x[i])+
					(atom_velocity_in->y[i]*atom_velocity_in->y[i])+
					(atom_velocity_in->z[i]*atom_velocity_in->z[i])*sysp.mol_mass_template[u];
	}

	result=0.5*ekin_tmp/(itep.dt*itep.dt*sysp.molecule_count_aligned);
	return result; //WARNING possible error in unit conversion

}
