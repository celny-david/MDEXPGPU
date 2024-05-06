/*
 * system_property.h
 *
 *  Created on: Jan 18, 2020
 *      Author: klimam
 */

#ifndef SYSTEM_PROPERTY_H_
#define SYSTEM_PROPERTY_H_

#include "value.h"


// == ========================== ==
// == system relevant parameters ==


/* Structures: System_prop
 * enveloping of the basic system properties
 *
 * useful mainly for setup on CPU
 *
 * atom_per_molecule	  - the number of atoms per molecule i.e. 3 for SCPE water model
 * atom_count			  - the requested number of atoms
 * atom_count_aligned	  - the used number of atoms aligned by to safely fit to GPU (power of 2)
 * molecule_count		  - the requested number of molecules
 * molecule_count_aligned - the used number of molecules aligned by to safely fit to GPU (power of 2)
 * cutoff				  - the cutoff distance for standard interaction
 * cutoff_verlet		  - the cutoff distance for verlet construction
 * temperature			  - requested temperature
 * mol_mass_template	  - the [atom_per_molecule+1] array with the atom m. masses the final element the m. mass of molecule
 * lx					  - the x dimension of computation box
 * ly					  - the y dimension of computation box
 * lz					  - the z dimension of computation box
 * name					  - the name of simulated molecule
 */
struct System_prop
{
	unsigned int atom_per_molecule;
	unsigned int atom_count;
	unsigned int atom_count_aligned;
	unsigned int molecule_count;
	unsigned int molecule_count_aligned;
	char name[32];

	value_t cutoff;
	value_t cutoff_elstat;
	value_t cutoff_verlet;
	value_t temperature;
	value_t density;
	value_t system_real_mass; //[kg]
	value_t mol_mass_template[SUBSTANCE+1];

	value_t lx;
	value_t ly;
	value_t lz;
};
extern System_prop sysp;

// === set functions simulation ===
void set_molecule_count(unsigned int molecule_count_in,
		                unsigned int atom_per_molecule_in);
void set_molecule_count(unsigned int molecule_count_in,
		                unsigned int atom_per_molecule_in,
		                unsigned int align_factor);
void set_substance_name();
void set_cutoff(value_t cutoff_in,value_t cutoff_verlet_in);
void set_cutoff(value_t cutoff_in,value_t cutoff_elstat_in,value_t cutoff_verlet_in);
void set_temperature(value_t temperature_in);
void set_mol_mass();
void set_system_real_mass();
void set_l(value_t lx_in, value_t ly_in, value_t lz_in);
void set_l_from_density(float p_density);
void set_density_from_l(value_t l);
void set_density(void);


#endif /* SYSTEM_PROPERTY_H_ */
