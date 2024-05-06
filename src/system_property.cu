/*
 * system_property.cu
 *
 *  Created on: Jan 18, 2020
 *      Author: klimam
 */

#include "system_property.h"
#include "constants.h"

// ===== ===== ===== =====
// ===== SYSTEM PROP =====
// ===== ===== ===== ===== 

// === set functions system ===

/* Function: set_molecule_count
 * Set up the molecule count overloaded function
 *
 * Notes:
 * >the count is set to be aligned to the align factor
 *
 * Parameter:
 * molecule_count_in - the number of molecules
 * align_factor		 - align factor of molecules (usually power of 2 number i.e. 1024)
 */
void set_molecule_count(unsigned int molecule_count_in, unsigned int align_factor)
{
	// beware this works fine only with limited number of substances as argon is uniatomal, nitrogen diatomic, spce triatomic and tip4p quatroatomic
	set_molecule_count(molecule_count_in,SUBSTANCE, align_factor);
	return;
}

/* Function: set_molecule_count
 * Set up the molecule count the original function
 *
 * Notes:
 * >the count is set to be aligned to the align factor
 * > takes the atom per molecule into account
 * > !!! checking criterion
 * >   atom_count%align_factor==0 && atom_count%sysp.molecule_count==0
 *
 * Parameter:
 * molecule_count_in 	- the number of molecules [#]
 * atom_per_molecule_in - the number of atoms per molecule (i.e. 2 for N2)
 * align_factor		 	- align factor of molecules (usually power of 2 number i.e. 1024)
 */
void set_molecule_count(unsigned int molecule_count_in, unsigned int atom_per_molecule_in, unsigned int align_factor)
{
	sysp.molecule_count = molecule_count_in;
	sysp.atom_per_molecule = atom_per_molecule_in;
	sysp.atom_count = sysp.molecule_count*sysp.atom_per_molecule;
	if (sysp.atom_count%align_factor==0 && sysp.atom_count%sysp.molecule_count==0) // the input already aligned
	{
		sysp.molecule_count_aligned = sysp.molecule_count;
		sysp.atom_count_aligned = sysp.molecule_count*sysp.atom_per_molecule;
	}else // else add what is left to be aligned
	{
		sysp.atom_count_aligned = (sysp.atom_count/(align_factor*sysp.atom_per_molecule) +1 )*align_factor*sysp.atom_per_molecule;
		sysp.molecule_count_aligned = sysp.atom_count_aligned/sysp.atom_per_molecule;
	}
	return;
}

/* Function: set_substance_name
 * saves the substance name into System_prop structure sysp, exactly into char[32] sysp.name
 *
 * Notes:
 * > Enabled names are:
 * > argon, spce, nitrogen (wave 1)
 */
void set_substance_name()
{

#if	SUBSTANCE == ARGON

	strcpy(sysp.name,"argon");

#elif SUBSTANCE == SPCE

	strcpy(sysp.name,"spce");

#elif SUBSTANCE == NITROGEN

	strcpy(sysp.name,"nitrogen");

#endif

	return;
}

/* Function: set_cutoff
 * set method for the standard cutoff and verlet cutoff only
 *
 * electrostatic cutoff is set to be unused
 *
 * Notes:
 * > the cutoff_verlet should be larger then cutoff for any sensible results
 * >!no input correctness check is performed
 *
 * Parameter:
 * cutoff_in		- standard interaction cutoff distance [angstrom]
 * cutoff_verlet_in	- standard interaction cutoff distance for verlet list construction [angstrom]
 */
void set_cutoff(value_t cutoff_in, value_t cutoff_verlet_in)
{
	sysp.cutoff = cutoff_in;
	sysp.cutoff_verlet = cutoff_verlet_in;
	sysp.cutoff_elstat = -1;
	return;
}

/* Function: set_cutoff
 * set method for the standard cutoff and verlet cutoff and electrostatic cutoff
 *
 * Notes:
 * > the cutoff_verlet should be larger then both cutoffs for any sensible results
 * >!no input correctness check is performed
 *
 * Parameter:
 * cutoff_in		- standard interaction cutoff distance [angstrom]
 * cutoff_elstat_in	- electrostatic interaction cutoff distance [angstrom]
 * cutoff_verlet_in	- standard interaction cutoff distance for verlet list construction [angstrom]
 */
void set_cutoff(value_t cutoff_in, value_t cutoff_elstat_in, value_t cutoff_verlet_in)
{
	sysp.cutoff = cutoff_in;
	sysp.cutoff_elstat = cutoff_elstat_in;
	sysp.cutoff_verlet = cutoff_verlet_in;
	return;
}

/* Function: set_temperature
 * method for setting the temperature
 *
 * Notes:
 * >!no input correctness check is performed
 *
 * Parameter:
 * temperature_in - the input temperature in [Kelvin]
 */
void set_temperature(value_t temperature_in)
{
	sysp.temperature = temperature_in;
	return;
}

/* Function: set_mol_mass
 * method for setting the mass template for molecule
 *
 * Notes:
 * >the presently implemented models are hardcoded so no input parameter is required
 * >the mass template is for single molecule and -> homogeneous system in composition is expected
 * >the last element of the mass template is total mass of molecule
 */
void set_mol_mass()
{
	// = mass init section =
#if SUBSTANCE == ARGON
	sysp.mol_mass_template[0] = MASS_AR;
	sysp.mol_mass_template[1] = 1.0/MASS_AR;
#elif SUBSTANCE == NITROGEN
	sysp.mol_mass_template[0] = MASS_N;
	sysp.mol_mass_template[1] = MASS_N;
	sysp.mol_mass_template[2] = 1.0/(2*MASS_N);
#elif SUBSTANCE == SPCE
	sysp.mol_mass_template[0] = MASS_O;
	sysp.mol_mass_template[1] = MASS_H;
	sysp.mol_mass_template[2] = MASS_H;
	sysp.mol_mass_template[3] = 1.0/(MASS_O + 2*MASS_H);
#endif
}

void set_system_real_mass()
{
#if SUBSTANCE == ARGON
	sysp.system_real_mass=0.001*sysp.molecule_count_aligned*(REAL_MASS_AR)/N_Avogadro;
#elif SUBSTANCE == NITROGEN
	sysp.system_real_mass=0.001*sysp.molecule_count_aligned*2*(REAL_MASS_N)/N_Avogadro;
#elif SUBSTANCE == SPCE
	sysp.system_real_mass=0.001*sysp.molecule_count_aligned*(2*(REAL_MASS_H)+REAL_MASS_O)/N_Avogadro;
#endif

	return;
}

/* Function: set_l
 * the method for setting the dimension of the calculation volume
 *
 * Notes:
* >!no input correctness check is performed
 *
 * Parameter:
 * lx_in - the x size of the box
 * ly_in - the y size of the box
 * lz_in - the z size of the box
 */
void set_l(value_t lx_in, value_t ly_in, value_t lz_in)
{
	sysp.lx = (value_t) lx_in;
	sysp.ly = (value_t) ly_in;
	sysp.lz = (value_t) lz_in;
	return;
}

/* Function: set_l_from_density
 * calculates the box dimensions for given density of simulated system
 *
 *
 * Notes:
 * > Cubic simulation box is created
 * > Number of molecules in simulation is used
 * > Real molar mass is used
 *
 * Parameter:
 * p_density - float required density[kg/m^3]
 */
void set_l_from_density(float density_in)
{//density unit is kg*m^(-3)

	float l_from_density;
	// NOTE get the inputed density to the system property for later printout
	sysp.density = density_in;

#if	SUBSTANCE == ARGON

	l_from_density=cbrt(1e27*sysp.molecule_count_aligned*(REAL_MASS_AR)/(density_in*N_Avogadro)); // MASS_AR and no REAL_MASS_AR is correct here

#elif SUBSTANCE == SPCE

    l_from_density=cbrt(1e27*sysp.molecule_count_aligned*(REAL_MASS_O+2*REAL_MASS_H)/(density_in*N_Avogadro));

#elif SUBSTANCE == NITROGEN

    l_from_density=cbrt(1e27*sysp.molecule_count_aligned*(2*REAL_MASS_N)/(density_in*N_Avogadro));

#endif

	set_l(l_from_density,l_from_density,l_from_density);
	return;
}

void set_density_from_l(value_t l)//l [Angstrom]
{
	sysp.density=sysp.system_real_mass/(l*l*l*1e-30);
	return;
}

void set_density(void)
{
#if	SUBSTANCE == ARGON

	sysp.density =1e27*sysp.molecule_count_aligned*(REAL_MASS_AR)/(sysp.lx*sysp.ly*sysp.lz*N_Avogadro);

#elif SUBSTANCE == SPCE

	sysp.density =1e27*sysp.molecule_count_aligned*(REAL_MASS_O+2*REAL_MASS_H)/(sysp.lx*sysp.ly*sysp.lz*N_Avogadro);

#elif SUBSTANCE == NITROGEN

    sysp.density =1e27*sysp.molecule_count_aligned*(2*REAL_MASS_N)/(sysp.lx*sysp.ly*sysp.lz*N_Avogadro);

#endif


	return;
}
