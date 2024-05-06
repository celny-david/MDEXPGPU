/*
 ============================================================================
 Name        : TODO
 Author      : Dixiw
 Version     :
 Description : TODO

 ============================================================================
 */

#include <math.h>
#include <stdio.h>

#include "value.h"
#include "triplet.h"
#include "constants.h"
#include "cuda_definition.h"
#include "system_property.h"
#include "iterative_property.h"

// ===== ============== =====
// ==== helper functions ====
// ===== ============== =====

/* Function: fill_random_direction
 * procedure that generate random direction vector in unit sphere
 * method fills the given 3element array with the generated "vector"
 *
 * Notes:
 * >the allocated array of size 3 has to be given into the function
 * >use a static variable tmp_norm (reuses it for repeated generation)
 *
 * Beware:
 * > NO CONTROLL OF INPUT DIMENSION IS ENFORCED
 *
 * Parameter:
 * tmp_vect - value_t pointer to empty 3element array
 */
void fill_random_direction(value_t* tmp_vect)
{
	static value_t tmp_norm;

	tmp_vect[0] = ((value_t)rand()/(value_t)(RAND_MAX)) - 0.5;
	tmp_vect[1] = ((value_t)rand()/(value_t)(RAND_MAX)) - 0.5;
	tmp_vect[2] = ((value_t)rand()/(value_t)(RAND_MAX)) - 0.5;
	tmp_norm = sqrt( tmp_vect[0]*tmp_vect[0]
					+tmp_vect[1]*tmp_vect[1]
					+tmp_vect[2]*tmp_vect[2]);// euclidean norm
	tmp_vect[0] /= tmp_norm;
	tmp_vect[1] /= tmp_norm;
	tmp_vect[2] /= tmp_norm;
}

/* Function: fill_othogonal_direction
 * procedure that generate arbitrary orthogonal direction for given vector
 * method fills the given 3element array with the generated "vector"
 *
 * the result is normalised and in plane orthogonal to given direction
 *
 * Notes:
 * >the allocated result array has to be given into the function
 * >also requires the input nonempty 3element array as a base for orthogonal
 * >use a static variable tmp_norm (reuses it for repeated generation)
 * 
 * Beware:
 * > NO CONTROLL OF INPUT DIMENSION IS ENFORCED
 *
 * Parameter:
 * dir_vect - value_t pointer to 3element array containing direction
 * ort_vect - value_t pointer to empty 3element array (will be orthogonal)
 */
void fill_othogonal_direction(value_t* dir_vect, value_t* ort_vect)
{
	static value_t tmp_norm;

	if (dir_vect[0] !=0 )
	{// rare case of z axis vector
		ort_vect[0] = -(dir_vect[1]+dir_vect[2])/dir_vect[0];
		ort_vect[1] = 1.0;
		ort_vect[2] = 1.0;
		tmp_norm = sqrt(ort_vect[0]*ort_vect[0]+2);
	}
	else if ( dir_vect[1] !=0 )
	{
		ort_vect[0] = 1.0;
		ort_vect[1] = -(dir_vect[0]+dir_vect[2])/dir_vect[1];
		ort_vect[2] = 1.0;
		tmp_norm = sqrt(ort_vect[1]*ort_vect[1]+2);
	}
	else if ( dir_vect[2] !=0 )
	{
		ort_vect[0] = 1.0;
		ort_vect[1] = 1.0;
		ort_vect[2] = -(dir_vect[0]+dir_vect[1])/dir_vect[2];
		tmp_norm = sqrt(ort_vect[2]*ort_vect[2]+2);
	}
	ort_vect[0] /= tmp_norm;
	ort_vect[1] /= tmp_norm;
	ort_vect[2] /= tmp_norm;
	return;
}
// ===== ============ =====
// === placer functions ===
// ===== ============ =====
// placer function get the pointer to position place/ reference position place
// for the specified molecule it places it at designated position

/* Function: place_Ar
 * specific variant of placer function for Argon
 *
 * Notes:
 * >modifies the molecular and atomal triplets and adds one argon molecule
 *  at given position and index
 *
 * Parameter:
 * index 		- the index in the result array where to place the molecule
 * molecule_pos - triplet pointer of molecular(CoM) positions
 * atom_pos		- triplet pointer of atomal positions
 * where    	- the actual position where the molecule placement start from
 */
void place_Ar(unsigned int index,triplet* molecule_pos, triplet* atom_pos, value_t* where)
{
	molecule_pos->x[index] = where[0];
	molecule_pos->y[index] = where[1];
	molecule_pos->z[index] = where[2];
	atom_pos->x[index] = 0.0;
	atom_pos->y[index] = 0.0;
	atom_pos->z[index] = 0.0;
}

/* Function: place_N2
 * specific variant of placer function for Nitrogen
 *
 * Notes:
 * >modifies the molecular and atomal triplets and adds one nitrogen molecule
 *  at position calculated from given position
 * > the molecule is placed in such a way that where is the CoM(c)
 * >                     where
 * >                  N----c----N
 * >                    l2   l2
 *
 * Parameter:
 * index 		- the index in the result array where to place the molecule
 * molecule_pos - triplet pointer of molecular(CoM) positions
 * atom_pos		- triplet pointer of atomal positions
 * where    	- the actual position where the molecule placement start from
 */
void place_N2(unsigned int index,triplet* molecule_pos, triplet* atom_pos, value_t* where)
{
	// variant with the centre of mass in molecule position and differences from CoM in atom position
	static const value_t l2 = 1.08892776/2.0;
	value_t tmp_dir[3];

	fill_random_direction(tmp_dir);
	// leading nitrogen
	molecule_pos->x[index] = where[0];
	molecule_pos->y[index] = where[1];
	molecule_pos->z[index] = where[2];
	//reference position is numbered in atoms (for N2)
	atom_pos->x[2*index] = tmp_dir[0]*l2;
	atom_pos->y[2*index] = tmp_dir[1]*l2;
	atom_pos->z[2*index] = tmp_dir[2]*l2;
	//second nitrogen
	atom_pos->x[2*index+1] = - tmp_dir[0]*l2;
	atom_pos->y[2*index+1] = - tmp_dir[1]*l2;
	atom_pos->z[2*index+1] = - tmp_dir[2]*l2;
}

/* Function: place_SPCE
 * specific variant of placer function for SPCE water
 *
 * follows the geometrical constrains: H-O = 1, angle HOH = arccos(-1/3)
 *
 * Notes:
 * >modifies the molecular and atomal triplets and adds one argon molecule
 *  at given position and index
 * > the molecule is placed in such a way that where = CoM(c)
 * >                O
 * >   {          / | \
 * >   {       /    |    \
 * > l1{    /       c       \
 * >   { /       l3{|          \
 * >    H-----------v-----------H
 * >         l2           l2
 *
 * Parameter:
 * index 		- the index in the result array where to place the molecule
 * molecule_pos - triplet pointer of molecular(CoM) positions
 * atom_pos		- triplet pointer of atomal positions
 * where    	- the actual position where the molecule placement start from
 */
void place_SPCE(unsigned int index,triplet* molecule_pos, triplet* atom_pos, value_t* where)
{
	static const value_t angle = 0.5*acos(-1.0/3.0);
	static const value_t l1 = cos(angle); // the distance |vO|
	static const value_t l2 = sin(angle); // the distance |vH|
	static const value_t l3 = l1/(1.0+(2.0*MASS_H)/MASS_O); // the distance |vs| !!! wrong derivation - DC corrected
	// static const value_t l3 = l1/3.0; // geometric COM for testing purposes

	value_t rnd_dir[3]; // the second vector is {rnd_dir[1],rnd_dir[0],rnd_dir[1]}
	fill_random_direction(rnd_dir);

	value_t ort_dir[3];
	fill_othogonal_direction(rnd_dir,ort_dir);

	// the centre of mass c
	molecule_pos->x[index] = where[0];
	molecule_pos->y[index] = where[1];
	molecule_pos->z[index] = where[2];
	// the oxygen: c + (l1-l3) ortho
	atom_pos->x[3*index] = -ort_dir[0]*(l1-l3);
	atom_pos->y[3*index] = -ort_dir[1]*(l1-l3);
	atom_pos->z[3*index] = -ort_dir[2]*(l1-l3);
	//first hydrogen: c + l2 norm + l3 ortho
	atom_pos->x[3*index+1] = rnd_dir[0]*l2 + ort_dir[0]*l3;
	atom_pos->y[3*index+1] = rnd_dir[1]*l2 + ort_dir[1]*l3;
	atom_pos->z[3*index+1] = rnd_dir[2]*l2 + ort_dir[2]*l3;
	//second hydrogen: c - l2 norm + l3 ortho
	atom_pos->x[3*index+2] = -rnd_dir[0]*l2 + ort_dir[0]*l3;
	atom_pos->y[3*index+2] = -rnd_dir[1]*l2 + ort_dir[1]*l3;
	atom_pos->z[3*index+2] = -rnd_dir[2]*l2 + ort_dir[2]*l3;
}

// ===== ========== =====
// === pull functions ===
// ===== ========== =====

/* Function: fill_starting_properties
 * The unified way to obtain all properties for starting configuration
 *
 * sets up positions(mol,at), velocities (mol,at), forces(at)
 *
 * Notes:
 * >expect standard situation where number of molecules aligned
 *  and corresponding atom count aligned are used
 *
 * Parameter:
 * molecule_position - output for molecular position triplet
 * atom_position 	 - output for atomal position triplet
 * molecule_velocity - output for molecular velocity triplet
 * atom_velocity     - output for atomal velocity triplet
 * atom_force		 - output for atomal force triplet
 */
void fill_starting_properties(triplet *molecule_position, triplet *atom_position,
							  triplet *molecule_velocity, triplet *atom_velocity,
							  triplet *atom_force)
{
	unsigned int cnt, ii;
	unsigned short xx,yy,zz;
	value_t cell_edge[3], cell_edge2[3];
	value_t tmp_vec3[3]; // temporary place or direction
	value_t optimal_velocity;

	// = allocation section =
	triplet_alloc(molecule_position,sysp.molecule_count_aligned);
	triplet_alloc(atom_position,sysp.atom_count_aligned);

	triplet_alloc(molecule_velocity,sysp.molecule_count_aligned);
	triplet_alloc(atom_velocity,sysp.atom_count_aligned);

	triplet_alloc(atom_force,sysp.atom_count_aligned);

	// = grid preparation section =
	cell_edge[0] = pow(sysp.lx*sysp.ly*sysp.lz/(sysp.molecule_count_aligned),1.0/3.0); // no floor othervise when number exceeds lx*ly*lz then computation fails
	cell_edge2[0] = cell_edge[0]*0.5;
	xx = (unsigned short)floor(sysp.lx/ cell_edge[0]);
	yy = (unsigned short)floor(sysp.ly/ cell_edge[0]);
	zz = (unsigned short)floor(sysp.lz/ cell_edge[0]);
	cnt=0;

	// NOTE if cell is not suffitient due to the rounding expand it z then y then x direction untill it suffices
	while (xx*yy*zz<(sysp.molecule_count_aligned) && cnt<10)
	{
		if(cnt%3==0) zz++;
		if(cnt%3==1) yy++;
		if(cnt%3==2) xx++;
		cnt++;
	}

	// NOTE prevent the infinite expansion in case of some mishap
	if(cnt>10)
	{
		fprintf(stderr,"ERROR: expansion was not successful check the filling algorithm");
		exit(-1);
	}

	cell_edge[0] = sysp.lx/xx;
	cell_edge[1] = sysp.ly/yy;
	cell_edge[2] = sysp.lz/zz;

	cell_edge2[0] = cell_edge[0]*0.5;
	cell_edge2[1] = cell_edge[1]*0.5;
	cell_edge2[2] = cell_edge[2]*0.5;

	// NOTE place the molecule properties into the holding arrays
	for (cnt=0; cnt < sysp.molecule_count_aligned; cnt++)
	{
		// MOLECULE POSITION
		tmp_vec3[0] = cell_edge2[0] + cell_edge[0]*(cnt%xx);
		tmp_vec3[1] = cell_edge2[1] + cell_edge[1]*((cnt/xx)%yy);
		tmp_vec3[2] = cell_edge2[2] + cell_edge[2]*(cnt/(xx*yy));
#if SUBSTANCE == ARGON
		place_Ar(cnt,molecule_position,atom_position,tmp_vec3);
#elif SUBSTANCE == NITROGEN
		place_N2(cnt,molecule_position,atom_position,tmp_vec3);
#elif SUBSTANCE == SPCE
		place_SPCE(cnt,molecule_position,atom_position,tmp_vec3);
#endif
		// MOLECULE VELOCITY
		// NOTE the *itep.dt is mandatory for kolafa algorithm expects velocity dimension of dt*v
		// NOTE inverse of total molecule mass is the last element of mol_mass_template
		optimal_velocity = sqrt(2*kB*sysp.temperature*sysp.mol_mass_template[SUBSTANCE])*itep.dt;

		//printf("Optimal vel: %15f\n",optimal_velocity/itep.dt); // DEBUG write the assigned absolute value of velocity of molecule

		fill_random_direction(tmp_vec3);

		molecule_velocity->x[cnt] = tmp_vec3[0]*optimal_velocity;
		molecule_velocity->y[cnt] = tmp_vec3[1]*optimal_velocity;
		molecule_velocity->z[cnt] = tmp_vec3[2]*optimal_velocity;

		// NOTE set atoms in molecule as non-moving -> this will be smoothed out in first iteration of force calculation
		for (ii = cnt*SUBSTANCE; ii < (cnt+1)*SUBSTANCE; ii++)
		{
			// ATOM VELOCITY
			atom_velocity->x[ii] = 0.0;
			atom_velocity->y[ii] = 0.0;
			atom_velocity->z[ii] = 0.0;
			// ATOM FORCE
			atom_force->x[ii] = 0.0;
			atom_force->y[ii] = 0.0;
			atom_force->z[ii] = 0.0;
		}
	}
	return;

}
