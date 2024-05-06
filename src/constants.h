/*
 * constants.h
 *
 *  Created on: Dec 18, 2019
 *      Author: klimam
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// BUG how come this does not require #include "value.h"

// == HOST global/constant parameters ==
/* Constants: general Host constants

	m_pi 		- pi constant
	kB 			- reduced Boltzmann constant [ macsimus units]
	kB_real 	- Boltzmann constant [J/K]
	N_Avogadro 	- Avogadro constant [#]
	J_in_kcal 	- conversion constant for thermodynamic kcal to Joule [J]
*/
// == general constants section ==
static const value_t m_pi = 3.14159265358979323846;
static const value_t kB = 1; // kB=1 because of the reduced macsimus units used for mass
static const value_t kB_real = 1.38064852e-23;
static const value_t N_Avogadro = 6.02214076e23;
static const value_t J_in_kcal = 4184;

/* Constants: Output related parameters

	collumns_in_cp 	- number of columns in the cp files [#]
	NCP 			- TODO MK fill description []
	CPmark 			- TODO MK fill description []
	CPmarkU 		- TODO MK fill description []
	CPmarkV 		- TODO MK fill description []
	CPmarkW 		- TODO MK fill description []
*/
//static const unsigned int collumns_in_cp=8;
static const int collumns_in_cp=8;
static const int NCP=4;
static const float CPmark=((float)-1.755645e+37);
static const float CPmarkT=((float)-1.76398512e+37);
static const float CPmarkU=((float)-1.77232525e+37);
static const float CPmarkV=((float)-1.78066538e+37);
static const float CPmarkW=((float)-1.7890055e+37);

#endif /* CONSTANTS_H_ */
