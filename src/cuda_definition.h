/*
 * cuda_definition.h
 *
 *  Created on: 14 Jan 2020
 *      Author: dixiw
 */

#ifndef CUDA_DEFINITION_H_
#define CUDA_DEFINITION_H_

//====== DEBUG PURPOSES =======
//#define TESTING // uncomment if needed

//========= PRECISION =========
/* Defines: PREC
 *
 * FLOAT  = 1 - the sigle precision
 * DOUBLE = 2 - the double precision
 *
 * PREC switches the precision of all important calculated variables
 *
 * PREC also changes precision of the derived types as float3/double3
 */
#define FLOAT 1
#define DOUBLE 2

#define PREC FLOAT //markPrec
//#define EXPANSION //markExp


//========= SUBSTANCE =========
/* Defines: Substance
 *
 * SUBSTANCE = LJF 		= 0 	- Lennard Jones Fluid (sigma =1, epsilon =1)
 * SUBSTANCE = ARGON 	= 1 	- ARGON scaled LJF to (sigma =3.4 AA, epsilon =120K)
 * SUBSTANCE = NITROGEN	= 2 	- NITROGEN selfstanding (sigma =3,3078 AA, epsilon =36.6727K)
 * SUBSTANCE = SPCE		= 3 	- SPCE selfstanding (sigma =3,1655579 AA, epsilon ~=78.1974278134876 K)
 * SUBSTANCE = TIP4P	= 4 	- TIP4P !!! NOT IMPLEMENTED
 *
 * SUBSTANCE has to be changed manually in the source code 
 * TODO template this definition or generate it with build system
 */
#define LJF 0
#define ARGON 1
#define NITROGEN 2
#define SPCE 3
#define TIP4P 4 //BEWARE is not implemented

#define SUBSTANCE SPCE

//======== CUDA OPTION ========
/* Defines: CUDA option
 *
 * WORKLOAD 					- If workloading type of algorithm is implemented this mean the 3 of element thread evaluates
 * THREAD_P_BLOCK 				- determine # of threads in one block in forcefield related context - Beware the combination with WORKLOAD
 * THREAD_P_BLOCK_RECONSTRUCT 	- determine # of threads in one block in verlet related context
 *
 * CUDA specific parameter for kernel launch grid difvision
 * TODO calculate these dynamically depending on the GPU in use
 */
#define WORKLOAD 64 // <- set this for the number of element computed by thread change
					// influence occupancy of the GPU - too high result in computation bound case
					// too low result in the memory bound case
#define THREAD_P_BLOCK 16 //number of thread per block
					      // quite sufficient for most cases
#define THREAD_P_BLOCK_RECONSTRUCT 64 // <- set this for the number of thread per block during verlet reconstruct
									  // should be generally high as the algorithm is optimized for bigger threadblocks

//======= POTENTIAL TYPE ======
/* Defines: Potential types
 *
 * POT_TYPE = 6 	- LJF (sigma=1, epsilon=1) potential, For argon it has to be rescaled, type rr-3-32
 * POT_TYPE = 10 	- Nitrogen LJ rescaled potential, type rr-3-2.924633784943219 (grid/sigma^2)
 * POT_TYPE = 20 	- SPCE LJ O-O rescaled potential, type r-3-10.10880265341735 (grid/sigma)
 *
 * POT_TYPE is automatically swithced by the choice of SUBSTANCE
 */
#if SUBSTANCE==ARGON
	#define POT_TYPE 6 // 1-6
#elif SUBSTANCE== NITROGEN
	#define POT_TYPE 10
#elif SUBSTANCE== SPCE
	#define POT_TYPE 20
#endif // no other pot yet available

//======= ELECTROSTATIC POTENTIAL TYPE ======
/* Defines: Electrostatic potential types
 *
 * ELPOT_TYPE = 0 	- empty electrostatic potential - dummy that does not contain any interafction
 * ELPOT_TYPE = 20 	- SPCE Culomb H-H real cutoff potential, type qq4r-12.6875-32 - needs to be scaled to used units by K_E_INV
 *
 * ELPOT_TYPE is automatically swithced by the choice of SUBSTANCE
 */
#if SUBSTANCE==ARGON
	#define ELPOT_TYPE 0
#elif SUBSTANCE== NITROGEN
	#define ELPOT_TYPE 0
#elif SUBSTANCE== SPCE
	#define ELPOT_TYPE 20
#endif

//======= POTENTIAL FORCE ======
/* Defines: Force types
 *
 * FOR_TYPE = 0 - Potential require r without division valid only for {POT_TYPE=-1}
 * FOR_TYPE = 1 - Potential require r and division 1/r is performed
 * FOR_TYPE = 2 - Potential require r^2 without division
 *
 * FOR_TYPE switches the force calculation and should correspond with POT_TYPE
 *
 * the switching is done automatically according to the POT_TYPE
 */
#if POT_TYPE ==0
	#define FOR_TYPE 0
#elif POT_TYPE == 1 || POT_TYPE == 2 || POT_TYPE == 4 || POT_TYPE == 5 || POT_TYPE == 20
	#define FOR_TYPE 1
#elif POT_TYPE == -1 || POT_TYPE == 3 || POT_TYPE == 6 || POT_TYPE == 10
	#define FOR_TYPE 2
#endif

#if ELPOT_TYPE == 20
	#define ELFOR_TYPE 1
#endif


//==== SUBSTANCE PARAMETERS ====
/* Defines: Force types
 * 
 * MASS 		- define the mass of the substance in calculation macsimus units [molar_mass*Na/kB*1e4 = molar_mass*10/R]
 * REAL_MASS 	- define the real molar mass of the SUBSTANCE [g/mol]
 * K_E 			- Culoumb constant k_e [eV·Å·e−2]
 * K_E_INV 		- inverted Culoumb constant 1/k_e [1/(eV·Å·e−2)]
 */

#define MASS_LJF 1.0
#define MASS_AR 48.0464
#define MASS_N 16.8461937
#define MASS_O 19.2428617
#define MASS_H 1.21227359

#define REAL_MASS_LJF 1.0
#define REAL_MASS_AR 39.948
#define REAL_MASS_N 14.0067
#define REAL_MASS_O 15.9994
#define REAL_MASS_H 1.00794

#define K_E 14.3996 // culoumb constant [eV·Å·e−2]
#define K_E_INV 0.069446374 // inverse of culoumb constant [eV·Å·e−2]

#endif /* CUDA_DEFINITION_H_ */
