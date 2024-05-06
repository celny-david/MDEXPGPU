/*
 * err_check.cu
 *
 *  Created on: 25 Oct 2017
 *      Author: dixiw
 */

/* Function: cudaSafeCall
 * enables to check success of performed cuda operation
 *
 * Example:
 *   cudaSafeCall( cudaFree(variable) );
 *
 * Parameter:
 *    err - the cudaError builtin type
 */

/* Function: cudaSafeKernel
 * enables to check success of called kernel
 *
 * Example:
 * 	 yourkernell<<<1,1>>>(variable)
 *   cudaSafeKernel();
 */

#include <stdio.h>
#include "system_property.h"


void cutoff_box_check()
{
	if(sysp.cutoff >= sysp.lx/2 || sysp.cutoff >= sysp.ly/2 || sysp.cutoff >= sysp.lz/2 )
	{
		printf("!!! WARNING : The cutoff(%f) is higher then half of a box sizes(x/2: %f, y/2: %f, z/2: %f). !!!\n",sysp.cutoff,sysp.lx/2,sysp.ly/2,sysp.lz/2);
	}
}

void verletCutoff_box_check()
{
	if(sysp.cutoff_verlet >= sysp.lx/2 || sysp.cutoff_verlet >= sysp.ly/2 || sysp.cutoff_verlet >= sysp.lz/2 )
	{
		printf("!!! WARNING : The cutoff_verlet(%f) is higher then half of a box sizes(x/2: %f, y/2: %f, z/2: %f). !!!\n",sysp.cutoff_verlet,sysp.lx/2,sysp.ly/2,sysp.lz/2);
	}
}

void elstatCutoff_box_check()
{
	if(sysp.cutoff_elstat >= sysp.lx/2 || sysp.cutoff_elstat >= sysp.ly/2 || sysp.cutoff_elstat >= sysp.lz/2 )
	{
		printf("!!! WARNING : The cutoff_elstat(%f) is higher then half of a box sizes(x/2: %f, y/2: %f, z/2: %f). !!!\n",sysp.cutoff_elstat,sysp.lx/2,sysp.ly/2,sysp.lz/2);

	}
}

void check_settings()
{
	cutoff_box_check();
	verletCutoff_box_check();
	elstatCutoff_box_check();
}
