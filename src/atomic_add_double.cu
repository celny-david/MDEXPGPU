/*
 * atomic_add_double.cu
 *	This is copy of code from for cuda capable devices under 6:
 *  https://stackoverflow.com/questions/12626096/why-has-atomicadd-not-been-implemented-for-doubles
 *  Created on: Mar 13, 2018
 *      Author: celnyd
 */

/* Function: atomicAddOld
 * Atomic add for doubles using unsigned long long addressing
 *
 * Beware that this approach adds significant overhead compared to the float atomic add
 *
 * Necessary only for <6 cuda capable devices
 *
 * Parameters:
 *   address - the address of place where addition is performed
 *   val - value to be added to the address location
 */
__device__ double atomicAddOld(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do
    {
        assumed = old;
        old = atomicCAS(address_as_ull,
        		        assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

