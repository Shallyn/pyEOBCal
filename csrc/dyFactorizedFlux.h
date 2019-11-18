/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYFACTORIZEDFLUX__
#define __INCLUDE_DYFACTORIZEDFLUX__

#include "dyUtils.h"

REAL8
XLALInspiralSpinFactorizedFlux (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
				const REAL8 omega,	/**< orbital frequency */
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				INT const  lMax);


#endif

