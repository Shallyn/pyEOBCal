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
				INT const  lMax,
				INT allow_ecc);

REAL8 InspiralSpinFactorizedFlux_elip(REAL8Vector *values,
                                      const REAL8 rphivalues[],
                                      const REAL8 drphivalues[],
                                      EOBNonQCCoeffs *nqcCoeffs,
                                      const REAL8 omega,
                                      SpinEOBParams *ak,
                                      const REAL8 H,
                                      const UINT lMax);

INT InspiralSpinFactorizedFlux_New(REAL8Vector *values,
                                   EOBNonQCCoeffs *nqcCoeffs,
                                   const REAL8 omega,
                                   SpinEOBParams *ak,
                                   const REAL8 H,
                                   const UINT lMax,
                                   REAL8 *Eflux,
                                   REAL8 *Lflux);


#endif

