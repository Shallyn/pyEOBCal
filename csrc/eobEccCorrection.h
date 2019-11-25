/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_EOBECCCORRECTION__
#define __INCLUDE_EOBECCCORRECTION__

#include "dyUtils.h"
#include "dyHamiltonian.h"

INT EOBEccFactorizedWaveformCorrection(COMPLEX16 *hECC,
                                       const REAL8 rphivalues[],
                                       const REAL8 rdot,
                                       const REAL8 phidot,
                                       const REAL8 eta);

INT EOBCalculateEccCorrection(COMPLEX16 *hECC,
                              const REAL8Vector *values,
                              SpinEOBParams *seobParams);


#endif

