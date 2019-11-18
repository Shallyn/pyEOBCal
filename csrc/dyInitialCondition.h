/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYINITIALCONDITION__
#define __INCLUDE_DYINITIALCONDITION__

#include "dyUtils.h"
#include "dyHcapNumericalDerivative.h"

INT CalcInitialConditions(REAL8Vector *initConds,
                          const REAL8 mass1,
                          const REAL8 mass2,
                          const REAL8 fMin,
                          SpinEOBParams *params);


#endif

