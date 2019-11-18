/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_MYSPHERICAL__
#define __INCLUDE_MYSPHERICAL__

#include "myUtils.h"

INT SpinWeightedSphericalHarmonic(REAL8 theta,
                                  REAL8 phi,
                                  INT s,
                                  INT l,
                                  INT m,
                                  COMPLEX16 *ret);


#endif

