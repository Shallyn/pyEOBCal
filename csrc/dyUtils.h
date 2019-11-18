/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYUTILS__
#define __INCLUDE_DYUTILS__

#include "myUtils.h"
#include "dyDatatypes.h"
#include "dyConstants.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_odeiv.h>

void CalculateSigmaStar(REAL8Vector *sigmaStar,
                        REAL8 mass1,
                        REAL8 mass2,
                        REAL8Vector *s1,
                        REAL8Vector *s2);

void CalculateSigmaKerr(REAL8Vector *sigmaKerr,
                        REAL8 mass1,
                        REAL8 mass2,
                        REAL8Vector *s1,
                        REAL8Vector *s2);

void DestroySpinEOBDynamics(SpinEOBDynamics *dyEOB);
SpinEOBDynamics *SpinEOBDynamicsInit();


#endif

