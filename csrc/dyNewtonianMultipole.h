/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYNEWTONIANMULTIPOLE__
#define __INCLUDE_DYNEWTONIANMULTIPOLE__

#include "dyUtils.h"

INT ComputeNewtonMultipolePrefixes(NewtonMultipolePrefixes *prefix,
                                   const REAL8 m1 ,
                                   const REAL8 m2);

int
XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 REAL8 r,       /**<< Orbital separation (units of total mass M */
                 REAL8 phi,     /**<< Orbital phase (in radians) */
                 UINT  l,             /**<< Mode l */
                 INT  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );

int
XLALSimIMRSpinEOBCalculateNewtonianMultipole(
                COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                REAL8 r,       /**<< Orbital separation (units of total mass M */
                REAL8 phi,            /**<< Orbital phase (in radians) */
                UINT  l,             /**<< Mode l */
                INT  m,              /**<< Mode m */
                EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
);


#endif

