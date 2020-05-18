/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYHAMILTONIAN__
#define __INCLUDE_DYHAMILTONIAN__

#include "dyUtils.h"

REAL8 SpinEOBHamiltonian(const REAL8 eta,
                         REAL8Vector *x,
                         REAL8Vector *p,
                         REAL8Vector *s1Vec,
                         REAL8Vector *s2Vec,
                         REAL8Vector *sigmaKerr,
                         REAL8Vector *sigmaStar,
                         int tortoise,
                         SpinEOBHCoeffs * coeffs);

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff (const REAL8 values[],
					SpinEOBParams * funcParams
  );

REAL8
XLALSimIMRSpinAlignedEOBCalcOmega (const REAL8 values[],/**<< Dynamical variables */
				   SpinEOBParams * funcParams,
							/**<< EOB parameters */
                    REAL8 STEP_SIZE /**<< Step size for numerical derivation of H */
  );

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

double
GSLSpinAlignedHamiltonianWrapper (double x, void *params);

REAL8 auxCalculateCircularAngularMomentum(const REAL8 eta,
                         REAL8Vector *x,
                         REAL8Vector *s1Vec,
                         REAL8Vector *s2Vec,
                         REAL8Vector *sigmaKerr,
                         REAL8Vector *sigmaStar,
                         int tortoise,
                         SpinEOBHCoeffs * coeffs);

REAL8 PNCalcOrbitOmega(const REAL8 Hreal,
                       const REAL8 ecc,
                       const REAL8 eta);

#endif

