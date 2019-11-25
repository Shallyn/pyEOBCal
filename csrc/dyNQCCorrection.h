/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYNQCCORRECTION__
#define __INCLUDE_DYNQCCORRECTION__

#include "dyUtils.h"

int
XLALSimIMREOBNonQCCorrection (COMPLEX16 *nqc,	/**<< OUTPUT, The NQC correction */
			        REAL8Vector *values,
							/**<< Dynamics r, phi, pr, pphi */
			        const REAL8 omega,	/**<< Angular frequency */
			        EOBNonQCCoeffs *coeffs
							/**<< NQC coefficients */
  );

REAL8
XLALSimIMREOBGetNRSpinPeakDeltaTv4 (INT l,				/**<< Mode l */
				    INT m,				/**<< Mode m */
				    REAL8 m1,				/**<< mass 1 */
				    REAL8 m2,				/**<< mass 2 */
				    REAL8 chi1,				       /**<< Dimensionless spin1 */
				    REAL8 chi2				       /**<< Dimensionless spin2 */
  );

int
XLALSimIMRSpinEOBCalculateNQCCoefficientsV4 (REAL8Vector *amplitude,			   /**<< Waveform amplitude, func of time */
					     REAL8Vector *phase,			   /**<< Waveform phase(rad), func of time */
					     SpinEOBDynamics * dyHi,			   /**<< Position-vector, function of time */
					     REAL8Vector *orbOmegaVec,		   /**<< Orbital frequency, func of time */
					     INT modeL,						   /**<< Mode index l */
					     INT modeM,						   /**<< Mode index m */
					     REAL8 timePeak,					   /**<< Time of peak orbital frequency */
               REAL8 deltaT,
					     REAL8 m1,						   /**<< Component mass 1 */
					     REAL8 m2,						   /**<< Component mass 2 */
					     REAL8 a,						   /**<< Normalized spin of deformed-Kerr */
					     REAL8 chiA,					   /**<< Assymmetric dimensionless spin combination */
					     REAL8 chiS,					   /**<< Symmetric dimensionless spin combination */
					     EOBNonQCCoeffs *coeffs);

 int XLALSimIMRGetEOBCalibratedSpinNQC(EOBNonQCCoeffs *coeffs, 
                                    INT  l, 
                                    INT  m, 
                                    REAL8 eta, 
                                    REAL8 a );

#endif

