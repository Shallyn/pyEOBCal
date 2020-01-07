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

INT PNwaveformPRD100_044018_22mode_CalculateCoefficients(REAL8 eta,
                                                         PNEccCorrectCoeffs *EccCoeffs);

REAL8 Calculate3PNxi(const REAL8 E, /*Hreal - 1*/
                     const REAL8 eta,
                     const REAL8 v,
                     const REAL8 l);
REAL8 Calculate3PNEccentricity(const REAL8 eta,
                               REAL8Vector *values,
                               const REAL8 E /*Hreal - 1*/);


INT PNwaveformPRD100_044018_22mode_CalculateCoefficients(REAL8 eta,
                                                         PNEccCorrectCoeffs *EccCoeffs);

INT
SpinEOBCalculateFactorizedWaveform_ecc(COMPLEX16 *  hlm,
						      /**< OUTPUT, hlm waveforms */
					    REAL8Vector *  values,
						      /**< dyanmical variables */
					    const REAL8 v,
						      /**< velocity */
					    const REAL8 Hreal,
						      /**< real Hamiltonian */
					    const INT l,
						      /**< l mode index */
					    const INT m,
						      /**< m mode index */
                        const REAL8 ecc ,
                              /**< eccentricity */
                        const REAL8 xi,
					    SpinEOBParams *  params
						       /**< Spin EOB parameters */
  );

int XLALSpinAlignedHcapDerivative_ecc(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

#endif

