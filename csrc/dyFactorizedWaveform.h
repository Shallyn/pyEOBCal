/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYFACTORIZEDWAVEFORM__
#define __INCLUDE_DYFACTORIZEDWAVEFORM__

#include "dyUtils.h"

INT CalculateSpinFactorizedWaveformCoefficients(FacWaveformCoeffs *const coeffs,
                                                SpinEOBParams * params,
                                                const REAL8 m1,
                                                const REAL8 m2,
                                                const REAL8 eta,
                                                REAL8 a,
                                                const REAL8 chiS,
                                                const REAL8 chiA);

INT
XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform (COMPLEX16 *  hlm,
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
						SpinEOBParams *  params);

INT
XLALSimIMRSpinEOBGetSpinFactorizedWaveform (COMPLEX16 *  hlm,
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
					    SpinEOBParams *  params
						       /**< Spin EOB parameters */
  );

#endif

