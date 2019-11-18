/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyFactorizedFlux.h"
#include "dyNewtonianMultipole.h"
#include "dyFactorizedWaveform.h"
#include "dyNQCCorrection.h"

/**
 * This function calculates the spin factorized-resummed GW energy flux
 * for given dynamical variables.
 */

REAL8
XLALInspiralSpinFactorizedFlux (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
				const REAL8 omega,	/**< orbital frequency */
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				INT const  lMax)
{

  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  COMPLEX16 hLM;
  INT l, m;

  /* Omegs is the derivative of phi */
  omegaSq = omega * omega;

  v = cbrt (omega);

  for (l = 2; l <= (INT) lMax; l++)
    {
      for (m = 1; m <= l; m++)
	{


        if (XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform
            (&hLM, values, v, H, l, m, ak) == CEV_FAILURE)
        {
            return CEV_FAILURE;
        }

	  /* For the 2,2 mode, we apply NQC correction to the flux */
	  if (l == 2 && m == 2)
	    {
	      COMPLEX16 hNQC;
	      XLALSimIMREOBNonQCCorrection (&hNQC, values, omega, nqcCoeffs);
	      /* Eq. 16 */
	      hLM *= hNQC;
	    }
	  //printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
	  /* Eq. 13 */
	  flux +=
	    (REAL8) (m * m) * omegaSq * (creal (hLM) * creal (hLM) +
					 cimag (hLM) * cimag (hLM));
	}
    }
  return flux * CST_1_PI / 8.0;
}

