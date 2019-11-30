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
#include "eobEccCorrection.h"

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
				INT const  lMax,
				INT allow_ecc)
{

  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  COMPLEX16 hLM, hNQC, hECC;
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
	      
	      	XLALSimIMREOBNonQCCorrection (&hNQC, values, omega, nqcCoeffs);
		  	
#if ALLOW_ECC
			if(allow_ecc)
			{
				if (EOBCalculateEccCorrection(&hECC, values, ak) == CEV_FAILURE)
				{
					print_warning("Error _%s: Failed to apply ECC Factorized waveform correction.\n", __func__);
					return REAL8_FAIL_NAN;
				}
				//print_err("DEBUG: hECC = %.5e + i%.5e, hLM = %.5e + i%.5e\n", C_REAL(hECC), C_IMAG(hECC), C_REAL(hLM), C_IMAG(hLM));
				hLM += hECC;
				//print_err("DEBUG: hLM_re / hLM = %.5e\n", C_ABS(hLM_re) / C_ABS(hLM));
			}
#endif	
            hLM *= hNQC;
	      /* Eq. 16 */
	      
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


REAL8 InspiralSpinFactorizedFlux_elip(REAL8Vector *values,
                                      const REAL8 rphivalues[],
                                      const REAL8 drphivalues[],
                                      EOBNonQCCoeffs *nqcCoeffs,
                                      const REAL8 omega,
                                      SpinEOBParams *ak,
                                      const REAL8 H,
                                      const UINT lMax)
{
    if ( nqcCoeffs==NULL )
    {
        return REAL8_FAIL_NAN;
    }
    REAL8 flux = 0.0;// flux_re = 0.0;
    COMPLEX16 hECC, hNQC;
    REAL8 v;
    REAL8 omegaSq;
    COMPLEX16 hLM, hLM_re;
    INT l, m;
    
    if (lMax < 2)
    {
        return REAL8_FAIL_NAN;
    }
    
    /* Omegs is the derivative of phi */
    omegaSq = omega * omega;
    
    v = GET_CBRT (omega);
    
    //  printf( "v = %.16e\n", v );
    for (l = 2; l <= (INT) lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {
			
            if (XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform(&hLM, values, v, H, l, m, ak) == CEV_FAILURE)
                {
                    return REAL8_FAIL_NAN;
                }
			//print_debug("values = (%e, %e, %e, %e), v = %e, H = %e, hLM = %e + i%e\n", 
				//values->data[0], values->data[1], values->data[2], values->data[3], v, H, creal(hLM), cimag(hLM));
            /* For the 2,2 mode, we apply NQC correction to the flux */
            if (l == 2 && m == 2)
            {
                XLALSimIMREOBNonQCCorrection (&hNQC, values, omega, nqcCoeffs);
                /* Eq. 16 */
                
#if ALLOW_ECC
                if (EOBEccFactorizedWaveformCorrection(&hECC, rphivalues, drphivalues[0], drphivalues[3], ak->eobParams->eta) == CEV_FAILURE)
                {
                    print_warning("Error _%s: Failed to apply ECC Factorized waveform correction.\n", __func__);
                    return REAL8_FAIL_NAN;
                }
                //print_debug("DEBUG: hECC = %.5e + i%.5e, hLM = %.5e + i%.5e\n", C_REAL(hECC), C_IMAG(hECC), C_REAL(hLM), C_IMAG(hLM));
                hLM += hECC;
                //print_err("DEBUG: hLM_re / hLM = %.5e\n", C_ABS(hLM_re) / C_ABS(hLM));
#endif
                hLM *= hNQC;
            }
            
            //printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
            /* Eq. 13 */
            flux +=
            (REAL8) (m * m) * omegaSq * (C_REAL (hLM) * C_REAL(hLM) + C_IMAG (hLM) * C_IMAG (hLM));
            /* DEBUG:
            if (isnan(flux))
            {
                print_warning("Flux is nan: flux = %e\n", flux);
                print_err("\thECC = %f + i%f\n", C_REAL(hECC), C_IMAG(hECC));
                print_err("\tomegaSq = %f\n", omegaSq);
            }*/
            //flux_re += (REAL8) (m * m) * omegaSq * (C_REAL (hLM_re) * C_REAL(hLM_re) + C_IMAG (hLM_re) * C_IMAG (hLM_re));
        }
    }
    //print_err("DEBUG: flux = %.5e, flux_re = %.5e\n", flux, flux_re);
    return flux * CST_1_PI / 8.0;
}

