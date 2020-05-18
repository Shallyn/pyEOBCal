/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "eobEccCorrection.h"
#include "dyHcapNumericalDerivative.h"
#include "dyNewtonianMultipole.h"
#include "dyNQCCorrection.h"
#include "dyFactorizedWaveform.h"
#include "dyFactorizedFlux.h"
#include <gsl/gsl_deriv.h>

#define STEP_SIZE 1.0e-4

#define RRtype 1

REAL8 PNCalcOrbitOmega(const REAL8 Hreal,
                       const REAL8 ecc,
                       const REAL8 eta);


REAL8 CalculatePNNonKeplerCoeff(const REAL8 values[],
                                const REAL8 Hreal,
                                const REAL8 ecc,
                                const REAL8 eta);

int XLALSpinAlignedHcapDerivative(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

REAL8 Calculate3PNEccentricity(const REAL8 eta,
                               REAL8Vector *values,
                               const REAL8 E /*Hreal - 1*/)
{
    REAL8 e2, twoEh2, sqrttwoEh2, twoE, twoE2, twoE3;
    REAL8 PN0, PN1, PN2, PN3;
    REAL8 eta2, eta3, PIsq;
    REAL8 pphi = values->data[3];
    twoEh2 = -2*E*pphi*pphi/eta;
    sqrttwoEh2 = sqrt(twoEh2);
    twoE = -2*E/eta;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    twoE2 = twoE*twoE;
    twoE3 = twoE2*twoE;
    PIsq = CST_PI * CST_PI;
    PN0 = 1 - twoEh2;
    PN1 = twoE * (24 - 4*eta +5*twoEh2*(-3+eta))/4;
    PN2 = twoE2*( 52+2*eta+2*eta2 -
                     (80-55*eta+4*eta2)*twoEh2 -
                     8*(11*eta-17)/twoEh2)/8;
    PN3 = twoE3*( -768-6*eta*PIsq-344*eta - 216*eta2 + 
                3*twoEh2*(-1488+1556*eta-319*eta2+4*eta3)-\
                4*(588 - 8212*eta +177*eta*PIsq+480*eta2)/twoEh2 +
                192*(134-218*eta+5*PIsq*eta+16*eta2)/twoEh2/twoEh2 )/192;
    e2 = PN0 + PN1 + PN2 + PN3;
    if(e2<=1e-5)
        return 0;
    return sqrt(e2);
}

REAL8 Calculate3PNxi(const REAL8 E, /*Hreal - 1*/
                     const REAL8 eta,
                     const REAL8 v,
                     const REAL8 l)
{
    REAL8 x, M, xi, constant, n;
    x = v*v;
    //M = 1-x*eta/2;
    //n = pow(-2*E/eta, 3/2);
    constant = 17/18 - 2*CST_GAMMA/3 - 2*CST_LN2;
    //xi = l - 3*M*n*(log(x) - constant);
    xi = l - 3*(pow(x, 3./2.) - pow(x, 5./2.)*(3+eta/2))*(log(x) - constant);
    return xi;
}

INT PNwaveformPRD100_044018_22mode_CalculateCoefficients(REAL8 eta,
                                                         PNEccCorrectCoeffs *EccCoeffs)
{
    REAL8 eta2, eta3, PIsq, LN32;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    PIsq = CST_PI*CST_PI;
    LN32 = log(3/2);

    EccCoeffs->cv0em = 1/4;
    EccCoeffs->cv0ep = 5/4;

    EccCoeffs->cv2em = -257/168 + 169*eta/168;
    EccCoeffs->cv2ep = -31/24 + 35*eta/24;

    EccCoeffs->cv3em = 11*CST_PI/4 + 27*I*LN32;
    EccCoeffs->cv3ep = 13*CST_PI + 3*I*CST_LN2;

    EccCoeffs->cv4em = -4271/756 - 35131*eta/6048 + 421*eta2/864;
    EccCoeffs->cv4ep = -2155/252-1655*eta/672 + 371*eta2/288;

    EccCoeffs->cv5em = -27*I/2-1081*CST_PI/168+ 
        (-1013*I/140+137*CST_PI/42)*eta +(27*I/4+9*I*eta)*LN32;
    EccCoeffs->cv5ep = -9*I/2 + 229*CST_PI/168 + 
        (-436571*I/420 + 61*CST_PI/42)*eta + (473*I/28-3*I*eta/7)*CST_LN2;

    EccCoeffs->cv6em = 219775769/1663200 + 749*I*CST_PI/60 + 49*PIsq/24 - 749*CST_GAMMA/30 + 
        (-121717/20790-41*PIsq/192)*eta-86531*eta2/8316 - 33331*eta3/399168 + 
        (-2889/70 + 81*I*CST_PI/2)*LN32 - 81*LN32*LN32/2;
    EccCoeffs->cv6ep = 55608313/1058400 + 3103*I*CST_PI/420 + 29*PIsq/24 - 3103*CST_GAMMA/210 + 
        (-199855/3024 + 41*PIsq/48)*eta - 9967*eta2/1008 + 35579*eta3/36288 + 
        (-6527/210 + 3*I*CST_PI/2)*CST_LN2 + 3*CST_LN2*CST_LN2/2;
    return CEV_SUCCESS;
}

static int PNwaveformPRD100_044018_22modeCorrect(const REAL8 ecc,
                                                 const REAL8 eta,
                                                 const REAL8 xi,
                                                 const REAL8 v,
                                                 PNEccCorrectCoeffs *EccCoeffs,
                                                 COMPLEX16 *out)
{
    COMPLEX16 PN0e, PN1e, PN1_5e, PN2e, PN2_5e, PN3e;
    COMPLEX16 emxi, epxi;
    REAL8 v2, v3, v4, v5, v6;
    v2 = v*v;
    v3 = v2*v;
    v4 = v3*v;
    v5 = v4*v;
    v6 = v5*v;
    emxi = cexp(-I*xi);
    epxi = cexp(I*xi);

    PN0e =  ecc*(EccCoeffs->cv0em * emxi + EccCoeffs->cv0ep * epxi);
    
    PN1e = v2*ecc*(EccCoeffs->cv2em * emxi + EccCoeffs->cv2ep * epxi);
    
    PN1_5e = v3*ecc*(EccCoeffs->cv3em * emxi + EccCoeffs->cv3ep * epxi);
    
    PN2e = v4*ecc*(EccCoeffs->cv4em * emxi + EccCoeffs->cv4ep * epxi);
    
    PN2_5e = v5*ecc*(EccCoeffs->cv5em * emxi + EccCoeffs->cv5ep * epxi);
    
    PN3e = v6*ecc*(EccCoeffs->cv6em * emxi + EccCoeffs->cv6ep * epxi -
        749*log(16*v2)*emxi/60 -3103*log(v2)*epxi/420 );
    *out = PN0e + PN1e + PN1_5e + PN2e + PN2_5e + PN3e;
    return CEV_SUCCESS;
}


/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of the paper.
 */
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
  )
{
  /* Status of function calls */
  INT status;
  INT i;


  REAL8 eta;
  REAL8 r, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;	//pr
  REAL8 Slm, deltalm, rholm;
  COMPLEX16 auxflm = 0.0;
  COMPLEX16 Tlm, rholmPwrl;
  COMPLEX16 hNewton, hEcc = 0.;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 vPhi, vPhi2;

  /* Pre-computed coefficients */
  FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

  if (abs (m) > l)
    {
      return CEV_FAILURE;
    }
  if (m == 0)
    {
      return CEV_FAILURE;
    }
  eta = params->eobParams->eta;

  /* Check our eta was sensible */
  if (eta > 0.25 && eta < 0.25 +1e-4) {
      eta = 0.25;
  }
  if (eta > 0.25)
    {
      print_warning
	("XLAL Error - %s: Eta = %.16f seems to be > 0.25 - this isn't allowed!\n",
	 __func__,eta);
      return CEV_FAILURE;
    }
  /*else if ( eta == 0.25 && m % 2 )
     {
     // If m is odd and dM = 0, hLM will be zero
     memset( hlm, 0, sizeof( COMPLEX16 ) );
     return XLAL_SUCCESS;
     } */

  r = values->data[0];
  //pr    = values->data[2];
  pp = values->data[3];

  v2 = v * v;
  Omega = v2 * v;
  vh3 = Hreal * Omega;
  vh = cbrt (vh3);
  eulerlogxabs = CST_GAMMA + log (2.0 * (REAL8) m * v);

  /* Calculate the non-Keplerian velocity */
  //params->alignedSpins
  if (1)
    {
      // YP: !!!!! SEOBNRv3devel temporary change !!!!!

	  //vPhi = CalculatePNNonKeplerCoeff(values->data, Hreal, ecc, eta) ;
    vPhi = 1.0 / (Omega * Omega * r*r*r);
      // YP: !!!!! SEOBNRv3devel temporary change !!!!!

      if (IS_REAL8_FAIL_NAN (vPhi))
	{
	  return CEV_FAILURE;
	}

      vPhi = r * cbrt (vPhi);
      vPhi *= Omega;
      vPhi2 = vPhi * vPhi;
    }
  else
    {
      vPhi = v;
      vPhi2 = v2;
    }
  /* Calculate the newtonian multipole, 1st term in Eq. 17, given by Eq. A1 */
  // YP: !!!!! SEOBNRv3devel temporary change !!!!!
  status = XLALSimIMRSpinEOBCalculateNewtonianMultipole (&hNewton, vPhi2, r,
							 values->data[1],
							 (UINT) l, m,
							 params->eobParams);
//print_debug("hNewton = %.2e + i%.2e\n", hNewton);

  if (status == CEV_FAILURE)
    {
      return CEV_FAILURE;
    }

  /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
  if (((l + m) % 2) == 0)
    {
      Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
  else
    {
      Slm = v * pp;
    }
  //printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

  /* Calculate the Tail term, 3rd term in Eq. 17, given by Eq. A6 */
  k = m * Omega;
  hathatk = Hreal * k;
  status = gsl_sf_lngamma_complex_e (l + 1.0, -2.0 * hathatk, &lnr1,
					  &arg1);
  if (status != GSL_SUCCESS)
    {
      print_warning ("XLAL Error - %s: Error in GSL function: %s\n",
                      __func__, gsl_strerror (status));
      return CEV_FAILURE;
    }
  status = gsl_sf_fact_e (l, &z2);
  if (status != GSL_SUCCESS)
    {
      print_warning ("XLAL Error - %s: Error in GSL function: %s\n",
                      __func__, gsl_strerror (status));
      return CEV_FAILURE;
    }
  Tlm =
    cexp ((lnr1.val + CST_PI * hathatk) +
	  I * (arg1.val + 2.0 * hathatk * log (4.0 * k / sqrt (CST_E))));
  Tlm /= z2.val;


  /* Calculate the residue phase and amplitude terms */
  /* deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15, others  */
  /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
  /* auxflm is a special part of the 5th term in Eq. 17, given by Eq. A15 */
  /* Actual values of the coefficients are defined in the next function of this file */
  switch (l)
    {
    case 2:
      switch (abs (m))
	{
	case 2:
	  deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
							+
							vh * vh *
							(hCoeffs->delta22vh9 *
							 vh))) +
	    hCoeffs->delta22v5 * v * v2 * v2 +
	    hCoeffs->delta22v6 * v2 * v2 * v2 +
	    hCoeffs->delta22v8 * v2 * v2 * v2 * v2;
	  rholm =
	    1. + v2 * (hCoeffs->rho22v2 +
		       v * (hCoeffs->rho22v3 +
			    v * (hCoeffs->rho22v4 +
				 v * (hCoeffs->rho22v5 +
				      v * (hCoeffs->rho22v6 +
					   hCoeffs->rho22v6l * eulerlogxabs +
					   v * (hCoeffs->rho22v7 +
						v * (hCoeffs->rho22v8 +
						     hCoeffs->rho22v8l *
						     eulerlogxabs +
						     (hCoeffs->rho22v10 +
						      hCoeffs->rho22v10l *
						      eulerlogxabs) *
						     v2)))))));
        status = PNwaveformPRD100_044018_22modeCorrect(ecc, eta, xi, v, params->eccCoeffs, &hEcc);

	  break;
	case 1:
	  {
	    deltalm = vh3 * (hCoeffs->delta21vh3 + vh3 * (hCoeffs->delta21vh6
							  +
							  vh *
							  (hCoeffs->
							   delta21vh7 +
							   (hCoeffs->
							    delta21vh9) * vh *
							   vh))) +
	      hCoeffs->delta21v5 * v * v2 * v2 +
	      hCoeffs->delta21v7 * v2 * v2 * v2 * v;
	    rholm =
	      1. + v * (hCoeffs->rho21v1 +
			v * (hCoeffs->rho21v2 +
			     v * (hCoeffs->rho21v3 +
				  v * (hCoeffs->rho21v4 +
				       v * (hCoeffs->rho21v5 +
					    v * (hCoeffs->rho21v6 +
						 hCoeffs->rho21v6l *
						 eulerlogxabs +
						 v * (hCoeffs->rho21v7 +
						      hCoeffs->rho21v7l *
						      eulerlogxabs +
						      v * (hCoeffs->rho21v8 +
							   hCoeffs->rho21v8l *
							   eulerlogxabs +
							   (hCoeffs->
							    rho21v10 +
							    hCoeffs->
							    rho21v10l *
							    eulerlogxabs) *
							   v2))))))));
                   auxflm = v * hCoeffs->f21v1;
               }
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 3:
      switch (m)
	{
	case 3:
    deltalm =
        vh3 * (hCoeffs->delta33vh3 +
             vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
             hCoeffs->delta33v5 * v * v2 * v2;
            rholm =
                1. + v2 * (hCoeffs->rho33v2 +
                           v * (hCoeffs->rho33v3 +
                                v * (hCoeffs->rho33v4 +
                                     v * (hCoeffs->rho33v5 +
                                          v * (hCoeffs->rho33v6 +
                                               hCoeffs->rho33v6l * eulerlogxabs +
                                               v * (hCoeffs->rho33v7 +
                                                    v * (hCoeffs->rho33v8 +
                                                     hCoeffs->rho33v8l *
                                                     eulerlogxabs)))))));
            auxflm = v * v2 * hCoeffs->f33v3;
	  break;
	case 2:
	  deltalm =
	    vh3 * (hCoeffs->delta32vh3 +
		   vh * (hCoeffs->delta32vh4 +
			 vh * vh * (hCoeffs->delta32vh6 +
				    hCoeffs->delta32vh9 * vh3)));
	  rholm =
	    1. + v * (hCoeffs->rho32v +
		      v * (hCoeffs->rho32v2 +
			   v * (hCoeffs->rho32v3 +
				v * (hCoeffs->rho32v4 +
				     v * (hCoeffs->rho32v5 +
					  v * (hCoeffs->rho32v6 +
					       hCoeffs->rho32v6l *
					       eulerlogxabs +
					       (hCoeffs->rho32v8 +
						hCoeffs->rho32v8l *
						eulerlogxabs) * v2))))));
	  break;
	case 1:
	  deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
							+
							vh *
							(hCoeffs->delta31vh7 +
							 hCoeffs->delta31vh9 *
							 vh * vh))) +
	    hCoeffs->delta31v5 * v * v2 * v2;
	  rholm =
	    1. + v2 * (hCoeffs->rho31v2 +
		       v * (hCoeffs->rho31v3 +
			    v * (hCoeffs->rho31v4 +
				 v * (hCoeffs->rho31v5 +
				      v * (hCoeffs->rho31v6 +
					   hCoeffs->rho31v6l * eulerlogxabs +
					   v * (hCoeffs->rho31v7 +
						(hCoeffs->rho31v8 +
						 hCoeffs->rho31v8l *
						 eulerlogxabs) * v))))));
	  auxflm = v * v2 * hCoeffs->f31v3;
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 4:
      switch (m)
	{
	case 4:
	  deltalm = vh3 * (hCoeffs->delta44vh3 + hCoeffs->delta44vh6 * vh3)
	    + hCoeffs->delta44v5 * v2 * v2 * v;

    rholm = 1. + v2 * (hCoeffs->rho44v2
             + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                    +
                    v *
                    (hCoeffs->
                     rho44v5 +
                     (hCoeffs->
                      rho44v6 +
                      hCoeffs->
                      rho44v6l *
                      eulerlogxabs) *
                     v))));
	  break;
	case 3:
	  deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
						       +
						       hCoeffs->delta43vh6 *
						       vh * vh));
	  rholm =
	    1. + v * (hCoeffs->rho43v +
		      v * (hCoeffs->rho43v2 +
			   v2 * (hCoeffs->rho43v4 +
				 v * (hCoeffs->rho43v5 +
				      (hCoeffs->rho43v6 +
				       hCoeffs->rho43v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f43v;
	  break;
	case 2:
	  deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
	  rholm = 1. + v2 * (hCoeffs->rho42v2
			     + v * (hCoeffs->rho42v3 +
				    v * (hCoeffs->rho42v4 +
					 v * (hCoeffs->rho42v5 +
					      (hCoeffs->rho42v6 +
					       hCoeffs->rho42v6l *
					       eulerlogxabs) * v))));
	  break;
	case 1:
	  deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
						       +
						       hCoeffs->delta41vh6 *
						       vh * vh));
	  rholm =
	    1. + v * (hCoeffs->rho41v +
		      v * (hCoeffs->rho41v2 +
			   v2 * (hCoeffs->rho41v4 +
				 v * (hCoeffs->rho41v5 +
				      (hCoeffs->rho41v6 +
				       hCoeffs->rho41v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f41v;
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 5:
      switch (m)
	{
	case 5:
	  deltalm =
	    hCoeffs->delta55vh3 * vh3 + hCoeffs->delta55v5 * v2 * v2 * v;
	  rholm =
	    1. + v2 * (hCoeffs->rho55v2 +
		       v * (hCoeffs->rho55v3 +
			    v * (hCoeffs->rho55v4 +
				 v * (hCoeffs->rho55v5 +
				      hCoeffs->rho55v6 * v))));
	  break;
	case 4:
	  deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
	  rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
						     + hCoeffs->rho54v4 * v));
	  break;
	case 3:
	  deltalm = hCoeffs->delta53vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho53v2
			     + v * (hCoeffs->rho53v3 +
				    v * (hCoeffs->rho53v4 +
					 hCoeffs->rho53v5 * v)));
	  break;
	case 2:
	  deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
	  rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
						     + hCoeffs->rho52v4 * v));
	  break;
	case 1:
	  deltalm = hCoeffs->delta51vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho51v2
			     + v * (hCoeffs->rho51v3 +
				    v * (hCoeffs->rho51v4 +
					 hCoeffs->rho51v5 * v)));
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 6:
      switch (m)
	{
	case 6:
	  deltalm = hCoeffs->delta66vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
						     + hCoeffs->rho66v4 * v));
	  break;
	case 5:
	  deltalm = hCoeffs->delta65vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
	  break;
	case 4:
	  deltalm = hCoeffs->delta64vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
						     + hCoeffs->rho64v4 * v));
	  break;
	case 3:
	  deltalm = hCoeffs->delta63vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
	  break;
	case 2:
	  deltalm = hCoeffs->delta62vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
						     + hCoeffs->rho62v4 * v));
	  break;
	case 1:
	  deltalm = hCoeffs->delta61vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 7:
      switch (m)
	{
	case 7:
	  deltalm = hCoeffs->delta77vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
	  break;
	case 6:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho76v2 * v2;
	  break;
	case 5:
	  deltalm = hCoeffs->delta75vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
	  break;
	case 4:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho74v2 * v2;
	  break;
	case 3:
	  deltalm = hCoeffs->delta73vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
	  break;
	case 2:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho72v2 * v2;
	  break;
	case 1:
	  deltalm = hCoeffs->delta71vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 8:
      switch (m)
	{
	case 8:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho88v2 * v2;
	  break;
	case 7:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho87v2 * v2;
	  break;
	case 6:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho86v2 * v2;
	  break;
	case 5:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho85v2 * v2;
	  break;
	case 4:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho84v2 * v2;
	  break;
	case 3:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho83v2 * v2;
	  break;
	case 2:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho82v2 * v2;
	  break;
	case 1:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho81v2 * v2;
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    default:
      return CEV_FAILURE;
      break;
    }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  for(i = 0 ; i < l ; i++) rholmPwrl *= rholm;
  /* In the equal-mass odd m case, there is no contribution from nonspin terms,
   * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
   * In this case, we must ignore the nonspin terms directly, since the leading term defined by
   * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
   */
  if (eta == 0.25 && m % 2)
    {
      rholmPwrl = auxflm;
    }
  else
    {
      rholmPwrl += auxflm;
    }
  /* Put all factors in Eq. 17 together */
  *hlm = Tlm * cexp (I * deltalm) * Slm * rholmPwrl + hEcc;
  
  *hlm *= hNewton;
//print_debug("hNew = %.2e + i%.2e, correct = %.2e + i%.2e\n", hNewton, Slm * rholmPwrl);
  /*if (r > 8.5)
     {
     printf("YP::FullWave: Reh = %.16e, Imh = %.16e, hAmp = %.16e, hPhi = %.16e\n",creal(*hlm),cimag(*hlm),cabs(*hlm),carg(*hlm));
     } */
  return CEV_SUCCESS;
}

REAL8 InspiralSpinFactorizedFlux_withecc(REAL8Vector *values,
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
    COMPLEX16 hNQC;
    REAL8 v;
    REAL8 omegaSq;
    REAL8 ecc, xi;
    COMPLEX16 hLM;
    INT l, m;
    
    if (lMax < 2)
    {
        return REAL8_FAIL_NAN;
    }
    
    /* Omegs is the derivative of phi */
    omegaSq = omega * omega;
    
    v = GET_CBRT (omega);
    ecc = Calculate3PNEccentricity(ak->eobParams->eta, values, H-1);
    xi = Calculate3PNxi(H-1, ak->eobParams->eta, v, values->data[1]);
    //  printf( "v = %.16e\n", v );
    for (l = 2; l <= (INT) lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {
			
            if (SpinEOBCalculateFactorizedWaveform_ecc(&hLM, values, v, H, l, m, ecc, xi, ak) == CEV_FAILURE)
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

int XLALSpinAlignedHcapDerivative_ecc(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
  HcapDerivParams params;

  /* Since we take numerical derivatives wrt dynamical variables */
  /* but we want them wrt time, we use this temporary vector in  */
  /* the conversion */
  REAL8           tmpDValues[6];

  /* Cartesian values for calculating the Hamiltonian */
  REAL8           cartValues[6];

  REAL8           H; //Hamiltonian
  REAL8           flux;

  gsl_function F;
  INT4         gslStatus;
  UINT i;

  REAL8Vector rVec, pVec;
  REAL8 rData[3], pData[3];

  /* We need r, phi, pr, pPhi to calculate the flux */
  REAL8       r;
  REAL8Vector polarDynamics;
  REAL8       polData[4];

  REAL8 Mtotal, eta;

  /* Spins */
  REAL8Vector *s1Vec = NULL;
  REAL8Vector *s2Vec = NULL;
  REAL8Vector *sKerr = NULL;
  REAL8Vector *sStar = NULL;

  REAL8 a;

  REAL8 omega, omegaorb;

  /* EOB potential functions */
  REAL8 DeltaT, DeltaR;
  REAL8 csi;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Declare NQC coefficients */
  EOBNonQCCoeffs *nqcCoeffs = NULL;

  /* Set up pointers for GSL */ 
  params.values  = cartValues;
  params.params  = (SpinEOBParams *)funcParams;
  nqcCoeffs = params.params->nqcCoeffs;

  s1Vec = params.params->s1VecOverMtMt;
  s2Vec = params.params->s2VecOverMtMt;
  sKerr = params.params->sigmaKerr;
  sStar = params.params->sigmaStar;

  F.function = &GSLSpinAlignedHamiltonianWrapper;
  F.params   = &params;

  Mtotal = params.params->eobParams->Mtotal;
  eta   = params.params->eobParams->eta;

  r = values[0];

  /* Since this is spin aligned, I make the assumption */
  /* that the spin vector is along the z-axis.         */
  a  = sKerr->data[2];

  /* Calculate the potential functions and the tortoise coordinate factor csi,
     given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
  DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );

  DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );

  csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
  //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
  //printf( "csi in derivatives function = %.16e\n", csi );

  /* Populate the Cartesian values vector, using polar coordinate values */
  /* We can assume phi is zero wlog */
  memset( cartValues, 0, sizeof( cartValues ) );
  cartValues[0] = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate derivatives w.r.t. each Cartesian variable */
  for ( i = 0; i < 6; i++ )
  {
    params.varyParam = i;
    gslStatus = gsl_deriv_central( &F, cartValues[i], 
                    STEP_SIZE, &tmpDValues[i], &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
      print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
      return CEV_FAILURE;
    }
  }

  /* Calculate the Cartesian vectors rVec and pVec */
  polarDynamics.length = 4;
  polarDynamics.data   = polData;

  memcpy( polData, values, sizeof( polData ) );

  rVec.length = pVec.length = 3;
  rVec.data   = rData;
  pVec.data   = pData;

  memset( rData, 0, sizeof(rData) );
  memset( pData, 0, sizeof(pData) );

  rData[0] = values[0];
  pData[0] = values[2];
  pData[1] = values[3] / values[0];
  /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
  H =  SpinEOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );

  //printf( "csi = %.16e, ham = %.16e ( tortoise = %d)\n", csi, H, params.params->tortoise );
  //exit(1);
  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
  H = H * Mtotal;

  /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
      tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/
  /* Now calculate omega, and hence the flux */

  omega = tmpDValues[4] / r;
  //omegaorb = PNCalcOrbitOmega(H/Mtotal, values[3], eta);

  REAL8 Fr, Fphi;

    dvalues[0] = csi * tmpDValues[3];
    dvalues[1] = omega;
    dvalues[2] = -tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    dvalues[2] = dvalues[2] * csi;
    dvalues[3] = 0;
//print_debug("RRtype = %d\n", RRtype);
  switch(RRtype)
  {
    case 0:
    {
      /* Default */
      flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, 8);
      if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
      {
          print_warning("Failed to calculate flux.\n");
          //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
          return CEV_FAILURE;
      }
      flux = flux / eta;
      Fr = ( values[2] / values[3] ) * flux / omega;
      Fphi = flux / omega;
      break;
    }
    case 1:
    {
      REAL8 Fr_N, Fr_1PN, Fr_2PN, Fphi_N, Fphi_1PN, Fphi_2PN;
      REAL8 x2, x2_2, x2_3, x2_4, x3, x3_2, x3_3, x3_4; 
      REAL8 x4, x4_2, x4_3, x4_4;
      REAL8 eta2, eta3, eta4, eta5, eta6;
      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;
      eta5 = eta4*eta;
      eta6 = eta5*eta;

      x2 = pow(values[2] / csi,2);
      x2_2 = x2*x2;
      x2_3 = x2_2*x2;
      x2_4 = x2_3*x2;

      x3 = 1./r;
      x3_2 = x3*x3;
      x3_3 = x3_2*x3;
      x3_4 = x3_3*x3;

      x4 = tmpDValues[0];
      x4 *= r;
      x4_2 = x4*x4;
      x4_3 = x4_2*x4;
      x4_4 = x4_3*x4;

      Fr_N = (32*eta*x3)/3. - (56*eta*x4)/5.;
      Fr_1PN = (4.761904761904762 + (4*eta)/7.)*eta*x2*x3 + (-43.161904761904765 - 
              (3776*eta)/105.)*eta*x3_2 + (-2.2095238095238097 - 
              (76*eta)/105.)*eta*x2*x4 + eta*(9.504761904761905 + 
              (1172*eta)/35.)*x3*x4 + (17.571428571428573 - 
              (400*eta)/21.)*eta*x4_2;
      Fr_2PN = (-1124*x2_2*x3)/63. - (30152728*x2*x3_2)/2835. + 
              (3011492*x3_3)/2835. + (9194*x2_2*x4)/315. + 
              (21596*x2*x3*x4)/105. - (153847*x3_2*x4)/63. - 
              (772*x2*x4_2)/21. + (78208*x3*x4_2)/63. - 
              (7361*x4_3)/105.;
      Fr = (Fr_N + Fr_1PN + Fr_2PN) * values[2] * x3_3;

      /* Fphi N */
      Fphi_N = (8*eta*x2)/5. - (32*eta*x3)/5. + (16*eta*x4)/5.;
      Fphi_1PN = (-5.219047619047619 - (236*eta)/105.)*eta*x2_2 + 
          eta*(-12.495238095238095 + (722*eta)/105.)*x2*x3 + 
          eta*(34.304761904761904 + (56*eta)/5.)*x3_2 + 
          (1.2190476190476192 - (38*eta)/21.)*eta*x2*x4 + (-4.3619047619047615 
          - (352*eta)/35.)*eta*x3*x4 + eta*(-2.6476190476190475 + 
          (256*eta)/105.)*x4_2;

      Fphi_2PN = (355*x2_3)/21. + (48404*x2_2*x3)/945. - 
          (8416*x2*x3_2)/45. - (967774*x3_3)/2835. - 
          (4226*x2_2*x4)/105. + (300247*x2*x3*x4)/630. + 
          (35846*x3_2*x4)/63. - (850833*x2*x4_2)/3710. -   
          (15604*x3*x4_2)/105. - (7247*x4_3)/105.;
      Fphi = (Fphi_N + Fphi_1PN + Fphi_2PN) * values[3] * x3_3;
      break;
    }
    case 2:
    {
      REAL8 x2, x2_2, x2_3, x2_4, x3, x3_2, x3_3, x3_4; 
      REAL8 x4, x4_2, x4_3, x4_4;
      REAL8 eta2, eta3, eta4, eta5, eta6;
      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;
      eta5 = eta4*eta;
      eta6 = eta5*eta;

      x2 = pow(values[2] / csi,2);
      x2_2 = x2*x2;
      x2_3 = x2_2*x2;
      x2_4 = x2_3*x2;

      x3 = 1./r;
      x3_2 = x3*x3;
      x3_3 = x3_2*x3;
      x3_4 = x3_3*x3;

      x4 = tmpDValues[0];
      x4 *= r;
      x4_2 = x4*x4;
      x4_3 = x4_2*x4;
      x4_4 = x4_3*x4;

      REAL8 A0r, A0f, A1r, A1f, A2r, A2f;
      REAL8 FrCirc, FphiCirc;
      FphiCirc = (2*(15876*eta2*x3 - 483887*x3_2 + 27*eta*(-336 + 1801*x3)))/2835.;
      FrCirc = (4*(-25488*eta2*x3 + 752873*x3_2 + eta*(7560 - 30591*x3))) / 2835.;

      A0f=(-188150*x2_3 + 447956*x2_2*x4 + 2552499*x2*x4_2 + 
              768182*x4_3 + 212*eta2*(118*x2_2 + 95*x2*x4 - 128*x4_2) + 
              212*eta*(274*x2_2 - 4*x2*(21 + 16*x4) + x4*(-168 + 139*x4)))/(71232.*eta);
      A1f=(-28*(96808*x2_2 + 900741*x2*x4 - 280872*x4_2) + 
              1764*eta3*(118*x2_2 + 95*x2*x4 - 128*x4_2) + 
              3*eta*(112896 + 493474*x2_2 + x2*(69132 - 115264*x4) - 225624*x4 + 250339*x4_2) + 
              3*eta2*(373630*x2_2 + 4*(19656 - 37199*x4)*x4 + x2*(-170688 + 133463*x4)))/(338688.*eta);
      A2f=((172961936 + 13340007*eta - 28581372*eta2 - 10753344*eta3)*x2 + 
              2*(-126468496 - 21768687*eta + 478656*eta2 + 2476656*eta3)*x4)/(4.064256e6*eta);

      A0r = (-21*x4)/20. - (29*x2*x4)/140. - (19*eta*x2*x4)/280. + 
              (4597*x2_2*x4)/(1680.*eta) + (369*x4_2)/224. - 
              (25*eta*x4_2)/14. - (193*x2*x4_2)/(56.*eta) - 
              (7361*x4_3)/(1120.*eta);
      A1r = 1 + (25*x2)/56. + (3*eta*x2)/56. - (281*x2_2)/(168.*eta) - 
              (18803*x4)/5600. - (561*eta*x4)/1400. - (32857*x2*x4)/39200. + 
              (5399*x2*x4)/(280.*eta) - (10897*eta*x2*x4)/11200. - 
              (1121*eta2*x2*x4)/4900. + (418077*x4_2)/62720. + 
              (2444*x4_2)/(21.*eta) - (3277*eta*x4_2)/1960. - 
              (295*eta2*x4_2)/49.;
      A2r = (5665*x2)/3136. - (3769091*x2)/(3780.*eta) + (551*eta*x2)/320. + 
              (177*eta2*x2)/980. - (21303799*x4)/1.568e6 - 
              (3134207*x4)/(25200.*eta) - (5073121*eta*x4)/392000. - 
              (33099*eta2*x4)/24500.;
    Fr = FrCirc * (A0r + A1r*x3 + A2r*x3_2) * values[2] * x3_3;
    Fphi = FphiCirc * (A0f + A1f*x3 + A2f*x3_2) * values[3] * x3_3;
    break;
    }
    case 3:
    {
      REAL8 Jcirc, omegaCirc;
      cartValues[0] = values[0];
      cartValues[3] = 0;
      Jcirc = auxCalculateCircularAngularMomentum(eta, &rVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs);
      if(! isnan(Jcirc) )
      {
          cartValues[4] = Jcirc / r;
          params.varyParam = 4;
          gslStatus = gsl_deriv_central( &F, cartValues[4], 
                          STEP_SIZE, &omegaCirc, &absErr );

          if ( gslStatus != GSL_SUCCESS )
          {
              print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
              return CEV_FAILURE;
          }
          omegaCirc /= r;
          polData[2] = 0;
          polData[3] = Jcirc;
          flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omegaCirc, params.params, H, 8 ,0);
          //flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, lMax);
          if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
          {
              print_warning("Failed to calculate flux.\n");
              //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
              return CEV_FAILURE;
          }
          flux /= eta;
      }
      else
      {
        cartValues[4] = values[3] / r;
        omegaCirc = omega;
        flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omegaCirc, params.params, H, 8 ,0);
        if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
        {
            print_warning("Failed to calculate flux.\n");
            //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
            return CEV_FAILURE;
        }
      }
      REAL8 x2, x2_2, x2_3, x2_4, x3, x3_2, x3_3, x3_4; 
      REAL8 x4, x4_2, x4_3, x4_4;
      REAL8 eta2, eta3, eta4, eta5, eta6;
      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;
      eta5 = eta4*eta;
      eta6 = eta5*eta;

      x2 = pow(values[2] / csi,2);
      x2_2 = x2*x2;
      x2_3 = x2_2*x2;
      x2_4 = x2_3*x2;

      x3 = 1./r;
      x3_2 = x3*x3;
      x3_3 = x3_2*x3;
      x3_4 = x3_3*x3;

      x4 = tmpDValues[0];
      x4 *= r;
      x4_2 = x4*x4;
      x4_3 = x4_2*x4;
      x4_4 = x4_3*x4;


      REAL8 A0r, A0f, A1r, A1f, A2r, A2f;
      REAL8 FrCirc, FrNC, FphiNC;
      //FphiCirc = (2*(15876*eta2*x3 - 483887*x3_2 + 27*eta*(-336 + 1801*x3)))/2835.;
      FrCirc = (4*(-25488*eta2*x3 + 752873*x3_2 + eta*(7560 - 30591*x3))) / 2835.;

      A0f=(-188150*x2_3 + 447956*x2_2*x4 + 2552499*x2*x4_2 + 
              768182*x4_3 + 212*eta2*(118*x2_2 + 95*x2*x4 - 128*x4_2) + 
              212*eta*(274*x2_2 - 4*x2*(21 + 16*x4) + x4*(-168 + 139*x4)))/(71232.*eta);
      A1f=(-28*(96808*x2_2 + 900741*x2*x4 - 280872*x4_2) + 
              1764*eta3*(118*x2_2 + 95*x2*x4 - 128*x4_2) + 
              3*eta*(112896 + 493474*x2_2 + x2*(69132 - 115264*x4) - 225624*x4 + 250339*x4_2) + 
              3*eta2*(373630*x2_2 + 4*(19656 - 37199*x4)*x4 + x2*(-170688 + 133463*x4)))/(338688.*eta);
      A2f=((172961936 + 13340007*eta - 28581372*eta2 - 10753344*eta3)*x2 + 
              2*(-126468496 - 21768687*eta + 478656*eta2 + 2476656*eta3)*x4)/(4.064256e6*eta);

      A0r = (-21*x4)/20. - (29*x2*x4)/140. - (19*eta*x2*x4)/280. + 
              (4597*x2_2*x4)/(1680.*eta) + (369*x4_2)/224. - 
              (25*eta*x4_2)/14. - (193*x2*x4_2)/(56.*eta) - 
              (7361*x4_3)/(1120.*eta);
      A1r = 1 + (25*x2)/56. + (3*eta*x2)/56. - (281*x2_2)/(168.*eta) - 
              (18803*x4)/5600. - (561*eta*x4)/1400. - (32857*x2*x4)/39200. + 
              (5399*x2*x4)/(280.*eta) - (10897*eta*x2*x4)/11200. - 
              (1121*eta2*x2*x4)/4900. + (418077*x4_2)/62720. + 
              (2444*x4_2)/(21.*eta) - (3277*eta*x4_2)/1960. - 
              (295*eta2*x4_2)/49.;
      A2r = (5665*x2)/3136. - (3769091*x2)/(3780.*eta) + (551*eta*x2)/320. + 
              (177*eta2*x2)/980. - (21303799*x4)/1.568e6 - 
              (3134207*x4)/(25200.*eta) - (5073121*eta*x4)/392000. - 
              (33099*eta2*x4)/24500.;
    FphiNC = pow(A0f,3)/(pow(A0f,2) + pow(A1f,2)*x3_2 - A0f*x3*(A1f + A2f*x3));
    //FrNC = pow(A0r,3)/(pow(A0r,2) + pow(A1r,2)*x3_2 - A0r*x3*(A1r + A2r*x3));
    FrNC = A0r + A1r * x3 + A2r * x3_2;
    Fphi = FrNC * flux / omegaCirc;
    Fr = FrCirc * FrNC * values[2] * x3_3;
    break;
    }
    case 4:
    {
      /* Default */
      flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, 8);
      if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
      {
          print_warning("Failed to calculate flux.\n");
          //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
          return CEV_FAILURE;
      }
      flux = flux / eta;
      Fphi = flux / omega;

      REAL8 A0, A1, A2, Fr_circ, Fr_nc;
      REAL8 x2, x2_2, x2_3, x2_4, x3, x3_2, x3_3, x3_4; 
      REAL8 x4, x4_2, x4_3, x4_4;
      REAL8 eta2, eta3, eta4, eta5, eta6;
      eta2 = eta*eta;
      eta3 = eta2*eta;
      eta4 = eta3*eta;
      eta5 = eta4*eta;
      eta6 = eta5*eta;

      x2 = pow(values[2] / csi,2);
      x2_2 = x2*x2;
      x2_3 = x2_2*x2;
      x2_4 = x2_3*x2;

      x3 = 1./r;
      x3_2 = x3*x3;
      x3_3 = x3_2*x3;
      x3_4 = x3_3*x3;

      x4 = tmpDValues[0];
      x4 *= r;
      x4_2 = x4*x4;
      x4_3 = x4_2*x4;
      x4_4 = x4_3*x4;

      A0 = (-21*x4)/20. - (29*x2*x4)/140. - (19*eta*x2*x4)/280. + 
              (4597*x2_2*x4)/(1680.*eta) + (369*x4_2)/224. - 
              (25*eta*x4_2)/14. - (193*x2*x4_2)/(56.*eta) - 
              (7361*x4_3)/(1120.*eta);
      A1 = 1 + (25*x2)/56. + (3*eta*x2)/56. - (281*x2_2)/(168.*eta) - 
              (18803*x4)/5600. - (561*eta*x4)/1400. - (32857*x2*x4)/39200. + 
              (5399*x2*x4)/(280.*eta) - (10897*eta*x2*x4)/11200. - 
              (1121*eta2*x2*x4)/4900. + (418077*x4_2)/62720. + 
              (2444*x4_2)/(21.*eta) - (3277*eta*x4_2)/1960. - 
              (295*eta2*x4_2)/49.;
      A2 = (5665*x2)/3136. - (3769091*x2)/(3780.*eta) + (551*eta*x2)/320. + 
              (177*eta2*x2)/980. - (21303799*x4)/1.568e6 - 
              (3134207*x4)/(25200.*eta) - (5073121*eta*x4)/392000. - 
              (33099*eta2*x4)/24500.;

      Fr_circ = x3 * (4*(-25488*eta2*x3 + 752873*x3_2 + eta*(7560 - 30591*x3))) / 2835.;
      Fr_nc = pow(A0,3)/(pow(A0,2) + pow(A1,2)*x3_2 - A0*x3*(A1 + A2*x3)) / x3;
      Fr = Fr_circ * Fr_nc * values[2] * x3_3;
      //Fr = ( values[2] / values[3] ) * flux / omega;
      break;
    }
    case 5:
    {
      // New Mode Sum Formulation
      REAL8 rdot, phidot;
      rdot = tmpDValues[0];
      phidot = tmpDValues[4] / r;
      flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, 8);
      if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
      {
          print_warning("Failed to calculate flux.\n");
          //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
          return CEV_FAILURE;
      }
      flux = flux / eta;
      Fphi = flux / omega;
      Fr = (flux - phidot * Fphi ) / rdot;
      break;
    }
    default:
    {
      print_warning("_%s: Invalid RRtype option %d\n", __func__, RRtype);
    }
  }
//print_debug("RR = %d, Fr = %.2e, Fphi = %.2e, v = (%.2e, %.2e, %.2e, %.2e)\n", RRtype, Fr, Fphi, r, values[2], values[3], tmpDValues[0]*r);
  //flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omega, params.params, H/Mtotal, lMax ,1);
    //flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, lMax);
//flux = InspiralSpinFactorizedFlux_withecc(&polarDynamics, nqcCoeffs, omegaorb, params.params, H/Mtotal, 8);

  /* Looking at the non-spinning model, I think we need to divide the flux by eta */

  //printf( "Flux in derivatives function = %.16e\n", flux );

  /* Now we can calculate the final (spherical) derivatives */
  /* csi is needed because we use the tortoise co-ordinate */
  /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
  //dvalues[0] = csi * tmpDValues[3];
  //dvalues[1] = omega;
  /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
  dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
  dvalues[2] = dvalues[2] * csi - Fr;
  dvalues[3] = - Fphi;

  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Values:\n%f %f %f %f\n", values[0], values[1], values[2], values[3] );

  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Derivatives:\n%f %f %f %f\n", dvalues[0], r*dvalues[1], dvalues[2], dvalues[3] );

  if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
  {
    print_warning( "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
    return CEV_FAILURE;
  }
  return CEV_SUCCESS;
}

/* SEOBNRE */
static int PNwaveformPRD544813rdotc_22mode(REAL8 *hr,
                                           REAL8 *hi,
                                           const REAL8 x1,
                                           const REAL8 x2,
                                           const REAL8 x3, // test particle position
                                           const REAL8 v1,
                                           const REAL8 v2,
                                           const REAL8 v3, // test particle velocity
                                           const REAL8 eta) // symmetric mass ratio
{
    const REAL8 dm = GET_SQRT(1-4*eta);
    REAL8 r = GET_SQRT(x1*x1+x2*x2+x3*x3);
    REAL8 n1 = x1/r,n2 = x2/r,n3 = x3/r;
    REAL8 rdot = v1*n1+v2*n2+v3*n3;
    REAL8 vsqr = v1*v1+v2*v2+v3*v3;
    
    REAL8 vn1,vn2,vn3;
    vn1 = rdot*n1;
    vn2 = rdot*n2;
    vn3 = rdot*n3;
    REAL8 lambda1,lambda2,lambda3;
    lambda1 = v1-vn1;
    lambda2 = v2-vn2;
    lambda3 = v3-vn3;
    REAL8 ln = GET_SQRT(lambda1*lambda1+lambda2*lambda2+lambda3*lambda3);
    lambda1 = lambda1/ln;
    lambda2 = lambda2/ln;
    lambda3 = lambda3/ln;
    
    COMPLEX16 hm22;
    
    COMPLEX16 h11,h12,h13,h22,h23,h33;
    COMPLEX16 Q11,Q12,Q13,Q22,Q23,Q33;
    
    ///////////////////////////////////////00 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)/3;
    Q12 = -GET_SQRT(CST_PI/5)/3*I;
    Q13 = 0;
    Q22 = -GET_SQRT(CST_PI/5)/3;
    Q23 = 8/GET_SQRT(5*CST_PI)/3;
    Q33 = 0;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    h11 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v1+v1*n1)+(3*(1-3*eta)*rdot*rdot)/r*n1*n1);
    h12 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v2+v1*n2)+(3*(1-3*eta)*rdot*rdot)/r*n1*n2);
    h13 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v3+v1*n3)+(3*(1-3*eta)*rdot*rdot)/r*n1*n3);
    h22 += 1.0/3*(2/r*rdot*(5+3*eta)*(n2*v2+v2*n2)+(3*(1-3*eta)*rdot*rdot)/r*n2*n2);
    h23 += 1.0/3*(2/r*rdot*(5+3*eta)*(n2*v3+v2*n3)+(3*(1-3*eta)*rdot*rdot)/r*n2*n3);
    h33 += 1.0/3*(2/r*rdot*(5+3*eta)*(n3*v3+v3*n3)+(3*(1-3*eta)*rdot*rdot)/r*n3*n3);
    
    hm22 = Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////10 part///////////////////////////////////////
    Q11 = 0;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(x1 - I*x2)/(3.*r);
    Q22 = 0;
    Q23 = GET_SQRT(CST_PI/5)*(I*x1 + x2)/(3.*r);
    Q33 = 0;
    
    h11 = dm*3/r*(-rdot*n1*n1);
    h12 = dm*3/r*(-rdot*n1*n2);
    h13 = dm*3/r*(-rdot*n1*n2);
    h22 = dm*3/r*(-rdot*n2*n2);
    h23 = dm*3/r*(-rdot*n2*n3);
    h33 = dm*3/r*(-rdot*n3*n3);
    
    h11 += dm/12/r*((n1*v1+v1*n1)*(rdot*rdot*(63+54*eta))
                    +n1*n1*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v1*(186+24*eta));
    h12 += dm/12/r*((n1*v2+v1*n2)*(rdot*rdot*(63+54*eta))
                    +n1*n2*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v2*(186+24*eta));
    h13 += dm/12/r*((n1*v3+v1*n3)*(rdot*rdot*(63+54*eta))
                    +n1*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v3*(186+24*eta));
    h22 += dm/12/r*((n2*v2+v2*n2)*(rdot*rdot*(63+54*eta))
                    +n2*n2*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v2*v2*(186+24*eta));
    h23 += dm/12/r*((n2*v3+v2*n3)*(rdot*rdot*(63+54*eta))
                    +n2*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v2*v3*(186+24*eta));
    h33 += dm/12/r*((n3*v3+v3*n3)*(rdot*rdot*(63+54*eta))
                    +n3*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v3*v3*(186+24*eta));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////01 part///////////////////////////////////////
    Q11 = 0;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(v1 - I*v2)/3.;
    Q22 = 0;
    Q23 = GET_SQRT(CST_PI/5)*(I*v1 + v2)/3.;
    Q33 = 0;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    h11 += dm*(-(n1*v1+v1*n1)/2/r*rdot*(7+4*eta)
               -n1*n1/r*(0.75*(1-2*eta)*rdot*rdot));
    h12 += dm*(-(n1*v2+v1*n2)/2/r*rdot*(7+4*eta)
               -n1*n2/r*(0.75*(1-2*eta)*rdot*rdot));
    h13 += dm*(-(n1*v3+v1*n3)/2/r*rdot*(7+4*eta)
               -n1*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    h22 += dm*(-(n2*v2+v2*n2)/2/r*rdot*(7+4*eta)
               -n2*n2/r*(0.75*(1-2*eta)*rdot*rdot));
    h23 += dm*(-(n2*v3+v2*n3)/2/r*rdot*(7+4*eta)
               -n2*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    h33 += dm*(-(n3*v3+v3*n3)/2/r*rdot*(7+4*eta)
               -n3*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////20 part///////////////////////////////////////
    Q11 = -GET_SQRT(CST_PI/5)*(x1*x1-8.*I*x1*x2-7*x2*x2-x3*x3)/21./r/r;
    Q12 = -I*GET_SQRT(CST_PI/5)*(3*(x1*x1+x2*x2)+x3*x3)/21./r/r;
    Q13 = 2*GET_SQRT(CST_PI/5)*(x1-I*x2)*x3/21./r/r;
    Q22 = GET_SQRT(CST_PI/5)*(-7*x1*x1+8.*I*x1*x2+x2*x2-x3*x3)/21./r/r;
    Q23 = -2.*I*GET_SQRT(CST_PI/5)*(x1-I*x2)*x3/21./r/r;
    Q33 = 8*GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)/21./r/r;
    
    h11 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n1+15*rdot*(n1*v1+v1*n1));
    h12 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n2+15*rdot*(n1*v2+v1*n2));
    h13 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n3+15*rdot*(n1*v3+v1*n3));
    h22 = (1-3*eta)/3/r*((-15*rdot*rdot)*n2*n2+15*rdot*(n2*v2+v2*n2));
    h23 = (1-3*eta)/3/r*((-15*rdot*rdot)*n2*n3+15*rdot*(n2*v3+v2*n3));
    h33 = (1-3*eta)/3/r*((-15*rdot*rdot)*n3*n3+15*rdot*(n3*v3+v3*n3));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////11 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(-v1*x1+4.*I*v2*x1+4.*I*v1*x2+7*v2*x2+v3*x3)/21./r;
    Q12 = -I*GET_SQRT(CST_PI/5)*(3*v1*x1+3*v2*x2+v3*x3)/21./r;
    Q13 = GET_SQRT(CST_PI/5)*(v3*(x1-I*x2)+x3*(v1-I*v2))/21./r;
    Q22 = GET_SQRT(CST_PI/5)*(-7*v1*x1+4.*I*v2*x1+4.*I*v1*x2+v2*x2-v3*x3)/21./r;
    Q23 = -GET_SQRT(CST_PI/5)*(v3*(I*x1+x2)+x3*(I*v1+v2))/21./r;
    Q33 = 8*GET_SQRT(CST_PI/5)*(v1-I*v2)*(x1-I*x2)/21./r;
    
    h11 = (1-3*eta)/3/r*(12*rdot*n1*n1);
    h12 = (1-3*eta)/3/r*(12*rdot*n1*n2);
    h13 = (1-3*eta)/3/r*(12*rdot*n1*n3);
    h22 = (1-3*eta)/3/r*(12*rdot*n2*n2);
    h23 = (1-3*eta)/3/r*(12*rdot*n2*n3);
    h33 = (1-3*eta)/3/r*(12*rdot*n3*n3);
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////02 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(-v1*v1+8.*I*v1*v2+7*v2*v2+v3*v3)/21.;
    Q12 = -I*GET_SQRT(CST_PI/5)*(3*(v1*v1+v2*v2)+v3*v3)/21.;
    Q13 = 2*GET_SQRT(CST_PI/5)*(v1-I*v2)*v3/21.;
    Q22 = GET_SQRT(CST_PI/5)*(-7*v1*v1+8.*I*v1*v2+v2*v2-v3*v3)/21.;
    Q23 = -2.*I*GET_SQRT(CST_PI/5)*(v1-I*v2)*v3/21.;
    Q33 = 8*GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)/21.;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////30 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)*x3/7./r/r/r;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(x1-I*x2)*(2.*x1*x1+5.*I*x1*x2+7*x2*x2+3*x3*x3)/21./r/r/r;
    Q22 = GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)*x3/7./r/r/r;
    Q23 = GET_SQRT(CST_PI/5)*(I*x1+x2)*(7*x1*x1-5.*I*x1*x2+2*x2*x2+3*x3*x3)/21./r/r/r;
    Q33 = -2*GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)*x3/7./r/r/r;
    
    h11 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n1-8.5*rdot*v1*v1-(-105*rdot*rdot)/12*(n1*v1+v1*n1));
    h12 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n2-8.5*rdot*v1*v2-(-105*rdot*rdot)/12*(n1*v2+v1*n2));
    h13 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n3-8.5*rdot*v1*v3-(-105*rdot*rdot)/12*(n1*v3+v1*n3));
    h22 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n2*n2-8.5*rdot*v2*v2-(-105*rdot*rdot)/12*(n2*v2+v2*n2));
    h23 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n2*n3-8.5*rdot*v2*v3-(-105*rdot*rdot)/12*(n2*v3+v2*n3));
    h33 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n3*n3-8.5*rdot*v3*v3-(-105*rdot*rdot)/12*(n3*v3+v3*n3));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////21 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*((x1-I*x2)*(2*v1*x1+I*x1*v2+4.*I*v1*x2+7*x2*v2)+2*v3*(x1-I*x2)*x3+(v1-I*v2)*x3*x3)/21./r/r;
    Q22 = GET_SQRT(CST_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    Q23 = GET_SQRT(CST_PI/5)*((x1-I*x2)*(v2*(4*x1+2.*I*x2)+v1*(7.*I*x1+x2))+2*v3*(I*x1+x2)*x3+(I*v1+v2)*x3*x3)/21./r/r;
    Q33 = -2*GET_SQRT(CST_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    
    h11 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n1-54*rdot*(n1*v1+v1*n1));
    h12 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n2-54*rdot*(n1*v2+v1*n2));
    h13 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n3-54*rdot*(n1*v3+v1*n3));
    h22 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n2*n2-54*rdot*(n2*v2+v2*n2));
    h23 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n2*n3-54*rdot*(n2*v3+v2*n3));
    h33 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n3*n3-54*rdot*(n3*v3+v3*n3));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////12 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(v3*v3*(x1-I*x2)+v1*v1*(2*x1+I*x2)+v2*v2*(4*x1-7.*I*x2)-2.*I*v2*v3*x3+2*v1*(I*v2*x1+4*v2*x2+v3*x3))/21./r;
    Q22 = GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    Q23 = GET_SQRT(CST_PI/5)*(v3*v3*(I*x1+x2)+v2*v2*(-I*x1+2*x2)+v1*v1*(7.*I*x1+4*x2)+2*v2*v3*x3+v1*(8*v2*x1-2.*I*v2*x2+2.*I*v3*x3))/21./r;
    Q33 = -2*GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    
    h11 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n1);
    h12 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n2);
    h13 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n3);
    h22 = dm*(1-2*eta)*1.5/r*(-3*rdot*n2*n2);
    h23 = dm*(1-2*eta)*1.5/r*(-3*rdot*n2*n3);
    h33 = dm*(1-2*eta)*1.5/r*(-3*rdot*n3*n3);
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////03 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)*v3/7.;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v1*v1+5.*I*v1*v2+7*v2*v2+3*v3*v3)/21.;
    Q22 = GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)*v3/7.;
    Q23 = GET_SQRT(CST_PI/5)*(I*v1+v2)*(7*v1*v1-5.*I*v1*v2+2*v2*v2+3*v3*v3)/21.;
    Q33 = -2*GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)*v3/7.;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    //-----------------------------------complete-----------------------------------------
    
    hm22 = hm22*2.*eta;
    *hr = C_REAL(hm22);
    *hi = C_IMAG(hm22);
    if (isnan(*hr) || isnan(*hi))
    {
        print_warning("Error _%s: hECC is nan.\n",__func__);
        //print_err("\thm22 = %f + i%f", C_REAL(hm22), C_IMAG(hm22));
        //print_err("\tx1 = %f, x2 = %f, x3 = %f\n", x1, x2, x3);
        //print_err("\tv1 = %f, v2 = %f, v3 = %f\n", v1, v2, v3);
        //print_err("")
        return CEV_FAILURE;
    }

    return CEV_SUCCESS;
}

INT EOBEccFactorizedWaveformCorrection(COMPLEX16 *hECC,
                                       const REAL8 rphivalues[],
                                       const REAL8 rdot,
                                       const REAL8 phidot,
                                       const REAL8 eta)
{
    REAL8 r,phi,prt,pphi;
    r = rphivalues[0];
    prt = rphivalues[2];
    phi = rphivalues[1];
    pphi = rphivalues[3];

    REAL8 x1,x2,x3;
    x1 = r*GET_COS(phi);
    x2 = r*GET_SIN(phi);
    x3 = 0;
    
    REAL8 v1,v2,v3;
    v1 = rdot*GET_COS(phi)-r*phidot*GET_SIN(phi);
    v2 = rdot*GET_SIN(phi)+r*phidot*GET_COS(phi);
    v3 = 0;
    
    REAL8 hr, hi;
    if(PNwaveformPRD544813rdotc_22mode(&hr, &hi,x1,x2,x3,v1,v2,v3,eta) == CEV_FAILURE)
    {
        print_warning("Error _%s: hECC is nan.\n", __func__);
        //print_err("\tphidot = %f, pphi = %f, rdot = %f\n", phidot, pphi, rdot);
        //print_err("\tv1 = %f, v2 = %f, v3 = %f\n", v1, v2, v3);
        return CEV_FAILURE;
    }
    //print_debug("hECC = %e + i%e\n", hr, hi);
    *hECC = hr + I*hi;
    return CEV_SUCCESS;
}

INT EOBCalculateEccCorrection(COMPLEX16 *hECC,
                              const REAL8Vector *values,
                              SpinEOBParams *seobParams)
{
    REAL8 r,phi,pr,pphi,rdot,phidot;
    
    r = values->data[0];
    phi = values->data[1];
    pr = values->data[2];
    pphi = values->data[3];
    
    REAL8 dydt[4],ttmp=0;
    if (XLALSpinAlignedHcapDerivative(ttmp, values->data, dydt, seobParams) == GSL_FAILURE)
        return CEV_FAILURE;
    
    rdot = dydt[0];
    phidot = dydt[1];
    
    REAL8 x1,x2,x3;
    x1 = r*GET_COS(phi);
    x2 = r*GET_SIN(phi);
    x3 = 0;
    
    REAL8 v1,v2,v3;
    v1 = rdot*GET_COS(phi)-r*phidot*GET_SIN(phi);
    v2 = rdot*GET_SIN(phi)+r*phidot*GET_COS(phi);
    v3 = 0;
    
    REAL8 hr, hi;
    if (PNwaveformPRD544813rdotc_22mode(&hr,&hi,x1,x2,x3,v1,v2,v3,seobParams->eobParams->eta) == CEV_FAILURE)
    {
        print_warning("Error _%s: hECC is nan.\n", __func__);
        //print_err("\tphidot = %f, pphi = %f, rdot = %f\n", phidot, pphi, rdot);
        return CEV_FAILURE;
    }
    *hECC = hr + I*hi;
    return CEV_SUCCESS;
}

