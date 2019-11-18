/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyNQCCorrection.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int
XLALSimIMREOBNonQCCorrection (COMPLEX16 *nqc,	/**<< OUTPUT, The NQC correction */
			        REAL8Vector *values,
							/**<< Dynamics r, phi, pr, pphi */
			        const REAL8 omega,	/**<< Angular frequency */
			        EOBNonQCCoeffs *coeffs
							/**<< NQC coefficients */
  )
{

  REAL8 rOmega, rOmegaSq;
  REAL8 r, p, sqrtR;

  REAL8 mag, phase;


  r = values->data[0];
  p = values->data[2];

  sqrtR = sqrt (r);

  rOmega = r * omega;
  rOmegaSq = rOmega * rOmega;
/*printf("a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 = %.16e\n",coeffs->a1,coeffs->a2,coeffs->a3,coeffs->a3S, coeffs->a4,coeffs->a5);
printf("b1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n",coeffs->b1,coeffs->b2,coeffs->b3,coeffs->b4);*/
  /* In EOBNRv2, coeffs->a3S, coeffs->a4 and coeffs->a5 are set to zero */
  /* through XLALSimIMREOBGetCalibratedNQCCoeffs() */
  /* and XLALSimIMREOBCalculateNQCCoefficients() */
  mag = 1. + (p * p / rOmegaSq) * (coeffs->a1
				   + coeffs->a2 / r + (coeffs->a3 +
						       coeffs->a3S) / (r *
								       sqrtR)
				   + coeffs->a4 / (r * r) +
				   coeffs->a5 / (r * r * sqrtR));
//printf("NQC INFO mag = %.16e, r = %.16e, p = %.16e\n",mag,r,p);
  phase = coeffs->b1 * p / rOmega + p * p * p / rOmega * (coeffs->b2
							  +
							  coeffs->b3 / sqrtR +
							  coeffs->b4 / r);

  *nqc = mag * cos (phase);
  *nqc += I * mag * sin (phase);
/*printf("r = %.16e, pr = %.16e, omega = %.16e\n",r,p,omega);
printf("NQC mag = %.16e, arg = %.16e\n",mag,phase);*/
  return CEV_SUCCESS;
}

static inline REAL8 CombineTPLEQMFits (REAL8 eta, REAL8 A1, REAL8 fEQ, REAL8 fTPL)
{
    REAL8 A0, A2;
    REAL8 eta2 = eta * eta;
    // Impose that TPL and equal-mass limit are exactly recovered
    A0 = -0.00099601593625498 * A1 - 0.00001600025600409607 * fEQ + 1.000016000256004 * fTPL;
    A2 = -3.984063745019967 * A1 + 16.00025600409607 * fEQ - 16.0002560041612 * fTPL;
    // Final formula
    return A0 + A1 * eta + A2 * eta2;
}


REAL8
XLALSimIMREOBGetNRSpinPeakDeltaTv4 (INT l,				/**<< Mode l */
				    INT m,				/**<< Mode m */
				    REAL8 m1,				/**<< mass 1 */
				    REAL8 m2,				/**<< mass 2 */
				    REAL8 chi1,				       /**<< Dimensionless spin1 */
				    REAL8 chi2				       /**<< Dimensionless spin2 */
  )
{
  REAL8 eta = m1 * m2 / (m1 + m2) / (m1 + m2);
  REAL8 chi =
    0.5 * (chi1 + chi2) + 0.5 * (chi1 - chi2) * (m1 - m2) / (m1 + m2) / (1. -
									 2. *
									 eta);
  REAL8 eta2 = eta * eta, eta3 = eta2 * eta;
  REAL8 chiTo2 = chi * chi, chiTo3 = chiTo2 * chi;
  REAL8 coeff00, coeff01, coeff02, coeff03;
  REAL8 coeff10, coeff11, coeff12, coeff13;
  REAL8 coeff20, coeff21, coeff22, coeff23;
  REAL8 coeff30, coeff31, coeff32, coeff33;
  REAL8 res;
                  // Calibrationv21_Sep8a
                  coeff00 = 2.50499;
                  coeff01 = 13.0064;
                  coeff02 = 11.5435;
                  coeff03 = 0;
                  coeff10 = 45.8838;
                  coeff11 = -40.3183;
                  coeff12 = 0;
                  coeff13 = -19.0538;
                  coeff20 = 13.0879;
                  coeff21 = 0;
                  coeff22 = 0;
                  coeff23 = 0.192775;
                  coeff30 = -716.044;
                  coeff31 = 0;
                  coeff32 = 0;
                  coeff33 = 0;
                  res = coeff00 + coeff01 * chi + coeff02 * chiTo2 + coeff03 * chiTo3 +
                    coeff10 * eta + coeff11 * eta * chi + coeff12 * eta * chiTo2 +
                    coeff13 * eta * chiTo3 + coeff20 * eta2 + coeff21 * eta2 * chi +
                    coeff22 * eta2 * chiTo2 + coeff23 * eta2 * chiTo3 + coeff30 * eta3 +
                    coeff31 * eta3 * chi + coeff32 * eta3 * chiTo2 + coeff33 * eta3 * chiTo3;

	//RC: for the 55 mode the attachment is done at tpeak22 -10M, note that since here deltat22 is defined as -deltat22 with respect
	//to SEOBNRv4 paper, here I need to add 10 M
if((l == 5)&&(m == 5)){
	res = res + 10;
}
  return res;
}


/**
 * Peak amplitude predicted by fitting NR results
 * The coefficients for the mode (2,2) are in Eq.(A1) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.95.044028
 * The coefficients for the modes (2,1), (3,3), (4,4), (5,5) are in Eqs. (B7-B10) of https://arxiv.org/pdf/1803.10701.pdf
 */
static REAL8
XLALSimIMREOBGetNRSpinPeakAmplitudeV4 (INT modeL, INT modeM, REAL8 m1, REAL8 m2,
			  REAL8 chiS, REAL8 chiA)
{
	REAL8 eta = (m1 * m2) / ((m1 + m2) * (m1 + m2));
	REAL8 dM = sqrt(1.-4.*eta);
	REAL8 tempm1;
	if (m1 < m2)
		{
			//RC: The fits for the HMs are done under the assumption m1>m2, so if m2>m1 we just swap the two bodies
			tempm1 = m1;
			m1 = m2;
			m2 = tempm1;
			chiA = -chiA;
		}
	REAL8 eta2 = eta*eta;
  REAL8 chi = chiS + chiA * (m1 - m2) / (m1 + m2) / (1. -
					 2. *
					 eta);
	REAL8 chi21 = chiS*dM/(1.-1.3*eta) + chiA;
	REAL8 chi33 = chiS*dM + chiA;
	REAL8 chi44 = chiS*(1-5*eta) + chiA*dM;
	REAL8 chi2 = chi * chi, chi3 = chi * chi2;
  REAL8 res;
  REAL8 fTPL, fEQ, A1, e0, e1, e2, e3;
  switch (modeL) {
      case 2:
          switch (modeM) {
              case 2:
                  // TPL fit
                  fTPL = 1.4528573105413543 + 0.16613449160880395 * chi + 0.027355646661735258 * chi2 - 0.020072844926136438 * chi3;
                  // Equal-mass fit
                  fEQ = 1.577457498227 - 0.0076949474494639085 * chi +  0.02188705616693344 * chi2 + 0.023268366492696667 * chi3;
                  // Global fit coefficients
                  e0 = -0.03442402416125921;
                  e1 = -1.218066264419839;
                  e2 = -0.5683726304811634;
                  e3 = 0.4011143761465342;
                  A1 = e0 + e1 * chi + e2 * chi2 + e3 * chi3;
                  res = eta * CombineTPLEQMFits(eta, A1, fEQ, fTPL);
                  break;
							case 1:
									res = -((0.29256703361640224-0.19710255145276584*eta)*eta*chi21 + dM*eta*(-0.42817941710649793 + 0.11378918021042442*eta-0.7736772957051212*eta2
									+chi21*chi21*(0.047004057952214004 - eta*0.09326128322462478)) +dM*eta*chi21*(-0.010195081244587765 + 0.016876911550777807*chi21*chi21));
                  break;
              default:
                  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                  return REAL8_FAIL_NAN;
                  break;
          }
          break;

      case 3:
          switch (modeM) {
              case 3:
									res = (0.10109183988848384*eta - 0.4704095462146807*eta2 + 1.0735457779890183*eta2*eta)*chi33 +
									dM*(0.5636580081367962*eta - 0.054609013952480856*eta2 + 2.3093699480319234*eta2*eta + chi33*chi33*(0.029812986680919126*eta - 0.09688097244145283*eta2) );
                  break;
              default:
                  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                  return REAL8_FAIL_NAN;
                  break;
          }
          break;
			case 4:
						switch (modeM) {
							case 4:
								res = eta*(0.2646580063832686 + 0.067584186955327*chi44 +0.02925102905737779*chi44*chi44) + eta2 *(-0.5658246076387973 -0.8667455348964268*chi44 +0.005234192027729502*chi44*chi44)
								+ eta*eta2*(-2.5008294352355405 + 6.880772754797872*chi44 -1.0234651570264885*chi44*chi44) + eta2*eta2*(7.6974501716202735 -16.551524307203252*chi44);
								break;
							default:
									print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
									return REAL8_FAIL_NAN;
								break;
						}
					break;
			case 5:
						switch (modeM) {
							case 5:
								res = 0.128621*dM *eta -0.474201 *dM *eta*eta +1.0833 * dM *eta*eta*eta + 0.0322784 * eta * chi33 -0.134511 *chi33 *eta * eta +0.0990202 *chi33*eta*eta*eta;
                break;
							default:
									print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
									return REAL8_FAIL_NAN;
								break;
						}
					break;

      default:
          print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
          return REAL8_FAIL_NAN;
          break;
  }
//    printf("A %.16e\n", res);
  return res;
}

/**
 * Peak amplitude slope predicted by fitting NR results.
 * The coefficient for the (2,2) is 0 because the
 * attachment is done at the peak of the (2,2) mode see Eq.(2.7) in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.95.044028
 * The coefficients for the modes (2,1), (3,3), (4,4) and (5,5) are in Eqs.(B11-B14) of https://arxiv.org/pdf/1803.10701.pdf
 */
static  REAL8
XLALSimIMREOBGetNRSpinPeakADotV4 (INT modeL, INT modeM, REAL8  m1, REAL8 m2,
			  REAL8  chiS, REAL8  chiA)

{
	  REAL8 eta = (m1 * m2) / ((m1 + m2) * (m1 + m2));
		REAL8 dM = sqrt(1.-4.*eta);
		REAL8 dM2 = dM*dM;
		REAL8 tempm1;
		if (m1 < m2)
			{
				//RC: The fits for the HMs are done under the assumption m1>m2, so if m2>m1 we just swap the two bodies
				tempm1 = m1;
				m1 = m2;
				m2 = tempm1;
				chiA = -chiA;
			}
		REAL8 eta2 = eta*eta;
		REAL8 chi21 = chiS*dM/(1.-2.*eta) + chiA;
		REAL8 chi33 = chiS*dM + chiA;
		REAL8 chi44 = chiS*(1-7*eta) + chiA*dM;
		REAL8 res;
    switch (modeL) {
        case 2:
            switch (modeM) {
                case 2:
                    res = 0.;
                    break;
								case 1:
										res = dM*eta*(0.007147528020812309-eta*0.035644027582499495) + dM*eta*chi21*(-0.0087785131749995 + eta*0.03054672006241107) + eta*0.00801714459112299*fabs(-dM*(0.7875612917853588 + eta*1.161274164728927 + eta2*11.306060006923605)+chi21);
								break;
                default:
                    print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                    return REAL8_FAIL_NAN;
                    break;
            }
            break;
        case 3:
            switch (modeM) {
                case 3:
                    res = dM*eta*(-0.00309943555972098 + eta*0.010076527264663805)*chi33*chi33 +
										eta*0.0016309606446766923*sqrt(dM2*(8.811660714437027 + 104.47752236009688*eta) + dM*chi33*(-5.352043503655119 + eta*49.68621807460999) + chi33*chi33);
                    break;
                default:
                    print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                    return REAL8_FAIL_NAN;
                    break;
            }
            break;
				case 4:
						switch (modeM) {
							case 4:
										res = eta*(0.004347588211099233 -0.0014612210699052148*chi44 -0.002428047910361957*chi44*chi44) + eta2*(0.023320670701084355-0.02240684127113227*chi44+0.011427087840231389*chi44*chi44)+
										eta*eta2*(-0.46054477257132803 + 0.433526632115367*chi44) + eta2*eta2*(1.2796262150829425-1.2400051122897835*chi44);
										break;
							default:
										print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
										return REAL8_FAIL_NAN;
										break;
						}
						break;

				case 5:
						switch (modeM) {
							case 5:
										res = eta * (dM*(-0.008389798844109389 + 0.04678354680410954*eta) + dM*chi33*(-0.0013605616383929452 + 0.004302712487297126*eta) +dM*chi33*chi33*(-0.0011412109287400596 + 0.0018590391891716925*eta) +
									  0.0002944221308683548*fabs(dM*(37.11125499129578 - 157.79906814398277*eta) + chi33));
										break;
							default:
										print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
										return REAL8_FAIL_NAN;
										break;
						}
						break;
        default:
            print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
            return REAL8_FAIL_NAN;
            break;
    }
    //    printf("ddA %.16e\n", res);
    return res;
}

static REAL8
XLALSimIMREOBGetNRSpinPeakADDotV4 (INT modeL, INT modeM, REAL8 m1, REAL8 m2,
			  REAL8 chiS, REAL8 chiA)
{
	 REAL8 eta = (m1 * m2) / ((m1 + m2) * (m1 + m2));
	 REAL8 dM = sqrt(1.-4.*eta);
	 REAL8 tempm1;
	 if (m1 < m2)
		 {
			 //RC: The fits for the HMs are done under the assumption m1>m2, so if m2>m1 we just swap the two bodies
			 tempm1 = m1;
			 m1 = m2;
			 m2 = tempm1;
			 chiA = -chiA;
		 }
	 REAL8 eta2 = eta*eta;
	 REAL8 chi21 = chiS*dM/(1.-2.*eta) + chiA;
	 REAL8 chi = chiS + chiA * (m1 - m2) / (m1 + m2) / (1. -
					 2. *
					 eta);
		REAL8 chi33 = chiS*dM + chiA;
    REAL8 chiMinus1 = -1. + chi;
    REAL8 res;
    REAL8 fTPL, fEQ, A1, e0, e1;
    switch (modeL) {
        case 2:
            switch (modeM) {
                case 2:
                    // TPL fit
                    fTPL = 0.002395610769995033 * chiMinus1 -  0.00019273850675004356 * chiMinus1 * chiMinus1 - 0.00029666193167435337 * chiMinus1 * chiMinus1 * chiMinus1;
                    // Equal-mass fit
                    fEQ = -0.004126509071377509 + 0.002223999138735809 * chi;
                    // Global fit coefficients
                    e0 = -0.005776537350356959;
                    e1 = 0.001030857482885267;
                    A1 = e0 + e1 * chi;
                    res = eta * CombineTPLEQMFits(eta, A1, fEQ, fTPL);;
                    break;
								case 1:
										res = eta*dM*0.00037132201959950333 -
										fabs(dM*eta*(-0.0003650874948532221 - eta*0.003054168419880019)
										+dM*eta*chi21*chi21*(-0.0006306232037821514-eta*0.000868047918883389 + eta2*0.022306229435339213)+eta*chi21*chi21*chi21*0.0003402427901204342+dM*eta*chi21*0.00028398490492743);
										break;
                default:
                    print_warning("XLAL Error - %s: At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n", __func__);
                    return REAL8_FAIL_NAN;
                    break;
            }
            break;

        case 3:
            switch (modeM) {
                case 3:
                    res = dM*eta*(0.0009605689249339088 - 0.00019080678283595965*eta)*chi33 - 0.00015623760412359145*eta*fabs(dM*(4.676662024170895 + 79.20189790272218*eta - 1097.405480250759*eta2 + 6512.959044311574*eta*eta2 -13263.36920919937*eta2*eta2) + chi33);
                    break;
                default:
                    print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                    return REAL8_FAIL_NAN;
                    break;
            }
            break;
				case 4:
						switch (modeM) {
							  case 4:
										res = eta*(-0.000301722928925693 + 0.0003215952388023551*chi) + eta2*(0.006283048344165004 + 0.0011598784110553046*chi) + eta2*eta*(-0.08143521096050622 - 0.013819464720298994*chi)+
										eta2*eta2*(0.22684871200570564 + 0.03275749240408555*chi);
										break;
								default:
										print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
										return REAL8_FAIL_NAN;
										break;
						}
						break;
				case 5:
						switch (modeM) {
								case 5:
										res = eta * (dM *(0.00012727220842255978 + 0.0003211670856771251*eta) + dM*chi33*(-0.00006621677859895541 + 0.000328855327605536*eta) + chi33*chi33*(-0.00005824622885648688 + 0.00013944293760663706*eta));
										break;
								default:
									  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
										return REAL8_FAIL_NAN;
										break;
								}
								break;

        default:
            print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
            return REAL8_FAIL_NAN;
            break;
    }
    //    printf("ddA %.16e\n", res);
    return res;
}

/**
 * Peak frequency predicted by fitting NR results
 * The coefficients for the (2,2) are in Eq.(A6) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.95.044028
 * The coefficients for the modes (2,1), (3,3), (4,4) and (5,5) are in Eqs.(B19-B22) of https://arxiv.org/pdf/1803.10701.pdf
 */
static REAL8
XLALSimIMREOBGetNRSpinPeakOmegaV4 (INT modeL, INT modeM, REAL8 eta, REAL8 a)
{
  REAL8 chi = a, eta2 = eta*eta;
  REAL8 res;
  REAL8 c0, c1, c2, c3, c4, d2, d3, A3, A4;
  switch (modeL) {
      case 2:
          switch (modeM) {
              case 2:
                  // From TPL fit
                  c0 = 0.5626787200433265;
                  c1 = -0.08706198756945482;
                  c2 = 25.81979479453255;
                  c3 = 25.85037751197443;
                  // From equal-mass fit
                  d2 = 7.629921628648589;
                  d3 = 10.26207326082448;
                  // Combine TPL and equal-mass
                  A4 = d2 + 4 * (d2 - c2) * (eta - 0.25);
                  A3 = d3 + 4 * (d3 - c3) * (eta - 0.25);
                  c4 = 0.00174345193125868;
                  // Final formula
                  res = c0 + (c1 + c4 * chi) * log(A3 - A4 * chi);
                  break;
							case 1:
									res = (0.1743194440996283 + eta*0.1938944514123048 + 0.1670063050527942*eta2 + 0.053508705425291826 *chi - eta*chi*0.18460213023455802 + eta2*chi*0.2187305149636044
									+chi*chi*0.030228846150378793 -  eta*chi*chi*0.11222178038468673);
									break;
              default:
                    print_warning("XLAL Error - %s: At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n", __func__);
                    return REAL8_FAIL_NAN;
                    break;
            }
            break;

      case 3:
          switch (modeM) {
              case 3:
                  res = 0.3973947703114506 + 0.16419332207671075*chi + 0.1635531186118689*chi*chi + 0.06140164491786984*chi*chi*chi+
									eta*(0.6995063984915486-0.3626744855912085*chi -0.9775469868881651*chi*chi)+ eta2*(-0.3455328417046369+0.31952307610699876*chi + 1.9334166149686984*chi*chi);
                  break;
              default:
                  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                  return REAL8_FAIL_NAN;
                  break;
          }
          break;
			case 4:
					switch (modeM) {
						case 4:
									res = 0.5389359134370971 + 0.16635177426821202*chi + 0.2075386047689103*chi*chi + 0.15268115749910835*chi*chi*chi +
									eta*(0.7617423831337586 + 0.009587856087825369*chi - 1.302303785053009*chi*chi - 0.5562751887042064*chi*chi*chi)
									+ eta2*(0.9675153069365782 - 0.22059322127958586*chi + 2.678097398558074*chi*chi) - eta2*eta*4.895381222514275;
									break;
						default:
		              print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
		              return REAL8_FAIL_NAN;
		              break;
					}
					break;
			case 5:
					 switch (modeM) {
						case 5:
									res = 0.6437545281817488 + 0.22315530037902315*chi + 0.2956893357624277*chi*chi + 0.17327819169083758*chi*chi*chi +
									eta*(-0.47017798518175785 - 0.3929010618358481*chi - 2.2653368626130654*chi*chi - 0.5512998466154311*chi*chi*chi) +
									eta2*(2.311483807604238 + 0.8829339243493562*chi + 5.817595866020152*chi*chi);
									break;
						default:
									print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
									return REAL8_FAIL_NAN;
									break;
							}
							break;

        default:
          print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
          return REAL8_FAIL_NAN;
          break;
    }
//    printf("w %.16e\n", res);
  return res;
}

/**
 * Peak frequency slope predicted by fitting NR results
 * The coefficients for the (2,2) are in Eq.(A10) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.95.044028
 * The coefficients for the modes (2,1), (3,3), (4,4) and (5,5) are in Eqs.(B23-B26) of https://arxiv.org/pdf/1803.10701.pdf
 */
static REAL8
XLALSimIMREOBGetNRSpinPeakOmegaDotV4 (INT modeL, INT modeM, REAL8 eta,
			 REAL8 a)
{
  REAL8 chi = a, eta2 = eta*eta;
  REAL8 res;
  REAL8 fTPL, fEQ, A1, e0, e1;
  switch (modeL) {
      case 2:
          switch (modeM) {
              case 2:
                  // TPL fit
                  fTPL = -0.011209791668428353 +  (0.0040867958978563915 + 0.0006333925136134493 * chi) * log(68.47466578100956 - 58.301487557007206 * chi);
                  // Equal-mass fit
                  fEQ = 0.01128156666995859 + 0.0002869276768158971* chi;
                  // Global fit coefficients
                  e0 = 0.01574321112717377;
                  e1 = 0.02244178140869133;
                  A1 = e0 + e1 * chi;
                  res = CombineTPLEQMFits(eta, A1, fEQ, fTPL);;
                  break;
							case 1:
							res = (0.0070987396362959514 + eta*0.024816844694685373 -eta2*0.050428973182277494 + eta*eta2*0.03442040062259341-chi*0.0017751850002442097+eta*chi*0.004244058872768811
							-eta2*chi*0.031996494883796855-chi*chi*0.0035627260615894584+eta*chi*chi*0.01471807973618255 - chi*chi*chi*0.0019020967877681962);
							break;
              default:
                  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                  return REAL8_FAIL_NAN;
                  break;
          }
          break;

      case 3:
          switch (modeM) {
              case 3:
                  res = 0.010337157192240338 - 0.0053067782526697764*chi*chi - 0.005087932726777773*chi*chi*chi+
									eta*(0.027735564986787684 + 0.018864151181629343*chi + 0.021754491131531044*chi*chi + 0.01785477515931398*chi*chi*chi)+
									eta2*(0.018084233854540898 - 0.08204268775495138*chi);
                  break;

              default:
                  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
                  return REAL8_FAIL_NAN;
                  break;
          }
          break;

			case 4:
					switch (modeM) {
						case 4:
								  res = 0.013997911323773867 - 0.0051178205260273574*chi - 0.0073874256262988*chi*chi +
									eta*(0.0528489379269367 + 0.01632304766334543*chi + 0.02539072293029433*chi*chi)
									+eta2*(-0.06529992724396189 + 0.05782894076431308*chi);
						      break;

						  default:
								  print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
								  return REAL8_FAIL_NAN;
								  break;
					}
					break;
			case 5:
					switch (modeM) {
						case 5:
								  res = 0.01763430670755021 - 0.00024925743340389135*chi - 0.009240404217656968*chi*chi - 0.007907831334704586*chi*chi*chi+
									eta*(-0.1366002854361568 + 0.0561378177186783*chi + 0.16406275673019852*chi*chi + 0.07736232247880881*chi*chi*chi)+
									eta2*(0.9875890632901151 - 0.31392112794887855*chi - 0.5926145463423832*chi*chi) - 1.6943356548192614*eta2*eta;
									break;

							default:
									print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
									return REAL8_FAIL_NAN;
									break;
					}
					break;
      default:
          print_warning("XLAL Error - %s: (%d,%d) mode. At present only fits for the (2,2), (2,1), (3,3), (4,4) and (5,5) mode are available.\n",  __func__, modeL, modeM);
          return REAL8_FAIL_NAN;
          break;
  }
//    printf("dw %.16e\n", res);
  return res;
}



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
					     EOBNonQCCoeffs *coeffs)
{
  int signum, i;

  INT length = dyHi->length;
  /* Vectors which are used in the computation of the NQC coefficients */
  REAL8Vector *timeVec;
  REAL8Vector *q3 = NULL, *q4 = NULL, *q5 = NULL;
  REAL8Vector *p3 = NULL, *p4 = NULL;

  REAL8Vector *qNS = NULL, *pNS = NULL;

  /* Since the vectors we actually want are q etc * A, we will have to generate them here */
  REAL8Vector *q3LM = NULL;
  REAL8Vector *q4LM = NULL;
  REAL8Vector *q5LM = NULL;
  REAL8Vector *qNSLM = NULL;

  REAL8 eta = (m1 * m2) / ((m1 + m2) * (m1 + m2));
  REAL8 amp, aDot, aDDot;
  REAL8 omega, omegaDot;
  REAL8 rootR;
  REAL8 rOmega;
  //UNUSED REAL8 qNSLMPeak, qNSLMDot, qNSLMDDot;
  //UNUSED REAL8 pNSLMDot, pNSLMDDot;

  REAL8 nra, nraDot, nraDDot;
  REAL8 nromega, nromegaDot;

  REAL8 nrDeltaT, nrTimePeak;
  REAL8 chi1 = chiS + chiA;
  REAL8 chi2 = chiS - chiA;
  REAL8 q1, q2, p1, p2;
  /* Stuff for finding numerical derivatives */
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Matrix stuff for calculating coefficients */
  gsl_matrix *qMatrix = NULL, *pMatrix = NULL;
  gsl_vector *aCoeff = NULL, *bCoeff = NULL;

  gsl_vector *amps = NULL, *omegaVec = NULL;

  gsl_permutation *perm1 = NULL, *perm2 = NULL;

  memset (coeffs, 0, sizeof (EOBNonQCCoeffs));

  /* Populate the time vector */
  /* It is okay to assume initial t = 0 */
  timeVec = CreateREAL8Vector(length);
  q3 = CreateREAL8Vector (length);
  q4 = CreateREAL8Vector (length);
  q5 = CreateREAL8Vector (length);
  p3 = CreateREAL8Vector (length);
  p4 = CreateREAL8Vector (length);
  qNS = CreateREAL8Vector (length);
  pNS = CreateREAL8Vector (length);
  q3LM = CreateREAL8Vector (length);
  q4LM = CreateREAL8Vector (length);
  q5LM = CreateREAL8Vector (length);
  qNSLM = CreateREAL8Vector (length);

  if (!q3 || !q4 || !q5 || !p3 || !p4 || !qNS || !pNS || !q3LM
      || !q4LM || !q5LM || !qNSLM)
    {
      DestroyREAL8Vector (timeVec);
      DestroyREAL8Vector (q3);
      DestroyREAL8Vector (q4);
      DestroyREAL8Vector (q5);
      DestroyREAL8Vector (p3);
      DestroyREAL8Vector (p4);
      DestroyREAL8Vector (qNS);
      DestroyREAL8Vector (pNS);
      DestroyREAL8Vector (q3LM);
      DestroyREAL8Vector (q4LM);
      DestroyREAL8Vector (q5LM);
      DestroyREAL8Vector (qNSLM);
      return CEV_FAILURE;
    }

  /* Populate vectors as necessary. Eqs. 14 - 17 of the LIGO DCC document T1100433v2 */
//        FILE *out = fopen( "out.dat","w");
  for (i = 0; i < length; i++)
    {
      timeVec->data[i] = i*deltaT;
      rootR = sqrt (dyHi->rVec->data[i]);
      rOmega = dyHi->rVec->data[i] * orbOmegaVec->data[i];

      /* We don't need these as vectors as their coefficients are calibrated */

      q1 = dyHi->prVec->data[i] * dyHi->prVec->data[i] / (rOmega * rOmega);
      q2 = q1 / dyHi->rVec->data[i];
      q3->data[i] = q1;
      q4->data[i] = q2;
      q5->data[i] = q2 / rootR;

      p1 = dyHi->prVec->data[i] / rOmega;
      p2 = p1 * dyHi->prVec->data[i] * dyHi->prVec->data[i];
      p3->data[i] = p1;
      p4->data[i] = p2;

      qNS->data[i] = 0.;
      pNS->data[i] = 0.;
      q3LM->data[i] = q3->data[i] * amplitude->data[i];
      q4LM->data[i] = q4->data[i] * amplitude->data[i];
      q5LM->data[i] = q5->data[i] * amplitude->data[i];


      qNSLM->data[i] = 0.;
//            fprintf(out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", amplitude->data[i], prVec->data[i],rVec->data[i],orbOmegaVec->data[i], q3->data[i],q4->data[i],q5->data[i], p3->data[i], p4->data[i],phase->data[i]);
    }
//        fclose(out);

  /* Allocate all the memory we need */
		 /* a stuff */
		 qMatrix = gsl_matrix_alloc (3, 3);
		 aCoeff = gsl_vector_alloc (3);
		 amps = gsl_vector_alloc (3);
		 perm1 = gsl_permutation_alloc (3);
		 /* b stuff */
		 pMatrix = gsl_matrix_alloc (2, 2);
		 bCoeff = gsl_vector_alloc (2);
		 omegaVec = gsl_vector_alloc (2);
		 perm2 = gsl_permutation_alloc (2);

  if (!qMatrix || !aCoeff || !amps || !pMatrix || !bCoeff || !omegaVec)
    {
      DestroyREAL8Vector (timeVec);
      DestroyREAL8Vector (q3);
      DestroyREAL8Vector (q4);
      DestroyREAL8Vector (q5);
      DestroyREAL8Vector (p3);
      DestroyREAL8Vector (p4);
      DestroyREAL8Vector (qNS);
      DestroyREAL8Vector (pNS);
      DestroyREAL8Vector (q3LM);
      DestroyREAL8Vector (q4LM);
      DestroyREAL8Vector (q5LM);
      DestroyREAL8Vector (qNSLM);
      return CEV_FAILURE;
    }

  /* The time we want to take as the peak time depends on l and m */
  /* Calculate the adjustment we need to make here */
    nrDeltaT = XLALSimIMREOBGetNRSpinPeakDeltaTv4 (modeL, modeM, m1, m2, chi1, chi2);


  if (IS_REAL8_FAIL_NAN (nrDeltaT))
    {
      DestroyREAL8Vector (timeVec);
      DestroyREAL8Vector (q3);
      DestroyREAL8Vector (q4);
      DestroyREAL8Vector (q5);
      DestroyREAL8Vector (p3);
      DestroyREAL8Vector (p4);
      DestroyREAL8Vector (qNS);
      DestroyREAL8Vector (pNS);
      DestroyREAL8Vector (q3LM);
      DestroyREAL8Vector (q4LM);
      DestroyREAL8Vector (q5LM);
      DestroyREAL8Vector (qNSLM);
      return CEV_FAILURE;
    }

  /* nrDeltaT defined in XLALSimIMREOBGetNRSpinPeakDeltaT is a minus sign different from Eq. (33) of Taracchini et al.
   * Therefore, the plus sign in Eq. (21) of Taracchini et al and Eq. (18) of DCC document T1100433v2 is
   * changed to a minus sign here.
   */
  nrTimePeak = timePeak - nrDeltaT;
  /*if (debugAT)
    printf ("nrTimePeak, timePeak %.16e %.16e\n", nrTimePeak, timePeak);*/
  /* We are now in a position to use the interp stuff to calculate the derivatives we need */
  /* We will start with the quantities used in the calculation of the a coefficients */
  spline = gsl_spline_alloc (gsl_interp_cspline, amplitude->length);
  acc = gsl_interp_accel_alloc ();

  /* Populate the Q matrix in Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Q3 */
  gsl_spline_init (spline, timeVec->data, q3LM->data, q3LM->length);
  gsl_matrix_set (qMatrix, 0, 0, gsl_spline_eval (spline, nrTimePeak, acc));
  gsl_matrix_set (qMatrix, 1, 0,
		  gsl_spline_eval_deriv (spline, nrTimePeak, acc));
  gsl_matrix_set (qMatrix, 2, 0,
		  gsl_spline_eval_deriv2 (spline, nrTimePeak, acc));

  /* Q4 */
  gsl_spline_init (spline, timeVec->data, q4LM->data, q4LM->length);
  gsl_interp_accel_reset (acc);
  gsl_matrix_set (qMatrix, 0, 1, gsl_spline_eval (spline, nrTimePeak, acc));
  gsl_matrix_set (qMatrix, 1, 1,
		  gsl_spline_eval_deriv (spline, nrTimePeak, acc));
  gsl_matrix_set (qMatrix, 2, 1,
		  gsl_spline_eval_deriv2 (spline, nrTimePeak, acc));

  /* Q5 */
  gsl_spline_init (spline, timeVec->data, q5LM->data, q5LM->length);
  gsl_interp_accel_reset (acc);
  gsl_matrix_set (qMatrix, 0, 2, gsl_spline_eval (spline, nrTimePeak, acc));
  gsl_matrix_set (qMatrix, 1, 2,
		  gsl_spline_eval_deriv (spline, nrTimePeak, acc));
  gsl_matrix_set (qMatrix, 2, 2,
		  gsl_spline_eval_deriv2 (spline, nrTimePeak, acc));

  /* Populate the r.h.s vector of Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Amplitude */
  gsl_spline_init (spline, timeVec->data, amplitude->data, amplitude->length);
  gsl_interp_accel_reset (acc);
  amp = gsl_spline_eval (spline, nrTimePeak, acc);
  aDot = gsl_spline_eval_deriv (spline, nrTimePeak, acc);
  aDDot = gsl_spline_eval_deriv2 (spline, nrTimePeak, acc);

  /* qNSLM */
  gsl_spline_init (spline, timeVec->data, qNSLM->data, qNSLM->length);
  gsl_interp_accel_reset (acc);
  //qNSLMPeak = gsl_spline_eval (spline, nrTimePeak, acc);
  //qNSLMDot = gsl_spline_eval_deriv (spline, nrTimePeak, acc);
  //qNSLMDDot = gsl_spline_eval_deriv2 (spline, nrTimePeak, acc);

  nra =
    fabs(XLALSimIMREOBGetNRSpinPeakAmplitudeV4 (modeL, modeM, m1,m2,chiS,chiA));

  nraDot =
    XLALSimIMREOBGetNRSpinPeakADotV4 (modeL, modeM, m1,m2,chiS,chiA);
	//RC: In SEOBNRv4 nraDot is zero because the NQC are defining the peak of the 22 mode
	// which by definition has a first derivative	equal to 0.
	//For SEOBNRv4HM we are not computing the NQC at the peak fo the modes (see Eq.(4.3))
	//of https://arxiv.org/pdf/1803.10701.pdf, so first the derivative
	//is entering as a fitted coefficient.

  nraDDot =
    XLALSimIMREOBGetNRSpinPeakADDotV4 (modeL, modeM, m1,m2,chiS,chiA);
//    printf("eta, chiS, chiA, dM/M, chi = %.16e %.16e %.16e %.16e %.16e\n",eta,chiS,chiA, (m1 - m2)/(m1 + m2),chiS + chiA*(m1 - m2)/(m1 + m2)/(1. - 2.*eta));
  if (IS_REAL8_FAIL_NAN (nra) || IS_REAL8_FAIL_NAN (nraDot) || IS_REAL8_FAIL_NAN (nraDDot))
    {
      DestroyREAL8Vector (timeVec);
      DestroyREAL8Vector (q3);
      DestroyREAL8Vector (q4);
      DestroyREAL8Vector (q5);
      DestroyREAL8Vector (p3);
      DestroyREAL8Vector (p4);
      DestroyREAL8Vector (qNS);
      DestroyREAL8Vector (pNS);
      DestroyREAL8Vector (q3LM);
      DestroyREAL8Vector (q4LM);
      DestroyREAL8Vector (q5LM);
      DestroyREAL8Vector (qNSLM);
      return CEV_FAILURE;
    }

  gsl_vector_set (amps, 0, nra - amp);
  gsl_vector_set (amps, 1, nraDot -aDot);
  gsl_vector_set (amps, 2, nraDDot - aDDot);
//    printf("Amps %.16e %.16e %.16e\n", nra, amp, qNSLMPeak);
//    printf("dAmps %.16e %.16e\n", aDot, qNSLMDot);
//    printf("ddAmps %.16e %.16e %.16e\n", nraDDot, aDDot, qNSLMDDot);


  /* We have now set up all the stuff to calculate the a coefficients */
  /* So let us do it! */
  gsl_linalg_LU_decomp (qMatrix, perm1, &signum);
  gsl_linalg_LU_solve (qMatrix, perm1, amps, aCoeff);

/*
  if (debugAT)
    {
      printf ("Q MATRIX\n");
      for (unsigned int i = 0; i < 3; i++)
	{
	  for (unsigned int j = 0; j < 3; j++)
	    {
	      printf ("%.12e\t", gsl_matrix_get (qMatrix, i, j));
	    }
	  printf ("= %.12e\n", gsl_vector_get (amps, i));
	}
    }*/



  /* Now we (should) have calculated the a values. Now we can do the b values */

  /* Populate the P matrix in Eq. 18 of the LIGO DCC document T1100433v2 */
  /* P3 */
  gsl_spline_init (spline, timeVec->data, p3->data, p3->length);
  gsl_interp_accel_reset (acc);
  gsl_matrix_set (pMatrix, 0, 0,
		  -gsl_spline_eval_deriv (spline, nrTimePeak, acc));
  gsl_matrix_set (pMatrix, 1, 0,
		  -gsl_spline_eval_deriv2 (spline, nrTimePeak, acc));

  /* P4 */
  gsl_spline_init (spline, timeVec->data, p4->data, p4->length);
  gsl_interp_accel_reset (acc);
  gsl_matrix_set (pMatrix, 0, 1,
		  -gsl_spline_eval_deriv (spline, nrTimePeak, acc));
  gsl_matrix_set (pMatrix, 1, 1,
		  -gsl_spline_eval_deriv2 (spline, nrTimePeak, acc));

  /* Populate the r.h.s vector of Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Phase */
  gsl_spline_init (spline, timeVec->data, phase->data, phase->length);
  gsl_interp_accel_reset (acc);
  omega = gsl_spline_eval_deriv (spline, nrTimePeak, acc);
  omegaDot = gsl_spline_eval_deriv2 (spline, nrTimePeak, acc);

  /* pNSLM */
  gsl_spline_init (spline, timeVec->data, pNS->data, pNS->length);
  gsl_interp_accel_reset (acc);
  // = gsl_spline_eval_deriv (spline, nrTimePeak, acc);
  //pNSLMDDot = gsl_spline_eval_deriv2 (spline, nrTimePeak, acc);

  /* Since the phase can be decreasing, we need to take care not to have a -ve frequency */
  if (omega * omegaDot > 0.0)
    {
      omega = fabs (omega);
      omegaDot = fabs (omegaDot);
    }
  else
    {
      omega = fabs (omega);
      omegaDot = -fabs (omegaDot);
    }

  nromega =
    XLALSimIMREOBGetNRSpinPeakOmegaV4 (modeL, modeM, eta,
			  chiS + chiA * (m1 - m2) / (m1 + m2) / (1. -
								 2. * eta));
  nromegaDot =
    XLALSimIMREOBGetNRSpinPeakOmegaDotV4 (modeL, modeM, eta,
			     chiS + chiA * (m1 - m2) / (m1 + m2) / (1. -
								    2. *
								    eta));

/*
  if (debugAT)
    printf ("NR inputs: %.16e, %.16e, %.16e, %.16e\n", nra, nraDDot, nromega,
	    nromegaDot);*/

/*     printf("NR inputs: %.16e, %.16e, %.16e, %.16e\n",pNSLMDot, pNSLMDDot,omega,omegaDot);*/

  if (IS_REAL8_FAIL_NAN (nromega) || IS_REAL8_FAIL_NAN (nromegaDot))
    {
      DestroyREAL8Vector (timeVec);
      DestroyREAL8Vector (q3);
      DestroyREAL8Vector (q4);
      DestroyREAL8Vector (q5);
      DestroyREAL8Vector (p3);
      DestroyREAL8Vector (p4);
      DestroyREAL8Vector (qNS);
      DestroyREAL8Vector (pNS);
      DestroyREAL8Vector (q3LM);
      DestroyREAL8Vector (q4LM);
      DestroyREAL8Vector (q5LM);
      DestroyREAL8Vector (qNSLM);
      return CEV_FAILURE;
    }

  gsl_vector_set (omegaVec, 0, nromega - omega);
  gsl_vector_set (omegaVec, 1, nromegaDot - omegaDot);

  /*
  if (debugAT)
    {
      printf ("P MATRIX\n");
      for (unsigned int i = 0; i < 2; i++)
	{
	  for (unsigned int j = 0; j < 2; j++)
	    {
	      printf ("%.12e\t", gsl_matrix_get (pMatrix, i, j));
	    }
	  printf ("= %.12e\n", gsl_vector_get (omegaVec, i));
	}
    }*/

  /* And now solve for the b coefficients */
  gsl_linalg_LU_decomp (pMatrix, perm2, &signum);
  gsl_linalg_LU_solve (pMatrix, perm2, omegaVec, bCoeff);

  /* We can now populate the coefficients structure */
  coeffs->a3S = 0.;
  coeffs->a4 = 0.;
  coeffs->a5 = 0.;
  coeffs->b3 = 0.;
  coeffs->b4 = 0.;


  coeffs->a1 = gsl_vector_get (aCoeff, 0);
  coeffs->a2 = gsl_vector_get (aCoeff, 1);
  coeffs->a3 = gsl_vector_get (aCoeff, 2);

  coeffs->b1 = gsl_vector_get (bCoeff, 0);
  coeffs->b2 = gsl_vector_get (bCoeff, 1);



//    printf( "NQC coefficients:\n" );
//    printf( "{%f,%f,%f,%f,%f,%f}\n",  coeffs->a1, coeffs->a2, coeffs->a3, coeffs->a3S, coeffs->a4, coeffs->a5 );
//
//    printf( "{%f,%f,%f,%f}\n",  coeffs->b1, coeffs->b2, coeffs->b3, coeffs->b4 );


  /* Free memory and exit */
  gsl_matrix_free (qMatrix);
  gsl_vector_free (amps);
  gsl_vector_free (aCoeff);
  gsl_permutation_free (perm1);

  gsl_matrix_free (pMatrix);
  gsl_vector_free (omegaVec);
  gsl_vector_free (bCoeff);
  gsl_permutation_free (perm2);

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  DestroyREAL8Vector (timeVec);
  DestroyREAL8Vector (q3);
  DestroyREAL8Vector (q4);
  DestroyREAL8Vector (q5);
  DestroyREAL8Vector (p3);
  DestroyREAL8Vector (p4);
  DestroyREAL8Vector (qNS);
  DestroyREAL8Vector (pNS);
  DestroyREAL8Vector (q3LM);
  DestroyREAL8Vector (q4LM);
  DestroyREAL8Vector (q5LM);
  DestroyREAL8Vector (qNSLM);


  return CEV_SUCCESS;
}