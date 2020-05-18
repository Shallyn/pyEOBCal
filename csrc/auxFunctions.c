/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "auxFunctions.h"
#include "dyFactorizedWaveform.h"
#include "dyHamiltonian.h"
#include "pnWaveform.h"
#include "dyFactorizedFlux.h"
#include <gsl/gsl_deriv.h>

int
XLALSimIMREOBNonQCCorrection (COMPLEX16 *nqc,	/**<< OUTPUT, The NQC correction */
			        REAL8Vector *values,
							/**<< Dynamics r, phi, pr, pphi */
			        const REAL8 omega,	/**<< Angular frequency */
			        EOBNonQCCoeffs *coeffs
							/**<< NQC coefficients */
  );

REAL8
FactorizedAngularMomentumFlux (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
				const REAL8 omega,	/**< orbital frequency */
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				INT const  lMax);


#define STEP_SIZE 1.0e-3
#define Power(a,b) (a##b)
INT auxCalculateCompareFlux(const REAL8 values[],
                            SpinEOBParams *seobParams,
                            REAL8 *Eflux,
                            REAL8 *FphiEOB,
                            REAL8 *EfluxPN,
                            REAL8 *LfluxPN,
                            REAL8 *EfluxNagar,
                            REAL8 *LfluxNagar,
                            REAL8 *EfluxFac,
                            REAL8 *LfluxFac)
{
    HcapDerivParams params;

    REAL8 eta, Mtotal;
    Mtotal = seobParams->eobParams->Mtotal;
    eta   = seobParams->eobParams->eta;

    REAL8       cartValues[6];
    REAL8Vector polarDynamics;
    REAL8       polData[4];
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    memcpy( polData, values, sizeof( polData ) );

    REAL8 r, pr, prt, pphi;
    r = values[0];
    prt = values[2];
    pphi = values[3];

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;

    memset( rData, 0, sizeof(rData) );
    memset( pData, 0, sizeof(pData) );

    rData[0] = r;
    pData[0] = prt;
    pData[1] = pphi / r;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;
    params.values  = cartValues;
    params.params  = seobParams;
    s1Vec = params.params->s1VecOverMtMt;
    s2Vec = params.params->s2VecOverMtMt;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;
    REAL8  a  = sKerr->data[2];

    REAL8 H, omega, flux;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;
    nqcCoeffs = params.params->nqcCoeffs;

    gsl_function F;
    INT4         gslStatus;
    UINT i;
    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params   = &params;
    REAL8 tmpDValues[6];
    REAL8 absErr;

    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = r;
    cartValues[3] = prt;
    cartValues[4] = pphi / r;

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

    REAL8 csi, DeltaT, DeltaR;
    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );

    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );

    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    pr = prt / csi;


    H =  SpinEOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, seobParams->tortoise, seobParams->seobCoeffs );
    omega = tmpDValues[4] / r;
    flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omega, seobParams, H, 8 ,0);
        //print_debug("H = %.4f, flux = %.4f\n", flux);
    *Eflux = flux;
    *FphiEOB = flux / omega;

    // Calculate 2PN flux
    REAL8 r2, r3, r4, r5, r6, r7, r8;
    REAL8 eta2, eta3, eta4;
    REAL8 p2, p4, p6, p8, p10, pr2, pr4, pr6, pr8, phidot;
    REAL8 tmpEflux, tmpLflux;
    REAL8 FN, F1PN, F2PN, GN, G1PN, G2PN;
    r2 = r*r;
    r3 = r2*r;
    r4 = r3*r;
    r5 = r4*r;
    r6 = r3*r3;
    r7 = r6*r;
    r8 = r7*r;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    eta4 = eta3*eta;
    pr2 = pr*pr;
    pr4 = pr2*pr2;
    pr6 = pr4*pr2;
    pr8 = pr6*pr2;
    p2 = pr2 + pphi*pphi/r2;
    p4 = p2*p2;
    p6 = p4*p2;
    p8 = p6*p2;
    p10 = p8*p2;
    phidot = pphi / r2;
    FN = (8*(12*p2 - 11*pr2)*eta2)/(15.*r4);
    F1PN = (2*eta2*(16 + 2880*pr2*r + 449*p4*r2 + 
       2061*pr4*r2 - 2*p2*r*(1360 + 1333*pr2*r) + 
       14*eta2*(12*p2 - 11*pr2)*
        (2 + 3*pr2*r + p4*r2 + 8*pr4*r2 - 
          3*p2*r*(1 + 4*pr2*r)) + 
       2*eta*(-32 + pr2*(2926 - 256*r) + 84*p6*r2 - 
          pr4*r*(4235 + 314*r) + 
          p4*r*(84 + (-426 + 259*pr2)*r) + 
          p2*(-3192 + (248 + 4543*pr2)*r + 
             (720*pr2 - 308*pr4)*r2))))/(105.*r6);
    F2PN = ((eta2)*(12*eta*(-8208 + 36*(288 + 29114*(p2) - 31077*(pr2))*r - 
          6*(30153*(p4) + 6872*(pr2) - 155268*(pr4) + 
             (p2)*(-12308 + 95715*(pr2)))*(r2) + 
          (-23535*(p6) + 3*(p4)*(49888 + 39237*(pr2)) - 
             9*(p2)*(pr2)*(53704 + 131345*(pr2)) + 
             (pr4)*(315544 + 1020195*(pr2)))*(r3) + 
          3*(1347*(p8) + 31472*(pr6) + 
             2748*(p2)*(pr4)*(-26 + 9*(pr2)) - 
             10*(p6)*(436 + 261*(pr2)) + 
             (p4)*(45312*(pr2) - 25809*(pr4)))*(r4)) + 
       945*(eta4)*(12*(p2) - 11*(pr2))*
        pow(2 + 3*(pr2)*r + (p4)*(r2) + 
          8*(pr4)*(r2) - 3*(p2)*r*(1 + 4*(pr2)*r),2) - 
       8*r*(3480 + 243789*(pr2)*r + 703332*(pr4)*(r2) + 
          2862*(p6)*(r3) + 45018*(pr6)*(r3) - 
          162*(p4)*(r2)*(-1058 + 367*(pr2)*r) + 
          (p2)*r*(-198961 - 875988*(pr2)*r + 
             2322*(pr4)*(r2))) + 
       54*(eta3)*(-384 + 4*(pr4)*(385 - 1748*r)*r + 
          (pr6)*(53515 - 13056*r)*(r2) + 420*(p10)*(r4) - 
          140*(pr2)*(-209 + 24*r) - 
          40*(pr8)*(r3)*(4235 + 128*r) - 
          5*(p8)*(r3)*(168 + (408 + 749*(pr2))*r) - 
          2*(p6)*(r2)*
           (8190 - (3764 + 7525*(pr2))*r + 
             2*(pr2)*(-6749 + 3430*(pr2))*(r2)) + 
          (p4)*r*(48720 + 3*(-2832 + 47285*(pr2))*r - 
             2*(pr2)*(15978 + 140945*(pr2))*(r2) + 
             56*(pr4)*(-842 + 515*(pr2))*(r3)) - 
          2*(p2)*(15960 + 2*(-848 + 11585*(pr2))*r + 
             (-7868*(pr2) + 87325*(pr4))*(r2) - 
             98*(pr4)*(191 + 2200*(pr2))*(r3) + 
             176*(pr6)*(-79 + 35*(pr2))*(r4))) + 
       3*(eta2)*(64*(2079 - 272*r) + 3780*(p10)*(r4) + 
          593568*(pr8)*(r4) + 
          9*(p8)*(r3)*(840 + (-2284 + 2975*(pr2))*r) - 
          12*(pr2)*(416955 - 83556*r + 3248*(r2)) + 
          4*(pr4)*r*(3620925 - 106992*r + 9016*(r2)) - 
          3*(pr6)*(r2)*(3493875 - 180372*r + 29824*(r2)) + 
          6*(p6)*(r2)*
           (-47250 + 3*(-9454 + 24395*(pr2))*r + 
             4*(3244 - 16314*(pr2) + 1365*(pr4))*(r2)) - 
          3*(p4)*r*(95760 + 27*(-21992 + 5845*(pr2))*r + 
             (32384 + 39636*(pr2) - 418110*(pr4))*(r2) + 
             12*(pr2)*(5564 - 42357*(pr2) + 1540*(pr4))*(r3)
             ) - 2*(p2)*(-2729160 + 18*(33048 + 431585*(pr2))*r + 
             (-18352 + 628092*(pr2) - 6053355*(pr4))*(r2) + 
             18*(pr2)*(-1800 + 12279*(pr2) + 42350*(pr4))*
              (r3) + 72*(pr4)*(-1516 + 11835*(pr2))*(r4)))))
    /(11340.*(r8));
    tmpEflux = FN + F1PN + F2PN;
    *EfluxPN = tmpEflux;
    //GN = (8*(2 + 2*p2*r - 3*pr2*r)*eta2*omega)/(5.*r2);
    GN = (8*Power(eta,2)*phidot*(2 + 2*Power(p,2)*r - 3*Power(pr,2)*r))/(5.*Power(r,2));
    G1PN = (Power(eta,2)*phidot*(eta*(-6384 + (344 - 4620*Power(p,2) + 16422*Power(pr,2))*r + 
          (294*Power(p,4) + 5*Power(pr,2)*(40 - 2079*Power(pr,2)) + 
             Power(p,2)*(-1096 + 7413*Power(pr,2)))*Power(r,2) + 
          (126*Power(p,6) - 648*Power(pr,4) + Power(p,4)*(-424 + 315*Power(pr,2)) + 
             Power(p,2)*(1308*Power(pr,2) - 756*Power(pr,4)))*Power(r,3)) + 
       2*r*(-1322 + 1920*Power(pr,2)*r + 55*Power(p,4)*Power(r,2) + 
          285*Power(pr,4)*Power(r,2) - 2*Power(p,2)*r*(494 + 33*Power(pr,2)*r)) + 
       21*Power(eta,2)*(16 + 6*Power(pr,2)*r + 37*Power(pr,4)*Power(r,2) + 
          6*Power(p,6)*Power(r,3) - 72*Power(pr,6)*Power(r,3) - 
          Power(p,4)*Power(r,2)*(10 + 81*Power(pr,2)*r) + 
          3*Power(p,2)*r*(-4 - 17*Power(pr,2)*r + 52*Power(pr,4)*Power(r,2)))))/
   (105.*Power(r,4));
    G2PN = (Power(eta,2)*phidot*(9*eta*r*(696084 + 
          2*(32312 + 149073*Power(p,2) - 831381*Power(pr,2))*r + 
          (-45456*Power(p,4) - 24*Power(p,2)*(-9591 + 20479*Power(pr,2)) + 
             10*Power(pr,2)*(-30212 + 85095*Power(pr,2)))*Power(r,2) + 
          (-7833*Power(p,6) + Power(p,4)*(21424 + 10563*Power(pr,2)) + 
             3*Power(p,2)*Power(pr,2)*(-62980 + 12957*Power(pr,2)) + 
             Power(pr,4)*(179092 + 141075*Power(pr,2)))*Power(r,3) + 
          (495*Power(p,8) + 70760*Power(pr,6) + 
             3*Power(p,4)*Power(pr,2)*(11168 + 63*Power(pr,2)) + 
             2*Power(p,6)*(-2858 + 693*Power(pr,2)) + 
             12*Power(p,2)*Power(pr,4)*(-8729 + 855*Power(pr,2)))*Power(r,4)) + 
       54*Power(eta,3)*(-15960 + 4*(257 + 3696*Power(p,2) + 3381*Power(pr,2))*r + 
          (6426*Power(p,4) + 514*Power(pr,2) - 28434*Power(pr,4) + 
             Power(p,2)*(-3734 + 40992*Power(pr,2)))*Power(r,2) + 
          (-5334*Power(p,6) - 6*Power(p,2)*Power(pr,2)*(1013 + 34832*Power(pr,2)) + 
             Power(p,4)*(2914 + 52563*Power(pr,2)) + 
             Power(pr,4)*(2354 + 118671*Power(pr,2)))*Power(r,3) + 
          (-42*Power(p,8) + Power(p,6)*(236 + 2982*Power(pr,2)) - 
             4*Power(pr,6)*(1109 + 20790*Power(pr,2)) + 
             2*Power(p,2)*Power(pr,4)*(239 + 91560*Power(pr,2)) + 
             Power(p,4)*(8552*Power(pr,2) - 95466*Power(pr,4)))*Power(r,4) + 
          3*(42*Power(p,10) - 288*Power(pr,8) - 
             288*Power(p,2)*Power(pr,6)*(-10 + 7*Power(pr,2)) - 
             6*Power(p,6)*Power(pr,2)*(-347 + 196*Power(pr,2)) - 
             Power(p,8)*(148 + 399*Power(pr,2)) + 
             4*Power(p,4)*Power(pr,4)*(-1223 + 966*Power(pr,2)))*Power(r,5)) + 
       9*Power(eta,2)*(909720 + 36*(-4273 + 13832*Power(p,2) - 95893*Power(pr,2))*
           r + (1376 - 75978*Power(p,4) + 266346*Power(pr,2) + 
             4275810*Power(pr,4) - 42*Power(p,2)*(-6935 + 39504*Power(pr,2)))*
           Power(r,2) + (-27090*Power(p,6) + 
             3*Power(p,4)*(35132 + 14511*Power(pr,2)) + 
             4*Power(p,2)*(-8140 - 104955*Power(pr,2) + 382662*Power(pr,4)) - 
             3*Power(pr,2)*(-5352 + 78226*Power(pr,2) + 571725*Power(pr,4)))*
           Power(r,3) + (1386*Power(p,8) + 
             9*Power(p,6)*(-2117 + 5390*Power(pr,2)) + 
             Power(pr,4)*(41476 + 110007*Power(pr,2)) + 
             Power(p,4)*(33728 - 49953*Power(pr,2) + 109494*Power(pr,4)) + 
             Power(p,2)*(-51308*Power(pr,2) + 6927*Power(pr,4) - 249480*Power(pr,6))
             )*Power(r,4) + (378*Power(p,10) + 
             9*Power(p,8)*(-241 + 273*Power(pr,2)) + 
             8*Power(pr,6)*(-1606 + 2565*Power(pr,2)) + 
             6*Power(p,6)*(248 - 1947*Power(pr,2) + 252*Power(pr,4)) - 
             3*Power(p,4)*Power(pr,2)*
              (3520 - 11679*Power(pr,2) + 3024*Power(pr,4)) + 
             Power(p,2)*(16224*Power(pr,4) - 38124*Power(pr,6)))*Power(r,5)) + 
       567*Power(eta,4)*(10 + 6*Power(p,2)*r - 9*Power(pr,2)*r)*
        pow(2 + 3*Power(pr,2)*r + Power(p,4)*Power(r,2) + 
          8*Power(pr,4)*Power(r,2) - 3*Power(p,2)*r*(1 + 4*Power(pr,2)*r),2) - 
       4*Power(r,2)*(450*Power(p,6)*Power(r,3) + 
          27*Power(p,4)*Power(r,2)*(2019 + 971*Power(pr,2)*r) - 
          9*Power(p,2)*r*(-8624 + 26869*Power(pr,2)*r + 
             6450*Power(pr,4)*Power(r,2)) + 
          4*(-46085 + 26973*Power(pr,2)*r + 94590*Power(pr,4)*Power(r,2) + 
             12285*Power(pr,6)*Power(r,3)))))/(11340.*Power(r,6));
    tmpLflux = GN + G1PN + G2PN;
    *LfluxPN = tmpLflux;
        // print_debug("Eflux = %.4e, Lflux = %.4e, p2 = %.4e, pr2 = %.2e, r = %.2e, eta2 = %.2e, omega = %.2e\n", 
        //     tmpEflux, tmpLflux, p2, pr2, r, eta2, omega);


    // Calculate Nagar flux
    REAL8 x2, x2_2, x2_3, x2_4, x3, x3_2, x3_3, x3_4; 
    REAL8 x4, x4_2, x4_3, x4_4;
    REAL8 Fr_N, Fr_1PN, Fr_2PN, Fphi_N, Fphi_1PN, Fphi_2PN;
    params.varyParam = 0;
    gslStatus = gsl_deriv_central( &F, cartValues[0], 
                    STEP_SIZE, &x4, &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
        print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
        return CEV_FAILURE;
    }

    
    x2 = pr2;
    x2_2 = pr4;
    x2_3 = pr6;
    x2_4 = pr8;

    x3 = 1./r;
    x3_2 = x3*x3;
    x3_3 = x3_2*x3;
    x3_4 = x3_3*x3;

    x4 *= r;
    x4_2 = x4*x4;
    x4_3 = x4_2*x4;
    x4_4 = x4_3*x4;
    /* Fr N */
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
    REAL8 FrNagar, FphiNagar;
    FrNagar = (Fr_N + Fr_1PN + Fr_2PN) * pr * x3_3;
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
    FphiNagar = (Fphi_N + Fphi_1PN + Fphi_2PN) * pphi * x3_3;

    *LfluxNagar = FphiNagar;
    *EfluxNagar = tmpDValues[0] * FrNagar + omega * FphiNagar;

    // Calculate Factorized flux

    *EfluxFac = 0;
    REAL8 LfluxNew;
    LfluxNew = FactorizedAngularMomentumFlux(&polarDynamics, nqcCoeffs, omega, seobParams, H, 8);
    *LfluxFac = LfluxNew;

    return CEV_SUCCESS;
}

INT auxSpinEOBRadiationReactionForceAll(const REAL8 values[],
                                     SpinEOBParams *seobParams,
                                     REAL8 *Fr,
                                     REAL8 *Fphi,
                                     REAL8 *fluxout,
                                     REAL8 *fluxCircout,
                                     REAL8 *Fr2PN,
                                     REAL8 *Fphi2PN,
                                     REAL8 *FrFac,
                                     REAL8 *FphiFac,
                                     REAL8 *FrResum,
                                     REAL8 *FphiResum)
{
    //print_debug("values = (%.2e, %.2e, %.2e, %.2e)\n", values[0], values[1], values[2], values[3]);
    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[6];

    /* Cartesian values for calculating the Hamiltonian */
    REAL8           cartValues[6];

    REAL8           H; //Hamiltonian
    REAL8           flux, fluxE;

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r, prt, pr, pphi;
    REAL8Vector polarDynamics;
    REAL8       polData[4];

    REAL8 Mtotal, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;

    REAL8 omega;

    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = seobParams;
    nqcCoeffs = params.params->nqcCoeffs;

    s1Vec = params.params->s1VecOverMtMt;
    s2Vec = params.params->s2VecOverMtMt;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;

    gsl_function F;
    INT4         gslStatus;
    UINT i;


    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params   = &params;

    Mtotal = params.params->eobParams->Mtotal;
    eta   = params.params->eobParams->eta;

    r = values[0];
    prt = values[2];
    pphi = values[3];
    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );

    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );

    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    pr = prt / csi;
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    //print_debug( "csi in derivatives function = %.16e\n", csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    //print_debug("cartValues = (%.2e, %.2e, %.2e, %.2e, %.2e, %.2e)\n", cartValues[0], cartValues[1], cartValues[2], cartValues[3], cartValues[4], cartValues[5]);
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
    //print_debug("tmpDValues = (%.2e, %.2e, %.2e, %.2e, %.2e, %.2e)\n", tmpDValues[0], tmpDValues[1], tmpDValues[2], tmpDValues[3], tmpDValues[4], tmpDValues[5]);


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
    //print_debug("p = %e, %e, %e\n", pData[0], pData[1], pData[2]);
    /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */

    H =  SpinEOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, seobParams->tortoise, seobParams->seobCoeffs );

    //print_debug( "csi = %.16e, ham = %.16e ( tortoise = %d)\n", csi, H, params.params->tortoise );
    //return CEV_FAILURE;
    //exit(1);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );

    /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
        tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/
    /* Now calculate omega, and hence the flux */

    omega = tmpDValues[4] / r;
    //dvalues[0] = csi * tmpDValues[3];
    //dvalues[1] = omega;

    //dvalues[2] = -tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    //dvalues[2] = dvalues[2] * csi;
    //dvalues[3] = 0;

    flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omega, seobParams, H, 8 ,0);
    //fluxE  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, lMax);
    if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
    {
        print_warning("Failed to calculate flux.\n");
        //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
        return CEV_FAILURE;
    }

    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    flux = flux / eta;


    *Fr = ( values[2] / values[3] ) * flux / omega;
    *Fphi = flux / omega;
    *fluxout = flux;
    *fluxCircout = 0;
    REAL8 Lflux;
    Lflux = flux / omega;
    *Fr2PN = 0;
    *Fphi2PN = Lflux;
        COMPLEX16 h22Newt;
        REAL8 phi, rdot, rdot2, phi2, omega2, r2;
        phi = values[1];
        phi2 = phi*phi;
        rdot = GET_ABS(tmpDValues[0]);
        rdot2 = rdot*rdot;
        omega2 = omega*omega;
        r2 = r*r;
    h22Newt =  Power(rdot,2)*cos(2*phi) - sin(2*phi)/r - 
                    Power(r,2)*Power(omega,2)*cos(2*phi) - 
                    2*rdot*r*omega*sin(2*phi) + 
                    I*(-2*rdot*r*omega*cos(2*phi) - 
                    Power(rdot,2)*sin(2*phi) + sin(2*phi)/r + 
                    Power(r,2)*Power(omega,2)*sin(2*phi));
    h22Newt *=  12*eta*sqrt(CST_PI/5)/3;
    *FrFac = creal(h22Newt);
    *FphiFac = cimag(h22Newt);
        REAL8 FrModeSum;
        REAL8 pphi2, pphi4, r3, r6, pr2, prdot;
        pphi2 = pphi*pphi;
        pphi4 = pphi2*pphi2;
        r3 = r2*r;
        r6 = r3*r3;
        pr2 = pr*pr;
        prdot = tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    FrModeSum = pow(2*pr*prdot*Power(r,3)*
      (Power(pphi,2) + r*(-1 + Power(pr,2)*r)) - 
     (2*Power(pphi,4) + Power(r,2) - 
        Power(pr,2)*Power(r,3) + 
        Power(pphi,2)*r*(3 + 2*Power(pr,2)*r))*
      rdot,2)/
   (Power(r,6)*(Power(pphi,4) + 
       Power(r,2)*pow(-1 + Power(pr,2)*r,2) + 
       2*Power(pphi,2)*r*(1 + Power(pr,2)*r)));
    *FrResum = FrModeSum * 12*eta*sqrt(CST_PI/5)/3/8/CST_PI/rdot;
    *FphiResum = 0;

    return CEV_SUCCESS;
}
#undef Power


INT aux2PNEOBRadiationReactionForce(const REAL8 values[],
                                    SpinEOBParams *seobParams,
                                    REAL8 *Fr,
                                    REAL8 *Fphi)
{
    //print_debug("values = (%.2e, %.2e, %.2e, %.2e)\n", values[0], values[1], values[2], values[3]);
    HcapDerivParams params;

    REAL8           cartValues[6];

    gsl_function F;
    INT4         gslStatus;
    UINT i;

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r, prt, pr, pphi;
    REAL8Vector polarDynamics;
    REAL8       polData[4];

    REAL8 Mtotal, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;
    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = seobParams;
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
    prt = values[2];
    pphi = values[3];
    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );

    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );

    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);

    /* tortoise transformation */
    /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    pr = prt / csi;

    //pr = prt * csi;
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    //print_debug( "csi in derivatives function = %.16e\n", csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    //print_debug("cartValues = (%.2e, %.2e, %.2e, %.2e, %.2e, %.2e)\n", cartValues[0], cartValues[1], cartValues[2], cartValues[3], cartValues[4], cartValues[5]);
    REAL8 x4, x4_2, x4_3;
    params.varyParam = 0;
    gslStatus = gsl_deriv_central( &F, cartValues[0], 
                    STEP_SIZE, &x4, &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
        print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
        return CEV_FAILURE;
    }
    x4 *= r;
    x4_2 = x4*x4;
    x4_3 = x4_2*x4;

    REAL8 x2, x2_2, x2_3, x3, x3_2, x3_3;
    x2 = pr*pr;
    x2_2 = x2*x2;
    x2_3 = x2_2*x2;
    x3 = 1/r;
    x3_2 = x3*x3;
    x3_3 = x3_2*x3;

    REAL8 Fr_N, Fr_1PN, Fr_2PN;
    REAL8 Fphi_N, Fphi_1PN, Fphi_2PN;

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
    *Fr = (Fr_N + Fr_1PN + Fr_2PN) * prt * x3_3;

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
    *Fphi = (Fphi_N + Fphi_1PN + Fphi_2PN) * pphi * x3_3;

    return CEV_SUCCESS;
}


INT aux2PNFactorizedRadiationReactionForce(const REAL8 values[],
                                        SpinEOBParams *seobParams,
                                        REAL8 *Fr,
                                        REAL8 *Fphi)
{
    HcapDerivParams params;

    REAL8           cartValues[6];
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];

    params.values  = cartValues;
    params.params  = seobParams;

    gsl_function F;
    INT4         gslStatus;
    UINT i;

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r, prt, pr, pphi;
    REAL8Vector polarDynamics;
    REAL8       polData[4];
    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    REAL8 Mtotal, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a, omega, omegaCirc, Jcirc, Hreal, flux;
    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    nqcCoeffs = seobParams->nqcCoeffs;

    s1Vec = seobParams->s1VecOverMtMt;
    s2Vec = seobParams->s2VecOverMtMt;
    sKerr = seobParams->sigmaKerr;
    sStar = seobParams->sigmaStar;

    Mtotal = seobParams->eobParams->Mtotal;
    eta   = seobParams->eobParams->eta;

    r = values[0];
    prt = values[2];
    pphi = values[3];
    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( seobParams->seobCoeffs, r, eta, a );

    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( seobParams->seobCoeffs, r, eta, a );

    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);

    /* tortoise transformation */
    /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    pr = prt / csi;

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

    Hreal =  SpinEOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, seobParams->tortoise, seobParams->seobCoeffs );

    omega = PNCalcOrbitOmega(Hreal, pphi, eta);
    Jcirc = auxCalculateCircularAngularMomentum(eta, &rVec, s1Vec, s2Vec, sKerr, sStar, seobParams->tortoise, seobParams->seobCoeffs);
    cartValues[3] = values[2];
    cartValues[4] = pphi / r;
    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params   = &params;

    params.varyParam = 4;
    gslStatus = gsl_deriv_central( &F, cartValues[4], 
                    STEP_SIZE, &omegaCirc, &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
        print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
        return CEV_FAILURE;
    }

    omegaCirc /= r;


    flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omegaCirc, seobParams, Hreal, 8 ,0);
    //print_debug("flux = %.2e, omegaCirc = %.2e, r = %.2e, pr = %.2e, pphi = %.2e, Jcirc = %.2e\n", flux, omegaCirc, r, pr, pphi, Jcirc);
    if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
    {
        print_warning("Failed to calculate flux.\n");
        //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
        return CEV_FAILURE;
    }

    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    flux = flux / eta;

    /* Elliptic correction */
    REAL8 fN, f1PN, f2PN;
    REAL8 eta2, eta3, u, u2, u3, u4, pr2, pr4, pr6, pphi2, pphi4, pphi6, pphi8;
    u = 1./r;
    u2 = u*u;
    u3 = u2*u;
    u4 = u3*u;
    pr2 = pr*pr;
    pr4 = pr2*pr2;
    pr6 = pr4*pr2;
    pphi2 = pphi*pphi;
    pphi4 = pphi2*pphi2;
    pphi6 = pphi4*pphi2;
    pphi8 = pphi4*pphi4;
    fN = (pr2 - 6*u + 2 * pphi2 * u2 ) / ( 2.*u*(-3 + pphi2*u));
    f1PN = (-8*(79 + 40*eta)*pr4*(-3 + pphi2*u) + 
        pr2*u*(3353 - 4064*eta + (-1243 + 2980*eta)*pphi2*u + 
        9*(11 - 54*eta)*pphi4*u2)) / 
        (336.*u*pow(-3 + pphi2*u,2));
    f2PN = (40068*(3550 + 104*eta + 31*eta2 + 
        95*eta3)*pr6*pow(-3 + pphi2*u,2) + 
        636*pr4*u*(-9*(-1210132 + 423245*eta + 239556*eta2 + 
        42673*eta3) + 3*(-4017692 + 833396*eta + 523371*eta2 
        + 162667*eta3)*pphi2*u - (-4404988 + 315687*eta + 
        207558*eta2 + 157155*eta3)*pphi4*u2 + 
        3*(-177492 - 10496*eta - 8473*eta2 + 
        4129*eta3)*pphi6*u3) - 
        pr2*u2*(3*(26705157332 + 809392521*eta - 
        902131656*eta2 - 709549584*eta3) + 4*(-32097537395 + 
        107418492*eta + 207140589*eta2 + 
        533986713*eta3)*pphi2*u - 9*(-8236993436 + 
        212021253*eta - 82274338*eta2 + 
        77222272*eta3)*pphi4*u2 + 18*(-1018366510 + 
        40736065*eta - 22020228*eta2 + 
        3367832*eta3)*pphi6*u3 + 27*(60715424 - 
        2957347*eta + 1834330*eta2 + 
        100064*eta3)*pphi8*u4))/(2.6925696e7*eta*u*
        pow(-3 + pphi2*u,3));

    *Fr = fN + f1PN + f2PN;
    *Fphi = flux * (fN + f1PN + f2PN) / omegaCirc;

    return CEV_SUCCESS;
}




INT auxResummedRadiationReactionForce(const REAL8 values[],
                                    SpinEOBParams *seobParams,
                                    REAL8 *Fr,
                                    REAL8 *Fphi)
{
    //print_debug("values = (%.2e, %.2e, %.2e, %.2e)\n", values[0], values[1], values[2], values[3]);
    HcapDerivParams params;

    REAL8           cartValues[6];

    gsl_function F;
    INT4         gslStatus;
    UINT i;

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r, prt, pr, pphi;
    REAL8Vector polarDynamics;
    REAL8       polData[4];

    REAL8 Mtotal, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;
    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = seobParams;
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
    prt = values[2];
    pphi = values[3];
    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );

    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );

    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);

    /* tortoise transformation */
    /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    pr = prt / csi;

    //pr = prt * csi;
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    //print_debug( "csi in derivatives function = %.16e\n", csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    //print_debug("cartValues = (%.2e, %.2e, %.2e, %.2e, %.2e, %.2e)\n", cartValues[0], cartValues[1], cartValues[2], cartValues[3], cartValues[4], cartValues[5]);
    REAL8 x4, x4_2, x4_3;
    params.varyParam = 0;
    gslStatus = gsl_deriv_central( &F, cartValues[0], 
                    STEP_SIZE, &x4, &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
        print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
        return CEV_FAILURE;
    }
    x4 *= r;
    x4_2 = x4*x4;
    x4_3 = x4_2*x4;

    REAL8 x2, x2_2, x2_3, x3, x3_2, x3_3;
    x2 = pr*pr;
    x2_2 = x2*x2;
    x2_3 = x2_2*x2;
    x3 = 1/r;
    x3_2 = x3*x3;
    x3_3 = x3_2*x3;

    REAL8 Fr_N, Fr_1PN, Fr_2PN;
    REAL8 Fphi_N, Fphi_1PN, Fphi_2PN;


    return CEV_SUCCESS;
}



REAL8
FactorizedAngularMomentumFlux (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
				const REAL8 omega,	/**< orbital frequency */
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				INT const  lMax)
{

  REAL8 flux = 0.0;
  REAL8 v;
  COMPLEX16 hLM, hNQC, hECC;
  INT l, m;

  /* Omegs is the derivative of phi */

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
#if 0
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
	  //print_debug( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
	  /* Eq. 13 */
	  flux +=
	    (REAL8) (m * m) * omega * (creal (hLM) * creal (hLM) +
					 cimag (hLM) * cimag (hLM));
	}
    }
  return flux * CST_1_PI / 8.0;
}
