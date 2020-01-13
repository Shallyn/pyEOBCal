/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyEccEvolution.h"

#include "dyBHRingdown.h"
#include "dyNQCCorrection.h"
#include "dyHybridRingdown.h"
#include "eobEccCorrection.h"
#include "dyHcapNumericalDerivative.h"

#include <gsl/gsl_deriv.h>


#define EPS_ABS 1.0e-10
#define EPS_REL 1.0e-9
#define STEP_SIZE 1.0e-4


int CalculateSpinEOBHCoeffs (SpinEOBHCoeffs * coeffs,
					                const REAL8 eta,
					                const REAL8 a,
                                    AdjParams *adjParams);

REAL8 PNCalcOrbitOmega(const REAL8 Hreal,
                       const REAL8 ecc,
                       const REAL8 eta);


INT playDynamicEcc(REAL8 tend,
                REAL8 deltaT,
                REAL8 *init,
                int (* stop) (double t, const double y[], double dydt[], void *params),
                SpinEOBParams       *seobParams,
                SpinEOBDynamics     **dyout,
                CtrlParams          *ctrlParams);

REAL8 LocateEOBOmegaPeak(SpinEOBDynamics *dy,
                         REAL8Vector *omegaVec,
                         INT peakIdx);

INT CalculateNRPeakParams(REAL8Vector *tNR,
                          REAL8Vector *hNRreal,
                          REAL8Vector *hNRimag,
                          NRPeakParams *NRpeak);

INT ConstructFullWaveform(SpinEOBParams *seobParams,
                          SpinEOBDynamics *dy,
                          SpinEOBDynamics *dyHi,
                          REAL8Vector *sigReHi,
                          REAL8Vector *sigImHi,
                          REAL8Vector *omegaHi,
                          REAL8 nrTimePeak,
                          REAL8 deltaTHigh,
                          REAL8Vector *rdMatchPoint,
                          INT hiSRndx,
                          INT resampFac,
                          COMPLEX16TimeSeries **hLMAll_out);


/**
 * Note that here
 * values[0] = r
 * values[1] = phi
 * values[2] = pr
 * values[3] = pphi
 * dvalues[0] = dr/dt
 * dvalues[1] = dphi/dt
 * dvalues[2] = dpr/dt
 * dvalues[3] = dpphi/dt = omega
 */

int XLALSpinAlignedHiSRStopCondition(double t,  /**< UNUSED */
                           const double values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void *funcParams       /**< physical parameters */
                          );


static INT EOBHighSRStopConditionV4(double t,
                                    const double values[],
                                    double dvalues[],
                                    void *funcParams)
{
    REAL8 omega, r;
    UINT counter;
    SpinEOBParams *params = (SpinEOBParams *) funcParams;
    r = values[0];
    omega = dvalues[1];
    counter = params->eobParams->omegaPeaked;
    //printf("function 2: r = %.16e, omega = %.16e, pr = %.16e, dpr = %.16e, count = %.16u \n",values[0],dvalues[1],values[2],dvalues[2],counter);
    if (r < 6. && omega < params->eobParams->omega)
    {
        //        printf("Peak detection %.16e %.16e\n", omega, params->eobParams->omega);
        params->eobParams->omegaPeaked = counter + 1;
    }
    if (dvalues[2] >= 0. || params->eobParams->omegaPeaked == 5
        || isnan (dvalues[3]) || isnan (dvalues[2]) || isnan (dvalues[1])
        || isnan (dvalues[0]))
    {
        //        if ( dvalues[2] >= 0 ) printf("dvalues[2] >= 0\n");
        //        if ( params->eobParams->omegaPeaked == 5 ) printf("params->eobParams->omegaPeaked == 5\n");
        //        if ( isnan( dvalues[3] ) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) ) printf("%.16e %.16e %.16e %.16e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3]);
        return 1;
    }
    params->eobParams->omega = omega;
    return GSL_SUCCESS;
}

static int
EOBSpinAlignedStopCondition (double t, /**< UNUSED */
				 const double values[],
						  /**< dynamical variable values */
				 double dvalues[],/**< dynamical variable time derivative values */
				 void *funcParams /**< physical parameters */
  )
{

  REAL8 omega, r;
  SpinEOBParams *params = (SpinEOBParams *) funcParams;

  r = values[0];
  omega = dvalues[1];

//  printf("function 1: r = %.16e, omega = %.16e, pr = %.16e, dpr = %.16e, t = %.16e \n",values[0],dvalues[1],values[2],dvalues[2],t);
  //if ( omega < params->eobParams->omega )
  if (r < 6 && (omega < params->eobParams->omega || dvalues[2] >= 0))
    {
      params->eobParams->omegaPeaked = 0;
      return 1;
    }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/*-------------------------------------------------------------
*
*
*                       Full EOB Waveform
*
*
*--------------------------------------------------------------*/
INT EvolutionCoreEcc(const REAL8 m1,
                     const REAL8 m2,
                     const REAL8 fMin,
                     const REAL8 ecc,
                     const REAL8 deltaT,
                     const REAL8 spin1z,
                     const REAL8 spin2z,
                     COMPLEX16TimeSeries **hout,
                     AdjParams    *adjParams,
                     CtrlParams   *ctrlParams)
{
    INT i, status, failed = 0;
    REAL8 Mtotal = m1 + m2;
    REAL8 mTScaled = Mtotal * CST_MTSUN_SI;
    REAL8 eta = m1 * m2 / Mtotal / Mtotal;
#if DEBUG
print_debug("Mtotal = %.5f, eta = %.5f\n", Mtotal, eta);
#endif
    SpinEOBParams seobParams;
    EOBParams eobParams;
    FacWaveformCoeffs hCoeffs;
    NewtonMultipolePrefixes prefixes;
    EOBNonQCCoeffs nqcCoeffs;
    SpinEOBHCoeffs seobCoeffs;
    PNEccCorrectCoeffs eccCoeffs;
    memset (&seobParams, 0, sizeof (seobParams));
    memset (&seobCoeffs, 0, sizeof (seobCoeffs));
    memset (&eobParams, 0, sizeof (eobParams));
    memset (&hCoeffs, 0, sizeof (hCoeffs));
    memset (&prefixes, 0, sizeof (prefixes));
    memset (&eccCoeffs, 0, sizeof (eccCoeffs));

    /* Allocate spin parameters */
    REAL8Vector *sigmaStar = CreateREAL8Vector(3);
    REAL8Vector *sigmaKerr = CreateREAL8Vector(3);
    REAL8Vector *s1VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s2VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s1Vec = CreateREAL8Vector(3);
    REAL8Vector *s2Vec = CreateREAL8Vector(3);

    /* Dynamics */
    REAL8Vector *dyValues = NULL, *tmpValues = NULL;
    SpinEOBDynamics *dy = NULL,*dyHi = NULL;

    /* Waveform */
    REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
    REAL8Vector *ampNQC = NULL, *phaseNQC = NULL;
    REAL8Vector *omegaHi = NULL;
    COMPLEX16TimeSeries *hLMAll = NULL;

    REAL8 deltaTNU = deltaT / mTScaled;
    seobParams.deltaT = deltaTNU;
    seobParams.eccentricity = ecc;
    seobParams.s1Vec = s1Vec;
    seobParams.s2Vec = s2Vec;
    seobParams.s1VecOverMtMt = s1VecOverMtMt;
    seobParams.s2VecOverMtMt = s2VecOverMtMt;
    seobParams.sigmaStar = sigmaStar;
    seobParams.sigmaKerr = sigmaKerr;
    seobParams.eobParams = &eobParams;
    seobParams.nqcCoeffs = &nqcCoeffs;
    seobParams.seobCoeffs = &seobCoeffs;
    seobParams.eccCoeffs = &eccCoeffs;
    seobParams.tortoise = 1;

    eobParams.hCoeffs = &hCoeffs;
    eobParams.prefixes = &prefixes;
    eobParams.m1 = m1;
    eobParams.m2 = m2;
    eobParams.eta = eta;
    eobParams.Mtotal = Mtotal;

    REAL8 tplspin;
    REAL8   chiS, chiA;

    /*----------------------------------------------------*/
    /*                   Calculate spin                   */
    /*----------------------------------------------------*/
    s1Vec->data[0] = 0;
    s1Vec->data[1] = 0;
    s1Vec->data[2] = spin1z * m1 * m1;

    s2Vec->data[0] = 0;
    s2Vec->data[1] = 0;
    s2Vec->data[2] = spin2z * m2 * m2;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] / Mtotal / Mtotal;
    }

    if(!sigmaStar || !sigmaKerr || !s1VecOverMtMt || !s2VecOverMtMt)
    {
        failed = 1;
        goto END;
    }
    CalculateSigmaStar(sigmaStar, m1, m2, s1Vec, s2Vec);
    CalculateSigmaKerr(sigmaKerr, m1, m2, s1Vec, s2Vec);

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = spin1z;
    seobParams.chi2 = spin2z;


    chiS = 0.5 * (spin1z + spin2z);
    chiA = 0.5 * (spin1z - spin2z);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;
#if DEBUG
print_debug("chiS = %.5f, chiA = %.5f\n\tsigmaStar = %.5f, sigmaKerr = %.5f\n", chiS, chiA, sigmaStar->data[2], sigmaKerr->data[2]);
print_debug("a = %.5f\n", seobParams.a);
#endif
    /*----------------------------------------------------*/
    /*            Factorized  Waveform Coeffs             */
    /*----------------------------------------------------*/
    status = CalculateSpinFactorizedWaveformCoefficients(&hCoeffs, &seobParams, m1, m2, eta, tplspin, chiS, chiA);
    if (status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }

    /*----------------------------------------------------*/
    /*                Eccentric Coefficients              */
    /*----------------------------------------------------*/
    status = PNwaveformPRD100_044018_22mode_CalculateCoefficients(eta, &eccCoeffs);
    if (status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }

    /*----------------------------------------------------*/
    /*                 Newton Multipolars                 */
    /*----------------------------------------------------*/
    status = ComputeNewtonMultipolePrefixes(&prefixes, m1, m2);
    if (status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }

    /*----------------------------------------------------*/
    /*                   Hyper Coeffs                     */
    /*----------------------------------------------------*/
    status = CalculateSpinEOBHCoeffs(&seobCoeffs, eta, seobParams.a, adjParams);
    if (status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }
#if DEBUG
        print_debug("KK = %.5f, dSOv1 = %.5f, dSOv2 = %.5f\n\tdheffSS = %.5f, dheffSSv2 = %.5f\n",
                  seobCoeffs.KK, seobCoeffs.d1, seobCoeffs.d1v2,
                  seobCoeffs.dheffSS, seobCoeffs.dheffSSv2);

#endif
    /*----------------------------------------------------*/
    /*                     NQC Coeffs                     */
    /*----------------------------------------------------*/
    nqcCoeffs.a1 = 0.;
    nqcCoeffs.a2 = 0.;
    nqcCoeffs.a3 = 0.;
    nqcCoeffs.a3S = 0.;
    nqcCoeffs.a4 = 0.;
    nqcCoeffs.a5 = 0.;
    nqcCoeffs.b1 = 0.;
    nqcCoeffs.b2 = 0.;
    nqcCoeffs.b3 = 0.;
    nqcCoeffs.b4 = 0.;
#if NQCv1
    XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, seobParams.a );   
#endif

#if DEBUG
print_debug("Initial condition:\n");
#endif
    /*----------------------------------------------------*/
    /*                 Initial Condition                  */
    /*----------------------------------------------------*/
    REAL8Vector cartPosVec, cartMomVec;
    REAL8       cartPosData[3], cartMomData[3];
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;

    dyValues = CreateREAL8Vector(4);
    tmpValues = CreateREAL8Vector (14);
    if (!dyValues || !tmpValues)
    {
        failed = 1;
        goto END;
    }
    memset(dyValues->data, 0, dyValues->length * sizeof(REAL8));
    memset (tmpValues->data, 0, tmpValues->length * sizeof (REAL8));
    status = CalcInitialConditions(tmpValues, m1, m2, fMin, &seobParams);
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }
#if DEBUG
print_err("\t(x,y,z) = (%.3e, %.3e, %.3e), (px,py,pz) = (%.3e, %.3e, %.3e)\n", 
    tmpValues->data[0], tmpValues->data[1], tmpValues->data[2], tmpValues->data[3],
    tmpValues->data[4], tmpValues->data[5], tmpValues->data[6]);

print_debug("Evolve EOB LowSR.\n");
#endif
    /*----------------------------------------------------*/
    /*                     Evolve EOB                     */
    /*----------------------------------------------------*/
    INT retLen;
    dyValues->data[0] = tmpValues->data[0];
    dyValues->data[1] = 0.;
    dyValues->data[2] = tmpValues->data[3];
    dyValues->data[3] = tmpValues->data[0] * tmpValues->data[4];
#if DEBUG
print_debug("Spherical initial conditions: %e %e %e %e\n", dyValues->data[0], dyValues->data[1], dyValues->data[2], dyValues->data[3] );
#endif
    eobParams.rad = dyValues->data[0];
    eobParams.omegaPeaked = 0;

    status = playDynamicEcc(20 / mTScaled, 
            deltaTNU, dyValues->data, 
            EOBSpinAlignedStopCondition, &seobParams, 
            &dy, ctrlParams);

    if (status != CEV_SUCCESS)
    {
        print_warning("Dynamic evolution failed.\n");
        failed = 1;
        goto END;
    }
    retLen = dy->length;
    //*dyout = dy;

#if DEBUG
print_debug("retLen = %d\n", retLen);
print_debug("Evolve EOBHiSR.\n");
#endif
    /*----------------------------------------------------*/
    /*                   Evolve EOB HiSR                  */
    /*----------------------------------------------------*/
    INT retLenHi;
    INT resampFac = 1;
    INT resampPwr;
    REAL8 resampEstimate = 50 * deltaTNU;
    if (resampEstimate > 1.)
    {
        resampPwr = (UINT) ceil (log2 (resampEstimate));
        while (resampPwr--)
	    {
	        resampFac *= 2u;
	    }
    }
    REAL8 deltaTHigh = deltaT / (REAL8) resampFac;
    REAL8 deltaTHighNU = deltaTHigh / mTScaled;

    INT nStepBack;
    REAL8 tStepBack = 150 * mTScaled;
    if(tStepBack > retLen * deltaT)
        tStepBack = 0.5 * retLen * deltaT;
    nStepBack = ceil(tStepBack / deltaT);
    INT hiSRndx = retLen - nStepBack;

    dyValues->data[0] = dy->rVec->data[hiSRndx];
    dyValues->data[1] = dy->phiVec->data[hiSRndx];
    dyValues->data[2] = dy->prVec->data[hiSRndx];
    dyValues->data[3] = dy->pPhiVec->data[hiSRndx];
    eobParams.rad = dyValues->data[0];
    eobParams.omegaPeaked = 0;

    status = playDynamicEcc(20 / mTScaled, 
        deltaTHighNU, dyValues->data, 
        XLALSpinAlignedHiSRStopCondition, &seobParams, 
        &dyHi, ctrlParams);
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }
    retLenHi = dyHi->length;

#if DEBUG
print_debug("retLenHi = %df\n", retLenHi);
print_debug("Compute QNMfreq.\n");
#endif
    /*----------------------------------------------------*/
    /*                  Compute QNMfreq                   */
    /*----------------------------------------------------*/
    COMPLEX16Vector modefreqVec;
    COMPLEX16 modeFreq;
    modefreqVec.length = 1;
    modefreqVec.data = &modeFreq;
    status = XLALSimIMREOBGenerateQNMFreqV2(&modefreqVec, m1, m2, spin1z, spin2z, 2, 2, 1);
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }

#if DEBUG
print_debug("Populate HiSR wf.\n");
#endif
    /*----------------------------------------------------*/
    /*    Populate HiSR waveform and locate omegaPeak     */
    /*----------------------------------------------------*/
    REAL8 omega, v, omegaOld = 0.0;
    REAL8 eccentricity, xi;
    UINT rdLen = retLenHi + (UINT)ceil(20 / (cimag(modeFreq) * deltaTHigh));
    sigReHi = CreateREAL8Vector(rdLen);
    sigImHi = CreateREAL8Vector(rdLen);
    omegaHi = CreateREAL8Vector(rdLen);
    ampNQC = CreateREAL8Vector(retLenHi);
    phaseNQC = CreateREAL8Vector(retLenHi);
    if (!sigReHi || !sigImHi || 
        !omegaHi || !ampNQC || 
        !phaseNQC)
    {
        failed = 1;
        goto END;
    }
    INT phaseCounter = 0, peakIdx = 0, finalIdx = 0;
    COMPLEX16 hLM, hECC;
    REAL8 ham;
    for(i=0; i<retLenHi; i++)
    {
        dyValues->data[0] = dyHi->rVec->data[i];
        dyValues->data[1] = dyHi->phiVec->data[i];
        dyValues->data[2] = dyHi->prVec->data[i];
        dyValues->data[3] = dyHi->pPhiVec->data[i];

        cartPosVec.data[0] = dyValues->data[0];
        cartMomVec.data[0] = dyValues->data[2];
        cartMomVec.data[1] = dyValues->data[3] / dyValues->data[0];
        ham = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec,
                s1VecOverMtMt, s2VecOverMtMt,
                sigmaKerr, sigmaStar,
                seobParams.tortoise, &seobCoeffs);
        
        eccentricity = Calculate3PNEccentricity(eta, dyValues, ham-1);
        //omega = PNCalcOrbitOmega(ham, dyValues->data[3], eta);
        omega = XLALSimIMRSpinAlignedEOBCalcOmega(dyValues->data, &seobParams, STEP_SIZE);

        if(omega < 1.0e-15)
            omega = 1.0e-9;
        omegaHi->data[i] = omega;
        v = cbrt(omega);
        xi = Calculate3PNxi(ham-1, eta, v, dyHi->phiVec->data[i]);

        status = SpinEOBCalculateFactorizedWaveform_ecc(&hLM, dyValues, v, ham, 2, 2, eccentricity, xi, &seobParams);
        if(status != CEV_SUCCESS)
        {
            failed = 1;
            goto END;
        }


        ampNQC->data[i] = cabs(hLM);
        sigReHi->data[i] = creal(hLM);
        sigImHi->data[i] = cimag(hLM);
        phaseNQC->data[i] = carg(hLM) + phaseCounter * CST_2PI;
        if(i && phaseNQC->data[i] > phaseNQC->data[i-1])
        {
            phaseCounter--;
            phaseNQC->data[i] -= CST_2PI;
        }

        if (omega <= omegaOld && !peakIdx)
        {
            peakIdx = i;
        }
        omegaOld = omega;
    }
    finalIdx = retLenHi - 1;
    if(!peakIdx)
        peakIdx = finalIdx;

    /* Stuff to find the actual peak time */
    gsl_spline *spline = NULL;
    gsl_interp_accel *acc = NULL;
    REAL8 omegaDeriv1;		//, omegaDeriv2;
    REAL8 time1, time2;
    REAL8 timePeak, timewavePeak = 0., omegaDerivMid;
    REAL8 sigAmpSqHi = 0., oldsigAmpSqHi = 0.;
    INT peakCount = 0;

    //timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaTv4 (2, 2, m1, m2, s1Vec->data[2], s2Vec->data[2]);
    timewavePeak = adjParams->dtPeak;

    spline = gsl_spline_alloc (gsl_interp_cspline, retLenHi);
    acc = gsl_interp_accel_alloc ();

    time1 = dyHi->tVec->data[peakIdx];

    gsl_spline_init (spline, dyHi->tVec->data, omegaHi->data, retLenHi);
    omegaDeriv1 = gsl_spline_eval_deriv (spline, time1, acc);
    if(omegaDeriv1 > 0.)
    {
        time2 = dyHi->tVec->data[peakIdx + 1];
    }
    else
    {
        time2 = time1;
        time1 = dyHi->tVec->data[peakIdx-1];
        peakIdx--;
        omegaDeriv1 = gsl_spline_eval_deriv(spline, time1, acc);
    }

    do
    {
        timePeak = (time1 + time2) / 2.;
        omegaDerivMid = gsl_spline_eval_deriv(spline, timePeak, acc);
        if (omegaDerivMid * omegaDeriv1 < 0.0)
        {
            time2 = timePeak;
        }
        else
        {
            omegaDeriv1 = omegaDerivMid;
            time1 = timePeak;
        }
    }
    while(time2 - time1 > 1.0e-5);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

#if DEBUG
print_debug("timePeak = %.5f\n", timePeak);
print_debug("Calculate NQC.\n");
#endif
    /*----------------------------------------------------*/
    /*              Calculate NQC correction              */
    /*----------------------------------------------------*/
#if NQCv1

#else
    //if(timewavePeak < 1.0e-16 || peakCount == 0)
    timewavePeak = timePeak - timewavePeak;
    if(timewavePeak < 0 || timewavePeak > dyHi->tVec->data[retLenHi-1])
    {
        print_warning("Too much dtPeak stepback: %.3e\n", timewavePeak);
        failed = 1;
        goto END;
    }

    status = XLALSimIMRSpinEOBCalculateNQCCoefficientsV4
                (ampNQC, phaseNQC, dyHi, omegaHi, 2, 2, timewavePeak,
                 deltaTHighNU, m1, m2, sigmaKerr->data[2], chiA, chiS, &nqcCoeffs);
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }
#endif

#if DEBUG
print_debug("a1 = %.5f, a2 = %.5f a3 = %.5f\n", 
            nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3);
print_debug("a3S = %.5f, a4 = %.5f, a5 = %.5f\n", 
            nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5);
print_debug("b1 = %.5f, b2 = %.5f, b3 = %.5f, b4 = %.5f\n", 
            nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4);

print_debug("Apply high SR waveform.\n");
#endif
    /*----------------------------------------------------*/
    /*               Apply high SR waveform               */
    /*----------------------------------------------------*/
    COMPLEX16 hNQC;

    for(i = 0; i < retLenHi; i++)
    {
        dyValues->data[0] = dyHi->rVec->data[i];
        dyValues->data[1] = dyHi->phiVec->data[i];
        dyValues->data[2] = dyHi->prVec->data[i];
        dyValues->data[3] = dyHi->pPhiVec->data[i];
        status = XLALSimIMREOBNonQCCorrection(&hNQC, dyValues, omegaHi->data[i], &nqcCoeffs);
        if (status != CEV_SUCCESS)
        {
            failed = 1;
            goto END;
        }

        hLM = (sigReHi->data[i] + I*sigImHi->data[i]+ hECC) * hNQC;
        sigReHi->data[i] = (REAL8) creal(hLM);
        sigImHi->data[i] = (REAL8) cimag(hLM);
    }

#if DEBUG
print_debug("Calculate QNM excitation coeffs.\n");
#endif
    /*----------------------------------------------------*/
    /*           Calculate QNM excitation coeffs          */
    /*----------------------------------------------------*/

    REAL8 combSize;
    REAL8 chi = chiS + chiA * ((m1 - m2) / (m1 + m2)) / (1. - 2. * eta);
    
    REAL8 timeshiftPeak;
    timeshiftPeak = (timePeak - timewavePeak) > 0. ? (timePeak - timewavePeak) : 0.;
    combSize = (spin1z == 0.
		  && spin2z == 0.) ? 11. : 
          ((eta > 10. / 121. && chi >= 0.8) ? 8.5 : 12.);
    if ((eta > 30. / 31. / 31. && eta <= 10. / 121. && chi >= 0.8)
	  || (eta <= 30. / 31. / 31. && chi >= 0.8 && chi < 0.9))
	    combSize = 13.5;

    REAL8Vector *rdMatchPoint = CreateREAL8Vector(4);
    if(!rdMatchPoint)
    {
        failed = 1;
        goto END;
    }

    rdMatchPoint->data[0] =
        combSize <
        timePeak - timeshiftPeak ? timePeak - timeshiftPeak - combSize : 0;
    //rdMatchPoint->data[0] = timePeak + 2.0*deltaTHigh;
    rdMatchPoint->data[1] = timePeak - timeshiftPeak;
    rdMatchPoint->data[2] = dyHi->tVec->data[finalIdx];
    rdMatchPoint->data[0] -=
        fmod (rdMatchPoint->data[0], deltaTHighNU);
    rdMatchPoint->data[1] -=
        fmod (rdMatchPoint->data[1], deltaTHighNU);

    for(i = retLenHi; i < sigReHi->length; i++)
    {
        sigReHi->data[i] = 0;
        sigImHi->data[i] = 0;
    }
#if 0
status = XLALSimIMREOBHybridAttachRingdown( sigReHi, sigImHi, 2, 2,
                                           deltaTHigh, m1, m2, spin1z, spin2z,
                                           dyHi->tVec, rdMatchPoint);
#else
    status = XLALSimIMREOBAttachFitRingdown(sigReHi, sigImHi,deltaTHigh, m1, m2, spin1z, spin2z, dyHi->tVec, rdMatchPoint);
#endif
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }

#if DEBUG
print_debug("Get full IMRwaveform.\n");
#endif
    /*----------------------------------------------------*/
    /*               Get full IMRwaveform                 */
    /*----------------------------------------------------*/

    /* ...Dump... */
    gboolean debug_dump = FALSE;
    CHAR file_dump_dynamics[256];
    FILE *out_dynamics = NULL;
    if(strcmp(ctrlParams->dump, "None")!=0)
    {
        cmd_mkdir(ctrlParams->dump);
        sprintf(file_dump_dynamics, "%s/dynamics.dat", ctrlParams->dump);
        print_debug("Dump Energy flux to %s...\n", file_dump_dynamics);
        out_dynamics = fopen(file_dump_dynamics, "w");
        fprintf(out_dynamics, "#time #r #phi #pr #pPhi #omega #eccentricity #omegaold #Ham #hreal #himag\n");
        debug_dump = TRUE;
    }
    /* ...Dump... */
    hLMAll = CreateCOMPLEX16TimeSeries(0, deltaT, retLen + ceil (sigReHi->length / resampFac));
    if(!hLMAll)
    {
        failed = 1;
        goto END;
    }
    memset (hLMAll->data->data, 0, hLMAll->data->length * sizeof (COMPLEX16));
    REAL8 omega_old;
    //print_err("eta = %.16e\n", eta);
    for(i=0; i<retLen; i++)
    {
        dyValues->data[0] = dy->rVec->data[i];
        dyValues->data[1] = dy->phiVec->data[i];
        dyValues->data[2] = dy->prVec->data[i];
        dyValues->data[3] = dy->pPhiVec->data[i];
        cartPosVec.data[0] = dyValues->data[0];
        cartMomVec.data[0] = dyValues->data[2];
        cartMomVec.data[1] = dyValues->data[3] / dyValues->data[0];
        ham = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec, 
            s1VecOverMtMt, s2VecOverMtMt, sigmaKerr, sigmaStar, 
            seobParams.tortoise, &seobCoeffs);
        
        eccentricity = Calculate3PNEccentricity(eta, dyValues, ham-1);
        omega = PNCalcOrbitOmega(ham, dyValues->data[3], eta);
        
        v = cbrt(omega);
        xi = Calculate3PNxi(ham-1, eta, v, dyHi->phiVec->data[i]);
//print_err("%.16e %.16e %.16e\n", i*deltaTNU, omega, omega_old);
        status = SpinEOBCalculateFactorizedWaveform_ecc(&hLM, dyValues, v, ham, 2, 2, eccentricity, xi, &seobParams);
        if(status != CEV_SUCCESS)
        {
            failed = 1;
            goto END;
        }
        status = XLALSimIMREOBNonQCCorrection(&hNQC, dyValues, omega, &nqcCoeffs);
        if(status != CEV_SUCCESS)
        {
            failed = 1;
            goto END;
        }
    //print_err("%.16e %.16e %.16e\n", i*deltaTNU, creal(hLM), cimag(hLM));

        hLMAll->data->data[i] = hLM * hNQC;
        /* ...Dump... */
        if(debug_dump)
        {
            omega_old = XLALSimIMRSpinAlignedEOBCalcOmega(dyValues->data, &seobParams, STEP_SIZE);
            // #time #r #phi #pr #pPhi #omega #eccentricity #omegaold #Ham #hreal #himag
            fprintf(out_dynamics, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                i*deltaTNU, dyValues->data[0], dyValues->data[1], dyValues->data[2], dyValues->data[3],
                omega, eccentricity, omega_old, ham, creal(hLM), cimag(hLM));
        }
    }

    if(debug_dump)
    {
        fclose(out_dynamics);
    }

#if DEBUG
print_debug("Attach ringdown.\n");
print_debug("hiSRndx = %d\n", hiSRndx);
print_debug("hLMAllLen = %d, sigHiLen = %d\n", hLMAll->data->length, (INT)(retLenHi / resampFac));
#endif
    /* Attach ringdown */
    for(i = 0; i < (INT)(sigReHi->length / resampFac); i++ )
    {
        hLMAll->data->data[i + hiSRndx] = sigReHi->data[i*resampFac] + I*sigImHi->data[i*resampFac];
    }
    /*----------------------------------------------------*/
    /*                 Free memory, exit                  */
    /*----------------------------------------------------*/
#if DEBUG
print_debug("END, free memory.\n");
#endif
    END:
    if(sigmaStar)
        DestroyREAL8Vector(sigmaStar);
    if(sigmaKerr)
        DestroyREAL8Vector(sigmaKerr);
    if(s1Vec)
        DestroyREAL8Vector(s1Vec);
    if(s2Vec)
        DestroyREAL8Vector(s2Vec);
    if(s1VecOverMtMt)
        DestroyREAL8Vector(s1VecOverMtMt);
    if(s2VecOverMtMt)
        DestroyREAL8Vector(s2VecOverMtMt);
    if(tmpValues)
        DestroyREAL8Vector(tmpValues);
    if(dyValues)
        DestroyREAL8Vector(dyValues);
    if(dyHi)
        DestroySpinEOBDynamics(dyHi);
    if(sigReHi)
        DestroyREAL8Vector(sigReHi);
    if(sigImHi)
        DestroyREAL8Vector(sigImHi);
    if(ampNQC)
        DestroyREAL8Vector(ampNQC);
    if(phaseNQC)
        DestroyREAL8Vector(phaseNQC);
    if(omegaHi)
        DestroyREAL8Vector(omegaHi);
    if(dy)
        DestroySpinEOBDynamics(dy);

    if (failed)
    {
        if(hLMAll)
            DestroyCOMPLEX16TimeSeries(hLMAll);
        return CEV_FAILURE;
    }
    *hout = hLMAll;
    return CEV_SUCCESS;
}


/*
*   This Function Calculate EOB dynamics evolution.
*   from initial condition *init, evolve to tend.
*   The output dynamics will store in dyout
*/
INT playDynamicEcc(REAL8 tend,
                   REAL8 deltaT,
                   REAL8 *init,
                   int (* stop) (double t, const double y[], double dydt[], void *params),
                   SpinEOBParams       *seobParams,
                   SpinEOBDynamics     **dyout,
                   CtrlParams          *ctrlParams)
{    
    INT status;
    INT failed = 0;

    INT retLen;
    RK4GSLIntegrator        *integrator = NULL;
    REAL8Array              *dynamics   = NULL;
   
    REAL8 omega, v;
    SpinEOBDynamics *dyAll = NULL;

    integrator = InitRK4GSLIntegrator(4, XLALSpinAlignedHcapDerivative_ecc, stop, EPS_ABS, EPS_REL);
    if(!integrator || !init || !stop)
    {
        print_warning("Failed to allocate integrator.\n");
        failed = 1;
        goto EXIT;
    }
    integrator->stopontestonly = 1;
    integrator->retries = 1;

    retLen = playRungeKutta4(integrator, seobParams, init, 0, tend, deltaT, &dynamics);
    if (retLen == CEV_FAILURE || dynamics == NULL)
    {
        print_warning("Run RungeKutta4 failed, exit %d.\n", retLen);
        failed = 1;
        goto EXIT;
    } 

    dyAll = SpinEOBDynamicsInit(retLen);
    
    memcpy(dyAll->tVec->data, dynamics->data, retLen*sizeof(REAL8));
    memcpy(dyAll->rVec->data, dynamics->data + retLen, retLen*sizeof(REAL8));
    memcpy(dyAll->phiVec->data, dynamics->data + 2 * retLen, retLen*sizeof(REAL8));
    memcpy(dyAll->prVec->data, dynamics->data + 3 * retLen, retLen*sizeof(REAL8));
    memcpy(dyAll->pPhiVec->data, dynamics->data + 4 * retLen, retLen*sizeof(REAL8));

    memcpy(dyAll->drVec->data, dynamics->data + 5 * retLen, retLen*sizeof(REAL8));
    memcpy(dyAll->dphiVec->data, dynamics->data + 6 * retLen, retLen * sizeof(REAL8));
    memcpy(dyAll->dprVec->data, dynamics->data + 7 * retLen, retLen * sizeof(REAL8));
    memcpy(dyAll->dpPhiVec->data, dynamics->data + 8 * retLen, retLen * sizeof(REAL8));

    *dyout = dyAll;

    EXIT:
    if(integrator)
        DestroyIntegrator(integrator);
    if(dynamics)
        DestroyREAL8Array(dynamics);

    if(failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}
