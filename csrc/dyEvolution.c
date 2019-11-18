/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyEvolution.h"
#include "dyBHRingdown.h"
#include "dyNQCCorrection.h"
#include "dyHybridRingdown.h"
#include <gsl/gsl_deriv.h>

#define DEBUG 1

#define EPS_ABS 1.0e-10
#define EPS_REL 1.0e-9
#define STEP_SIZE 1.0e-4

static int CalculateSpinEOBHCoeffs (SpinEOBHCoeffs * coeffs,
					                const REAL8 eta,
					                const REAL8 a,
                                    AdjParams *adjParams);
static int XLALSpinAlignedHcapDerivative(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

INT playDynamic(REAL8 tend,
                REAL8 deltaT,
                REAL8 *init,
                int (* stop) (double t, const double y[], double dydt[], void *params),
                SpinEOBParams       *seobParams,
                SpinEOBDynamics     **dyout,
                CtrlParams          *ctrlParams);


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


INT EvolutionCore(const REAL8 m1,
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
    memset (&seobParams, 0, sizeof (seobParams));
    memset (&seobCoeffs, 0, sizeof (seobCoeffs));
    memset (&eobParams, 0, sizeof (eobParams));
    memset (&hCoeffs, 0, sizeof (hCoeffs));
    memset (&prefixes, 0, sizeof (prefixes));

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
print_err("\tr = %e, pr = %e\n\tphi = %e, pPhi = %e\n", dyValues->data[0], dyValues->data[1], dyValues->data[2], dyValues->data[3]);
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

    eobParams.rad = dyValues->data[0];
    eobParams.omegaPeaked = 0;

    status = playDynamic(20 / mTScaled, 
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

    status = playDynamic(20 / mTScaled, 
        deltaTHighNU, dyValues->data, 
        EOBHighSRStopConditionV4, &seobParams, 
        &dyHi, ctrlParams);
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }
    retLenHi = dyHi->length;

#if DEBUG
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
    COMPLEX16 hLM;
    REAL8 ham;
    for(i=0; i<retLenHi; i++)
    {
        dyValues->data[0] = dyHi->rVec->data[i];
        dyValues->data[1] = dyHi->phiVec->data[i];
        dyValues->data[2] = dyHi->prVec->data[i];
        dyValues->data[3] = dyHi->pPhiVec->data[i];

        omega = XLALSimIMRSpinAlignedEOBCalcOmega(dyValues->data, &seobParams, STEP_SIZE);
        if(omega < 1.0e-15)
            omega = 1.0e-9;
        omegaHi->data[i] = omega;
        v = cbrt(omega);
    
        cartPosVec.data[0] = dyValues->data[0];
        cartMomVec.data[0] = dyValues->data[2];
        cartMomVec.data[1] = dyValues->data[3] / dyValues->data[0];
        ham = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec,
                s1VecOverMtMt, s2VecOverMtMt,
                sigmaKerr, sigmaStar,
                seobParams.tortoise, &seobCoeffs);
        status = XLALSimIMRSpinEOBGetSpinFactorizedWaveform(&hLM, dyValues, v, ham, 2, 2, &seobParams);
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
print_debug("Calculate NQC.\n");
#endif
    /*----------------------------------------------------*/
    /*              Calculate NQC correction              */
    /*----------------------------------------------------*/
    status = XLALSimIMRSpinEOBCalculateNQCCoefficientsV4
                (ampNQC, phaseNQC, dyHi, omegaHi, 2, 2, timePeak,
                 deltaTHighNU, m1, m2, sigmaKerr->data[2], chiA, chiS, &nqcCoeffs);
    if(status != CEV_SUCCESS)
    {
        failed = 1;
        goto END;
    }

#if DEBUG
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
        hLM = (sigReHi->data[i] + I*sigImHi->data[i]) * hNQC;
        sigReHi->data[i] = (REAL8) creal(hLM);
        sigImHi->data[i] = (REAL8) cimag(hLM);
    }

#if DEBUG
print_debug("Calculate QNM excitation coeffs.\n");
#endif
    /*----------------------------------------------------*/
    /*           Calculate QNM excitation coeffs          */
    /*----------------------------------------------------*/
    //timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaTv4 (2, 2, m1, m2, s1Vec->data[2], s2Vec->data[2]);
    timewavePeak = adjParams->dtPeak;

    if(timewavePeak < 1.0e-16 || peakCount == 0)
        timewavePeak = timePeak - timewavePeak;
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

    status = XLALSimIMREOBAttachFitRingdown(sigReHi, sigImHi,deltaTHigh, m1, m2, spin1z, spin2z, dyHi->tVec, rdMatchPoint);
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
    hLMAll = CreateCOMPLEX16TimeSeries(0, deltaT, retLen + ceil (sigReHi->length / resampFac));
    if(!hLMAll)
    {
        failed = 1;
        goto END;
    }
    memset (hLMAll->data->data, 0, hLMAll->data->length * sizeof (COMPLEX16));
    for(i=0; i<retLen; i++)
    {
        dyValues->data[0] = dy->rVec->data[i];
        dyValues->data[1] = dy->phiVec->data[i];
        dyValues->data[2] = dy->prVec->data[i];
        dyValues->data[3] = dy->pPhiVec->data[i];
        omega = XLALSimIMRSpinAlignedEOBCalcOmega(dyValues->data, &seobParams, STEP_SIZE);
        v = cbrt(omega);
        cartPosVec.data[0] = dyValues->data[0];
        cartMomVec.data[0] = dyValues->data[2];
        cartMomVec.data[1] = dyValues->data[3] / dyValues->data[0];
        ham = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec, 
            s1VecOverMtMt, s2VecOverMtMt, sigmaKerr, sigmaStar, 
            seobParams.tortoise, &seobCoeffs);
        status = XLALSimIMRSpinEOBGetSpinFactorizedWaveform(&hLM, dyValues, v, ham, 2, 2, &seobParams);
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
        hLM *= hNQC;
        hLMAll->data->data[i] = hLM;
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


INT playDynamic(REAL8 tend,
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

    integrator = InitRK4GSLIntegrator(4, XLALSpinAlignedHcapDerivative, stop, EPS_ABS, EPS_REL);
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

INT applyDefaultAdjustableParameters(AdjParams *adjParams,
                                     const REAL8 m1,
                                     const REAL8 m2,
                                     const REAL8 spin1z,
                                     const REAL8 spin2z)
{
    REAL8 eta = m1*m2 / (m1+m2) / (m1+m2);
    REAL8 a = (spin1z*m1*m1 + spin2z*m2*m2) / (m1+m2) / (m1+m2);
    REAL8 chi = a / (1. - 2. * eta);
    REAL8 eta2 = eta * eta, eta3 = eta2 * eta;
    REAL8 chi2 = chi * chi, chi3 = chi2 * chi;
//print_debug("a = %.5f, chi = %.5f\n", a, chi);
    // K
    static const REAL8 coeff00K = 1.7336;
    static const REAL8 coeff01K = -1.62045;
    static const REAL8 coeff02K = -1.38086;
    static const REAL8 coeff03K = 1.43659;
    static const REAL8 coeff10K = 10.2573;
    static const REAL8 coeff11K = 2.26831;
    static const REAL8 coeff12K = 0;
    static const REAL8 coeff13K = -0.426958;
    static const REAL8 coeff20K = -126.687;
    static const REAL8 coeff21K = 17.3736;
    static const REAL8 coeff22K = 6.16466;
    static const REAL8 coeff23K = 0;
    static const REAL8 coeff30K = 267.788;
    static const REAL8 coeff31K = -27.5201;
    static const REAL8 coeff32K = 31.1746;
    static const REAL8 coeff33K = -59.1658;

    // dSO
    static const REAL8 coeff00dSO = -44.5324;
    static const REAL8 coeff01dSO = 0;
    static const REAL8 coeff02dSO = 0;
    static const REAL8 coeff03dSO = 66.1987;
    static const REAL8 coeff10dSO = 0;
    static const REAL8 coeff11dSO = 0;
    static const REAL8 coeff12dSO = -343.313;
    static const REAL8 coeff13dSO = -568.651;
    static const REAL8 coeff20dSO = 0;
    static const REAL8 coeff21dSO = 2495.29;
    static const REAL8 coeff22dSO = 0;
    static const REAL8 coeff23dSO = 147.481;
    static const REAL8 coeff30dSO = 0;
    static const REAL8 coeff31dSO = 0;
    static const REAL8 coeff32dSO = 0;
    static const REAL8 coeff33dSO = 0;

    // dSS
    static const REAL8 coeff00dSS = 6.06807;
    static const REAL8 coeff01dSS = 0;
    static const REAL8 coeff02dSS = 0;
    static const REAL8 coeff03dSS = 0;
    static const REAL8 coeff10dSS = -36.0272;
    static const REAL8 coeff11dSS = 37.1964;
    static const REAL8 coeff12dSS = 0;
    static const REAL8 coeff13dSS = -41.0003;
    static const REAL8 coeff20dSS = 0;
    static const REAL8 coeff21dSS = 0;
    static const REAL8 coeff22dSS = -326.325;
    static const REAL8 coeff23dSS = 528.511;
    static const REAL8 coeff30dSS = 706.958;
    static const REAL8 coeff31dSS = 0;
    static const REAL8 coeff32dSS = 1161.78;
    static const REAL8 coeff33dSS = 0.;

    adjParams->KK =
        coeff00K + coeff01K * chi + coeff02K * chi2 + coeff03K * chi3 +
        coeff10K * eta + coeff11K * eta * chi + coeff12K * eta * chi2 +
        coeff13K * eta * chi3 + coeff20K * eta2 + coeff21K * eta2 * chi +
        coeff22K * eta2 * chi2 + coeff23K * eta2 * chi3 + coeff30K * eta3 +
        coeff31K * eta3 * chi + coeff32K * eta3 * chi2 + coeff33K * eta3 * chi3;
 
      // dSO
    adjParams->dSO =
        coeff00dSO + coeff01dSO * chi + coeff02dSO * chi2 + coeff03dSO * chi3 +
        coeff10dSO * eta + coeff11dSO * eta * chi + coeff12dSO * eta * chi2 +
        coeff13dSO * eta * chi3 + coeff20dSO * eta2 + coeff21dSO * eta2 * chi +
        coeff22dSO * eta2 * chi2 + coeff23dSO * eta2 * chi3 + coeff30dSO * eta3 +
        coeff31dSO * eta3 * chi + coeff32dSO * eta3 * chi2 + coeff33dSO * eta3 * chi3;

    // dSS
    adjParams->dSS =
        coeff00dSS + coeff01dSS * chi + coeff02dSS * chi2 + coeff03dSS * chi3 +
        coeff10dSS * eta + coeff11dSS * eta * chi + coeff12dSS * eta * chi2 +
        coeff13dSS * eta * chi3 + coeff20dSS * eta2 + coeff21dSS * eta2 * chi +
        coeff22dSS * eta2 * chi2 + coeff23dSS * eta2 * chi3 + coeff30dSS * eta3 +
        coeff31dSS * eta3 * chi + coeff32dSS * eta3 * chi2 + coeff33dSS * eta3 * chi3;

    adjParams->dtPeak = XLALSimIMREOBGetNRSpinPeakDeltaTv4 (2, 2, m1, m2, spin1z, spin2z);
    return CEV_SUCCESS;
}

static int CalculateSpinEOBHCoeffs (SpinEOBHCoeffs * coeffs,
					                const REAL8 eta,
					                const REAL8 a,
                                    AdjParams *adjParams)
{

    REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
    REAL8 m1PlusEtaKK;

    /* Constants are fits taken from Eq. 37 */
    /* needed to get the correct self-force results */

    static const REAL8 c0 = 1.4467;	
    static const REAL8 c1 = -1.7152360250654402;
    static const REAL8 c2 = -3.246255899738242;

    static const REAL8 c20 = 1.712;
    static const REAL8 c21 = -1.803949138004582;
    static const REAL8 c22 = -39.77229225266885;
    static const REAL8 c23 = 103.16588921239249;


    if (!coeffs)
    {
        return CEV_FAILURE;
    }


    coeffs->b3 = 0.;
    coeffs->bb3 = 0.;
    // K
    //coeffs->KK = KK = c0 + c1 * eta + c2 * eta * eta;
    coeffs->KK = KK = adjParams->KK;

    REAL8 chi = a / (1. - 2. * eta);
    REAL8 eta2 = eta * eta, eta3 = eta2 * eta;
    REAL8 chi2 = chi * chi, chi3 = chi2 * chi;


    m1PlusEtaKK = -1. + eta * KK;
    /* Eqs. 5.77 - 5.81 of BB1 */
    coeffs->k0 = k0 = KK * (m1PlusEtaKK - 1.);
    coeffs->k1 = k1 = -2. * (k0 + KK) * m1PlusEtaKK;
    k1p2 = k1 * k1;
    k1p3 = k1 * k1p2;
    coeffs->k2 = k2 =
        (k1 * (k1 - 4. * m1PlusEtaKK)) / 2. -
        a * a * k0 * m1PlusEtaKK * m1PlusEtaKK;
    coeffs->k3 = k3 =
        -k1 * k1 * k1 / 3. + k1 * k2 + k1 * k1 * m1PlusEtaKK - 2. * (k2 -
                                    m1PlusEtaKK)
        * m1PlusEtaKK - a * a * k1 * m1PlusEtaKK * m1PlusEtaKK;
    coeffs->k4 = k4 =
        (24. * k1 * k1 * k1 * k1 - 96. * k1 * k1 * k2 + 48. * k2 * k2 -
        64. * k1 * k1 * k1 * m1PlusEtaKK + 48. * a * a * (k1 * k1 -
                                2. * k2) *
        m1PlusEtaKK * m1PlusEtaKK + 96. * k1 * (k3 + 2. * k2 * m1PlusEtaKK) -
        m1PlusEtaKK * (192. * k3 +
                m1PlusEtaKK * (-3008. + 123. * CST_PI * CST_PI))) / 96.;
    coeffs->k5 = k5 = m1PlusEtaKK * m1PlusEtaKK
        * (-4237. / 60. + 128. / 5. * CST_GAMMA +
        2275. * CST_PI * CST_PI / 512. - 1. / 3. * a * a * (k1p3 -
                                    3. * k1 * k2 +
                                    3. * k3) -
        (k1p3 * k1p2 - 5. * k1p3 * k2 + 5. * k1 * k2 * k2 +
            5. * k1p2 * k3 - 5. * k2 * k3 -
            5. * k1 * k4) / 5. / m1PlusEtaKK / m1PlusEtaKK + (k1p2 * k1p2 -
                                    4. * k1p2 * k2 +
                                    2. * k2 * k2 +
                                    4. * k1 * k3 -
                                    4. * k4) / 2. /
        m1PlusEtaKK + 256. / 5. * log (2.) + (41. * CST_PI * CST_PI / 32. -
                            221. / 6.) * eta);
    coeffs->k5l = k5l = m1PlusEtaKK * m1PlusEtaKK * 64. / 5.;
    /*printf( "a = %.16e, k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e, k4 = %.16e, b3 = %.16e, bb3 = %.16e, KK = %.16e\n",
        a, coeffs->k0, coeffs->k1, coeffs->k2, coeffs->k3, coeffs->k4, coeffs->b3, coeffs->bb3, coeffs->KK );
    */

    /* Now calibrated parameters for spin models */
    // dSO
    coeffs->d1 = coeffs->d1v2 = 0.0;
    //coeffs->d1 = -69.5;
    coeffs->d1 = adjParams->dSO;
    // dSS
    coeffs->dheffSS = coeffs->dheffSSv2 = 0.0;
    //coeffs->dheffSS = 2.75;
    coeffs->dheffSS = adjParams->dSS;
    /*
    // dSO
    coeffs->d1v2 =
            coeff00dSO + coeff01dSO * chi + coeff02dSO * chi2 + coeff03dSO * chi3 +
            coeff10dSO * eta + coeff11dSO * eta * chi + coeff12dSO * eta * chi2 +
            coeff13dSO * eta * chi3 + coeff20dSO * eta2 + coeff21dSO * eta2 * chi +
            coeff22dSO * eta2 * chi2 + coeff23dSO * eta2 * chi3 + coeff30dSO * eta3 +
            coeff31dSO * eta3 * chi + coeff32dSO * eta3 * chi2 + coeff33dSO * eta3 * chi3;

    // dSS
    coeffs->dheffSSv2 =
            coeff00dSS + coeff01dSS * chi + coeff02dSS * chi2 + coeff03dSS * chi3 +
            coeff10dSS * eta + coeff11dSS * eta * chi + coeff12dSS * eta * chi2 +
            coeff13dSS * eta * chi3 + coeff20dSS * eta2 + coeff21dSS * eta2 * chi +
            coeff22dSS * eta2 * chi2 + coeff23dSS * eta2 * chi3 + coeff30dSS * eta3 +
            coeff31dSS * eta3 * chi + coeff32dSS * eta3 * chi2 + coeff33dSS * eta3 * chi3;
//          printf("dSO %.16e, dSS %.16e\n", coeffs->d1v2,coeffs->dheffSSv2);
    */
    return CEV_SUCCESS;
}

#define lMax 8
static int XLALSpinAlignedHcapDerivative(
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
  flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omega, params.params, H/Mtotal, lMax );

  /* Looking at the non-spinning model, I think we need to divide the flux by eta */
  flux = flux / eta;

  //printf( "Flux in derivatives function = %.16e\n", flux );

  /* Now we can calculate the final (spherical) derivatives */
  /* csi is needed because we use the tortoise co-ordinate */
  /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
  dvalues[0] = csi * tmpDValues[3];
  dvalues[1] = omega;
  /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
  dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
  dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux / omega;
  dvalues[3] = - flux / omega;

  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Values:\n%f %f %f %f\n", values[0], values[1], values[2], values[3] );

  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Derivatives:\n%f %f %f %f\n", dvalues[0], r*dvalues[1], dvalues[2], dvalues[3] );

  if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
  {
    //printf( "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
    return 1;
  }
  return CEV_SUCCESS;
}

#undef STEP_SIZE
#undef lMax