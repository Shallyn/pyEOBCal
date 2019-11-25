/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyHybridRingdown.h"
#include "dyBHRingdown.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

static INT XLALSimFindIndexMaxAmpli( UINT * indAmax, 
            REAL8Vector * timeVec, 
            REAL8Vector * ampWave, 
            REAL8 * valAmax, 
            REAL8 tofAmax ) {
    if ( indAmax == NULL || timeVec == NULL || ampWave == NULL || valAmax == NULL )
    {
        return CEV_FAILURE;
    }
    *indAmax = 0;
    INT found = 0;
    UINT i;
    for (i = 0; i < timeVec->length - 1; i++) 
    {
        if (timeVec->data[i] == tofAmax) 
        {
            found = 1;
            *indAmax = i;
            *valAmax = ampWave->data[i];
        }
    }
    if (found == 0) {
        return CEV_FAILURE;
    }
    else {
        return CEV_SUCCESS;
    }
}


/**
 * The main  function for performing the ringdown attachment for SEOBNRv4 (and beyond)
 * This is the function which gets called by the code generating the full IMR waveform once
 * generation of the inspiral part has been completed.
 * The ringdown is attached by factoring the less damped harmonics and apply tanh fit to the rest.
 * STEP 1) Get mass and spin of the final black hole and the complex ringdown frequencies
 * STEP 2) Apply the fit function from the attachment time
 * STEP 3) Construct the full RD by applying (factor) of less damped 220 mode
 * STEP 4) Constructing the RD stitched to the inspiral-merger
 * See Section 5.2 of https://dcc.ligo.org/T1600383
 */
 INT XLALSimIMREOBAttachFitRingdown(
    REAL8Vector * signal1, /**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
    REAL8Vector * signal2, /**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
    const REAL8 dt,        /**<< Sample time step (in seconds) */
    const REAL8 mass1,     /**<< First component mass (in Solar masses) */
    const REAL8 mass2,     /**<< Second component mass (in Solar masses) */
    REAL8 spin1z,    /**<<The spin of the first object;  */
    REAL8 spin2z,    /**<<The spin of the first object;  */
    REAL8Vector *timeVec, /**<< Vector containing the time values */
    REAL8Vector *matchrange /**<< Time values chosen as points for performing comb matching */
    ) 
{
    //INT4 debugSB = 0;
    UINT i;
    REAL8Vector *ampWave;
    REAL8Vector *phWave;
    REAL8Vector *rdtime;
    COMPLEX16Vector *modefreqs;
    //COMPLEX16Vector *modefreqs22; RC: we use this variable to compute the damping time of the 22 mode which will be used to set the lenght of the ringdown for all the modes

    ampWave = CreateREAL8Vector(signal1->length);
    phWave = CreateREAL8Vector(signal1->length);

    REAL8 mtot = mass1 + mass2;
    REAL8 eta = mass1 * mass2 / (mtot * mtot);

    /* Here I assume that the spins were properly projected (is precessing) and only spin1z, spin2z
     * are relevant, if not we need to add extra function to define what we call chi1, chi2 */
    REAL8 chi1 = spin1z;
    REAL8 chi2 = spin2z;

    gsl_spline *spline0 = NULL;
    gsl_interp_accel *acc0 = NULL;
/*
    if (debugSB) {
        printf("RDfit: we use spin1 %f, spin2 %f, and it should be dimensionless [-1,1] \n", chi1, chi2);
        printf("We use approximant = %d \n", appr);
    }*/
    REAL8 chi = 0.5 * (chi1 + chi2) + 0.5 * (chi1 - chi2) * (mass1 - mass2)/(mass1 + mass2) / (1.0 - 2.0 * eta);



    /*********************************************************************************************/
    /*                                          QNMs of the remnant                                        */
    /*********************************************************************************************/
    /* Getting  QNMs */
    modefreqs = CreateCOMPLEX16Vector(1);
    if (XLALSimIMREOBGenerateQNMFreqV2(modefreqs, mass1, mass2, chi1, chi2, 2, 2, 1) == CEV_FAILURE) {
        DestroyCOMPLEX16Vector(modefreqs);
        return CEV_FAILURE;
    }


    //RC: we use this variable to compute the damping time of the 22 mode which will be used to set the lenght of the ringdown for all the modes
    /*
    modefreqs22 = XLALCreateCOMPLEX16Vector(1);
    if (XLALSimIMREOBGenerateQNMFreqV2(modefreqs22, mass1, mass2, spin1, spin2, 2, 2, 1) == CEV_FAILURE) {
        DestroyCOMPLEX16Vector(modefreqs);
        return CEV_FAILURE;
    }*/
    /* Least-damped QNM */
    COMPLEX16 sigmalm0;
    sigmalm0 = (-cimag(modefreqs->data[0]) - I * creal(modefreqs->data[0])) * (mtot * CST_MTSUN_SI);
    // sigmalm0 = -0.08119703120243188 - I*0.5040714995216017;

/*
    if (debugSB) {
      printf("Mode %d %d \n",l, m);
        printf("matchpoints are: %.16f,  %.16f,   %.16f, %.16f\n", matchrange->data[0], matchrange->data[1], matchrange->data[2], matchrange->data[3]);
        printf("the 0-overtone is: %.16f + i %.16f \n", creal(sigmalm0), cimag(sigmalm0));
    }
*/


    /*********************************************************************************************/
    /*                                         Prepare inspiral signal                                         */
    /*                   This is needed to guarantee continuity with RD signal               */
    /*********************************************************************************************/
    /*  Compute amplitude and the unwrapped phase of the inspiral */
    for (i = 0; i < signal1->length; i++) 
    {
        ampWave->data[i] = 0.0;
        phWave->data[i] = 0.0;
    }
    REAL8 prev, now, corph;
    prev = atan2(signal2->data[0], signal1->data[0]);
    phWave->data[0] = prev;
    for (i = 0; i < timeVec->length; i++) 
    {
        ampWave->data[i] = sqrt(signal1->data[i] * signal1->data[i] + signal2->data[i] * signal2->data[i]);
        now = atan2(signal2->data[i], signal1->data[i]);
        if (i > 0) 
        {
            /* Unwrapping */
            corph = now - prev;
            corph = corph > CST_PI ? corph - CST_2PI : (corph < -CST_PI ? corph + CST_2PI : corph);

            phWave->data[i] = phWave->data[i - 1] + corph;
            prev = now;
        }
        //phWave->data[i] = now;
    }
    /*
    FILE *fout = NULL;
    if (debugSB) {
        char filename2[sizeof "CheckStasAmplPhaseXX.dat"];
        sprintf(filename2,"CheckStasAmplPhase%01d%01d.dat",l,m);
        fout = fopen(filename2, "w");
        printf("Check the length: time %d, signal %d \n", timeVec->length, signal1->length);
        for (i = 0; i < timeVec->length; i++) {
            fprintf(fout, "%f   %.16e   %.16e   %.16e   %.16e  \n", timeVec->data[i], ampWave->data[i], phWave->data[i], signal1->data[i], signal2->data[i]);
        }

        fclose(fout);
    }*/
    /* For the modes 22, 33, 21, 44 search for index at which the maximum of the amplitude of the 22 mode occurs */
    /* For the mode 55 search for the index corresponding at the time tpeak22-10M, where tpeak22 is the time of the maximum of the amplitude of the 22 mode*/
    /* See Eq.(4.3) of https://arxiv.org/pdf/1803.10701.pdf */

    REAL8 valAmax = ampWave->data[0];
    REAL8 valdAmax = 0;
    REAL8 tofAmax = matchrange->data[1];
    UINT indAmax;
    /*
    if (l == 5 && m==5) {
      tofAmax = matchrange->data[3];
    }*/
    //RC we should find indAmax only for the mode 22, for the other one the attaching point is the same, with the exception of the 55
    //For the 55 mode we need to find tpeak22 -10M, which in that case is fed as matchrange->data[3]
    if ( XLALSimFindIndexMaxAmpli( &indAmax, timeVec, ampWave, &valAmax, tofAmax ) == CEV_FAILURE )
    {
        print_warning("Time of maximum amplitude is not found .\n");
        DestroyCOMPLEX16Vector(modefreqs);
        return CEV_FAILURE;
    }
 /*
    if((l == 2 && m == 2) || (l == 5 && m == 5)){
      if ( XLALSimFindIndexMaxAmpli( &indAmax, timeVec, ampWave, &valAmax, tofAmax ) == CEV_FAILURE )
      {
          print_warning("Time of maximum amplitude is not found .\n");
          DestroyCOMPLEX16Vector(modefreqs);
          return CEV_FAILURE;
        }
       
      if (l==2 && m==2) {
        *indAmpMax = indAmax;
      }
      //RC for the 55 mode we also need to find the derivative at the attachment point, since we are not making the attachment
      //at the peak of the mode where the derivative is zero, see Eq.(4.3) of https://arxiv.org/pdf/1803.10701.pdf */
      /*
      if(l == 5 && m == 5){
        spline0 = gsl_spline_alloc (gsl_interp_cspline, timeVec->length);
        acc0 = gsl_interp_accel_alloc ();
        gsl_spline_init (spline0, timeVec->data, ampWave->data, timeVec->length);
        gsl_interp_accel_reset (acc0);
        valdAmax = gsl_spline_eval_deriv (spline0, tofAmax, acc0);
        }
      }
      else{
        indAmax = *indAmpMax;
        valAmax = ampWave->data[indAmax];
        //For the modes different from 22, this IS NOT the value of the maximum of the amplitude of the mode, but is the value of the amplitude at the peak of the 22 mode
        //For the SEOBNRv4HM model we need to compute also the derivative at the attaching point. This is not zero because we are attaching the HMs not at the paek of the lm mode, but at the peak of the 22 mode
        //See Eqs. (4.22-4.23) of https://arxiv.org/pdf/1803.10701.pdf 
        spline0 = gsl_spline_alloc (gsl_interp_cspline, timeVec->length);
        acc0 = gsl_interp_accel_alloc ();
        gsl_spline_init (spline0, timeVec->data, ampWave->data, timeVec->length);
        gsl_interp_accel_reset (acc0);
        valdAmax = gsl_spline_eval_deriv (spline0, tofAmax, acc0);
        // valAmax = 0.019835031225903733;
        // valdAmax = 0.00020163171965464758;
      }*/

/*
    if (debugSB) {
        printf("Check: The maximum of amplitude is %.16e found at t=%f, index = %d (out of %d) \n", valAmax, tofAmax, indAmax, timeVec->length);
        printf("Amp@Attachment = %.16f \n", valAmax);
        printf("dAmp@Attachment = %.16f \n", valdAmax);
        FILE *indMax = NULL;
        char filename2[sizeof "indexMaxXXHi.dat"];
        sprintf(filename2,"indexMax%01d%01dHi.dat",l,m);
        indMax = fopen (filename2, "w");
        fprintf(indMax, "%d \n",indAmax);
        fclose(indMax);
    }
*/


    /*********************************************************************************************/
    /*                  Constant coefficients entering RD fitting formulas                      */
    /*********************************************************************************************/
    /* Computing fit coefficients that enter amplitude and phase of the RD signal */
    /* Numeircal constants are defined at the top of this file */
    /* Eq. (54) of https://dcc.ligo.org/T1600383 for SEOBNRv4 and 22 mode SEOBNRv4HM*/
    /* Eqs. (C1, C3, C5, C7) of https://arxiv.org/pdf/1803.10701.pdf for SEOBNRv4HM*/

    REAL8 ampcf1;
    // ampcf1 = A1coeff00[l][m] + A1coeff01[l][m] * chi + A1coeff02[l][m] * chi * chi + A1coeff10[l][m] * eta + A1coeff11[l][m] * eta * chi + A1coeff20[l][m] * eta * eta;
    // ampcf1 = XLALCalculateRDAmplitudeCoefficient1(2, 2, eta, chi);
    ampcf1 = 0.0830664 - 0.0196758 * chi - 0.0136459 * chi * chi + 0.0612892 * eta +
        0.00146142 * eta * chi  - 0.0893454 * eta * eta;
    //ampcf1 = 0.08953902308302486;
    /*
    if(debugSB){
      printf("ampcf1 = %.16f \n", ampcf1);
    }*/

    /* Eq. (55) of https://dcc.ligo.org/T1600383 for SEOBNRv4 and 22 mode SEOBNRv4HM*/
    /* Eqs. (C2, C4, C6, C8) of https://arxiv.org/pdf/1803.10701.pdf for SEOBNRv4HM*/
    REAL8 ampcf2;
    // ampcf2 = A2coeff00[l][m] + A2coeff01[l][m] * chi + A2coeff10[l][m] * eta + A2coeff11[l][m] * eta * chi + A2coeff20[l][m] * eta * eta + A2coeff21[l][m] * eta * eta * chi;
    //ampcf2 = XLALCalculateRDAmplitudeCoefficient2(2, 2, eta, chi);
    ampcf2 = -0.623953  -0.371365 * chi + 1.39777 * eta +
        2.40203 * eta * chi -1.82173 * eta * eta -5.25339 * eta * eta * chi;
    //ampcf2 =  -0.506368869538912;
    /*
    if(debugSB){
      printf("ampcf2 = %.16f \n", ampcf2);
    }*/
//    printf("creal(sigma220), 2.*ampcf1*tanh(ampcf2) = %.16e %.16e\n",1./creal(sigma220), 2.*ampcf1*tanh(ampcf2));
    /* Eqs. (57)-(58) of https://dcc.ligo.org/T1600383 */
    if (creal(sigmalm0) > 2. * ampcf1 * tanh(ampcf2)) 
    {
        ampcf1 = creal(sigmalm0) / (2. * tanh(ampcf2));
    }
    /* Eq. (62) of https://dcc.ligo.org/T1600383 for SEOBNRv4 and 22 mode SEOBNRv4HM*/
    /* Eqs. (C9, C11, C13, C15) of https://arxiv.org/pdf/1803.10701.pdf for SEOBNRv4HM*/
    REAL8 phasecf1;
    // phasecf1 = P1coeff00[l][m] + P1coeff01[l][m] * chi + P1coeff02[l][m] * chi * chi + P1coeff10[l][m] * eta + P1coeff11[l][m] * eta * chi + P1coeff20[l][m] * eta * eta;
    //phasecf1 = XLALCalculateRDPhaseCoefficient1(2, 2, eta, chi);
    phasecf1 = 0.147584 + 0.00779176 * chi -0.0244358 * chi * chi + 0.263456 * eta -
        0.120853 * eta * chi -0.808987 * eta * eta; 
    /*
    if(debugSB){
      printf("phasecf1 = %.16f \n", phasecf1);
    }*/

    /* Eq. (63) of https://dcc.ligo.org/T1600383 for SEOBNRv4 and 22 mode SEOBNRv4HM*/
    /* Eqs. (C10, C12, C14, C16) of https://arxiv.org/pdf/1803.10701.pdf for SEOBNRv4HM*/
    REAL8 phasecf2;
    // phasecf2 = P2coeff00[l][m] + P2coeff01[l][m] * chi + P2coeff02[l][m] * chi * chi + P2coeff10[l][m] * eta + P2coeff11[l][m] * eta * chi + P2coeff20[l][m] * eta * eta;
    //phasecf2 = XLALCalculateRDPhaseCoefficient2(2, 2, eta, chi);
    phasecf2 = 2.46654 + 3.13067 * chi + 0.581626 * chi * chi -6.99396 * eta +
        9.61861 * eta * chi + 17.5646 * eta * eta;
    /*
    if(debugSB){
      printf("phasecf2 = %.16f \n", phasecf2);
    }*/

    /*********************************************************************************************/
    /*                                            RD fitting formulas                                           */
    /*********************************************************************************************/
    /* Ringdown signal length: 10 times the decay time of the n=0 mode */
    UINT Nrdwave = (INT) (10.0 / cimag(modefreqs->data[0]) / dt);
    //printf("Stas Nrdwave %d,  dt = %f", Nrdwave, dt);
    REAL8 dtM = dt / (mtot * CST_MTSUN_SI);     // go to geometric units
    rdtime = CreateREAL8Vector(Nrdwave);
    for (i = 0; i < Nrdwave; i++) 
    {
            rdtime->data[i] = i * dtM;      // this time for RD and it starts with 0
    }
    REAL8 tcons = 0.;

    /* Rescalings to guarantee continuity with inspiral */
    REAL8 Arescaledtcons = valAmax * exp(-creal(sigmalm0) * tcons) / eta;
    REAL8 dtArescaledtcons = (valdAmax - creal(sigmalm0) * valAmax) * exp(-creal(sigmalm0) * tcons) / eta;    // valdAmax = 0 - assumes extermum (max) of peak amplitude //RC: valdAmax = 0 for 22 mode
    /*Eq (4.22) of https://arxiv.org/pdf/1803.10701.pdf*/
    REAL8 ampcc1 = dtArescaledtcons * pow(cosh(ampcf2), 2) / ampcf1;
    /*Eq (4.23) of https://arxiv.org/pdf/1803.10701.pdf*/
    REAL8 ampcc2 = (Arescaledtcons * ampcf1 - dtArescaledtcons * cosh(ampcf2) * sinh(ampcf2)) / ampcf1;


    REAL8Vector *ampRD;
    REAL8Vector *phRD;
    ampRD = CreateREAL8Vector(Nrdwave);
    phRD = CreateREAL8Vector(Nrdwave);

    /* Construct RD amplitude */
    /* Eqs. (48)-(49) of https://dcc.ligo.org/T1600383 */
    for (i = 0; i < Nrdwave; i++) 
    {
        ampRD->data[i] = eta * exp(creal(sigmalm0) * rdtime->data[i]) * (ampcc1 * tanh(ampcf1 * rdtime->data[i] + ampcf2) + ampcc2);
    }
    /*
    if (debugSB) {
        char filename[sizeof "StasAmpRD_fullXX.dat"];
        sprintf(filename,"StasAmpRD_full%01d%01d.dat",l,m);
        fout = fopen (filename, "w");
        // for (i = 0; i < indAmax; i++) {
        //     fprintf(fout, "%.16e    %.16e \n", timeVec->data[i], ampWave->data[i]);
        // }
        for (i = 1; i < Nrdwave; i++) {
            fprintf(fout, "%.16e    %.16e \n", rdtime->data[i] + tofAmax*0, ampRD->data[i]);
        }
        fclose(fout);
    }*/

    /* Construct RD phase */
    REAL8 omegarescaledtcons = (phWave->data[indAmax] - phWave->data[indAmax - 1]) / dtM - cimag(sigmalm0);
    // omegarescaledtcons = -0.2354501785436996 - cimag(sigmalm0); 21 mode
    // omegarescaledtcons = -0.6076660197512114 - cimag(sigmalm0); 33 mode
    // omegarescaledtcons = -0.8063013198270685 - cimag(sigmalm0); 44 mode
    // omegarescaledtcons = -0.7813933465833801 - cimag(sigmalm0); 55 mode
    /*
    if(debugSB){
      printf("omega0 = %.16f\n", phWave->data[0]);
      printf("omega0fit = %.16f, omega0numder = %.16f \n", XLALSimIMREOBGetNRSpinPeakOmegaV4 (l, m, eta, chi), (phWave->data[indAmax] - phWave->data[indAmax - 1]) / dtM);
    }*/
    REAL8 phasecc1 = omegarescaledtcons * (phasecf2 + 1.0) / phasecf2 / phasecf1;
    REAL8 phi0 = phWave->data[indAmax];
    // phi0 = -165.37581501842033; 21 mode
    // phi0 = -486.98433862793974; 33 mode
    // phi0 = -651.2335864724689; 44 mode
    // phi0 = -805.2637568966843; 55 mode
    /*
    if(debugSB){
      printf("phi_0 = %.16f \n", phi0);
    }*/
    REAL8 logargnum, logargden;
    /* Eq. (59)-(60) of https://dcc.ligo.org/T1600383 */
    for (i = 0; i < Nrdwave; i++) 
    {
        logargnum = 1. + phasecf2 * exp(-phasecf1 * rdtime->data[i]);
        logargden = 1. + phasecf2;
        phRD->data[i] = phi0 - phasecc1 * log(logargnum / logargden) + cimag(sigmalm0) * rdtime->data[i];
    }

/*
    if (debugSB) {
        char filename[sizeof "StasPhRD_fullXX.dat"];
        sprintf(filename,"StasPhRD_full%01d%01d.dat",l,m);
        fout = fopen (filename, "w");
        // for (i = 0; i < indAmax; i++) {
        //     fprintf(fout, "%.16e    %.16e \n", timeVec->data[i], phWave->data[i]);
        // }
        for (i = 1; i < Nrdwave; i++) {
            fprintf(fout, "%.16e    %.16e \n", rdtime->data[i] + tofAmax*0, phRD->data[i]);
        }
        fclose(fout);

        // Compute the frequency of the ful signal
        UINT4 totSz = indAmax + Nrdwave;
        REAL8Vector *PhFull;
        PhFull = XLALCreateREAL8Vector(totSz);
        REAL8Vector *tFull;
        tFull = XLALCreateREAL8Vector(totSz);
        REAL8Vector *frFull;
        frFull = XLALCreateREAL8Vector(totSz);

        for (i = 0; i < indAmax; i++) {
            tFull->data[i] = timeVec->data[i];
            PhFull->data[i] = phWave->data[i];
        }
        for (i = 0; i < Nrdwave; i++) {
            tFull->data[i + indAmax] = rdtime->data[i] + tofAmax;
            PhFull->data[i + indAmax] = phRD->data[i];
        }
        fout = fopen("StasPhRD_full2.dat", "w");
        for (i = 0; i < totSz; i++) {
            fprintf(fout, "%.16e   %.16e \n", tFull->data[i], PhFull->data[i]);
        }
        fclose(fout);

        gsl_spline *spline = NULL;
        gsl_interp_accel *acc = NULL;
        spline = gsl_spline_alloc(gsl_interp_cspline, totSz);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, tFull->data, PhFull->data, totSz);
        for (i = 0; i < totSz; i++) {
            frFull->data[i] = gsl_spline_eval_deriv(spline, tFull->data[i], acc);
        }
        fout = fopen("StasFrRD_full.dat", "w");
        for (i = 0; i < totSz; i++) {
            fprintf(fout, "%.16e   %.16e \n", tFull->data[i], frFull->data[i]);
        }
        fclose(fout);

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

        XLALDestroyREAL8Vector(PhFull);
        XLALDestroyREAL8Vector(tFull);
        XLALDestroyREAL8Vector(frFull);

    }*/

    /* Construct RD signal */
    for (i = 0; i < Nrdwave; i++) {
        signal1->data[i + indAmax] = ampRD->data[i] * cos(phRD->data[i]);
        signal2->data[i + indAmax] = ampRD->data[i] * sin(phRD->data[i]);
    }

    gsl_spline_free (spline0);
    gsl_interp_accel_free (acc0);

    DestroyREAL8Vector(ampWave);
    DestroyREAL8Vector(phWave);
    DestroyCOMPLEX16Vector(modefreqs);
    DestroyREAL8Vector(rdtime);
    DestroyREAL8Vector(ampRD);
    DestroyREAL8Vector(phRD);

    return CEV_SUCCESS;

}

INT XLALSimIMREOBFinalMassSpin(REAL8 * finalMass,
                              /**<< OUTPUT, the final mass (scaled by original total mass) */
    REAL8 * finalSpin,        /**<< OUTPUT, the final spin (scaled by final mass) */
    const REAL8 mass1,        /**<< The mass of the 1st component of the system */
    const REAL8 mass2,        /**<< The mass of the 2nd component of the system */
    const REAL8 spin1z,     /**<< The spin of the 1st object; only needed for spin waveforms */
    const REAL8 spin2z     /**<< The spin of the 2nd object; only needed for spin waveforms */
    );

 static INT XLALGenerateHybridWaveDerivatives (
	REAL8Vector	*rwave,      /**<< OUTPUT, values of the waveform at comb points */
	REAL8Vector	*dwave,      /**<< OUTPUT, 1st deriv of the waveform at comb points */
	REAL8Vector	*ddwave,     /**<< OUTPUT, 2nd deriv of the waveform at comb points */
        REAL8Vector	*timeVec,    /**<< Vector containing the time */
	REAL8Vector	*wave,       /**<< Last part of inspiral waveform */
	REAL8Vector	*matchrange, /**<< Times which determine the size of the comb */
        REAL8           dt,          /**<< Sample time step */
        REAL8           mass1,       /**<< First component mass (in Solar masses) */
        REAL8           mass2        /**<< Second component mass (in Solar masses) */
	);

static INT XLALSimIMREOBHybridRingdownWave(
                                     REAL8Vector          *rdwave1,   /**<< OUTPUT, Real part of ringdown waveform */
                                     REAL8Vector          *rdwave2,   /**<< OUTPUT, Imag part of ringdown waveform */
                                     const REAL8           dt,        /**<< Sampling interval */
                                     const REAL8           mass1,     /**<< First component mass (in Solar masses) */
                                     const REAL8           mass2,     /**<< Second component mass (in Solar masses) */
                                     REAL8VectorSequence  *inspwave1, /**<< Values and derivs of real part inspiral waveform */
                                     REAL8VectorSequence  *inspwave2, /**<< Values and derivs of imag part inspiral waveform */
                                     COMPLEX16Vector      *modefreqs, /**<< Complex freqs of ringdown (scaled by total mass) */
                                     REAL8Vector          *matchrange /**<< Times which determine the comb of ringdown attachment */
);

static REAL8 GetNRSpinPeakOmega( INT l, INT  m, REAL8  eta, REAL8 a )
{
  /* Fit for HOMs missing */
  return 0.27581190323955274 + 0.19347381066059993*eta
       - 0.08898338208573725*log(1.0 - a/(1.0-2.0*eta))
       + eta*eta*(1.78832*(0.2690779744133912 + a/(2.0-4.0*eta))*(1.2056469070395925
       + a/(2.0-4.0*eta)) + 1.423734113371796*log(1.0 - a/(1.0-2.0*eta)));
}

/**
 * The main workhorse function for performing the ringdown attachment for EOB
 * models EOBNRv2 and SEOBNRv1. This is the function which gets called by the
 * code generating the full IMR waveform once generation of the inspiral part
 * has been completed.
 * The ringdown is attached using the hybrid comb matching detailed in
 * The method is describe in Sec. II C of Pan et al. PRD 84, 124052 (2011),
 * specifically Eqs. 30 - 32.. Further details of the
 * implementation of the found in the DCC document T1100433.
 * In SEOBNRv1, the last physical overtone is replace by a pseudoQNM. See
 * Taracchini et al. PRD 86, 024011 (2012) for details.
 * STEP 1) Get mass and spin of the final black hole and the complex ringdown frequencies
 * STEP 2) Based on least-damped-mode decay time, allocate memory for rigndown waveform
 * STEP 3) Get values and derivatives of inspiral waveforms at matching comb points
 * STEP 4) Solve QNM coefficients and generate ringdown waveforms
 * STEP 5) Stitch inspiral and ringdown waveoforms
 */
INT XLALSimIMREOBHybridAttachRingdown(
                                       REAL8Vector *signal1,    /**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
                                       REAL8Vector *signal2,    /**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
                                       const INT   l,          /**<< Current mode l */
                                       const INT   m,          /**<< Current mode m */
                                       const REAL8  dt,         /**<< Sample time step (in seconds) */
                                       const REAL8  mass1,      /**<< First component mass (in Solar masses) */
                                       const REAL8  mass2,      /**<< Second component mass (in Solar masses) */
                                       const REAL8  spin1z,     /**<<The spin of the first object; only needed for spin waveforms */
                                       const REAL8  spin2z,     /**<<The spin of the second object; only needed for spin waveforms */
                                       REAL8Vector *timeVec,    /**<< Vector containing the time values */
                                       REAL8Vector *matchrange /**<< Time values chosen as points for performing comb matching */
)
{
    
    COMPLEX16Vector *modefreqs;
    UINT Nrdwave;
    UINT j;
    INT failed = 0;
    
    UINT nmodes;
    REAL8Vector        *rdwave1;
    REAL8Vector        *rdwave2;
    REAL8Vector        *rinspwave;
    REAL8Vector        *dinspwave;
    REAL8Vector        *ddinspwave;
    REAL8VectorSequence    *inspwaves1;
    REAL8VectorSequence    *inspwaves2;
    REAL8 eta, a, NRPeakOmega22; /* To generate pQNM frequency */
    REAL8 mTot; /* In geometric units */
    REAL8 spin1[3] = { 0, 0, spin1z };
    REAL8 spin2[3] = { 0, 0, spin2z };
    REAL8 finalMass, finalSpin;
    
    mTot  = (mass1 + mass2) * CST_MTSUN_SI;
    eta       = mass1 * mass2 / ( (mass1 + mass2) * (mass1 + mass2) );
    
    /**
     * STEP 1) Get mass and spin of the final black hole and the complex ringdown frequencies
     */
    
    /* Create memory for the QNM frequencies */
    nmodes = 8;
    modefreqs = CreateCOMPLEX16Vector( nmodes );
    
    if ( XLALSimIMREOBGenerateQNMFreqV2( modefreqs, mass1, mass2, spin1z, spin2z, l, m, nmodes ) == CEV_FAILURE )
    {
        DestroyCOMPLEX16Vector( modefreqs );
        return CEV_FAILURE;
    }
    //printf("mode0 = %f + i%f\n", creal(modefreqs->data[0]), cimag(modefreqs->data[0]));
    
    /* Call XLALSimIMREOBFinalMassSpin_re() to get mass and spin of the final black hole */
    if ( XLALSimIMREOBFinalMassSpin(&finalMass, &finalSpin, mass1, mass2, spin1z, spin2z) == CEV_FAILURE )
    {
        DestroyCOMPLEX16Vector( modefreqs );
        return CEV_FAILURE;
    }
    
    //if ( approximant == SEOBNRv1 )
    //{
    /* Replace the last QNM with pQNM */
    /* We assume aligned/antialigned spins here */
    a  = (spin1[2] + spin2[2]) / 2. * (1.0 - 2.0 * eta) + (spin1[2] - spin2[2]) / 2. * (mass1 - mass2) / (mass1 + mass2);
    NRPeakOmega22 = GetNRSpinPeakOmega( l, m, eta, a ) / mTot;
    /*printf("a and NRomega in QNM freq: %.16e %.16e %.16e %.16e %.16e\n",spin1[2],spin2[2],
     mTot/LAL_MTSUN_SI,a,NRPeakOmega22*mTot);*/
    modefreqs->data[7] = (NRPeakOmega22/finalMass + creal(modefreqs->data[0])) / 2.;
    modefreqs->data[7] += I * 10./3. * cimag(modefreqs->data[0]);
    //}
    
    /*for (j = 0; j < nmodes; j++)
     {
     printf("QNM frequencies: %d %d %d %e %e\n",l,m,j,modefreqs->data[j].re*mTot,1./modefreqs->data[j].im/mTot);
     }*/
    
    /* Ringdown signal length: 10 times the decay time of the n=0 mode */
    Nrdwave = (INT4) (10.0 / cimag(modefreqs->data[0]) / dt);
    
    /* Check the value of attpos, to prevent memory access problems later */
    if ( matchrange->data[0] * mTot / dt < 5 || matchrange->data[1]*mTot/dt > matchrange->data[2] *mTot/dt - 2 )
    {
        print_warning("More inspiral points needed for ringdown matching.\n" );
        //printf("%.16e,%.16e,%.16e\n",matchrange->data[0] * mTot / dt, matchrange->data[1]*mTot/dt, matchrange->data[2] *mTot/dt - 2);
        DestroyCOMPLEX16Vector( modefreqs );
        return CEV_FAILURE;
    }
    
    /**
     * STEP 2) Based on least-damped-mode decay time, allocate memory for rigndown waveform
     */
    
    /* Create memory for the ring-down and full waveforms, and derivatives of inspirals */
    
    rdwave1 = CreateREAL8Vector( Nrdwave );
    rdwave2 = CreateREAL8Vector( Nrdwave );
    rinspwave = CreateREAL8Vector( 6 );
    dinspwave = CreateREAL8Vector( 6 );
    ddinspwave = CreateREAL8Vector( 6 );
    inspwaves1 = CreateREAL8VectorSequence( 3, 6 );
    inspwaves2 = CreateREAL8VectorSequence( 3, 6 );
    
    /* Check memory was allocated */
    if ( !rdwave1 || !rdwave2 || !rinspwave || !dinspwave
        || !ddinspwave || !inspwaves1 || !inspwaves2 )
    {
        DestroyCOMPLEX16Vector( modefreqs );
        if (rdwave1)    DestroyREAL8Vector( rdwave1 );
        if (rdwave2)    DestroyREAL8Vector( rdwave2 );
        if (rinspwave)  DestroyREAL8Vector( rinspwave );
        if (dinspwave)  DestroyREAL8Vector( dinspwave );
        if (ddinspwave) DestroyREAL8Vector( ddinspwave );
        if (inspwaves1) DestroyREAL8VectorSequence( inspwaves1 );
        if (inspwaves2) DestroyREAL8VectorSequence( inspwaves2 );
        return CEV_FAILURE;
    }
    
    memset( rdwave1->data, 0, rdwave1->length * sizeof( REAL8 ) );
    memset( rdwave2->data, 0, rdwave2->length * sizeof( REAL8 ) );
    
    /**
     * STEP 3) Get values and derivatives of inspiral waveforms at matching comb points
     */
    
    /* Generate derivatives of the last part of inspiral waves */
    /* Get derivatives of signal1 */
    if ( XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, timeVec, signal1,
                                           matchrange, dt, mass1, mass2 ) == CEV_FAILURE )
    {
        DestroyCOMPLEX16Vector( modefreqs );
        DestroyREAL8Vector( rdwave1 );
        DestroyREAL8Vector( rdwave2 );
        DestroyREAL8Vector( rinspwave );
        DestroyREAL8Vector( dinspwave );
        DestroyREAL8Vector( ddinspwave );
        DestroyREAL8VectorSequence( inspwaves1 );
        DestroyREAL8VectorSequence( inspwaves2 );
        return CEV_FAILURE;
    }
    for (j = 0; j < 6; j++)
    {
        inspwaves1->data[j] = rinspwave->data[j];
        inspwaves1->data[j + 6] = dinspwave->data[j];
        inspwaves1->data[j + 12] = ddinspwave->data[j];
    }
    
    /* Get derivatives of signal2 */
    if ( XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, timeVec, signal2,
                                           matchrange, dt, mass1, mass2 ) == CEV_FAILURE )
    {
        DestroyCOMPLEX16Vector( modefreqs );
        DestroyREAL8Vector( rdwave1 );
        DestroyREAL8Vector( rdwave2 );
        DestroyREAL8Vector( rinspwave );
        DestroyREAL8Vector( dinspwave );
        DestroyREAL8Vector( ddinspwave );
        DestroyREAL8VectorSequence( inspwaves1 );
        DestroyREAL8VectorSequence( inspwaves2 );
        return CEV_FAILURE;
    }
    for (j = 0; j < 6; j++)
    {
        inspwaves2->data[j] = rinspwave->data[j];
        inspwaves2->data[j + 6] = dinspwave->data[j];
        inspwaves2->data[j + 12] = ddinspwave->data[j];
    }
    
    
    /**
     * STEP 4) Solve QNM coefficients and generate ringdown waveforms
     */
    
    /* Generate ring-down waveforms */
    //XLALSimIMREOBHybridRingdownWave locates in the current file
    if ( XLALSimIMREOBHybridRingdownWave( rdwave1, rdwave2, dt, mass1, mass2, inspwaves1, inspwaves2,
                                         modefreqs, matchrange ) == CEV_FAILURE )
    {
        DestroyCOMPLEX16Vector( modefreqs );
        DestroyREAL8Vector( rdwave1 );
        DestroyREAL8Vector( rdwave2 );
        DestroyREAL8Vector( rinspwave );
        DestroyREAL8Vector( dinspwave );
        DestroyREAL8Vector( ddinspwave );
        DestroyREAL8VectorSequence( inspwaves1 );
        DestroyREAL8VectorSequence( inspwaves2 );
        return CEV_FAILURE;
    }
    
    /**
     * STEP 5) Stitch inspiral and ringdown waveoforms
     */
    
    /* Generate full waveforms, by stitching inspiral and ring-down waveforms */
    UINT attachIdx = (UINT)(matchrange->data[1] * mTot / dt);
    //printf("attachIdx = %d\n", attachIdx);
    //printf("len(signal1) = %d\n", signal1->length);
    //printf("Nrdwave = %d\n", Nrdwave);

    for (j = 1; j < Nrdwave; ++j)
    {
        signal1->data[j + attachIdx] = rdwave1->data[j];
        signal2->data[j + attachIdx] = rdwave2->data[j];
    }
    
    memset( signal1->data+Nrdwave+attachIdx, 0, (signal1->length - Nrdwave - attachIdx)*sizeof(REAL8) );
    memset( signal2->data+Nrdwave+attachIdx, 0, (signal2->length - Nrdwave - attachIdx)*sizeof(REAL8) );
    
    /* Free memory */
    DestroyCOMPLEX16Vector( modefreqs );
    DestroyREAL8Vector( rdwave1 );
    DestroyREAL8Vector( rdwave2 );
    DestroyREAL8Vector( rinspwave );
    DestroyREAL8Vector( dinspwave );
    DestroyREAL8Vector( ddinspwave );
    DestroyREAL8VectorSequence( inspwaves1 );
    DestroyREAL8VectorSequence( inspwaves2 );
    
    return CEV_SUCCESS;
}




/**
 * Function which calculates the value of the waveform, plus its
 * first and second derivatives, for the points which will be required
 * in the hybrid comb attachment of the ringdown.
 */
 static INT XLALGenerateHybridWaveDerivatives (
	REAL8Vector	*rwave,      /**<< OUTPUT, values of the waveform at comb points */
	REAL8Vector	*dwave,      /**<< OUTPUT, 1st deriv of the waveform at comb points */
	REAL8Vector	*ddwave,     /**<< OUTPUT, 2nd deriv of the waveform at comb points */
        REAL8Vector	*timeVec,    /**<< Vector containing the time */
	REAL8Vector	*wave,       /**<< Last part of inspiral waveform */
	REAL8Vector	*matchrange, /**<< Times which determine the size of the comb */
        REAL8           dt,          /**<< Sample time step */
        REAL8           mass1,       /**<< First component mass (in Solar masses) */
        REAL8           mass2        /**<< Second component mass (in Solar masses) */
	)
{

  /* XLAL error handling */
  INT errcode = CEV_SUCCESS;

  /* For checking GSL return codes */
  INT gslStatus;

  UINT j;
  UINT vecLength;
  REAL8 m;
  double *y;
  double ry, dy, dy2;
  double rt;
  double *tlist;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  /* Total mass in geometric units */
  m  = (mass1 + mass2) * CST_MTSUN_SI;

  tlist = (double *) malloc(6 * sizeof(double));
  rt = (matchrange->data[1] - matchrange->data[0]) / 5.;
  tlist[0] = matchrange->data[0];
  tlist[1] = tlist[0] + rt;
  tlist[2] = tlist[1] + rt;
  tlist[3] = tlist[2] + rt;
  tlist[4] = tlist[3] + rt;
  tlist[5] = matchrange->data[1];

  /* Set the length of the interpolation vectors */
  vecLength = (UINT)( m * matchrange->data[2] / dt ) + 1;

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initiate arrays and supporting variables for gsl */
  y = (double *) malloc(vecLength * sizeof(double));

  if ( !y )
  {
    return CEV_FAILURE;
  }

  for (j = 0; j < vecLength; ++j)
  {
	y[j] = wave->data[j];
  }


  //XLAL_CALLGSL( acc = (gsl_interp_accel*) gsl_interp_accel_alloc() );
  acc = (gsl_interp_accel*) gsl_interp_accel_alloc();
  //XLAL_CALLGSL( spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, vecLength) );
  spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, vecLength);
  if ( !acc || !spline )
  {
    if ( acc )    gsl_interp_accel_free(acc);
    if ( spline ) gsl_spline_free(spline);
    free( y );
    return CEV_FAILURE;
  }

  /* Gall gsl spline interpolation */
  gslStatus = gsl_spline_init(spline, timeVec->data, y, vecLength);
  if ( gslStatus != GSL_SUCCESS )
  { 
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free( y );
	
    return CEV_FAILURE;
  }

  /* Getting first and second order time derivatives from gsl interpolations */
  for (j = 0; j < 6; ++j)
  {
    gslStatus = gsl_spline_eval_e( spline, tlist[j], acc, &ry );
    if ( gslStatus == GSL_SUCCESS )
    {
      gslStatus = gsl_spline_eval_deriv_e(spline, tlist[j], acc, &dy );
      gslStatus = gsl_spline_eval_deriv2_e(spline, tlist[j], acc, &dy2 );
    }
    if (gslStatus != GSL_SUCCESS )
    {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      free( y );
        return CEV_FAILURE;
    }
    rwave->data[j]  = (REAL8)(ry);
    dwave->data[j]  = (REAL8)(dy/m);
    ddwave->data[j] = (REAL8)(dy2/m/m);

  }
  
  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  free( tlist );
  free(y);

  return errcode;
}


/**
 * Generates the ringdown wave associated with the given real
 * and imaginary parts of the inspiral waveform. The parameters of
 * the ringdown, such as amplitude and phase offsets, are determined
 * by solving the linear equations defined in the DCC document T1100433.
 * In the linear equations Ax=y,
 * A is a 16-by-16 matrix depending on QNM (complex) frequencies,
 * x is a 16-d vector of the 8 unknown complex QNM amplitudes,
 * y is a 16-d vector depending on inspiral-plunge waveforms and their derivatives near merger.
 */
static INT XLALSimIMREOBHybridRingdownWave(
                                     REAL8Vector          *rdwave1,   /**<< OUTPUT, Real part of ringdown waveform */
                                     REAL8Vector          *rdwave2,   /**<< OUTPUT, Imag part of ringdown waveform */
                                     const REAL8           dt,        /**<< Sampling interval */
                                     const REAL8           mass1,     /**<< First component mass (in Solar masses) */
                                     const REAL8           mass2,     /**<< Second component mass (in Solar masses) */
                                     REAL8VectorSequence  *inspwave1, /**<< Values and derivs of real part inspiral waveform */
                                     REAL8VectorSequence  *inspwave2, /**<< Values and derivs of imag part inspiral waveform */
                                     COMPLEX16Vector      *modefreqs, /**<< Complex freqs of ringdown (scaled by total mass) */
                                     REAL8Vector          *matchrange /**<< Times which determine the comb of ringdown attachment */
)
{
#if debugOutput
    fprintf(stderr,"matchrange = %f, %f, %f\n",matchrange->data[0],matchrange->data[1],matchrange->data[2]);
#endif
    /* XLAL error handling */
    INT errcode = CEV_SUCCESS;
    
    /* For checking GSL return codes */
    INT gslStatus;
    
    UINT i, j, k, nmodes = 8;
    
    /* Sampling rate from input */
    REAL8 t1, t2, t3, t4, t5, rt;
    gsl_matrix *coef;
    gsl_vector *hderivs;
    gsl_vector *x;
    gsl_permutation *p;
    REAL8Vector *modeamps;
    int s;
    REAL8 tj;
    REAL8 m;
    
    /* mass in geometric units */
    m  = (mass1 + mass2) * CST_MTSUN_SI;
    t5 = (matchrange->data[0] - matchrange->data[1]) * m;
    rt = -t5 / 5.;
    
    t4 = t5 + rt;
    t3 = t4 + rt;
    t2 = t3 + rt;
    t1 = t2 + rt;
    
    if ( inspwave1->length != 3 || inspwave2->length != 3 ||
        modefreqs->length != nmodes )
    {
        return CEV_FAILURE;
    }
    
    /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
    /* Initiate matrices and supporting variables */
    //XLAL_CALLGSL( coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes) );
    coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes);
    //XLAL_CALLGSL( hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
    hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes);
    //XLAL_CALLGSL( x = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
    x = (gsl_vector *) gsl_vector_alloc(2 * nmodes);
    //XLAL_CALLGSL( p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes) );
    p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes);
    
    /* Check all matrices and variables were allocated */
    if ( !coef || !hderivs || !x || !p )
    {
        if (coef)    gsl_matrix_free(coef);
        if (hderivs) gsl_vector_free(hderivs);
        if (x)       gsl_vector_free(x);
        if (p)       gsl_permutation_free(p);
        return CEV_FAILURE;
    }
    
    /* Define the linear system Ax=y */
    /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
    /* The half of A defined here corresponds to matrices M1 and -M2 in the DCC document T1100433 */
    /* Define y here as "hderivs" */
    for (i = 0; i < nmodes; ++i)
    {
        gsl_matrix_set(coef, 0, i, 1);
        gsl_matrix_set(coef, 1, i, - cimag(modefreqs->data[i]));
        gsl_matrix_set(coef, 2, i, exp(-cimag(modefreqs->data[i])*t1) * cos(creal(modefreqs->data[i])*t1));
        gsl_matrix_set(coef, 3, i, exp(-cimag(modefreqs->data[i])*t2) * cos(creal(modefreqs->data[i])*t2));
        gsl_matrix_set(coef, 4, i, exp(-cimag(modefreqs->data[i])*t3) * cos(creal(modefreqs->data[i])*t3));
        gsl_matrix_set(coef, 5, i, exp(-cimag(modefreqs->data[i])*t4) * cos(creal(modefreqs->data[i])*t4));
        gsl_matrix_set(coef, 6, i, exp(-cimag(modefreqs->data[i])*t5) * cos(creal(modefreqs->data[i])*t5));
        gsl_matrix_set(coef, 7, i, exp(-cimag(modefreqs->data[i])*t5) *
                       (-cimag(modefreqs->data[i]) * cos(creal(modefreqs->data[i])*t5)
                        -creal(modefreqs->data[i]) * sin(creal(modefreqs->data[i])*t5)));
        gsl_matrix_set(coef, 8, i, 0);
        gsl_matrix_set(coef, 9, i, - creal(modefreqs->data[i]));
        gsl_matrix_set(coef, 10, i, -exp(-cimag(modefreqs->data[i])*t1) * sin(creal(modefreqs->data[i])*t1));
        gsl_matrix_set(coef, 11, i, -exp(-cimag(modefreqs->data[i])*t2) * sin(creal(modefreqs->data[i])*t2));
        gsl_matrix_set(coef, 12, i, -exp(-cimag(modefreqs->data[i])*t3) * sin(creal(modefreqs->data[i])*t3));
        gsl_matrix_set(coef, 13, i, -exp(-cimag(modefreqs->data[i])*t4) * sin(creal(modefreqs->data[i])*t4));
        gsl_matrix_set(coef, 14, i, -exp(-cimag(modefreqs->data[i])*t5) * sin(creal(modefreqs->data[i])*t5));
        gsl_matrix_set(coef, 15, i, exp(-cimag(modefreqs->data[i])*t5) *
                       (cimag(modefreqs->data[i]) * sin(creal(modefreqs->data[i])*t5)
                        -creal(modefreqs->data[i]) * cos(creal(modefreqs->data[i])*t5)));
    }
    for (i = 0; i < 2; ++i)
    {
        gsl_vector_set(hderivs, i, inspwave1->data[(i + 1) * inspwave1->vectorLength - 1]);
        gsl_vector_set(hderivs, i + nmodes, inspwave2->data[(i + 1) * inspwave2->vectorLength - 1]);
        gsl_vector_set(hderivs, i + 6, inspwave1->data[i * inspwave1->vectorLength]);
        gsl_vector_set(hderivs, i + 6 + nmodes, inspwave2->data[i * inspwave2->vectorLength]);
    }
    gsl_vector_set(hderivs, 2, inspwave1->data[4]);
    gsl_vector_set(hderivs, 2 + nmodes, inspwave2->data[4]);
    gsl_vector_set(hderivs, 3, inspwave1->data[3]);
    gsl_vector_set(hderivs, 3 + nmodes, inspwave2->data[3]);
    gsl_vector_set(hderivs, 4, inspwave1->data[2]);
    gsl_vector_set(hderivs, 4 + nmodes, inspwave2->data[2]);
    gsl_vector_set(hderivs, 5, inspwave1->data[1]);
    gsl_vector_set(hderivs, 5 + nmodes, inspwave2->data[1]);
    
    /* Complete the definition for the rest half of A */
    for (i = 0; i < nmodes; ++i)
    {
        for (k = 0; k < nmodes; ++k)
        {
            gsl_matrix_set(coef, i, k + nmodes, - gsl_matrix_get(coef, i + nmodes, k));
            gsl_matrix_set(coef, i + nmodes, k + nmodes, gsl_matrix_get(coef, i, k));
        }
    }
    
#if 0
    /* print ringdown-matching linear system: coefficient matrix and RHS vector */
    fprintf(stderr,"\nRingdown matching matrix:\n");
    for (i = 0; i < 16; ++i)
    {
        for (j = 0; j < 16; ++j)
        {
            fprintf(stderr,"%.12e ",gsl_matrix_get(coef,i,j));
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"RHS:  ");
    for (i = 0; i < 16; ++i)
    {
        fprintf(stderr,"%.12e   ",gsl_vector_get(hderivs,i));
    }
    fprintf(stderr,"\n");
#endif
    
    /* Call gsl LU decomposition to solve the linear system */
    //XLAL_CALLGSL( gslStatus = gsl_linalg_LU_decomp(coef, p, &s) );
    gslStatus = gsl_linalg_LU_decomp(coef, p, &s) ;
    if ( gslStatus == GSL_SUCCESS )
    {
        //XLAL_CALLGSL( gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x) );
        gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x);
    }
    if ( gslStatus != GSL_SUCCESS )
    {
        gsl_matrix_free(coef);
        gsl_vector_free(hderivs);
        gsl_vector_free(x);
        gsl_permutation_free(p);
        
        return CEV_FAILURE;
    }
    
    /* Putting solution to an XLAL vector */
    modeamps = CreateREAL8Vector(2 * nmodes);
    
    if ( !modeamps )
    {
        gsl_matrix_free(coef);
        gsl_vector_free(hderivs);
        gsl_vector_free(x);
        gsl_permutation_free(p);
        
        return CEV_FAILURE;
    }
    
    for (i = 0; i < nmodes; ++i)
    {
        modeamps->data[i] = gsl_vector_get(x, i);
        modeamps->data[i + nmodes] = gsl_vector_get(x, i + nmodes);
    }
#if debugOutput
    fprintf(stderr,"using 20Msun to do scaling.\n");

    for (i = 0; i < nmodes; ++i)
    {
        fprintf(stderr,"%d-th amps = %e+i %e, freq = %e+i %e\n",i,modeamps->data[i],
               modeamps->data[i + nmodes],creal(modefreqs->data[i])*20*4.92549095e-6,cimag(modefreqs->data[i])*20*4.92549095e-6);
    }
    
    fprintf(stderr,"using Mtotal to do scaling.\n");
    for (i = 0; i < nmodes; ++i)
    {
        fprintf(stderr,"%d-th amps = %e+i %e, freq = %e+i %e\n",i,modeamps->data[i],
               modeamps->data[i + nmodes],creal(modefreqs->data[i])*m,cimag(modefreqs->data[i])*m);
    }
#endif
    /* Free all gsl linear algebra objects */
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    
    /* Build ring-down waveforms */
    
    REAL8 timeOffset = fmod( matchrange->data[1], dt/m) * dt;
    
    for (j = 0; j < rdwave1->length; ++j)
    {
        tj = j * dt - timeOffset;
        rdwave1->data[j] = 0;
        rdwave2->data[j] = 0;
        for (i = 0; i < nmodes; ++i)
        {
            rdwave1->data[j] += exp(- tj * cimag(modefreqs->data[i]))
            * ( modeamps->data[i] * cos(tj * creal(modefreqs->data[i]))
               +   modeamps->data[i + nmodes] * sin(tj * creal(modefreqs->data[i])) );
            rdwave2->data[j] += exp(- tj * cimag(modefreqs->data[i]))
            * (- modeamps->data[i] * sin(tj * creal(modefreqs->data[i]))
               +   modeamps->data[i + nmodes] * cos(tj * creal(modefreqs->data[i])) );
        }
    }
    
    DestroyREAL8Vector(modeamps);
    return errcode;
}
