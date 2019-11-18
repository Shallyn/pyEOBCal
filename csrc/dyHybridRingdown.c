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

