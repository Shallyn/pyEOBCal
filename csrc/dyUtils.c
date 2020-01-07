/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyUtils.h"

INT CreateSpinEOBParams(REAL8 m1,
                        REAL8 m2,
                        REAL8 spin1x,
                        REAL8 spin1y,
                        REAL8 spin1z,
                        REAL8 spin2x,
                        REAL8 spin2y,
                        REAL8 spin2z,
                        REAL8 eccentricity,
                        INT tortoise,
                        SpinEOBParams *out)
{
    INT i;
    REAL8 Mtotal = m1 + m2;
    REAL8 eta = m1 * m2 / Mtotal / Mtotal;

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

    seobParams.deltaT = 1.;
    seobParams.eccentricity = eccentricity;
    seobParams.s1Vec = s1Vec;
    seobParams.s2Vec = s2Vec;
    seobParams.s1VecOverMtMt = s1VecOverMtMt;
    seobParams.s2VecOverMtMt = s2VecOverMtMt;
    seobParams.sigmaStar = sigmaStar;
    seobParams.sigmaKerr = sigmaKerr;
    seobParams.eobParams = &eobParams;
    seobParams.nqcCoeffs = &nqcCoeffs;
    seobParams.seobCoeffs = &seobCoeffs;
    seobParams.tortoise = tortoise;

    eobParams.hCoeffs = &hCoeffs;
    eobParams.prefixes = &prefixes;
    eobParams.m1 = m1;
    eobParams.m2 = m2;
    eobParams.eta = eta;
    eobParams.Mtotal = Mtotal;

    s1Vec->data[0] = spin1z;
    s1Vec->data[1] = spin1z;
    s1Vec->data[2] = spin1z;

    s2Vec->data[0] = spin2x;
    s2Vec->data[1] = spin2y;
    s2Vec->data[2] = spin2z;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] * m1 * m1 / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] * m2 * m2 / Mtotal / Mtotal;
        sigmaStar->data[i] = ( s1Vec->data[i] + s2Vec->data[i]) * eta;
        sigmaKerr->data[i] = ( m1*m1*s1Vec->data[i] + m2*m2*s2Vec->data[i]) / (Mtotal * Mtotal);
    }

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = spin1z;
    seobParams.chi2 = spin2z;
    *out = seobParams;

   return CEV_SUCCESS; 
}

void DestroySpinEOBParams(SpinEOBParams *seobParams)
{
    DestroyREAL8Vector(seobParams->s1Vec);
    seobParams->s1Vec = NULL;
    DestroyREAL8Vector(seobParams->s2Vec);
    seobParams->s2Vec = NULL;
    DestroyREAL8Vector(seobParams->sigmaStar);
    seobParams->sigmaStar = NULL;
    DestroyREAL8Vector(seobParams->sigmaKerr);
    seobParams->sigmaKerr = NULL;
    return;
}

void CalculateSigmaStar(REAL8Vector *sigmaStar,
                        REAL8 mass1,
                        REAL8 mass2,
                        REAL8Vector *s1,
                        REAL8Vector *s2)
{
    UINT i;
    REAL8 totalMass = mass1 + mass2;
    for (i = 0; i < 3; i++)
    {
        sigmaStar->data[i] = (mass2/mass1 * s1->data[i] + mass1/mass2 * s2->data[i]) / (totalMass * totalMass);
    }
}

void CalculateSigmaKerr(REAL8Vector *sigmaKerr,
                        REAL8 mass1,
                        REAL8 mass2,
                        REAL8Vector *s1,
                        REAL8Vector *s2)
{
    UINT i;
    REAL8 totalMass = mass1 + mass2;
    
    for (i=0; i<3; i++)
    {
        sigmaKerr->data[i] = (s1->data[i] + s2->data[i]) / (totalMass * totalMass);
    }
}

SpinEOBDynamics *SpinEOBDynamicsInit(UINT length)
{
    SpinEOBDynamics *dyout;
    dyout = (SpinEOBDynamics*)malloc(sizeof(*dyout));
    dyout->phiVec = CreateREAL8Vector(length);
    dyout->dphiVec = CreateREAL8Vector(length);
    dyout->pPhiVec = CreateREAL8Vector(length);
    dyout->dpPhiVec = CreateREAL8Vector(length);
    dyout->prVec = CreateREAL8Vector(length);
    dyout->dprVec = CreateREAL8Vector(length);
    dyout->rVec = CreateREAL8Vector(length);
    dyout->drVec = CreateREAL8Vector(length);
    dyout->tVec = CreateREAL8Vector(length);
    dyout->length = length;
    return dyout;
}

void DestroySpinEOBDynamics(SpinEOBDynamics *dyEOB)
{
    if (!dyEOB)
        return;
    
    if(dyEOB->tVec)
        DestroyREAL8Vector(dyEOB->tVec);
    if(dyEOB->rVec)
        DestroyREAL8Vector(dyEOB->rVec);
    if(dyEOB->prVec)
        DestroyREAL8Vector(dyEOB->prVec);
    if(dyEOB->phiVec)
        DestroyREAL8Vector(dyEOB->phiVec);
    if(dyEOB->pPhiVec)
        DestroyREAL8Vector(dyEOB->pPhiVec);
    if(dyEOB->drVec)
        DestroyREAL8Vector(dyEOB->drVec);
    if(dyEOB->dprVec)
        DestroyREAL8Vector(dyEOB->dprVec);
    if(dyEOB->dphiVec)
        DestroyREAL8Vector(dyEOB->dphiVec);
    if(dyEOB->dpPhiVec)
        DestroyREAL8Vector(dyEOB->dpPhiVec);

    if(dyEOB)
    {
        dyEOB->length = 0;
        free(dyEOB);
    }
    return;
}
