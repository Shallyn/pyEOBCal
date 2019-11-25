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
        dyEOB = NULL;
    }
    return;
}
