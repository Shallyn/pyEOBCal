/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pyUtils.h"
#include "../csrc/dyFactorizedWaveform.h"
#include "../csrc/dyHamiltonian.h"
#include "../csrc/pnWaveform.h"
#include "../csrc/auxFunctions.h"

#define STEP_SIZE 1.0e-4

int CalculateSpinEOBHCoeffs (SpinEOBHCoeffs * coeffs,
					                const REAL8 eta,
					                const REAL8 a,
                                    AdjParams *adjParams);

INT ComputeNewtonMultipolePrefixes(NewtonMultipolePrefixes *prefix,
                                   const REAL8 m1 ,
                                   const REAL8 m2);

REAL8
XLALInspiralSpinFactorizedFlux (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
				const REAL8 omega,	/**< orbital frequency */
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				INT const  lMax,
				INT allow_ecc);


static PyObject *Func_SpinEOBDynamicsGeneration(PyObject *self, PyObject *args)
{
    INT status,i;
    REAL8 m1, m2, s1z, s2z;
    REAL8 KK, dSS, dSO;
    PyObject *py_r, *py_pr, *py_phi, *py_pphi;
    if(!PyArg_ParseTuple(args, "dddddddOOOO", &m1, &m2, &s1z, &s2z, &KK, &dSS, &dSO, &py_r, &py_pr, &py_phi, &py_pphi))
    {
        return NULL;
    }
    if(fabs(s1z)>0.999 || fabs(s2z)>0.999)
    {
        return NULL;
    }
    if(m1 < m2)
    {
        SWAP(m1, m2);
    }
    
    MACRO_PyArray2REAL8Vector(py_r, rVec);
    MACRO_PyArray2REAL8Vector(py_pr, prVec);
    MACRO_PyArray2REAL8Vector(py_phi, phiVec);
    MACRO_PyArray2REAL8Vector(py_pphi, pphiVec);

    INT length = GET_MIN(GET_MIN(rVec->length, prVec->length), GET_MIN(phiVec->length, pphiVec->length));
    if(length < 1)
    {
        DestroyREAL8Vector(rVec);
        DestroyREAL8Vector(prVec);
        DestroyREAL8Vector(phiVec);
        DestroyREAL8Vector(pphiVec);
        return NULL;
    }
    AdjParams adjParams;
    adjParams.dSO = dSO;
    adjParams.dSS = dSS;
    adjParams.KK = KK;
    adjParams.dtPeak = 0;

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
    memset (&nqcCoeffs, 0, sizeof(nqcCoeffs));

    /* Allocate spin parameters */
    REAL8Vector *sigmaStar = CreateREAL8Vector(3);
    REAL8Vector *sigmaKerr = CreateREAL8Vector(3);
    REAL8Vector *s1VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s2VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s1Vec = CreateREAL8Vector(3);
    REAL8Vector *s2Vec = CreateREAL8Vector(3);

    seobParams.deltaT = 1.;
    seobParams.eccentricity = 0;
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

    s1Vec->data[0] = 0;
    s1Vec->data[1] = 0;
    s1Vec->data[2] = s1z;

    s2Vec->data[0] = 0;
    s2Vec->data[1] = 0;
    s2Vec->data[2] = s2z;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] * m1 * m1 / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] * m2 * m2 / Mtotal / Mtotal;
        sigmaStar->data[i] = ( s1Vec->data[i] + s2Vec->data[i]) * eta;
        sigmaKerr->data[i] = ( m1*m1*s1Vec->data[i] + m2*m2*s2Vec->data[i]) / (Mtotal * Mtotal);
    }

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = s1z;
    seobParams.chi2 = s2z;

    REAL8 chiS, chiA, tplspin;
    chiS = 0.5 * (s1z + s2z);
    chiA = 0.5 * (s1z - s2z);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;

    status = CalculateSpinEOBHCoeffs(&seobCoeffs, eta, sigmaKerr->data[2], &adjParams);
    status = CalculateSpinFactorizedWaveformCoefficients(&hCoeffs, &seobParams, m1, m2, eta, tplspin, chiS, chiA);
    status = ComputeNewtonMultipolePrefixes(&prefixes, m1, m2);

    REAL8 r, pr, phi, pphi;
    REAL8Vector *dyValues = CreateREAL8Vector(4);
    REAL8Vector *hamVec = CreateREAL8Vector(length);
    REAL8Vector *omegaVec = CreateREAL8Vector(length);
    REAL8Vector *hrealVec = CreateREAL8Vector(length);
    REAL8Vector *himagVec = CreateREAL8Vector(length);
    REAL8Vector *fluxVec = CreateREAL8Vector(length);
    
    REAL8Vector cartPosVec, cartMomVec;
    REAL8 cartPosData[3] = {0,0,0};
    REAL8 cartMomData[3] = {0,0,0};
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;
    COMPLEX16 hLM;
    REAL8 ham, omega, v;

    for(i=0; i<length; i++)
    {
        r = rVec->data[i];
        pr = prVec->data[i];
        phi = phiVec->data[i];
        pphi = pphiVec->data[i];

        cartPosData[0] = r;
        cartMomData[0] = pr;
        cartMomData[1] = pphi/r;

        ham = hamVec->data[i] = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec, 
            s1VecOverMtMt, s2VecOverMtMt, sigmaKerr, sigmaStar, 
            1, &seobCoeffs);
        dyValues->data[0] = r;
        dyValues->data[1] = phi;
        dyValues->data[2] = pr;
        dyValues->data[3] = pphi;
        
        omega = omegaVec->data[i] = XLALSimIMRSpinAlignedEOBCalcOmega(dyValues->data, &seobParams, STEP_SIZE);

        fluxVec->data[i]  = XLALInspiralSpinFactorizedFlux( dyValues, &nqcCoeffs, omega, &seobParams, ham, 8 ,0);
        v = cbrt(omega);
        status = XLALSimIMRSpinEOBGetSpinFactorizedWaveform(&hLM, dyValues, v, ham, 2, 2, &seobParams);
        hrealVec->data[i] = creal(hLM);
        himagVec->data[i] = cimag(hLM);

        if(status!=CEV_SUCCESS)
        {
            return NULL;
        }
    }

    DestroyREAL8Vector(s1Vec);
    DestroyREAL8Vector(s2Vec);
    DestroyREAL8Vector(s1VecOverMtMt);
    DestroyREAL8Vector(s2VecOverMtMt);
    DestroyREAL8Vector(sigmaStar);
    DestroyREAL8Vector(sigmaKerr);
    DestroyREAL8Vector(dyValues);

    DestroyREAL8Vector(rVec);
    DestroyREAL8Vector(prVec);
    DestroyREAL8Vector(phiVec);
    DestroyREAL8Vector(pphiVec);

    MACRO_REAL8Vector2PyArray(hamVec, py_ham);
    MACRO_REAL8Vector2PyArray(omegaVec, py_omega);
    MACRO_REAL8Vector2PyArray(hrealVec, py_hreal);
    MACRO_REAL8Vector2PyArray(himagVec, py_himag);
    MACRO_REAL8Vector2PyArray(fluxVec, py_flux);

    PyObject *out;
    out = Py_BuildValue("OOOOO", py_ham, py_omega, py_hreal, py_himag, py_flux);
    return out;
}


static PyObject *Func_PNwaveformCQG_25_165003(PyObject *self, PyObject *args)
{
    INT status,i;
    REAL8 m1, m2, s1z, s2z;
    REAL8 KK, dSS, dSO;
    PyObject *py_r, *py_pr, *py_phi, *py_pphi;
    if(!PyArg_ParseTuple(args, "dddddddOOOO", &m1, &m2, &s1z, &s2z, &KK, &dSS, &dSO, &py_r, &py_pr, &py_phi, &py_pphi))
    {
        return NULL;
    }
    if(fabs(s1z)>0.999 || fabs(s2z)>0.999)
    {
        return NULL;
    }
    if(m1 < m2)
    {
        SWAP(m1, m2);
    }
    MACRO_PyArray2REAL8Vector(py_r, rVec);
    MACRO_PyArray2REAL8Vector(py_pr, prVec);
    MACRO_PyArray2REAL8Vector(py_phi, phiVec);
    MACRO_PyArray2REAL8Vector(py_pphi, pphiVec);
    INT length = GET_MIN(GET_MIN(rVec->length, prVec->length), GET_MIN(phiVec->length, pphiVec->length));
    if(length < 1)
    {
        DestroyREAL8Vector(rVec);
        DestroyREAL8Vector(prVec);
        DestroyREAL8Vector(phiVec);
        DestroyREAL8Vector(pphiVec);
        return NULL;
    }
    AdjParams adjParams;
    adjParams.dSO = dSO;
    adjParams.dSS = dSS;
    adjParams.KK = KK;
    adjParams.dtPeak = 0;

    REAL8 Mtotal = m1 + m2;
    REAL8 eta = m1 * m2 / Mtotal / Mtotal;

    SpinEOBParams seobParams;
    EOBParams eobParams;
    EOBNonQCCoeffs nqcCoeffs;
    SpinEOBHCoeffs seobCoeffs;

    memset (&seobParams, 0, sizeof (seobParams));
    memset (&seobCoeffs, 0, sizeof (seobCoeffs));
    memset (&eobParams, 0, sizeof (eobParams));
    memset (&nqcCoeffs, 0, sizeof(nqcCoeffs));

    WaveformCoeffsInit(hCoeffs);
    /* Allocate spin parameters */
    REAL8Vector *sigmaStar = CreateREAL8Vector(3);
    REAL8Vector *sigmaKerr = CreateREAL8Vector(3);
    REAL8Vector *s1VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s2VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s1Vec = CreateREAL8Vector(3);
    REAL8Vector *s2Vec = CreateREAL8Vector(3);

    seobParams.deltaT = 1.;
    seobParams.eccentricity = 0;
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

    eobParams.hCoeffs = NULL;
    eobParams.prefixes = NULL;
    eobParams.m1 = m1;
    eobParams.m2 = m2;
    eobParams.eta = eta;
    eobParams.Mtotal = Mtotal;

    s1Vec->data[0] = 0;
    s1Vec->data[1] = 0;
    s1Vec->data[2] = s1z;

    s2Vec->data[0] = 0;
    s2Vec->data[1] = 0;
    s2Vec->data[2] = s2z;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] * m1 * m1 / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] * m2 * m2 / Mtotal / Mtotal;
        sigmaStar->data[i] = ( s1Vec->data[i] + s2Vec->data[i]) * eta;
        sigmaKerr->data[i] = ( m1*m1*s1Vec->data[i] + m2*m2*s2Vec->data[i]) / (Mtotal * Mtotal);
    }

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = s1z;
    seobParams.chi2 = s2z;

    REAL8 chiS, chiA, tplspin;
    chiS = 0.5 * (s1z + s2z);
    chiA = 0.5 * (s1z - s2z);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;

    status = CalculateSpinEOBHCoeffs(&seobCoeffs, eta, sigmaKerr->data[2], &adjParams);
    status = CalculateWaveformCoeffs_CQG_25_165003(eta, &hCoeffs);
    REAL8 r, pr, phi, pphi;
    REAL8Vector *dyValues = CreateREAL8Vector(4);
    REAL8Vector *hamVec = CreateREAL8Vector(length);
    REAL8Vector *omegaVec = CreateREAL8Vector(length);
    REAL8Vector *hrealVec = CreateREAL8Vector(length);
    REAL8Vector *himagVec = CreateREAL8Vector(length);
    REAL8Vector *fluxVec = CreateREAL8Vector(length);
    
    REAL8Vector cartPosVec, cartMomVec;
    REAL8 cartPosData[3] = {0,0,0};
    REAL8 cartMomData[3] = {0,0,0};
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;
    COMPLEX16 hLM;
    REAL8 ham, omega, v;

    for(i=0; i<length; i++)
    {
        r = rVec->data[i];
        pr = prVec->data[i];
        phi = phiVec->data[i];
        pphi = pphiVec->data[i];

        cartPosData[0] = r;
        cartMomData[0] = pr;
        cartMomData[1] = pphi/r;

        ham = hamVec->data[i] = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec, 
            s1VecOverMtMt, s2VecOverMtMt, sigmaKerr, sigmaStar, 
            1, &seobCoeffs);
        dyValues->data[0] = r;
        dyValues->data[1] = phi;
        dyValues->data[2] = pr;
        dyValues->data[3] = pphi;
        
        omega = omegaVec->data[i] = XLALSimIMRSpinAlignedEOBCalcOmega(dyValues->data, &seobParams, STEP_SIZE);

        fluxVec->data[i]  = 0;
        v = cbrt(omega);

        status = Waveform_CQG_25_165003(v, phi, eta, 2, 2, &hCoeffs, &hLM);
        hrealVec->data[i] = creal(hLM);
        himagVec->data[i] = cimag(hLM);

        if(status!=CEV_SUCCESS)
        {
            return NULL;
        }
    }

    DestroyREAL8Vector(s1Vec);
    DestroyREAL8Vector(s2Vec);
    DestroyREAL8Vector(s1VecOverMtMt);
    DestroyREAL8Vector(s2VecOverMtMt);
    DestroyREAL8Vector(sigmaStar);
    DestroyREAL8Vector(sigmaKerr);
    DestroyREAL8Vector(dyValues);

    DestroyREAL8Vector(rVec);
    DestroyREAL8Vector(prVec);
    DestroyREAL8Vector(phiVec);
    DestroyREAL8Vector(pphiVec);

    MACRO_REAL8Vector2PyArray(hamVec, py_ham);
    MACRO_REAL8Vector2PyArray(omegaVec, py_omega);
    MACRO_REAL8Vector2PyArray(hrealVec, py_hreal);
    MACRO_REAL8Vector2PyArray(himagVec, py_himag);
    MACRO_REAL8Vector2PyArray(fluxVec, py_flux);

    PyObject *out;
    out = Py_BuildValue("OOOOO", py_ham, py_omega, py_hreal, py_himag, py_flux);
    return out;

}

static PyObject *Func_SEOBNR_Hamiltonian(PyObject *self, PyObject *args)
{
    INT status,i;
    REAL8 m1, m2, s1z, s2z;
    REAL8 KK, dSS, dSO;
    REAL8 py_r, py_pr, py_phi, py_pphi;
    if(!PyArg_ParseTuple(args, "ddddddddddd", &m1, &m2, &s1z, &s2z, &KK, &dSS, &dSO, &py_r, &py_pr, &py_phi, &py_pphi))
    {
        return NULL;
    }
    if(fabs(s1z)>0.999 || fabs(s2z)>0.999)
    {
        return NULL;
    }
    if(m1 < m2)
    {
        SWAP(m1, m2);
    }

    AdjParams adjParams;
    adjParams.dSO = dSO;
    adjParams.dSS = dSS;
    adjParams.KK = KK;
    adjParams.dtPeak = 0;

    REAL8 Mtotal = m1 + m2;
    REAL8 eta = m1 * m2 / Mtotal / Mtotal;

    SpinEOBParams seobParams;
    EOBParams eobParams;
    EOBNonQCCoeffs nqcCoeffs;
    SpinEOBHCoeffs seobCoeffs;

    memset (&seobParams, 0, sizeof (seobParams));
    memset (&seobCoeffs, 0, sizeof (seobCoeffs));
    memset (&eobParams, 0, sizeof (eobParams));
    memset (&nqcCoeffs, 0, sizeof(nqcCoeffs));

    WaveformCoeffsInit(hCoeffs);
    /* Allocate spin parameters */
    REAL8Vector *sigmaStar = CreateREAL8Vector(3);
    REAL8Vector *sigmaKerr = CreateREAL8Vector(3);
    REAL8Vector *s1VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s2VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s1Vec = CreateREAL8Vector(3);
    REAL8Vector *s2Vec = CreateREAL8Vector(3);

    seobParams.deltaT = 1.;
    seobParams.eccentricity = 0;
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

    eobParams.hCoeffs = NULL;
    eobParams.prefixes = NULL;
    eobParams.m1 = m1;
    eobParams.m2 = m2;
    eobParams.eta = eta;
    eobParams.Mtotal = Mtotal;

    s1Vec->data[0] = 0;
    s1Vec->data[1] = 0;
    s1Vec->data[2] = s1z;

    s2Vec->data[0] = 0;
    s2Vec->data[1] = 0;
    s2Vec->data[2] = s2z;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] * m1 * m1 / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] * m2 * m2 / Mtotal / Mtotal;
        sigmaStar->data[i] = ( s1Vec->data[i] + s2Vec->data[i]) * eta;
        sigmaKerr->data[i] = ( m1*m1*s1Vec->data[i] + m2*m2*s2Vec->data[i]) / (Mtotal * Mtotal);
    }

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = s1z;
    seobParams.chi2 = s2z;

    REAL8 chiS, chiA, tplspin;
    chiS = 0.5 * (s1z + s2z);
    chiA = 0.5 * (s1z - s2z);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;

    status = CalculateSpinEOBHCoeffs(&seobCoeffs, eta, sigmaKerr->data[2], &adjParams);

    REAL8Vector *dyValues = CreateREAL8Vector(4);

    REAL8Vector cartPosVec, cartMomVec;
    REAL8 cartPosData[3] = {0,0,0};
    REAL8 cartMomData[3] = {0,0,0};
    REAL8 ham;
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;

    cartPosData[0] = py_r;
    cartMomData[0] = py_pr;
    cartMomData[1] = py_pphi/py_r;

    ham = SpinEOBHamiltonian(eta, &cartPosVec, &cartMomVec, 
        s1VecOverMtMt, s2VecOverMtMt, sigmaKerr, sigmaStar, 
        1, &seobCoeffs);


    DestroyREAL8Vector(s1Vec);
    DestroyREAL8Vector(s2Vec);
    DestroyREAL8Vector(s1VecOverMtMt);
    DestroyREAL8Vector(s2VecOverMtMt);
    DestroyREAL8Vector(sigmaStar);
    DestroyREAL8Vector(sigmaKerr);
    DestroyREAL8Vector(dyValues);

    PyObject *out;
    out = Py_BuildValue("d", ham);
    return out;
}

static PyObject *Func_SpinEOBRadiationReactionForce(PyObject *self, PyObject *args)
{
    INT status,i;
    REAL8 m1, m2, s1z, s2z;
    REAL8 KK, dSS, dSO;
    PyObject *py_r, *py_pr, *py_phi, *py_pphi;
    if(!PyArg_ParseTuple(args, "dddddddOOOO", &m1, &m2, &s1z, &s2z, &KK, &dSS, &dSO, &py_r, &py_pr, &py_phi, &py_pphi))
    {
        print_warning("Failed to Parse function arguments.\n");
        return NULL;
    }
    if(fabs(s1z)>0.999 || fabs(s2z)>0.999)
    {
        print_warning("The spin is out of range.\n");
        return NULL;
    }
    if(m1 < m2)
    {
        SWAP(m1, m2);
    }
    
    MACRO_PyArray2REAL8Vector(py_r, rVec);
    MACRO_PyArray2REAL8Vector(py_pr, prVec);
    MACRO_PyArray2REAL8Vector(py_phi, phiVec);
    MACRO_PyArray2REAL8Vector(py_pphi, pphiVec);

    INT length = GET_MIN(GET_MIN(rVec->length, prVec->length), GET_MIN(phiVec->length, pphiVec->length));
    if(length < 1)
    {
        print_warning("The length of dynamic vector is less than 1\n");
        DestroyREAL8Vector(rVec);
        DestroyREAL8Vector(prVec);
        DestroyREAL8Vector(phiVec);
        DestroyREAL8Vector(pphiVec);
        return NULL;
    }
    AdjParams adjParams;
    adjParams.dSO = dSO;
    adjParams.dSS = dSS;
    adjParams.KK = KK;
    adjParams.dtPeak = 0;

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
    memset (&nqcCoeffs, 0, sizeof(nqcCoeffs));

    /* Allocate spin parameters */
    REAL8Vector *sigmaStar = CreateREAL8Vector(3);
    REAL8Vector *sigmaKerr = CreateREAL8Vector(3);
    REAL8Vector *s1VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s2VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s1Vec = CreateREAL8Vector(3);
    REAL8Vector *s2Vec = CreateREAL8Vector(3);

    seobParams.deltaT = 1.;
    seobParams.eccentricity = 0;
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

    s1Vec->data[0] = 0;
    s1Vec->data[1] = 0;
    s1Vec->data[2] = s1z;

    s2Vec->data[0] = 0;
    s2Vec->data[1] = 0;
    s2Vec->data[2] = s2z;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] * m1 * m1 / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] * m2 * m2 / Mtotal / Mtotal;
        sigmaStar->data[i] = ( s1Vec->data[i] + s2Vec->data[i]) * eta;
        sigmaKerr->data[i] = ( m1*m1*s1Vec->data[i] + m2*m2*s2Vec->data[i]) / (Mtotal * Mtotal);
    }

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = s1z;
    seobParams.chi2 = s2z;

    REAL8 chiS, chiA, tplspin;
    chiS = 0.5 * (s1z + s2z);
    chiA = 0.5 * (s1z - s2z);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;

    status = CalculateSpinEOBHCoeffs(&seobCoeffs, eta, sigmaKerr->data[2], &adjParams);
    status = CalculateSpinFactorizedWaveformCoefficients(&hCoeffs, &seobParams, m1, m2, eta, tplspin, chiS, chiA);
    status = ComputeNewtonMultipolePrefixes(&prefixes, m1, m2);

    REAL8 values[4];
    REAL8 Fr, Fphi, flux, fluxCirc, Fr2PN, Fphi2PN, FrFac, FphiFac, FrResum, FphiResum;
    REAL8Vector *FrVec = CreateREAL8Vector(length);
    REAL8Vector *FphiVec = CreateREAL8Vector(length);
    REAL8Vector *fluxVec = CreateREAL8Vector(length);
    REAL8Vector *FrPNVec = CreateREAL8Vector(length);
    REAL8Vector *FphiPNVec = CreateREAL8Vector(length);
    REAL8Vector *fluxCircVec = CreateREAL8Vector(length);
    REAL8Vector *FrFacVec = CreateREAL8Vector(length);
    REAL8Vector *FphiFacVec = CreateREAL8Vector(length);
    REAL8Vector *FrResVec = CreateREAL8Vector(length);
    REAL8Vector *FphiResVec = CreateREAL8Vector(length);

    for(i=0; i<length; i++)
    {
        values[0] = rVec->data[i];
        values[1] = phiVec->data[i];
        values[2] = prVec->data[i];
        values[3] = pphiVec->data[i];

        status = 
            auxSpinEOBRadiationReactionForceAll(values, &seobParams, 
                            &Fr, &Fphi, &flux, &fluxCirc, 
                            &Fr2PN, &Fphi2PN, 
                            &FrFac, &FphiFac, 
                            &FrResum, &FphiResum);
        FrVec->data[i] = Fr;
        FphiVec->data[i] = Fphi;
        fluxVec->data[i] = flux;

        FrPNVec->data[i] = Fr2PN;
        FphiPNVec->data[i] = Fphi2PN;
        fluxCircVec->data[i] = fluxCirc;

        FrFacVec->data[i] = FrFac;
        FphiFacVec->data[i] = FphiFac;

        FrResVec->data[i] = FrResum;
        FphiResVec->data[i] = FphiResum;
        if(status!=CEV_SUCCESS)
        {
            print_warning("Failed to calculate radiation reaction force.\n");
            return NULL;
        }
    }

    DestroyREAL8Vector(s1Vec);
    DestroyREAL8Vector(s2Vec);
    DestroyREAL8Vector(s1VecOverMtMt);
    DestroyREAL8Vector(s2VecOverMtMt);
    DestroyREAL8Vector(sigmaStar);
    DestroyREAL8Vector(sigmaKerr);

    DestroyREAL8Vector(rVec);
    DestroyREAL8Vector(prVec);
    DestroyREAL8Vector(phiVec);
    DestroyREAL8Vector(pphiVec);

    MACRO_REAL8Vector2PyArray(FrVec, py_Fr);
    MACRO_REAL8Vector2PyArray(FphiVec, py_Fphi);
    MACRO_REAL8Vector2PyArray(fluxVec, py_flux);
    MACRO_REAL8Vector2PyArray(FrPNVec, py_FrPN);
    MACRO_REAL8Vector2PyArray(FphiPNVec, py_FphiPN);
    MACRO_REAL8Vector2PyArray(fluxCircVec, py_fluxCircPN);
    MACRO_REAL8Vector2PyArray(FrFacVec, py_FrFac);
    MACRO_REAL8Vector2PyArray(FphiFacVec, py_FphiFac);
    MACRO_REAL8Vector2PyArray(FrResVec, py_FrRes);
    MACRO_REAL8Vector2PyArray(FphiResVec, py_FphiRes);

    PyObject *out;
    out = Py_BuildValue("OOOOOOOOOO", py_Fr, py_Fphi, py_FrPN, py_FphiPN, py_FrRes, py_FphiRes, py_FrFac, py_FphiFac, py_flux, py_fluxCircPN);
    return out;
}


static PyObject *Func_CalculateCompareFlux(PyObject *self, PyObject *args)
{
    INT status,i;
    REAL8 m1, m2, s1z, s2z;
    REAL8 KK, dSS, dSO;
    PyObject *py_r, *py_pr, *py_phi, *py_pphi;
    if(!PyArg_ParseTuple(args, "dddddddOOOO", &m1, &m2, &s1z, &s2z, &KK, &dSS, &dSO, &py_r, &py_pr, &py_phi, &py_pphi))
    {
        print_warning("Failed to Parse function arguments.\n");
        return NULL;
    }
    if(fabs(s1z)>0.999 || fabs(s2z)>0.999)
    {
        print_warning("The spin is out of range.\n");
        return NULL;
    }
    if(m1 < m2)
    {
        SWAP(m1, m2);
    }
    
    MACRO_PyArray2REAL8Vector(py_r, rVec);
    MACRO_PyArray2REAL8Vector(py_pr, prVec);
    MACRO_PyArray2REAL8Vector(py_phi, phiVec);
    MACRO_PyArray2REAL8Vector(py_pphi, pphiVec);

    INT length = GET_MIN(GET_MIN(rVec->length, prVec->length), GET_MIN(phiVec->length, pphiVec->length));
    if(length < 1)
    {
        print_warning("The length of dynamic vector is less than 1\n");
        DestroyREAL8Vector(rVec);
        DestroyREAL8Vector(prVec);
        DestroyREAL8Vector(phiVec);
        DestroyREAL8Vector(pphiVec);
        return NULL;
    }
    AdjParams adjParams;
    adjParams.dSO = dSO;
    adjParams.dSS = dSS;
    adjParams.KK = KK;
    adjParams.dtPeak = 0;

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
    memset (&nqcCoeffs, 0, sizeof(nqcCoeffs));

    /* Allocate spin parameters */
    REAL8Vector *sigmaStar = CreateREAL8Vector(3);
    REAL8Vector *sigmaKerr = CreateREAL8Vector(3);
    REAL8Vector *s1VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s2VecOverMtMt = CreateREAL8Vector(3);
    REAL8Vector *s1Vec = CreateREAL8Vector(3);
    REAL8Vector *s2Vec = CreateREAL8Vector(3);

    seobParams.deltaT = 1.;
    seobParams.eccentricity = 0;
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

    s1Vec->data[0] = 0;
    s1Vec->data[1] = 0;
    s1Vec->data[2] = s1z;

    s2Vec->data[0] = 0;
    s2Vec->data[1] = 0;
    s2Vec->data[2] = s2z;

    for(i=0; i<3; i++)
    {
        s1VecOverMtMt->data[i] = s1Vec->data[i] * m1 * m1 / Mtotal / Mtotal;
        s2VecOverMtMt->data[i] = s2Vec->data[i] * m2 * m2 / Mtotal / Mtotal;
        sigmaStar->data[i] = ( s1Vec->data[i] + s2Vec->data[i]) * eta;
        sigmaKerr->data[i] = ( m1*m1*s1Vec->data[i] + m2*m2*s2Vec->data[i]) / (Mtotal * Mtotal);
    }

    seobParams.a = sigmaKerr->data[2];
    seobParams.chi1 = s1z;
    seobParams.chi2 = s2z;

    REAL8 chiS, chiA, tplspin;
    chiS = 0.5 * (s1z + s2z);
    chiA = 0.5 * (s1z - s2z);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;

    status = CalculateSpinEOBHCoeffs(&seobCoeffs, eta, sigmaKerr->data[2], &adjParams);
    status = CalculateSpinFactorizedWaveformCoefficients(&hCoeffs, &seobParams, m1, m2, eta, tplspin, chiS, chiA);
    status = ComputeNewtonMultipolePrefixes(&prefixes, m1, m2);

    REAL8 values[4];
    REAL8 Eflux, FphiEOB, EfluxPN, LfluxPN, EfluxNagar, LfluxNagar, EfluxFac, LfluxFac;

    REAL8Vector *EfluxVec = CreateREAL8Vector(length);
    REAL8Vector *FphiEOBVec = CreateREAL8Vector(length);
    REAL8Vector *EfluxPNVec = CreateREAL8Vector(length);
    REAL8Vector *LfluxPNVec = CreateREAL8Vector(length);
    REAL8Vector *EfluxNagarVec = CreateREAL8Vector(length);
    REAL8Vector *LfluxNagarVec = CreateREAL8Vector(length);
    REAL8Vector *EfluxFacVec = CreateREAL8Vector(length);
    REAL8Vector *LfluxFacVec = CreateREAL8Vector(length);

    for(i=0; i<length; i++)
    {
        values[0] = rVec->data[i];
        values[1] = phiVec->data[i];
        values[2] = prVec->data[i];
        values[3] = pphiVec->data[i];

        status = 
            auxCalculateCompareFlux(values, &seobParams, 
                &Eflux, &FphiEOB, 
                &EfluxPN, &LfluxPN,
                &EfluxNagar, &LfluxNagar,
                &EfluxFac, &LfluxFac);
        if(status!=CEV_SUCCESS)
        {
            print_warning("Failed to calculate radiation reaction force.\n");
            return NULL;
        }
        //print_debug("flux = %.4f\n", Eflux);
        EfluxVec->data[i] = Eflux;
        FphiEOBVec->data[i] = FphiEOB;

        EfluxPNVec->data[i] = EfluxPN;
        LfluxPNVec->data[i] = LfluxPN;

        EfluxNagarVec->data[i] = EfluxNagar;
        LfluxNagarVec->data[i] = LfluxNagar;

        EfluxFacVec->data[i] = EfluxFac;
        LfluxFacVec->data[i] = LfluxFac;
    }

    DestroyREAL8Vector(s1Vec);
    DestroyREAL8Vector(s2Vec);
    DestroyREAL8Vector(s1VecOverMtMt);
    DestroyREAL8Vector(s2VecOverMtMt);
    DestroyREAL8Vector(sigmaStar);
    DestroyREAL8Vector(sigmaKerr);

    DestroyREAL8Vector(rVec);
    DestroyREAL8Vector(prVec);
    DestroyREAL8Vector(phiVec);
    DestroyREAL8Vector(pphiVec);

    MACRO_REAL8Vector2PyArray(EfluxVec, py_Eflux);
    MACRO_REAL8Vector2PyArray(FphiEOBVec, py_FphiEOB);
    MACRO_REAL8Vector2PyArray(EfluxPNVec, py_EfluxPN);
    MACRO_REAL8Vector2PyArray(LfluxPNVec, py_LfluxPN);
    MACRO_REAL8Vector2PyArray(EfluxNagarVec, py_EfluxNagar);
    MACRO_REAL8Vector2PyArray(LfluxNagarVec, py_LfluxNagar);
    MACRO_REAL8Vector2PyArray(EfluxFacVec, py_EfluxFac);
    MACRO_REAL8Vector2PyArray(LfluxFacVec, py_LfluxFac);


    PyObject *out;
    out = Py_BuildValue("OOOOOOOO", py_Eflux, py_FphiEOB, 
                                py_EfluxPN, py_LfluxPN,
                                py_EfluxNagar, py_LfluxNagar,
                                py_EfluxFac, py_LfluxFac);
    return out;
}


static PyMethodDef PyUtils_Methods[] = {
    {"SpinEOBDynamicsGeneration", Func_SpinEOBDynamicsGeneration, METH_VARARGS},
    {"PNwaveformCQG_25_165003", Func_PNwaveformCQG_25_165003, METH_VARARGS},
    {"SpinEOBHamiltonian", Func_SEOBNR_Hamiltonian, METH_VARARGS},
    {"SpinEOBRadiationReactionForce", Func_SpinEOBRadiationReactionForce, METH_VARARGS},
    {"CalculateCompareFlux", Func_CalculateCompareFlux, METH_VARARGS},
    {NULL, NULL},
};

static struct PyModuleDef PyUtilsModule = {
    PyModuleDef_HEAD_INIT,
    "PyUtils",
    NULL,
    -1,
    PyUtils_Methods
};

PyMODINIT_FUNC PyInit_PyUtils(void)
{
    PyObject *m = PyModule_Create(&PyUtilsModule);
    import_array();
    return m;
}

