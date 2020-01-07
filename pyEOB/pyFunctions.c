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

static PyMethodDef PyUtils_Methods[] = {
    {"SpinEOBDynamicsGeneration", Func_SpinEOBDynamicsGeneration, METH_VARARGS},
    {"PNwaveformCQG_25_165003", Func_PNwaveformCQG_25_165003, METH_VARARGS},
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

