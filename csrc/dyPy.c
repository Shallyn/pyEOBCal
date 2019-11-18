/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyPy.h"


static PyObject *Func_SimIMRSpinAlignedEOBFullWaveformWithAdj(PyObject *self, PyObject *args)
{
    INT status;
    REAL8 m1, m2;
    REAL8 s1x, s1y, s1z, s2x, s2y, s2z, ecc;
    REAL8 fMin, srate;
    AdjParams adjParams;
    CtrlParams ctrlParams;
    REAL8 KK, dSS, dSO, dtPeak;
    INT adjDefault;
    if(!PyArg_ParseTuple(args, "dddddddddddddddi", &m1, &m2, 
            &s1x, &s1y, &s1z, &s2x, &s2y, &s2z,
            &ecc, &fMin, &srate, 
            &KK, &dSS, &dSO, &dtPeak, &adjDefault))
            return NULL;
    
    COMPLEX16TimeSeries *h22 = NULL;
    if(adjDefault)
    {
        applyDefaultAdjustableParameters(&adjParams, m1, m2, s1z, s2z);
    }
    else
    {
        adjParams.dSO = dSO;
        adjParams.KK = KK;
        adjParams.dSS = dSS;
        adjParams.dtPeak = dtPeak;
    }
    status = EvolutionCore(m1, m2, fMin, ecc, 1./srate, 
        s1z, s2z, &h22, &adjParams, &ctrlParams);
    if( (status != CEV_SUCCESS) || !h22)
    {
        print_warning("Failed! return code = %d", status);
        DestroyCOMPLEX16TimeSeries(h22);
        return NULL;
    }
print_debug("Return...\n");
    UINT length, i;
    REAL8 t0, dt;
    t0 = h22->epoch;
    dt = h22->deltaT;
    length = h22->data->length;
print_debug("Return...\n");
    REAL8 *ret_time = NULL;
    REAL8 *ret_hreal = NULL;
    REAL8 *ret_himag = NULL;
    ret_time = (REAL8*) malloc(length*sizeof(REAL8));
    ret_hreal = (REAL8*) malloc(length*sizeof(REAL8));
    ret_himag = (REAL8*) malloc(length*sizeof(REAL8));
print_debug("Return...\n");
    for (i=0;i<length;i++)
    {
        ret_time[i] = t0 + i*dt;
        ret_hreal[i] = creal(h22->data->data[i]);
        ret_himag[i] = cimag(h22->data->data[i]);
    }
    DestroyCOMPLEX16TimeSeries(h22);
print_debug("Return...\n");
    PyObject *py_time, *py_hreal, *py_himag;
    PyObject *out;
    npy_intp h_length[] = {0};
    h_length[0] = length;
    py_time = PyArray_SimpleNewFromData(1, h_length, NPY_DOUBLE, (void*)ret_time);
    py_hreal = PyArray_SimpleNewFromData(1, h_length, NPY_DOUBLE, (void*)ret_hreal);
    py_himag = PyArray_SimpleNewFromData(1, h_length, NPY_DOUBLE, (void*)ret_himag);
    out = Py_BuildValue("OOO", py_time, py_hreal, py_himag);
print_debug("Return...\n");
    return out;
}

static PyMethodDef PyEOBCal_Methods[] = {
    {"SimIMRSpinAlignedEOBFullWaveformBEvolveWithAdj", Func_SimIMRSpinAlignedEOBFullWaveformWithAdj, METH_VARARGS},
    {NULL, NULL},
};

static struct PyModuleDef PyEOBCalModule = {
    PyModuleDef_HEAD_INIT,
    "PyEOBCal",
    NULL,
    -1,
    PyEOBCal_Methods
};

PyMODINIT_FUNC PyInit_PyEOBCal(void)
{
    PyObject *m = PyModule_Create(&PyEOBCalModule);
    import_array();
    return m;
}

