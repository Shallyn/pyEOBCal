#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 00:37:51 2019

@author: drizl
"""

import numpy as np
import sys, os, json
from pathlib import Path
from WTestLib.SXS import SXSparameters, DEFAULT_SRCLOC, DEFAULT_TABLE, loadSXStxtdata
from WTestLib.h22datatype import get_Mtotal, h22base, dim_t, h22_alignment
from .HyperCalibrator import SEOBHCoeffsCalibrator
from . import playEOB_withAdj, playEOB

class SXSAdjustor(SXSparameters):
    def __init__(self, SXSnum, f_min_dimless = 0.002, f_min = 10, D = 100, srate = 16384,
                 srcloc = DEFAULT_SRCLOC, 
                 table = DEFAULT_TABLE):
        Mtotal = get_Mtotal(f_min_dimless, f_min)
        self._tprod = dim_t(Mtotal)
        super(SXSAdjustor, self).__init__(SXSnum, table = table, f_ini = f_min_dimless, Mtotal = Mtotal, D = D, verbose = False, ishertz = False)
        self._f_min_dimless = f_min_dimless
        self._f_min = f_min
        _cal = SEOBHCoeffsCalibrator('SEOBNRv4')
        self._initSHC = _cal(self.m1, self.m2, self.s1z, self.s2z)
        t, hr, hi = loadSXStxtdata(SXSnum, srcloc)
        self._SXSh22 = h22base(t/self._tprod, hr, hi, srate)

    @property
    def adjParamsV4(self):
        return self._initSHC

    @property
    def srate(self):
        return self._SXSh22.srate

    def get_waveform(self, pms, ecc = 0):
        KK, dSS, dSO, dtPeak = pms[0], pms[1], pms[2], pms[3]        
        ret = playEOB_withAdj(m1 = self.m1, m2 = self.m2,
                            spin1z = self.s1z, spin2z=self.s2z, eccentricity = ecc,
                            fMin = self._f_min, fs = self.srate, 
                            KK = KK, dSS = dSS, dSO = dSO, dtPeak = dtPeak)
        if ret is not None:
            t,hr,hi = ret
            ret = h22base(t, hr, hi, self.srate)
            return ret
        else:
            return None

    def get_lnprob(self, pms, ecc = 0):
        wf = self.get_waveform(pms, ecc)
        if wf is None:
            return -np.inf
        FF, dephase = calculate_FF_dephase(self._SXSh22, wf)
        return -(pow(FF/0.01,2) + pow(dephase/5/self._tprod,2 ))/2

def calculate_FF_dephase(wf1, wf2):
    wf_1, wf_2, _ = h22_alignment(wf1, wf2)
    dephase = abs(wf_1.tpeak - wf_2.tpeak)
    fs = wf_1.srate
    NFFT = len(wf_1)
    df = fs/NFFT
    #freqs = np.abs(np.fft.fftfreq(NFFT, 1./fs))
    #power_vec = psdfunc(freqs)
    htilde_1 = wf_1.h22f
    htilde_2 = wf_2.h22f
    O11 = np.sum(htilde_1 * htilde_1.conjugate()).real * df
    O22 = np.sum(htilde_2 * htilde_2.conjugate()).real * df
    Ox = htilde_1 * htilde_2.conjugate()
    Oxt = np.fft.ifft(Ox) * fs / np.sqrt(O11 * O22)
    FF = np.abs(Oxt).max()
    return (FF, dephase)