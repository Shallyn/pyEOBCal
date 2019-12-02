#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 00:37:51 2019

@author: drizl
"""

import numpy as np
import sys, os, json
from pathlib import Path
from WTestLib.SXS import SXSparameters, DEFAULT_SRCLOC, DEFAULT_TABLE, loadSXStxtdata, plot_fit, parse_ecc, CompResults
from WTestLib.h22datatype import get_Mtotal, h22base, dim_t, h22_alignment, get_fmin
from WTestLib.generator import self_adaptivor
from .HyperCalibrator import SEOBHCoeffsCalibrator
from . import playEOB_withAdj, playEOB

class SXSAdjustor(SXSparameters):
    def __init__(self, SXSnum, f_min_dimless = 0.002, Mtotal = 30, f_min = -1, D = 100, srate = 16384,
                 srcloc = DEFAULT_SRCLOC, 
                 table = DEFAULT_TABLE):
        if f_min > 0:
            Mtotal = get_Mtotal(f_min_dimless, f_min)
        else:
            f_min = get_fmin(f_min_dimless, Mtotal)
        print(f_min)
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

    def plot_fit(self, pms, ecc = 0, fname = 'save.png', fit = True, **kwargs):
        wf = self.get_waveform(pms, ecc)
        if fit:
            plot_fit(self._SXSh22, wf, fname, 
                     name1 = f'SXS:BBH:{self._SXSnum}',
                     name2 = 'SEOBNRE', **kwargs)
        else:
            wf.plot(fname)

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
        Eps, dephase = calculate_FF_dephase(self._SXSh22, wf)
        return -(pow(Eps/0.01,2) + pow(dephase/5/self._tprod,2 ))/2
    
    def get_FF(self, pms, ecc = 0):
        wf = self.get_waveform(pms, ecc)
        if wf is None:
            return 0,0,-1
        FF, tc, phic = calculate_FF(self._SXSh22, wf)
        return tc, phic, FF

        
    def sa_find_ecc(self, pms, estep = 0.02, eccrange = None,
                    maxitr = None, prec_x = 1e-6, prec_y = 1e-6, verbose = True):
        if eccrange is None:
            eccrange = parse_ecc(self.ecc, 0.8)
        def func(ecc):
            return self.get_FF(pms, ecc = ecc)
        SA = self_adaptivor(func, eccrange, estep, outindex = 2)
        wrapper =  SA.run(maxitr = maxitr, 
                      verbose = verbose,
                      prec_x = prec_x, 
                      prec_y = prec_y)
        return CompResultsNew(self, wrapper, pms)

class CompResultsNew(object):
    def __init__(self, adjustor, results, pms):
        self._core = adjustor
        self._results = results
        self._parse_results()
        self._pms = pms
        
    def _parse_results(self):
        ecc, olpout = self._results
        self._ecc = ecc
        self._tc = olpout[:,0]
        self._phic = olpout[:,1]
        self._FF = olpout[:,2]
        fitarg = self._FF.argmax()
        self._max_FF = self._FF[fitarg]
        self._fit_tc = self._tc[fitarg]
        self._fit_phic = self._phic[fitarg]
        self._fit_ecc = self._ecc[fitarg]
        
    def plot_fit(self, fname, **kwargs):
        return self._core.plot_fit(pms = self._pms, ecc = self._fit_ecc, fname = fname,
                                   fit = True, **kwargs)
        

def calculate_FF(wf1, wf2):
    wf_1, wf_2, _ = h22_alignment(wf1, wf2)
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
    Oxt_abs = np.abs(Oxt)
    idxmax = Oxt_abs.argmax()
    lth = len(Oxt_abs)
    if idxmax > lth / 2:
        tc = (idxmax - lth) / fs
    else:
        tc = idxmax / fs
    return Oxt_abs[idxmax], tc, np.angle(Oxt[idxmax])

def calculate_FF_dephase(wf1, wf2):
    wf_1, wf_2, _ = h22_alignment(wf1, wf2)
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
    Oxt = np.abs(np.fft.ifft(Ox) * fs / np.sqrt(O11 * O22))
    idxmax = Oxt.argmax()
    lth = len(Oxt)
    if idxmax > lth / 2:
        tc = (idxmax - lth) / fs
    else:
        tc = idxmax / fs

    return (1 - Oxt[idxmax], tc)