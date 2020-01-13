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
from . import playEOB_withAdj, playEOB, playEOB_iterNQC, playEOB_withecc

class SXSAdjustor(SXSparameters):
    def __init__(self, SXSnum, f_min_dimless = 0.002, Mtotal = 30, f_min = -1, D = 100, srate = 16384,
                 srcloc = DEFAULT_SRCLOC, 
                 table = DEFAULT_TABLE):
        if f_min > 0:
            Mtotal = get_Mtotal(f_min_dimless, f_min)
            ishertz = False
        else:
            f_min_dimless = get_fmin(f_min_dimless, Mtotal)
            ishertz = True
        self._tprod = dim_t(Mtotal)
        super(SXSAdjustor, self).__init__(SXSnum, table = table, f_ini = f_min_dimless, Mtotal = Mtotal, D = D, verbose = False, ishertz = ishertz)
        self._f_min_dimless = self.f_ini_dimless
        self._f_min = self.f_ini
        _cal = SEOBHCoeffsCalibrator('SEOBNRv4')
        self._initSHC = _cal(self.m1, self.m2, self.s1z, self.s2z)
        self._initSHC_V1 = _cal.myCalibratorV1(self.m1, self.m2, self.s1z, self.s2z)
        self._initSHC_V2 = _cal.myCalibratorV2(self.m1, self.m2, self.s1z, self.s2z)
        t, hr, hi = loadSXStxtdata(SXSnum, srcloc)
        self._SXSh22 = h22base(t/self._tprod, hr, hi, srate)
        self._fileSXS = Path(srcloc) / f'BBH_{self._SXSnum}.txt'

    @property
    def adjParamsV4(self):
        return self._initSHC

    @property
    def adjParamsNewV1(self):
        return self._initSHC_V1

    @property
    def adjParamsNewV2(self):
        return self._initSHC_V2

    @property
    def srate(self):
        return self._SXSh22.srate

    def plot_fit(self, pms, ecc = 0, fname = 'save.png', fit = True, **kwargs):
        wf = self.get_waveform_withecc(pms, ecc)
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
    
    def get_waveform_iterNQC(self, pms, ecc = 0,
                            eps = 1e-1, maxiterstep = 100):
        KK, dSS, dSO, dtPeak = pms[0], pms[1], pms[2], pms[3]        
        ret = playEOB_iterNQC(m1 = self.m1, m2 = self.m2,
                            spin1z = self.s1z, spin2z=self.s2z, eccentricity = ecc,
                            fMin = self._f_min, fs = self.srate, 
                            KK = KK, dSS = dSS, dSO = dSO, dtPeak = dtPeak, 
                            eps = eps, maxiterstep = maxiterstep,
                            SXSnum = self.SXSnum, srcloc = self._fileSXS.parent,
                            ret_waveform = True)
        if ret is not None:
            t,hr,hi = ret
            ret = h22base(t, hr, hi, self.srate)
            return ret
        else:
            return None

    def get_waveform_withecc(self, pms, ecc = 0):
        KK, dSS, dSO, dtPeak = pms[0], pms[1], pms[2], pms[3]        
        ret = playEOB_withecc(m1 = self.m1, m2 = self.m2,
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
        print(f'FF = {1-Eps}')
        return -(pow(Eps/0.01,2) + pow(dephase/5/self._tprod,2 ))/2

    def get_lnprob_nospin(self, pms, ecc = 0):
        KK, dtPeak = pms[0], pms[1]
        wf = self.get_waveform((KK, self.adjParamsV4.dSS, self.adjParamsV4.dSO, dtPeak), ecc)
        if wf is None:
            return -np.inf
        Eps, dephase = calculate_FF_dephase(self._SXSh22, wf)
        print(f'FF = {1-Eps}')
        return -(pow(Eps/0.01,2) + pow(dephase/5/self._tprod,2 ))/2

    def get_lnprob_nospin_iterNQC(self, pms, ecc = 0,
                                    eps = 1e-1, maxiterstep = 100):
        KK, dtPeak = pms[0], pms[1]
        wf = self.get_waveform_iterNQC((KK, self.adjParamsV4.dSS, self.adjParamsV4.dSO, dtPeak), 
                                        ecc, eps = eps, maxiterstep = maxiterstep)
        if wf is None:
            return -np.inf
        Eps, dephase = calculate_FF_dephase(self._SXSh22, wf)
        return -(pow(Eps/0.01,2) + pow(dephase/5/self._tprod,2 ))/2
    
    def get_lnprob_withecc(self, pms, ecc = 0):
        wf = self.get_waveform_withecc(pms, ecc)
        if wf is None:
            return -np.inf
        Eps, dephase = calculate_FF_dephase(self._SXSh22, wf)
        print(f'FF = {1-Eps}')
        return -(pow(Eps/0.01,2) + pow(dephase/5/self._tprod,2 ))/2


    def get_lnprob_nospin_withecc(self, pms, ecc = 0):
        KK, dtPeak = pms[0], pms[1]
        return self.get_lnprob_withecc([KK, self.adjParamsV4.dSS, self.adjParamsV4.dSO, dtPeak], ecc)

    
    def get_FF(self, pms, ecc = 0):
        wf = self.get_waveform(pms, ecc)
        if wf is None:
            return -1,0,-1
        FF, tc, phic = calculate_FF(self._SXSh22, wf)
        return tc, phic, FF

    def iterNQC(self, pms,
                ecc = 0,
                eps = 1e-3,
                maxiterstep = 100,
                srcloc = DEFAULT_SRCLOC):
        KK, dSS, dSO, dtPeak = pms[0], pms[1], pms[2], pms[3] 
        ret = playEOB_iterNQC(m1 = self.m1,
                              m2 = self.m2,
                              spin1z = self.s1z,
                              spin2z = self.s2z,
                              eccentricity = ecc,
                              fMin = self._f_min, fs = self.srate,
                              KK = KK, dSS = dSS, dSO = dSO, dtPeak = dtPeak,
                              eps = eps, maxiterstep = maxiterstep,
                              srcloc = srcloc, SXSnum = self.SXSnum)
        return ret
        
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

from WTestLib.SXS import save_namecol, add_csv
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

    def save(self, prefix):
        prefix = Path(prefix)
        if not prefix.exists():
            prefix.mkdir(parents = True)
        fsave_all = prefix / f'eccentricity_FF_{self._core.SXSnum}.txt'
        fsave_csv = prefix / f'All.csv'
        if not fsave_csv.exists():
            save_namecol(fsave_csv, data = [ ['#SXSid', '#mratio', '#spin1x', '#spin1y', '#spin1z', '#spin2x', '#spin2y', '#spin2z', '#ecc', '#ecc_fit', '#FF'] ])
        np.savetxt(fsave_all, np.stack([self._ecc, self._FF]).T)
        add_csv(fsave_csv, [ [self._core.SXSnum, \
                              str(self._core.q), \
                              str(self._core.s1x), \
                              str(self._core.s1y), \
                              str(self._core.s1z), \
                              str(self._core.s2x), \
                              str(self._core.s2y), \
                              str(self._core.s2z), \
                              self._core.ecc, \
                              str(self._fit_ecc), \
                              str(self._max_FF)] ])

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