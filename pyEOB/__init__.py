"""
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
"""

from pathlib import Path
import time
import numpy as np
from . import PyUtils as pu

__all__ = ['SpinEOBFactorizedWaveform', 'PNwaveformCQG_25_165003', 'SpinEOBRadiationReactionForce']


def SpinEOBFactorizedWaveform(m1 = 10, 
                m2 = 10, 
                spin1z = 0, 
                spin2z = 0,
                r = 50,
                pr = -1e-6,
                phi = 0,
                pphi = 6.78604,
                KK = 0,
                dSO = 0,
                dSS = 0):
    if not hasattr(r, '__len__'):
        r = [r]
    if not hasattr(pr, '__len__'):
        pr = [pr]
    if not hasattr(phi, '__len__'):
        phi = [phi]
    if not hasattr(pphi, '__len__'):
        pphi = [pphi]
    r = np.asarray(r)
    pr = np.asarray(pr)
    phi = np.asarray(phi)
    pphi = np.asarray(pphi)
    try:
        return pu.SpinEOBDynamicsGeneration(m1, m2, spin1z, spin2z, KK, dSS, dSO, r, pr, phi, pphi)
    except:
        return None

def PNwaveformCQG_25_165003(m1 = 10, 
                m2 = 10, 
                spin1z = 0, 
                spin2z = 0,
                r = 50,
                pr = -1e-6,
                phi = 0,
                pphi = 6.78604,
                KK = 0,
                dSO = 0,
                dSS = 0):
    if not hasattr(r, '__len__'):
        r = [r]
    if not hasattr(pr, '__len__'):
        pr = [pr]
    if not hasattr(phi, '__len__'):
        phi = [phi]
    if not hasattr(pphi, '__len__'):
        pphi = [pphi]
    r = np.asarray(r)
    pr = np.asarray(pr)
    phi = np.asarray(phi)
    pphi = np.asarray(pphi)
    try:
        return pu.PNwaveformCQG_25_165003(m1, m2, spin1z, spin2z, KK, dSS, dSO, r, pr, phi, pphi)
    except:
        return None

def SpinEOBRadiationReactionForce(m1 = 10, 
                m2 = 10, 
                spin1z = 0, 
                spin2z = 0,
                r = 50,
                pr = -1e-6,
                phi = 0,
                pphi = 6.78604,
                KK = 0,
                dSO = 0,
                dSS = 0):
    if not hasattr(r, '__len__'):
        r = [r]
    if not hasattr(pr, '__len__'):
        pr = [pr]
    if not hasattr(phi, '__len__'):
        phi = [phi]
    if not hasattr(pphi, '__len__'):
        pphi = [pphi]
    r = np.asarray(r)
    pr = np.asarray(pr)
    phi = np.asarray(phi)
    pphi = np.asarray(pphi)
    try:
        return pu.SpinEOBRadiationReactionForce(m1, m2, spin1z, spin2z, KK, dSS, dSO, r, pr, phi, pphi)
    except:
        return None

def CalculateCompareFlux(m1 = 10, 
                m2 = 10, 
                spin1z = 0, 
                spin2z = 0,
                r = 50,
                pr = -1e-6,
                phi = 0,
                pphi = 6.78604,
                KK = 0,
                dSO = 0,
                dSS = 0):
    if not hasattr(r, '__len__'):
        r = [r]
    if not hasattr(pr, '__len__'):
        pr = [pr]
    if not hasattr(phi, '__len__'):
        phi = [phi]
    if not hasattr(pphi, '__len__'):
        pphi = [pphi]
    r = np.asarray(r)
    pr = np.asarray(pr)
    phi = np.asarray(phi)
    pphi = np.asarray(pphi)
    try:
        return pu.CalculateCompareFlux(m1, m2, spin1z, spin2z, KK, dSS, dSO, r, pr, phi, pphi)
    except:
        return None