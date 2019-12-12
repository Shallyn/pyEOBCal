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
from numpy import random as rd
from WTestLib.SXS import DEFAULT_SRCLOC

"""
from .PyEOBCal import SimIMRSpinAlignedEOBFullWaveformBEvolveWithAdj
__all__ = ['playEOB', 'playEOB_withAdj']

def playEOB(m1 = 10, 
            m2 = 10,
            spin1x = 0,
            spin1y = 0,
            spin1z = 0,
            spin2x = 0,
            spin2y = 0,
            spin2z = 0,
            eccentricity = 0,
            fMin = 40,
            fs = 16384):

    KK = dSS = dSO = dtPeak = 0
    try:
        ret = SimIMRSpinAlignedEOBFullWaveformBEvolveWithAdj(m1, m2, 
                spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                eccentricity, fMin, fs, 
                KK, dSS, dSO, dtPeak, 1)
    except:
        ret = None
    return ret

def playEOB_withAdj(m1 = 10, 
                    m2 = 10,
                    spin1x = 0,
                    spin1y = 0,
                    spin1z = 0,
                    spin2x = 0,
                    spin2y = 0,
                    spin2z = 0,
                    eccentricity = 0,
                    fMin = 40,
                    fs = 16384,
                    KK = 0,
                    dSS = 0,
                    dSO = 0,
                    dtPeak = 0):
    try:
        ret = SimIMRSpinAlignedEOBFullWaveformBEvolveWithAdj(m1, m2, 
                spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                eccentricity, fMin, fs, 
                KK, dSS, dSO, dtPeak, 0)
    except:
        ret = None
    return ret


"""
from WTestLib.Utils import cmd_stdout_cev, CEV 

LOC = Path(__file__).parent
EXE = LOC / 'MAIN'
ITERNQC = LOC / 'ITERNQC'

def ConstructCMD(exe = EXE, **kwargs):
    CMD = [str(exe)]
    if kwargs:
        for (attr, val) in kwargs.items():
            opt = attr.replace('_','-')
            if val is None:
                CMD.append(f'--{opt}')
            else:
                CMD.append(f'--{opt}={val}')
    return ' '.join(CMD)

def get_random_jobtag():
    ttag = int(time.time() / 1000)
    rdtag = int(rd.randn()*ttag)
    return f'_job_{ttag}_{rdtag}.job'

def playEOB(m1 = 10, 
            m2 = 10,
            spin1x = 0,
            spin1y = 0,
            spin1z = 0,
            spin2x = 0,
            spin2y = 0,
            spin2z = 0,
            eccentricity = 0,
            fMin = 40,
            fs = 16384):

    CMD = ConstructCMD(exe = EXE, m1 = m1, m2 = m2, 
                    spin1z = spin1z,
                    spin2z = spin2z,
                    eccentricity = eccentricity, 
                    f_min = fMin, sample_rate = fs)
    status, data = cmd_stdout_cev(CMD, get_random_jobtag())
    if status is CEV.SUCCESS and len(data) != 0:
        return (data[:,0], data[:,1], data[:,2])
    else:
        return None

def playEOB_withAdj(m1 = 10, 
                    m2 = 10,
                    spin1x = 0,
                    spin1y = 0,
                    spin1z = 0,
                    spin2x = 0,
                    spin2y = 0,
                    spin2z = 0,
                    eccentricity = 0,
                    fMin = 40,
                    fs = 16384,
                    KK = 0,
                    dSS = 0,
                    dSO = 0,
                    dtPeak = 0):
    
    CMD = ConstructCMD(exe = EXE, m1 = m1, m2 = m2, 
                    spin1z = spin1z,
                    spin2z = spin2z,
                    eccentricity = eccentricity, 
                    f_min = fMin, sample_rate = fs,
                    KK = KK, dSS = dSS, dSO = dSO, dtPeak = dtPeak)
    status, data = cmd_stdout_cev(CMD, get_random_jobtag())
    if status is CEV.SUCCESS and len(data) != 0:
        return (data[:,0], data[:,1], data[:,2])
    else:
        return None


def playEOB_iterNQC(m1 = 10, 
                    m2 = 10,
                    spin1x = 0,
                    spin1y = 0,
                    spin1z = 0,
                    spin2x = 0,
                    spin2y = 0,
                    spin2z = 0,
                    eccentricity = 0,
                    fMin = 40,
                    fs = 16384,
                    KK = 0,
                    dSS = 0,
                    dSO = 0,
                    dtPeak = 0,
                    eps = 1e-2,
                    maxiterstep = 100,
                    SXSnum = '0001',
                    srcloc = DEFAULT_SRCLOC,
                    ret_waveform = False):
    
    fileNR = DEFAULT_SRCLOC / f'BBH_{SXSnum}.txt'
    if ret_waveform:
        CMD = ConstructCMD(exe = ITERNQC, m1 = m1, m2 = m2, 
                            spin1z = spin1z, spin2z = spin2z,
                            eccentricity = eccentricity,
                            f_min = fMin, sample_rate = fs,
                            KK = KK, dSS = dSS, dSO = dSO, dtPeak = dtPeak, eps = eps,
                            file_SXS = fileNR, max_iterstep = maxiterstep,
                            get_waveform = None)
    else:
        CMD = ConstructCMD(exe = ITERNQC, m1 = m1, m2 = m2, 
                            spin1z = spin1z, spin2z = spin2z,
                            eccentricity = eccentricity,
                            f_min = fMin, sample_rate = fs,
                            KK = KK, dSS = dSS, dSO = dSO, dtPeak = dtPeak,eps = eps,
                            file_SXS = fileNR, max_iterstep = maxiterstep)
    print(CMD)
    status, data = cmd_stdout_cev(CMD, get_random_jobtag())
    if status is CEV.SUCCESS and len(data) != 0:
        return (data[:,0], data[:,1], data[:,2])
    else:
        return None
