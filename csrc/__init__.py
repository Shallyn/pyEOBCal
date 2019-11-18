"""
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
"""

from WTestLib.Utils import cmd_stdout_cev, CEV 
from pathlib import Path
import time
from numpy import random as rd

LOC = Path(__file__).parent
EXE = LOC / 'MAIN'

__all__ = ['playEOB', 'playEOB_withAdj']
def ConstructCMD(exe = EXE, **kwargs):
    CMD = [str(exe)]
    if kwargs:
        for (attr, val) in kwargs.items():
            opt = attr.replace('_','-')
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
