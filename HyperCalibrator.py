#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 00:36:42 2019

@author: drizl
"""

from WTestLib.Utils import switch


def SEOBHyperCoefficients_v1(eta, a):
    c0 = 1.4467
    c1 = -1.7152360250654402
    c2 = -3.246255899738242
    KK = c0 + c1 * eta + c2 * eta * eta
    dSO = -69.5
    dSS = 2.75
    dtPeak = 0
    dic = {}
    dic['KK'] = KK
    dic['dSO'] = dSO
    dic['dSS'] = dSS
    dic['dtPeak'] = dtPeak
    return SEOBHCoeffs(dic)

def SEOBHyperCoefficients_v2(eta, a):
    c20 = 1.712
    c21 = -1.803949138004582
    c22 = -39.77229225266885
    c23 = 103.16588921239249
    KK = c20 + c21 * eta + c22 * eta * eta + c23 * eta * eta * eta
    dSO = -74.71 - 156. * eta + 627.5 * eta * eta
    dSS = 8.127 - 154.2 * eta + 830.8 * eta * eta
    dtPeak = 0
    dic = {}
    dic['KK'] = KK
    dic['dSO'] = dSO
    dic['dSS'] = dSS
    dic['dtPeak'] = dtPeak
    return SEOBHCoeffs(dic)

def SEOBHyperCoefficients_v4(eta, a):
    coeff00K = 1.7336
    coeff01K = -1.62045
    coeff02K = -1.38086
    coeff03K = 1.43659
    coeff10K = 10.2573
    coeff11K = 2.26831
    coeff12K = 0
    coeff13K = -0.426958
    coeff20K = -126.687
    coeff21K = 17.3736
    coeff22K = 6.16466
    coeff23K = 0
    coeff30K = 267.788
    coeff31K = -27.5201
    coeff32K = 31.1746
    coeff33K = -59.1658

    coeff00dSO = -44.5324
    coeff01dSO = 0
    coeff02dSO = 0
    coeff03dSO = 66.1987
    coeff10dSO = 0
    coeff11dSO = 0
    coeff12dSO = -343.313
    coeff13dSO = -568.651
    coeff20dSO = 0
    coeff21dSO = 2495.29
    coeff22dSO = 0
    coeff23dSO = 147.481
    coeff30dSO = 0
    coeff31dSO = 0
    coeff32dSO = 0
    coeff33dSO = 0

    coeff00dSS = 6.06807
    coeff01dSS = 0
    coeff02dSS = 0
    coeff03dSS = 0
    coeff10dSS = -36.0272
    coeff11dSS = 37.1964
    coeff12dSS = 0
    coeff13dSS = -41.0003
    coeff20dSS = 0
    coeff21dSS = 0
    coeff22dSS = -326.325
    coeff23dSS = 528.511
    coeff30dSS = 706.958
    coeff31dSS = 0
    coeff32dSS = 1161.78
    coeff33dSS = 0.

    coeff00DT = 2.50499
    coeff01DT = 13.0064
    coeff02DT = 11.5435
    coeff03DT = 0
    coeff10DT = 45.8838
    coeff11DT = -40.3183
    coeff12DT = 0
    coeff13DT = -19.0538
    coeff20DT = 13.0879
    coeff21DT = 0
    coeff22DT = 0
    coeff23DT = 0.192775
    coeff30DT = -716.044
    coeff31DT = 0
    coeff32DT = 0
    coeff33DT = 0


    chi = a / (1. - 2. * eta)
    eta2 = eta * eta
    eta3 = eta2 * eta
    chi2 = chi * chi
    chi3 = chi2 * chi
    KK = \
        coeff00K + coeff01K * chi + coeff02K * chi2 + coeff03K * chi3 + \
        coeff10K * eta + coeff11K * eta * chi + coeff12K * eta * chi2 + \
        coeff13K * eta * chi3 + coeff20K * eta2 + coeff21K * eta2 * chi + \
        coeff22K * eta2 * chi2 + coeff23K * eta2 * chi3 + coeff30K * eta3 + \
        coeff31K * eta3 * chi + coeff32K * eta3 * chi2 + coeff33K * eta3 * chi3
    
    dSO = \
            coeff00dSO + coeff01dSO * chi + coeff02dSO * chi2 + coeff03dSO * chi3 + \
            coeff10dSO * eta + coeff11dSO * eta * chi + coeff12dSO * eta * chi2 + \
            coeff13dSO * eta * chi3 + coeff20dSO * eta2 + coeff21dSO * eta2 * chi + \
            coeff22dSO * eta2 * chi2 + coeff23dSO * eta2 * chi3 + coeff30dSO * eta3 + \
            coeff31dSO * eta3 * chi + coeff32dSO * eta3 * chi2 + coeff33dSO * eta3 * chi3
    
    dSS = \
            coeff00dSS + coeff01dSS * chi + coeff02dSS * chi2 + coeff03dSS * chi3 + \
            coeff10dSS * eta + coeff11dSS * eta * chi + coeff12dSS * eta * chi2 + \
            coeff13dSS * eta * chi3 + coeff20dSS * eta2 + coeff21dSS * eta2 * chi + \
            coeff22dSS * eta2 * chi2 + coeff23dSS * eta2 * chi3 + coeff30dSS * eta3 + \
            coeff31dSS * eta3 * chi + coeff32dSS * eta3 * chi2 + coeff33dSS * eta3 * chi3
    
    dtPeak = \
        coeff00DT + coeff01DT * chi + coeff02DT * chi2 + coeff03DT * chi3 + \
        coeff10DT * eta + coeff11DT * eta * chi + coeff12DT * eta * chi2 + \
        coeff13DT * eta * chi3 + coeff20DT * eta2 + coeff21DT * eta2 * chi + \
        coeff22DT * eta2 * chi2 + coeff23DT * eta2 * chi3 + coeff30DT * eta3 + \
        coeff31DT * eta3 * chi + coeff32DT * eta3 * chi2 + coeff33DT * eta3 * chi3
    dic = {}
    dic['KK'] = KK
    dic['dSO'] = dSO
    dic['dSS'] = dSS
    dic['dtPeak'] = dtPeak
    return SEOBHCoeffs(dic)

class SEOBHCoeffs(object):
    def __init__(self, valdic = None):
        if valdic:
            self._space = valdic
            for (attr, val) in valdic.items():
                setattr(self, attr, val)
        else:
            self._space = {}

    def append(self, attr, val):
        setattr(self, attr, val)
        self._space[attr] = val

    def get_values(self):
        return [val for val in self._space.values()]

    def get_values_with_name(self):
        return [(attr, val) for (attr, val) in self._space.items()]


class SEOBHCoeffsCalibrator(object):
    def __init__(self, approx):
        for case in switch(approx):
            if case(1) or case('v1') or case('SEOBNRv1'):
                self._core = SEOBHyperCoefficients_v1
                break
            if case(2) or case('v2') or case('SEOBNRv2'):
                self._core = SEOBHyperCoefficients_v2
                break
            if case(4) or case('v4') or case('SEOBNRv4'):
                self._core = SEOBHyperCoefficients_v4
                break
            raise Exception(f'Invalid version tag {approx}')
    
    def __call__(self, m1, m2, spin1z, spin2z):
        eta = m1 * m2 / (m1+m2) / (m1+m2)
        a = (spin1z*m1*m1 + spin2z*m2*m2) / (m1+m2) / (m1+m2)
        return self._core(eta, a)