/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYDATATYPES__
#define __INCLUDE_DYDATATYPES__

#include "myDatatypes.h"

/* ECC */
typedef struct tagPNEccCorrectCoeffs
{
    COMPLEX16 cv0em;
    COMPLEX16 cv0ep;

    COMPLEX16 cv2em;
    COMPLEX16 cv2ep;

    COMPLEX16 cv3em;
    COMPLEX16 cv3ep;

    COMPLEX16 cv4em;
    COMPLEX16 cv4ep;

    COMPLEX16 cv5em;
    COMPLEX16 cv5ep;

    COMPLEX16 cv6em;
    COMPLEX16 cv6ep;
}PNEccCorrectCoeffs;


/* NQC */
typedef struct
tagEOBNonQCCoeffs
{
    REAL8 a1;
    REAL8 a2;
    REAL8 a3;
    REAL8 a3S;
    REAL8 a4;
    REAL8 a5;
    REAL8 b1;
    REAL8 b2;
    REAL8 b3;
    REAL8 b4;
} EOBNonQCCoeffs;

/* Spin EOBH Coeffs */
typedef struct
tagSpinEOBHCoeffs
{
    double KK;
    double k0;
    double k1;
    double k2;
    double k3;
    double k4;
    double k5;
    double k5l;
    double b3;
    double bb3;
    double d1;
    double d1v2;
    double dheffSS;
    double dheffSSv2;
}
SpinEOBHCoeffs;

/* Factorized waveform Coeffs */
typedef struct
tagFacWaveformCoeffs
{
    REAL8 delta22vh3;
    REAL8 delta22vh6;
    REAL8 delta22vh6S;
    REAL8 delta22v8;
    REAL8 delta22v8S;
    REAL8 delta22vh9;
    REAL8 delta22v5;
    REAL8 delta22v6;
    REAL8 delta22v6S;
    
    REAL8 rho22v2;
    REAL8 rho22v3;
    REAL8 rho22v3S;
    REAL8 rho22v4;
    REAL8 rho22v4S;
    REAL8 rho22v5;
    REAL8 rho22v5S;
    REAL8 rho22v6;
    REAL8 rho22v6S;
    REAL8 rho22v6l;
    REAL8 rho22v7;
    REAL8 rho22v7S;
    REAL8 rho22v8;
    REAL8 rho22v8S;
    REAL8 rho22v8l;
    REAL8 rho22v10;
    REAL8 rho22v10l;
    
    REAL8 delta21vh3;
    REAL8 delta21vh6;
    REAL8 delta21vh6S;
    REAL8 delta21vh7;
    REAL8 delta21vh7S;
    REAL8 delta21vh9;
    REAL8 delta21v5;
    REAL8 delta21v7;
    
    REAL8 rho21v1;
    REAL8 rho21v2;
    REAL8 rho21v2S;
    REAL8 rho21v3;
    REAL8 rho21v3S;
    REAL8 rho21v4;
    REAL8 rho21v4S;
    REAL8 rho21v5;
    REAL8 rho21v5S;
    REAL8 rho21v6;
    REAL8 rho21v6S;
    REAL8 rho21v6l;
    REAL8 rho21v7;
    REAL8 rho21v7S;
    REAL8 rho21v7l;
    REAL8 rho21v7lS;
    REAL8 rho21v8;
    REAL8 rho21v8l;
    REAL8 rho21v10;
    REAL8 rho21v10l;
    
    REAL8 f21v1;
    REAL8 f21v1S;
    REAL8 f21v3;
    REAL8 f21v3S;
    REAL8 f21v4;
    REAL8 f21v5;
    REAL8 f21v6;
    REAL8 f21v7c;

    REAL8 delta33vh3;
    REAL8 delta33vh6;
    REAL8 delta33vh6S;
    REAL8 delta33vh9;
    REAL8 delta33v5;
    REAL8 delta33v7;
    
    REAL8 rho33v2;
    REAL8 rho33v3;
    REAL8 rho33v4;
    REAL8 rho33v4S;
    REAL8 rho33v5;
    REAL8 rho33v5S;
    REAL8 rho33v6;
    REAL8 rho33v6S;
    REAL8 rho33v6l;
    REAL8 rho33v7;
    REAL8 rho33v7S;
    REAL8 rho33v8;
    REAL8 rho33v8l;
    REAL8 rho33v10;
    REAL8 rho33v10l;

    REAL8 f33v3;
    REAL8 f33v4;
    REAL8 f33v5;
    REAL8 f33v6;
    REAL8 f33v3S;
    REAL8 f33vh6;

    REAL8 delta32vh3;
    REAL8 delta32vh4;
    REAL8 delta32vh4S;
    REAL8 delta32vh6;
    REAL8 delta32vh6S;
    REAL8 delta32vh9;
    
    REAL8 rho32v;
    REAL8 rho32vS;
    REAL8 rho32v2;
    REAL8 rho32v2S;
    REAL8 rho32v3;
    REAL8 rho32v3S;
    REAL8 rho32v4;
    REAL8 rho32v4S;
    REAL8 rho32v5;
    REAL8 rho32v5S;
    REAL8 rho32v6;
    REAL8 rho32v6S;
    REAL8 rho32v6l;
    REAL8 rho32v8;
    REAL8 rho32v8l;
    
    REAL8 delta31vh3;
    REAL8 delta31vh6;
    REAL8 delta31vh6S;
    REAL8 delta31vh7;
    REAL8 delta31vh7S;
    REAL8 delta31vh9;
    REAL8 delta31v5;
    
    REAL8 rho31v2;
    REAL8 rho31v3;
    REAL8 rho31v4;
    REAL8 rho31v4S;
    REAL8 rho31v5;
    REAL8 rho31v5S;
    REAL8 rho31v6;
    REAL8 rho31v6S;
    REAL8 rho31v6l;
    REAL8 rho31v7;
    REAL8 rho31v7S;
    REAL8 rho31v8;
    REAL8 rho31v8l;
    
    REAL8 f31v3;
    REAL8 f31v3S;
    
    REAL8 delta44vh3;
    REAL8 delta44vh6;
    REAL8 delta44vh9;
    REAL8 delta44vh6S;
    REAL8 delta44v5;

    REAL8 rho44v2;
    REAL8 rho44v3;
    REAL8 rho44v3S;
    REAL8 rho44v4;
    REAL8 rho44v4S;
    REAL8 rho44v5;
    REAL8 rho44v5S;
    REAL8 rho44v6;
    REAL8 rho44v6S;
    REAL8 rho44v6l;
    REAL8 rho44v8;
    REAL8 rho44v8l;
    REAL8 rho44v10;
    REAL8 rho44v10l;

    REAL8 delta43vh3;
    REAL8 delta43vh4;
    REAL8 delta43vh4S;
    REAL8 delta43vh6;
    
    REAL8 rho43v;
    REAL8 rho43v2;
    REAL8 rho43v4;
    REAL8 rho43v4S;
    REAL8 rho43v5;
    REAL8 rho43v5S;
    REAL8 rho43v6;
    REAL8 rho43v6l;
    
    REAL8 f43v;
    REAL8 f43vS;
    
    REAL8 delta42vh3;
    REAL8 delta42vh6;
    REAL8 delta42vh6S;
    
    REAL8 rho42v2;
    REAL8 rho42v3;
    REAL8 rho42v3S;
    REAL8 rho42v4;
    REAL8 rho42v4S;
    REAL8 rho42v5;
    REAL8 rho42v5S;
    REAL8 rho42v6;
    REAL8 rho42v6S;
    REAL8 rho42v6l;
    
    REAL8 delta41vh3;
    REAL8 delta41vh4;
    REAL8 delta41vh4S;
    REAL8 delta41vh6;
    
    REAL8 rho41v;
    REAL8 rho41v2;
    REAL8 rho41v4;
    REAL8 rho41v4S;
    REAL8 rho41v5;
    REAL8 rho41v5S;
    REAL8 rho41v6;
    REAL8 rho41v6l;
    
    REAL8 f41v;
    REAL8 f41vS;
    
    REAL8 delta55vh3;
    REAL8 delta55vh6;
    REAL8 delta55vh9;
    REAL8 delta55v5;
    REAL8 rho55v2;
    REAL8 rho55v3;
    REAL8 rho55v3S;
    REAL8 rho55v4;
    REAL8 rho55v4S;
    REAL8 rho55v5;
    REAL8 rho55v5S;
    REAL8 rho55v6;
    REAL8 rho55v6l;
    REAL8 rho55v8;
    REAL8 rho55v8l;
    REAL8 rho55v10;
    REAL8 rho55v10l;
    REAL8 f55v3;
    REAL8 f55v4;
    REAL8 f55v5c;

    REAL8 delta54vh3;
    REAL8 delta54vh4;
    REAL8 delta54vh4S;
    REAL8 rho54v2;
    REAL8 rho54v3;
    REAL8 rho54v3S;
    REAL8 rho54v4;
    REAL8 rho54v4S;
    
    REAL8 delta53vh3;
    REAL8 rho53v2;
    REAL8 rho53v3;
    REAL8 rho53v3S;
    REAL8 rho53v4;
    REAL8 rho53v4S;
    REAL8 rho53v5;
    REAL8 rho53v5S;
    
    REAL8 delta52vh3;
    REAL8 delta52vh4;
    REAL8 delta52vh4S;
    REAL8 rho52v2;
    REAL8 rho52v3;
    REAL8 rho52v3S;
    REAL8 rho52v4;
    REAL8 rho52v4S;
    
    REAL8 delta51vh3;
    REAL8 rho51v2;
    REAL8 rho51v3;
    REAL8 rho51v3S;
    REAL8 rho51v4;
    REAL8 rho51v4S;
    REAL8 rho51v5;
    REAL8 rho51v5S;
    
    REAL8 delta66vh3;
    REAL8 rho66v2;
    REAL8 rho66v3;
    REAL8 rho66v3S;
    REAL8 rho66v4;
    REAL8 rho66v4S;
    
    REAL8 delta65vh3;
    REAL8 rho65v2;
    REAL8 rho65v3;
    REAL8 rho65v3S;
    
    REAL8 delta64vh3;
    REAL8 rho64v2;
    REAL8 rho64v3;
    REAL8 rho64v3S;
    REAL8 rho64v4;
    REAL8 rho64v4S;
    
    REAL8 delta63vh3;
    REAL8 rho63v2;
    REAL8 rho63v3;
    REAL8 rho63v3S;
    
    REAL8 delta62vh3;
    REAL8 rho62v2;
    REAL8 rho62v3;
    REAL8 rho62v3S;
    REAL8 rho62v4;
    REAL8 rho62v4S;
    
    REAL8 delta61vh3;
    REAL8 rho61v2;
    REAL8 rho61v3;
    REAL8 rho61v3S;
    
    REAL8 delta77vh3;
    REAL8 rho77v2;
    REAL8 rho77v3;
    REAL8 rho77v3S;
    
    REAL8 rho76v2;
    
    REAL8 delta75vh3;
    REAL8 rho75v2;
    REAL8 rho75v3;
    REAL8 rho75v3S;
    
    REAL8 rho74v2;
    
    REAL8 delta73vh3;
    REAL8 rho73v2;
    REAL8 rho73v3;
    REAL8 rho73v3S;
    
    REAL8 rho72v2;
    
    REAL8 delta71vh3;
    REAL8 rho71v2;
    REAL8 rho71v3;
    REAL8 rho71v3S;
    
    REAL8 rho88v2;
    REAL8 rho87v2;
    REAL8 rho86v2;
    REAL8 rho85v2;
    REAL8 rho84v2;
    REAL8 rho83v2;
    REAL8 rho82v2;
    REAL8 rho81v2;
}
FacWaveformCoeffs;

/* Newton Multipole Prefix */
typedef struct tagNewtonMultipolePrefixes
{
    COMPLEX16 values[8+1][8+1];
}
NewtonMultipolePrefixes;


typedef struct tagEOBParams
{
    REAL8 m1;
    REAL8 m2;
    REAL8 eta; // mass ratio
    REAL8 Mtotal; // total mass
    REAL8 rad;
    REAL8 omega;
    INT omegaPeaked;
    FacWaveformCoeffs       *hCoeffs;
    NewtonMultipolePrefixes *prefixes;
} EOBParams;

/* Spin EOB basic */
typedef struct
tagSpinEOBParams
{
    EOBParams               *eobParams;
    EOBNonQCCoeffs          *nqcCoeffs;
    SpinEOBHCoeffs          *seobCoeffs;
    PNEccCorrectCoeffs      *eccCoeffs;
    REAL8Vector             *s1Vec; // s1VecOverMtMt
    REAL8Vector             *s2Vec; // s2VecOverMtMt
    REAL8Vector             *s1VecOverMtMt;
    REAL8Vector             *s2VecOverMtMt;
    REAL8Vector             *sigmaStar;
    REAL8Vector             *sigmaKerr;
    REAL8                   a;
    REAL8                   chi1; // s1z
    REAL8                   chi2; // s2z
    REAL8                   deltaT;
    REAL8                   eccentricity;
    int                     tortoise;
}
SpinEOBParams;


typedef struct tagCtrlParams
{
    gboolean verbose;
    CHAR dump[256];
} CtrlParams;

typedef struct tagAdjParams
{
    REAL8 KK;
    REAL8 dSO;
    REAL8 dSS;
    REAL8 dtPeak;
} AdjParams;

typedef
struct tagHcapDerivParams
{
    const REAL8   *values;
    SpinEOBParams *params;
    UINT         varyParam;
}
HcapDerivParams;

typedef
struct tagHcapSphDeriv2Params
{
  const REAL8     *sphValues;
  SpinEOBParams   *params;
  UINT           varyParam1;
  UINT           varyParam2;
}
HcapSphDeriv2Params;

typedef
struct tagSpinEOBDynamics
{
    UINT length;
    REAL8Vector *tVec;
    REAL8Vector *rVec;
    REAL8Vector *drVec;
    REAL8Vector *phiVec;
    REAL8Vector *dphiVec;
    REAL8Vector *prVec;
    REAL8Vector *dprVec;
    REAL8Vector *pPhiVec;
    REAL8Vector *dpPhiVec;
}SpinEOBDynamics;

typedef
struct tagNRPeakParams
{
    REAL8 ampPeak;
    REAL8 ampPeakDot;
    REAL8 ampPeakDDot;
    REAL8 omegaPeak;
    REAL8 omegaPeakDot;
}NRPeakParams;

#endif

