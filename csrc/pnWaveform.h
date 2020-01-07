/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNWAVEFORM__
#define __INCLUDE_PNWAVEFORM__

#include "dyUtils.h"

typedef struct
tagPNEccentricWaveformCoeffs
{
    REAL8 hCoeffv00;
    REAL8 hCoeffv0em;
    REAL8 hCoeffv0ep;

    REAL8 hCoeffv20;
    REAL8 hCoeffv2em;
    REAL8 hCoeffv2ep;

    REAL8 hCoeffv30;
    REAL8 hCoeffv3em;
    REAL8 hCoeffv3ep;

    REAL8 hCoeffv40;
    REAL8 hCoeffv4em;
    REAL8 hCoeffv4ep;

    REAL8 hCoeffv50;
    REAL8 hCoeffv5em;
    REAL8 hCoeffv5ep;

    REAL8 hCoeffv60;
    REAL8 hCoeffv6em;
    REAL8 hCoeffv6ep;

}PNEccentricWaveformCoeffs;

typedef struct
tagPNwaveformCoeffs
{
    COMPLEX16 Pre;
    COMPLEX16 Coeff0;
    COMPLEX16 Coeff1;
    COMPLEX16 Coeff2;
    COMPLEX16 Coeff3;
    COMPLEX16 Coeff4;
    COMPLEX16 Coeff5;
    COMPLEX16 Coeff6;
}PNwaveformCoeffs;

typedef struct
tagCommhCoeffs
{
    // l=2
    PNwaveformCoeffs *coeffs20;
    PNwaveformCoeffs *coeffs21;
    PNwaveformCoeffs *coeffs22;
    // l=3
    PNwaveformCoeffs *coeffs30;
    PNwaveformCoeffs *coeffs31;
    PNwaveformCoeffs *coeffs32;
    PNwaveformCoeffs *coeffs33;
    // l=4
    PNwaveformCoeffs *coeffs40;
    PNwaveformCoeffs *coeffs41;
    PNwaveformCoeffs *coeffs42;
    PNwaveformCoeffs *coeffs43;
    PNwaveformCoeffs *coeffs44;
    // l=5
    PNwaveformCoeffs *coeffs50;
    PNwaveformCoeffs *coeffs51;
    PNwaveformCoeffs *coeffs52;
    PNwaveformCoeffs *coeffs53;
    PNwaveformCoeffs *coeffs54;
    PNwaveformCoeffs *coeffs55;
    // l=6
    PNwaveformCoeffs *coeffs60;
    PNwaveformCoeffs *coeffs61;
    PNwaveformCoeffs *coeffs62;
    PNwaveformCoeffs *coeffs63;
    PNwaveformCoeffs *coeffs64;
    PNwaveformCoeffs *coeffs65;
    PNwaveformCoeffs *coeffs66;
    // l=7
    PNwaveformCoeffs *coeffs70;
    PNwaveformCoeffs *coeffs71;
    PNwaveformCoeffs *coeffs72;
    PNwaveformCoeffs *coeffs73;
    PNwaveformCoeffs *coeffs74;
    PNwaveformCoeffs *coeffs75;
    PNwaveformCoeffs *coeffs76;
    PNwaveformCoeffs *coeffs77;
    // l=8
    PNwaveformCoeffs *coeffs80;
    PNwaveformCoeffs *coeffs81;
    PNwaveformCoeffs *coeffs82;
    PNwaveformCoeffs *coeffs83;
    PNwaveformCoeffs *coeffs84;
    PNwaveformCoeffs *coeffs85;
    PNwaveformCoeffs *coeffs86;
    PNwaveformCoeffs *coeffs87;
    PNwaveformCoeffs *coeffs88;

}CommhCoeffs;

#define WaveformCoeffsInit(hCoeffs) _WaveformCoeffisInit(hCoeffs, __UNIQUE_ID(_id))
#define _WaveformCoeffisInit(hCoeffs, x) \
CommhCoeffs hCoeffs; \
PNwaveformCoeffs __PASTE(c20, x); \
PNwaveformCoeffs __PASTE(c21, x); \
PNwaveformCoeffs __PASTE(c22, x); \
memset (&__PASTE(c20, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c21, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c22, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs20 = &__PASTE(c20, x); \
hCoeffs.coeffs21 = &__PASTE(c21, x); \
hCoeffs.coeffs22 = &__PASTE(c22, x); \
PNwaveformCoeffs __PASTE(c30, x); \
PNwaveformCoeffs __PASTE(c31, x); \
PNwaveformCoeffs __PASTE(c32, x); \
PNwaveformCoeffs __PASTE(c33, x); \
memset (&__PASTE(c30, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c31, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c32, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c33, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs30 = &__PASTE(c30, x); \
hCoeffs.coeffs31 = &__PASTE(c31, x); \
hCoeffs.coeffs32 = &__PASTE(c32, x); \
hCoeffs.coeffs33 = &__PASTE(c33, x); \
PNwaveformCoeffs __PASTE(c40, x); \
PNwaveformCoeffs __PASTE(c41, x); \
PNwaveformCoeffs __PASTE(c42, x); \
PNwaveformCoeffs __PASTE(c43, x); \
PNwaveformCoeffs __PASTE(c44, x); \
memset (&__PASTE(c40, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c41, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c42, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c43, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c44, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs40 = &__PASTE(c40, x); \
hCoeffs.coeffs41 = &__PASTE(c41, x); \
hCoeffs.coeffs42 = &__PASTE(c42, x); \
hCoeffs.coeffs43 = &__PASTE(c43, x); \
hCoeffs.coeffs44 = &__PASTE(c44, x); \
PNwaveformCoeffs __PASTE(c50, x); \
PNwaveformCoeffs __PASTE(c51, x); \
PNwaveformCoeffs __PASTE(c52, x); \
PNwaveformCoeffs __PASTE(c53, x); \
PNwaveformCoeffs __PASTE(c54, x); \
PNwaveformCoeffs __PASTE(c55, x); \
memset (&__PASTE(c50, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c51, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c52, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c53, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c54, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c55, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs50 = &__PASTE(c50, x); \
hCoeffs.coeffs51 = &__PASTE(c51, x); \
hCoeffs.coeffs52 = &__PASTE(c52, x); \
hCoeffs.coeffs53 = &__PASTE(c53, x); \
hCoeffs.coeffs54 = &__PASTE(c54, x); \
hCoeffs.coeffs55 = &__PASTE(c55, x); \
PNwaveformCoeffs __PASTE(c60, x); \
PNwaveformCoeffs __PASTE(c61, x); \
PNwaveformCoeffs __PASTE(c62, x); \
PNwaveformCoeffs __PASTE(c63, x); \
PNwaveformCoeffs __PASTE(c64, x); \
PNwaveformCoeffs __PASTE(c65, x); \
PNwaveformCoeffs __PASTE(c66, x); \
memset (&__PASTE(c60, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c61, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c62, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c63, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c64, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c65, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c66, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs60 = &__PASTE(c60, x); \
hCoeffs.coeffs61 = &__PASTE(c61, x); \
hCoeffs.coeffs62 = &__PASTE(c62, x); \
hCoeffs.coeffs63 = &__PASTE(c63, x); \
hCoeffs.coeffs64 = &__PASTE(c64, x); \
hCoeffs.coeffs65 = &__PASTE(c65, x); \
hCoeffs.coeffs66 = &__PASTE(c66, x); \
PNwaveformCoeffs __PASTE(c70, x); \
PNwaveformCoeffs __PASTE(c71, x); \
PNwaveformCoeffs __PASTE(c72, x); \
PNwaveformCoeffs __PASTE(c73, x); \
PNwaveformCoeffs __PASTE(c74, x); \
PNwaveformCoeffs __PASTE(c75, x); \
PNwaveformCoeffs __PASTE(c76, x); \
PNwaveformCoeffs __PASTE(c77, x); \
memset (&__PASTE(c70, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c71, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c72, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c73, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c74, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c75, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c76, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c77, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs70 = &__PASTE(c70, x); \
hCoeffs.coeffs71 = &__PASTE(c71, x); \
hCoeffs.coeffs72 = &__PASTE(c72, x); \
hCoeffs.coeffs73 = &__PASTE(c73, x); \
hCoeffs.coeffs74 = &__PASTE(c74, x); \
hCoeffs.coeffs75 = &__PASTE(c75, x); \
hCoeffs.coeffs76 = &__PASTE(c76, x); \
hCoeffs.coeffs77 = &__PASTE(c77, x); \
PNwaveformCoeffs __PASTE(c80, x); \
PNwaveformCoeffs __PASTE(c81, x); \
PNwaveformCoeffs __PASTE(c82, x); \
PNwaveformCoeffs __PASTE(c83, x); \
PNwaveformCoeffs __PASTE(c84, x); \
PNwaveformCoeffs __PASTE(c85, x); \
PNwaveformCoeffs __PASTE(c86, x); \
PNwaveformCoeffs __PASTE(c87, x); \
PNwaveformCoeffs __PASTE(c88, x); \
memset (&__PASTE(c80, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c81, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c82, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c83, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c84, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c85, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c86, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c87, x), 0, sizeof(PNwaveformCoeffs)); \
memset (&__PASTE(c88, x), 0, sizeof(PNwaveformCoeffs)); \
hCoeffs.coeffs80 = &__PASTE(c80, x); \
hCoeffs.coeffs81 = &__PASTE(c81, x); \
hCoeffs.coeffs82 = &__PASTE(c82, x); \
hCoeffs.coeffs83 = &__PASTE(c83, x); \
hCoeffs.coeffs84 = &__PASTE(c84, x); \
hCoeffs.coeffs85 = &__PASTE(c85, x); \
hCoeffs.coeffs86 = &__PASTE(c86, x); \
hCoeffs.coeffs87 = &__PASTE(c87, x); \
hCoeffs.coeffs88 = &__PASTE(c88, x); 


INT Waveform_CQG_25_165003(REAL8 v, /*(omega)^1/3*/
                            REAL8 phi,
                            REAL8 eta,
                            INT l,
                            INT m,
                            CommhCoeffs *hCoeffs,
                            COMPLEX16 *out);

INT CalculateWaveformCoeffs_CQG_25_165003(REAL8 eta,
                                          CommhCoeffs *hCoeffs);

#endif

