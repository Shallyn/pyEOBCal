/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pnWaveform.h"

/* PN Flux */

/* QuasiCircular */
INT CalculateWaveformCoeffs_CQG_25_165003(REAL8 eta,
                                          CommhCoeffs *hCoeffs)
{
    INT l,m;
    if (eta < 0 || eta > 0.25)
    {
        return CEV_FAILURE;
    }
    REAL8 dM = sqrt(1-4*eta);
    REAL8 eta2, eta3;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    REAL8 PIsq ;
    PIsq = CST_PI * CST_PI;
    REAL8 LN32 = log(3/2), LN52 = log(5/2);
    for(l=2; l<=8; l++)
    {
        for(m=0; m<=l; m++)
        {
            switch(l)
            {
                case 2:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            // mode20 (9.4c)
                            hCoeffs->coeffs20->Pre = 1;
                            hCoeffs->coeffs20->Coeff0 = -5/14/sqrt(6);
                            break;
                        }
                        case 1:
                        {
                            // mode21 (9.4b)
                            hCoeffs->coeffs21->Pre = I * dM / 3;
                            hCoeffs->coeffs21->Coeff1 = 1;
                            hCoeffs->coeffs21->Coeff3 = -17/28 + 5*eta / 7;
                            hCoeffs->coeffs21->Coeff4 = CST_PI + I * (-1/2 - 2*CST_LN2);
                            hCoeffs->coeffs21->Coeff5 = -43/126 - 509*eta/126 + 79*eta2/168;
                            hCoeffs->coeffs21->Coeff6 = -17*CST_PI/28 + 3*CST_PI*eta/14 + 
                                I*(17/56 + eta*(-353/28-3*CST_LN2/7) + 17*CST_LN2/14);
                            break;
                        }
                        case 2:
                        {
                            // mode22 (9.4a)
                            hCoeffs->coeffs22->Pre = 1;
                            hCoeffs->coeffs22->Coeff0 = 1;
                            hCoeffs->coeffs22->Coeff2 = -107/42 + 55*eta/42;
                            hCoeffs->coeffs22->Coeff3 = 2*CST_PI;
                            hCoeffs->coeffs22->Coeff4 = -2173/1512 -1069*eta/216 +2047*eta2/1512;
                            hCoeffs->coeffs22->Coeff5 = -107*CST_PI/21 - 24*I*eta + 34*CST_PI*eta/21;
                            hCoeffs->coeffs22->Coeff6 = 27027409/646800 - 856*CST_GAMMA/105 + 428*I*CST_PI/105 + 2*PIsq/3+
                                (-278185/33264+41*PIsq/96)*eta - 20261*eta2/2772 + 114635*eta3/99792; // - 428*ln(16x)/105
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                case 3:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            hCoeffs->coeffs30->Pre = 1;
                            hCoeffs->coeffs30->Coeff5 = -2*I*eta*sqrt(6/7);
                            break;
                        }
                        case 1:
                        {
                            hCoeffs->coeffs31->Pre = I*dM/12/sqrt(14);
                            hCoeffs->coeffs31->Coeff1 = 1;
                            hCoeffs->coeffs31->Coeff3 = -8/3 - 2*eta/3;
                            hCoeffs->coeffs31->Coeff4 = CST_PI + I * (-7/5 - 2*CST_LN2);
                            hCoeffs->coeffs31->Coeff5 = 607/198 - 136*eta/99 - 247*eta2/198;
                            hCoeffs->coeffs31->Coeff6 = -8*CST_PI/3 - 7*CST_PI*eta/6 +
                                I*(56/15 + 16*CST_LN2/3 + eta*(-1/15 + 7*CST_LN2/2) );
                            break;
                        }
                        case 2:
                        {
                            hCoeffs->coeffs32->Pre = (1/3)*sqrt(5/7);
                            hCoeffs->coeffs32->Coeff2 = 1-3*eta;
                            hCoeffs->coeffs32->Coeff4 = (-193/90 + 145*eta/18 -73*eta2/18);
                            hCoeffs->coeffs32->Coeff5 = 2*CST_PI - 6*CST_PI*eta + I*(-3+66*eta/5);
                            hCoeffs->coeffs32->Coeff6 = -1451/3960 -17387*eta/3960 +5557*eta2/220 -5341*eta3/1320;
                            break;
                        }
                        case 3:
                        {
                            hCoeffs->coeffs33->Pre = -3*I*sqrt(15/14)*dM/4;
                            hCoeffs->coeffs33->Coeff1 = 1;
                            hCoeffs->coeffs33->Coeff3 = -4 + 2*eta;
                            hCoeffs->coeffs33->Coeff4 = 3*CST_PI + I*(-21/5 + 6*log(3/2));
                            hCoeffs->coeffs33->Coeff5 = 132/110 - 1838*eta/165 + 887*eta2/330;
                            hCoeffs->coeffs33->Coeff6 = -12*CST_PI + 9*CST_PI*eta +
                                I*(84/5 - 24*log(3/2) + eta*(-48103/1215 + 9*LN32) );
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                case 4:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            hCoeffs->coeffs40->Pre = 1/504/sqrt(2);
                            break;
                        }
                        case 1:
                        {
                            hCoeffs->coeffs41->Pre = I*dM/84/sqrt(10);
                            hCoeffs->coeffs41->Coeff3 = 1 - 2*eta;
                            hCoeffs->coeffs41->Coeff5 = -101/33 + 337*eta/44 - 83*eta2/33;
                            hCoeffs->coeffs41->Coeff6 = CST_PI - CST_2PI*eta +
                                I*(-32/15 - 2*CST_LN2 + eta*(1661/30 + 4*CST_LN2));
                            break;
                        }
                        case 2:
                        {
                            hCoeffs->coeffs42->Pre = sqrt(5) / 63;
                            hCoeffs->coeffs42->Coeff2 = 1-3*eta;
                            hCoeffs->coeffs42->Coeff4 = -437/110 + 805*eta/66 - 19*eta2/22;
                            hCoeffs->coeffs42->Coeff5 = CST_2PI - 3*CST_2PI*eta + I*(-21/5 + 84*eta/5);
                            hCoeffs->coeffs42->Coeff6 = 1038039/200200 - 606751*eta/28600 + 400457*eta2/25740 + 25783*eta3/17160;
                            break;
                        }
                        case 3:
                        {
                            hCoeffs->coeffs43->Pre = -9*I*dM/4/sqrt(70);
                            hCoeffs->coeffs43->Coeff3 = 1-2*eta;
                            hCoeffs->coeffs43->Coeff5 = -39/11 + 1267*eta/132 - 131*eta2/33;
                            hCoeffs->coeffs43->Coeff6 = 3*CST_PI - 3*CST_2PI*eta + 
                                I*(-32/5 + eta*(16301/810 - 12*LN32) + 6*LN32);
                            break;
                        }
                        case 4:
                        {
                            hCoeffs->coeffs44->Pre = -8*sqrt(5/7)/9;
                            hCoeffs->coeffs44->Coeff2 = 1-3*eta;
                            hCoeffs->coeffs44->Coeff4 = -593/110 + 1273*eta/66 -175*eta2/22;
                            hCoeffs->coeffs44->Coeff5 = 4*CST_PI - 12*CST_PI*eta + 
                                I*(-42/5 + eta*(1193/40 - 24*CST_LN2) +8*CST_LN2);
                            hCoeffs->coeffs44->Coeff6 = 1068671/200200 - 1088119*eta/28600 + 146879*eta2/2340 - 226097*eta3/17160;
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                case 5:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            break;
                        }
                        case 1:
                        {
                            hCoeffs->coeffs51->Pre = I*dM/288/sqrt(385);
                            hCoeffs->coeffs51->Coeff3 = 1-2*eta;
                            hCoeffs->coeffs51->Coeff5 = -179/39 + 352*eta/39 -4*eta2/39;
                            hCoeffs->coeffs51->Coeff6 = CST_PI-CST_2PI*eta + 
                                I*(-181/70 -2*CST_LN2 + eta*(626/5+4*CST_LN2) );
                            break;
                        }
                        case 2:
                        {
                            hCoeffs->coeffs52->Pre = 2/27/sqrt(55);
                            hCoeffs->coeffs52->Coeff4 = 1-5*eta+5*eta2;
                            hCoeffs->coeffs52->Coeff6 = -3911/910 + 3079*eta/130 -413*eta2/13+231*eta3/26;
                            break;
                        }
                        case 3:
                        {
                            hCoeffs->coeffs53->Pre = -I*9*sqrt(3/110)*dM/32;
                            hCoeffs->coeffs53->Coeff3 = 1-2*eta;
                            hCoeffs->coeffs53->Coeff5 = -69/13+464*eta/39 -88*eta2/39;
                            hCoeffs->coeffs53->Coeff6 = 3*CST_PI-6*CST_PI*eta +
                                I*(-543/70+eta*(83702/3645 -12*LN32) +6*LN32 );
                            break;
                        }
                        case 4:
                        {
                            hCoeffs->coeffs54->Pre = -32/9/sqrt(165);
                            hCoeffs->coeffs54->Coeff4 = 1-5*eta+5*eta2;
                            hCoeffs->coeffs54->Coeff6 = -4451/910+3691*eta/130-521*eta2/13+339*eta3/26;
                            break;
                        }
                        case 5:
                        {
                            hCoeffs->coeffs55->Pre = 625*I*dM/96/sqrt(66);
                            hCoeffs->coeffs55->Coeff3 = 1-2*eta;
                            hCoeffs->coeffs55->Coeff5 = -263/39 + 688*eta/39 -256*eta2/39;
                            hCoeffs->coeffs55->Coeff6 = 5*CST_PI -10*CST_PI*eta +
                                I*(-181/14 + eta*(105834/3125 -20*LN52) +10*LN52 );
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                case 6:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            break;
                        }
                        case 1:
                        {
                            hCoeffs->coeffs61->Pre = 1;
                            hCoeffs->coeffs61->Coeff5 = I*dM*(1-4*eta+3*eta2)/8316/sqrt(26);
                            break;
                        }
                        case 2:
                        {
                            hCoeffs->coeffs62->Pre = 2/297/sqrt(65);
                            hCoeffs->coeffs62->Coeff4 = 1-5*eta+5*eta2;
                            hCoeffs->coeffs62->Coeff6 = -81/14 +59*eta/2 -32*eta2 +7*eta3/2;
                            break;
                        }
                        case 3:
                        {
                            hCoeffs->coeffs63->Pre = 1;
                            hCoeffs->coeffs63->Coeff5 = -81*dM*(1-4*eta+3*eta2)/616/sqrt(65);
                            break;
                        }
                        case 4:
                        {
                            hCoeffs->coeffs64->Pre = -128*sqrt(2/39)/495;
                            hCoeffs->coeffs64->Coeff4 = 1-5*eta+5*eta2;
                            hCoeffs->coeffs64->Coeff6 = -93/14 +71*eta/2 -44*eta2 +19*eta3/2;
                            break;
                        }
                        case 5:
                        {
                            hCoeffs->coeffs65->Pre = 1;
                            hCoeffs->coeffs65->Coeff5 = 3125*I*dM*(1-4*eta+3*eta2)/504/sqrt(429);
                            break;
                        }
                        case 6:
                        {
                            hCoeffs->coeffs66->Pre = 54/5/sqrt(143);
                            hCoeffs->coeffs66->Coeff4 = 1-5*eta+5*eta2;
                            hCoeffs->coeffs66->Coeff6 = -113/14 + 91*eta/2 -64*eta2 + 39*eta3/2;
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                case 7:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            break;
                        }
                        case 1:
                        {
                            hCoeffs->coeffs71->Pre = 1;
                            hCoeffs->coeffs71->Coeff5 = I*dM*(1-4*eta+3*eta2)/864864/sqrt(2);
                            break;
                        }
                        case 2:
                        {
                            hCoeffs->coeffs72->Pre = 1;
                            hCoeffs->coeffs72->Coeff6 = (1-7*eta+14*eta2-7*eta3)/3003/sqrt(3);
                            break;
                        }
                        case 3:
                        {
                            hCoeffs->coeffs73->Pre = 1;
                            hCoeffs->coeffs73->Coeff5 = -243*I*dM*sqrt(3/2)*(1-4*eta+3*eta2)/160160;
                            break;
                        }
                        case 4:
                        {
                            hCoeffs->coeffs74->Pre = 1;
                            hCoeffs->coeffs74->Coeff6 = -128*sqrt(2/33)*(1-7*eta+14*eta2-7*eta3)/1365;
                            break;
                        }
                        case 5:
                        {
                            hCoeffs->coeffs75->Pre = 1;
                            hCoeffs->coeffs75->Coeff5 = 15625*I*dM*(1-4*eta+3*eta2)/26208/sqrt(66);
                            break;
                        }
                        case 6:
                        {
                            hCoeffs->coeffs76->Pre = 1;
                            hCoeffs->coeffs76->Coeff6 = 81*(1-7*eta+14*eta2-7*eta3)*sqrt(3/143)/35;
                            break;
                        }
                        case 7:
                        {
                            hCoeffs->coeffs77->Pre = 1;
                            hCoeffs->coeffs77->Coeff5 = -16807*I*dM*sqrt(7/858)*(1-4*eta+3*eta2)/1440;
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                case 8:
                {
                    switch(m)
                    {
                        case 0:
                        {
                            break;
                        }
                        case 1:
                        {
                            break;
                        }
                        case 2:
                        {
                            hCoeffs->coeffs82->Pre = 1;
                            hCoeffs->coeffs82->Coeff6 = (1-7*eta+14*eta2-7*eta3)/9009/sqrt(85);
                            break;
                        }
                        case 3:
                        {
                            break;
                        }
                        case 4:
                        {
                            hCoeffs->coeffs84->Pre = 1;
                            hCoeffs->coeffs84->Coeff6 = -128*sqrt(2/187)*(1-7*eta+14*eta2-7*eta3)/4095;
                            break;
                        }
                        case 5:
                        {
                            break;
                        }
                        case 6:
                        {
                            hCoeffs->coeffs86->Pre = 1;
                            hCoeffs->coeffs86->Coeff6 = 243 * sqrt(3/17017)*(1-7*eta+14*eta2-7*eta3)/35;
                            break;
                        }
                        case 7:
                        {
                            break;
                        }
                        case 8:
                        {
                            hCoeffs->coeffs88->Pre = 1;
                            hCoeffs->coeffs88->Coeff6 = -16384*sqrt(2/85085)*(1-7*eta+14*eta2-7*eta3)/63;
                            break;
                        }
                        default:
                            return CEV_FAILURE;
                    }
                    break;
                }
                default:
                    return CEV_FAILURE;
            }
        }
    }
    return CEV_SUCCESS;
}

/* QuasiCircular */
INT Waveform_CQG_25_165003(REAL8 v, /*(omega)^1/3*/
                            REAL8 phi,
                            REAL8 eta,
                            INT l,
                            INT m,
                            CommhCoeffs *hCoeffs,
                            COMPLEX16 *out)
{
    INT failed = 0, absm = abs(m);
    REAL8 v1, v2, v3, v4, v5, v6;
    v1 = v;
    v2 = v1*v1;
    v3 = v2*v1;
    v4 = v3*v1;
    v5 = v4*v1;
    v6 = v5*v1;
    COMPLEX16 Pre, hLM, extra=0;
    PNwaveformCoeffs *coeffs;
    COMPLEX16 correct;
    Pre = -2*eta*v2*sqrt(16*CST_PI/5)*cexp(-I*m*phi);
//print_debug("hPre = %.2e + i%.2e\n", Pre);
    switch(l)
    {
        case 2:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs20;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs21;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs22;
                    extra = -428*log(16*v2)*v6/105;
                    break;
                default:
                    return CEV_FAILURE;
            }
            break;
        }
        case 3:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs30;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs31;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs32;
                    break;
                case 3:
                    coeffs = hCoeffs->coeffs33;
                    break;
                default:
                    return CEV_FAILURE;
            }
            break;
        }
        case 4:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs40;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs41;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs42;
                    break;
                case 3:
                    coeffs = hCoeffs->coeffs43;
                    break;
                case 4:
                    coeffs = hCoeffs->coeffs44;
                    break;
                default:
                    return CEV_FAILURE;
            }
        }
        case 5:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs50;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs51;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs52;
                    break;
                case 3:
                    coeffs = hCoeffs->coeffs53;
                    break;
                case 4:
                    coeffs = hCoeffs->coeffs54;
                    break;
                case 5:
                    coeffs = hCoeffs->coeffs55;
                    break;
                default:
                    return CEV_FAILURE;
            }
        }
        case 6:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs60;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs61;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs62;
                    break;
                case 3:
                    coeffs = hCoeffs->coeffs63;
                    break;
                case 4:
                    coeffs = hCoeffs->coeffs64;
                    break;
                case 5:
                    coeffs = hCoeffs->coeffs65;
                    break;
                case 6:
                    coeffs = hCoeffs->coeffs66;
                    break;
                default:
                    return CEV_FAILURE;
            }
        }
        case 7:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs70;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs71;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs72;
                    break;
                case 3:
                    coeffs = hCoeffs->coeffs73;
                    break;
                case 4:
                    coeffs = hCoeffs->coeffs74;
                    break;
                case 5:
                    coeffs = hCoeffs->coeffs75;
                    break;
                case 6:
                    coeffs = hCoeffs->coeffs76;
                    break;
                case 7:
                    coeffs = hCoeffs->coeffs77;
                    break;
                default:
                    return CEV_FAILURE;
            }
        }
        case 8:
        {
            switch(absm)
            {
                case 0:
                    coeffs = hCoeffs->coeffs80;
                    break;
                case 1:
                    coeffs = hCoeffs->coeffs81;
                    break;
                case 2:
                    coeffs = hCoeffs->coeffs82;
                    break;
                case 3:
                    coeffs = hCoeffs->coeffs83;
                    break;
                case 4:
                    coeffs = hCoeffs->coeffs84;
                    break;
                case 5:
                    coeffs = hCoeffs->coeffs85;
                    break;
                case 6:
                    coeffs = hCoeffs->coeffs86;
                    break;
                case 7:
                    coeffs = hCoeffs->coeffs87;
                    break;
                case 8:
                    coeffs = hCoeffs->coeffs88;
                    break;
                default:
                    return CEV_FAILURE;
            }
        }
        default:
            return CEV_FAILURE;
    }
    correct = coeffs->Pre*( coeffs->Coeff0 +
                            v1*coeffs->Coeff1 +
                            v2*coeffs->Coeff2 +
                            v3*coeffs->Coeff3 +
                            v4*coeffs->Coeff4 +
                            v5*coeffs->Coeff5 +
                            v6*coeffs->Coeff6 ) + extra;
//print_debug("hNewt = %.2e + i%.2e, correct = %.2e + i%.2e\n", Pre, correct);
    hLM = Pre*correct;
    *out = hLM;
    return CEV_SUCCESS;
}


/*
*
*       Eccentric PN Waveform
*
*/
INT Calculate3PNEccentricWaveformCoefficients(REAL8 eta,
                                              PNEccentricWaveformCoeffs *eCoeffs)
{
    REAL8 eta2, eta3, PIsq, LN32;
    LN32 = log(3/2);
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    PIsq = CST_PI*CST_PI;

    eCoeffs->hCoeffv00 = 1;
    eCoeffs->hCoeffv0em = 1/4;
    eCoeffs->hCoeffv0ep = 5/4;

    eCoeffs->hCoeffv20 = -107/42 + 55*eta/42;
    eCoeffs->hCoeffv2em = -257/168 + 169*eta/168;
    eCoeffs->hCoeffv2ep = -31/24 + 35*eta/24;

    eCoeffs->hCoeffv30 = 2*CST_PI;
    eCoeffs->hCoeffv3em = 11*CST_PI/4 + 27*I*LN32;
    eCoeffs->hCoeffv3ep = 13*CST_PI + 3*I*CST_LN2;

    eCoeffs->hCoeffv40 = -2173/1512 -1069*eta/216 +2047*eta2/1512;
    eCoeffs->hCoeffv4em = -4271/756 - 35131*eta/6048 + 421*eta2/864;
    eCoeffs->hCoeffv4ep = -2155/252-1655*eta/672 + 371*eta2/288;

    eCoeffs->hCoeffv50 = -107*CST_PI/21 - 24*I*eta + 34*CST_PI*eta/21;
    eCoeffs->hCoeffv4em = -27*I/2-1081*CST_PI/168+
        (-1013*I/140+137*CST_PI/42)*eta +(27*I/4+9*I*eta)*LN32;
    eCoeffs->hCoeffv4ep = -9*I/2 + 229*CST_PI/168 + 
        (-436571*I/420 + 61*CST_PI/42)*eta + (473*I/28-3*I*eta/7)*CST_LN2;

    eCoeffs->hCoeffv60 = 27027409/646800 - 856*CST_GAMMA/105 + 428*I*CST_PI/105 + 2*PIsq/3+
        (-278185/33264+41*PIsq/96)*eta - 20261*eta2/2772 + 114635*eta3/99792; // - 428*ln(16x)/105
    eCoeffs->hCoeffv4em = 219775769/1663200 + 749*I*CST_PI/60 + 49*PIsq/24 - 749*CST_GAMMA/30 +
        (-121717/20790-41*PIsq/192)*eta-86531*eta2/8316 - 33331*eta3/399168 +
        (-2889/70 + 81*I*CST_PI/2)*LN32 - 81*LN32*LN32/2; // - 749*ln(16x)/60
    eCoeffs->hCoeffv4ep = 55608313/1058400 + 3103*I*CST_PI/420 + 29*PIsq/24 - 3103*CST_GAMMA/210 +
        (-199855/3024 + 41*PIsq/48)*eta - 9967*eta2/1008 + 35579*eta3/36288 +
        (-6527/210 + 3*I*CST_PI/2)*CST_LN2 + 3*CST_LN2*CST_LN2/2; //-3103*ln(x)/420

    return CEV_SUCCESS;
}


INT CalculateInspiral3PNEccentricWaveform22mode(REAL8 v,
                                                REAL8 e,
                                                REAL8Vector *eccVec,
                                                REAL8 xi,
                                                PNEccentricWaveformCoeffs *eCoeffs,
                                                COMPLEX16 *out)
{
    COMPLEX16 hLM, PN0, PN1, PN1_5, PN2, PN2_5, PN3;
    COMPLEX16 emxi, epxi;
    COMPLEX16 extra;
    REAL8 v2, v3, v4, v5, v6;
    v2 = v*v;
    v3 = v2*v;
    v4 = v3*v;
    v5 = v4*v;
    emxi = cexp(-I*xi);
    epxi = cexp(I*xi);
    PN0 = eCoeffs->hCoeffv00 + e*(eCoeffs->hCoeffv0em * emxi + eCoeffs->hCoeffv0ep * epxi);
    PN1 = v2*(eCoeffs->hCoeffv20 + e*(eCoeffs->hCoeffv2em * emxi + eCoeffs->hCoeffv2ep * epxi));
    PN1_5 = v3*(eCoeffs->hCoeffv30 + e*(eCoeffs->hCoeffv3em * emxi + eCoeffs->hCoeffv3ep * epxi));
    PN2 = v4*(eCoeffs->hCoeffv40 + e*(eCoeffs->hCoeffv4em * emxi + eCoeffs->hCoeffv4ep * epxi));
    PN2_5 = v5*(eCoeffs->hCoeffv50 + e*(eCoeffs->hCoeffv5em * emxi + eCoeffs->hCoeffv5ep * epxi));
    PN3 = v6*(eCoeffs->hCoeffv60 + e*(eCoeffs->hCoeffv6em * emxi + eCoeffs->hCoeffv6ep * epxi));
    extra = v6*(-428*log(16*v2)/105 - e*(-749*log(16*v2)*emxi/60 -3103*log(v2)*epxi/420 ));
    hLM = PN0 + PN1 + PN1_5 + PN2 + PN2_5 + PN3 + extra;
    *out = hLM;
    return CEV_SUCCESS;
}


