/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "eobEccCorrection.h"
#include "dyHcapNumericalDerivative.h"

int XLALSpinAlignedHcapDerivative(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );


static int PNwaveformPRD544813rdotc_22mode(REAL8 *hr,
                                           REAL8 *hi,
                                           const REAL8 x1,
                                           const REAL8 x2,
                                           const REAL8 x3, // test particle position
                                           const REAL8 v1,
                                           const REAL8 v2,
                                           const REAL8 v3, // test particle velocity
                                           const REAL8 eta) // symmetric mass ratio
{
    const REAL8 dm = GET_SQRT(1-4*eta);
    REAL8 r = GET_SQRT(x1*x1+x2*x2+x3*x3);
    REAL8 n1 = x1/r,n2 = x2/r,n3 = x3/r;
    REAL8 rdot = v1*n1+v2*n2+v3*n3;
    REAL8 vsqr = v1*v1+v2*v2+v3*v3;
    
    REAL8 vn1,vn2,vn3;
    vn1 = rdot*n1;
    vn2 = rdot*n2;
    vn3 = rdot*n3;
    REAL8 lambda1,lambda2,lambda3;
    lambda1 = v1-vn1;
    lambda2 = v2-vn2;
    lambda3 = v3-vn3;
    REAL8 ln = GET_SQRT(lambda1*lambda1+lambda2*lambda2+lambda3*lambda3);
    lambda1 = lambda1/ln;
    lambda2 = lambda2/ln;
    lambda3 = lambda3/ln;
    
    COMPLEX16 hm22;
    
    COMPLEX16 h11,h12,h13,h22,h23,h33;
    COMPLEX16 Q11,Q12,Q13,Q22,Q23,Q33;
    
    ///////////////////////////////////////00 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)/3;
    Q12 = -GET_SQRT(CST_PI/5)/3*I;
    Q13 = 0;
    Q22 = -GET_SQRT(CST_PI/5)/3;
    Q23 = 8/GET_SQRT(5*CST_PI)/3;
    Q33 = 0;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    h11 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v1+v1*n1)+(3*(1-3*eta)*rdot*rdot)/r*n1*n1);
    h12 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v2+v1*n2)+(3*(1-3*eta)*rdot*rdot)/r*n1*n2);
    h13 += 1.0/3*(2/r*rdot*(5+3*eta)*(n1*v3+v1*n3)+(3*(1-3*eta)*rdot*rdot)/r*n1*n3);
    h22 += 1.0/3*(2/r*rdot*(5+3*eta)*(n2*v2+v2*n2)+(3*(1-3*eta)*rdot*rdot)/r*n2*n2);
    h23 += 1.0/3*(2/r*rdot*(5+3*eta)*(n2*v3+v2*n3)+(3*(1-3*eta)*rdot*rdot)/r*n2*n3);
    h33 += 1.0/3*(2/r*rdot*(5+3*eta)*(n3*v3+v3*n3)+(3*(1-3*eta)*rdot*rdot)/r*n3*n3);
    
    hm22 = Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////10 part///////////////////////////////////////
    Q11 = 0;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(x1 - I*x2)/(3.*r);
    Q22 = 0;
    Q23 = GET_SQRT(CST_PI/5)*(I*x1 + x2)/(3.*r);
    Q33 = 0;
    
    h11 = dm*3/r*(-rdot*n1*n1);
    h12 = dm*3/r*(-rdot*n1*n2);
    h13 = dm*3/r*(-rdot*n1*n2);
    h22 = dm*3/r*(-rdot*n2*n2);
    h23 = dm*3/r*(-rdot*n2*n3);
    h33 = dm*3/r*(-rdot*n3*n3);
    
    h11 += dm/12/r*((n1*v1+v1*n1)*(rdot*rdot*(63+54*eta))
                    +n1*n1*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v1*(186+24*eta));
    h12 += dm/12/r*((n1*v2+v1*n2)*(rdot*rdot*(63+54*eta))
                    +n1*n2*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v2*(186+24*eta));
    h13 += dm/12/r*((n1*v3+v1*n3)*(rdot*rdot*(63+54*eta))
                    +n1*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v1*v3*(186+24*eta));
    h22 += dm/12/r*((n2*v2+v2*n2)*(rdot*rdot*(63+54*eta))
                    +n2*n2*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v2*v2*(186+24*eta));
    h23 += dm/12/r*((n2*v3+v2*n3)*(rdot*rdot*(63+54*eta))
                    +n2*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v2*v3*(186+24*eta));
    h33 += dm/12/r*((n3*v3+v3*n3)*(rdot*rdot*(63+54*eta))
                    +n3*n3*rdot*(rdot*rdot*(15-90*eta)-vsqr*(63-54*eta)+(242-24*eta)/r)-rdot*v3*v3*(186+24*eta));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////01 part///////////////////////////////////////
    Q11 = 0;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(v1 - I*v2)/3.;
    Q22 = 0;
    Q23 = GET_SQRT(CST_PI/5)*(I*v1 + v2)/3.;
    Q33 = 0;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    h11 += dm*(-(n1*v1+v1*n1)/2/r*rdot*(7+4*eta)
               -n1*n1/r*(0.75*(1-2*eta)*rdot*rdot));
    h12 += dm*(-(n1*v2+v1*n2)/2/r*rdot*(7+4*eta)
               -n1*n2/r*(0.75*(1-2*eta)*rdot*rdot));
    h13 += dm*(-(n1*v3+v1*n3)/2/r*rdot*(7+4*eta)
               -n1*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    h22 += dm*(-(n2*v2+v2*n2)/2/r*rdot*(7+4*eta)
               -n2*n2/r*(0.75*(1-2*eta)*rdot*rdot));
    h23 += dm*(-(n2*v3+v2*n3)/2/r*rdot*(7+4*eta)
               -n2*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    h33 += dm*(-(n3*v3+v3*n3)/2/r*rdot*(7+4*eta)
               -n3*n3/r*(0.75*(1-2*eta)*rdot*rdot));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////20 part///////////////////////////////////////
    Q11 = -GET_SQRT(CST_PI/5)*(x1*x1-8.*I*x1*x2-7*x2*x2-x3*x3)/21./r/r;
    Q12 = -I*GET_SQRT(CST_PI/5)*(3*(x1*x1+x2*x2)+x3*x3)/21./r/r;
    Q13 = 2*GET_SQRT(CST_PI/5)*(x1-I*x2)*x3/21./r/r;
    Q22 = GET_SQRT(CST_PI/5)*(-7*x1*x1+8.*I*x1*x2+x2*x2-x3*x3)/21./r/r;
    Q23 = -2.*I*GET_SQRT(CST_PI/5)*(x1-I*x2)*x3/21./r/r;
    Q33 = 8*GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)/21./r/r;
    
    h11 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n1+15*rdot*(n1*v1+v1*n1));
    h12 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n2+15*rdot*(n1*v2+v1*n2));
    h13 = (1-3*eta)/3/r*((-15*rdot*rdot)*n1*n3+15*rdot*(n1*v3+v1*n3));
    h22 = (1-3*eta)/3/r*((-15*rdot*rdot)*n2*n2+15*rdot*(n2*v2+v2*n2));
    h23 = (1-3*eta)/3/r*((-15*rdot*rdot)*n2*n3+15*rdot*(n2*v3+v2*n3));
    h33 = (1-3*eta)/3/r*((-15*rdot*rdot)*n3*n3+15*rdot*(n3*v3+v3*n3));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////11 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(-v1*x1+4.*I*v2*x1+4.*I*v1*x2+7*v2*x2+v3*x3)/21./r;
    Q12 = -I*GET_SQRT(CST_PI/5)*(3*v1*x1+3*v2*x2+v3*x3)/21./r;
    Q13 = GET_SQRT(CST_PI/5)*(v3*(x1-I*x2)+x3*(v1-I*v2))/21./r;
    Q22 = GET_SQRT(CST_PI/5)*(-7*v1*x1+4.*I*v2*x1+4.*I*v1*x2+v2*x2-v3*x3)/21./r;
    Q23 = -GET_SQRT(CST_PI/5)*(v3*(I*x1+x2)+x3*(I*v1+v2))/21./r;
    Q33 = 8*GET_SQRT(CST_PI/5)*(v1-I*v2)*(x1-I*x2)/21./r;
    
    h11 = (1-3*eta)/3/r*(12*rdot*n1*n1);
    h12 = (1-3*eta)/3/r*(12*rdot*n1*n2);
    h13 = (1-3*eta)/3/r*(12*rdot*n1*n3);
    h22 = (1-3*eta)/3/r*(12*rdot*n2*n2);
    h23 = (1-3*eta)/3/r*(12*rdot*n2*n3);
    h33 = (1-3*eta)/3/r*(12*rdot*n3*n3);
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////02 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(-v1*v1+8.*I*v1*v2+7*v2*v2+v3*v3)/21.;
    Q12 = -I*GET_SQRT(CST_PI/5)*(3*(v1*v1+v2*v2)+v3*v3)/21.;
    Q13 = 2*GET_SQRT(CST_PI/5)*(v1-I*v2)*v3/21.;
    Q22 = GET_SQRT(CST_PI/5)*(-7*v1*v1+8.*I*v1*v2+v2*v2-v3*v3)/21.;
    Q23 = -2.*I*GET_SQRT(CST_PI/5)*(v1-I*v2)*v3/21.;
    Q33 = 8*GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)/21.;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////30 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)*x3/7./r/r/r;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(x1-I*x2)*(2.*x1*x1+5.*I*x1*x2+7*x2*x2+3*x3*x3)/21./r/r/r;
    Q22 = GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)*x3/7./r/r/r;
    Q23 = GET_SQRT(CST_PI/5)*(I*x1+x2)*(7*x1*x1-5.*I*x1*x2+2*x2*x2+3*x3*x3)/21./r/r/r;
    Q33 = -2*GET_SQRT(CST_PI/5)*C_POW(x1-I*x2,2)*x3/7./r/r/r;
    
    h11 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n1-8.5*rdot*v1*v1-(-105*rdot*rdot)/12*(n1*v1+v1*n1));
    h12 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n2-8.5*rdot*v1*v2-(-105*rdot*rdot)/12*(n1*v2+v1*n2));
    h13 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n1*n3-8.5*rdot*v1*v3-(-105*rdot*rdot)/12*(n1*v3+v1*n3));
    h22 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n2*n2-8.5*rdot*v2*v2-(-105*rdot*rdot)/12*(n2*v2+v2*n2));
    h23 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n2*n3-8.5*rdot*v2*v3-(-105*rdot*rdot)/12*(n2*v3+v2*n3));
    h33 = dm*(1-2*eta)/r*(1.25*(3*vsqr-7*rdot*rdot+6/r)*rdot*n3*n3-8.5*rdot*v3*v3-(-105*rdot*rdot)/12*(n3*v3+v3*n3));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////21 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*((x1-I*x2)*(2*v1*x1+I*x1*v2+4.*I*v1*x2+7*x2*v2)+2*v3*(x1-I*x2)*x3+(v1-I*v2)*x3*x3)/21./r/r;
    Q22 = GET_SQRT(CST_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    Q23 = GET_SQRT(CST_PI/5)*((x1-I*x2)*(v2*(4*x1+2.*I*x2)+v1*(7.*I*x1+x2))+2*v3*(I*x1+x2)*x3+(I*v1+v2)*x3*x3)/21./r/r;
    Q33 = -2*GET_SQRT(CST_PI/5)*(x1-I*x2)*(v3*(x1-I*x2)+2.*(v1-I*v2)*x3)/21./r/r;
    
    h11 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n1-54*rdot*(n1*v1+v1*n1));
    h12 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n2-54*rdot*(n1*v2+v1*n2));
    h13 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n1*n3-54*rdot*(n1*v3+v1*n3));
    h22 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n2*n2-54*rdot*(n2*v2+v2*n2));
    h23 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n2*n3-54*rdot*(n2*v3+v2*n3));
    h33 = dm*(1-2*eta)/4/r*((45*rdot*rdot)*n3*n3-54*rdot*(n3*v3+v3*n3));
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////12 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(v3*v3*(x1-I*x2)+v1*v1*(2*x1+I*x2)+v2*v2*(4*x1-7.*I*x2)-2.*I*v2*v3*x3+2*v1*(I*v2*x1+4*v2*x2+v3*x3))/21./r;
    Q22 = GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    Q23 = GET_SQRT(CST_PI/5)*(v3*v3*(I*x1+x2)+v2*v2*(-I*x1+2*x2)+v1*v1*(7.*I*x1+4*x2)+2*v2*v3*x3+v1*(8*v2*x1-2.*I*v2*x2+2.*I*v3*x3))/21./r;
    Q33 = -2*GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v3*(x1-I*x2)+(v1-I*v2)*x3)/21./r;
    
    h11 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n1);
    h12 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n2);
    h13 = dm*(1-2*eta)*1.5/r*(-3*rdot*n1*n3);
    h22 = dm*(1-2*eta)*1.5/r*(-3*rdot*n2*n2);
    h23 = dm*(1-2*eta)*1.5/r*(-3*rdot*n2*n3);
    h33 = dm*(1-2*eta)*1.5/r*(-3*rdot*n3*n3);
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    
    ///////////////////////////////////////03 part///////////////////////////////////////
    Q11 = GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)*v3/7.;
    Q12 = 0;
    Q13 = -GET_SQRT(CST_PI/5)*(v1-I*v2)*(2*v1*v1+5.*I*v1*v2+7*v2*v2+3*v3*v3)/21.;
    Q22 = GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)*v3/7.;
    Q23 = GET_SQRT(CST_PI/5)*(I*v1+v2)*(7*v1*v1-5.*I*v1*v2+2*v2*v2+3*v3*v3)/21.;
    Q33 = -2*GET_SQRT(CST_PI/5)*C_POW(v1-I*v2,2)*v3/7.;
    
    h11 = 0;
    h12 = 0;
    h13 = 0;
    h22 = 0;
    h23 = 0;
    h33 = 0;
    
    hm22 += Q11*h11+Q22*h22+Q33*h33+2.*(Q12*h12+Q13*h13+Q23*h23);
    //-----------------------------------complete-----------------------------------------
    
    hm22 = hm22*2.*eta;
    *hr = C_REAL(hm22);
    *hi = C_IMAG(hm22);
    if (isnan(*hr) || isnan(*hi))
    {
        print_warning("Error _%s: hECC is nan.\n",__func__);
        //print_err("\thm22 = %f + i%f", C_REAL(hm22), C_IMAG(hm22));
        //print_err("\tx1 = %f, x2 = %f, x3 = %f\n", x1, x2, x3);
        //print_err("\tv1 = %f, v2 = %f, v3 = %f\n", v1, v2, v3);
        //print_err("")
        return CEV_FAILURE;
    }

    return CEV_SUCCESS;
}

INT EOBEccFactorizedWaveformCorrection(COMPLEX16 *hECC,
                                       const REAL8 rphivalues[],
                                       const REAL8 rdot,
                                       const REAL8 phidot,
                                       const REAL8 eta)
{
    REAL8 r,phi,prt,pphi;
    r = rphivalues[0];
    prt = rphivalues[2];
    phi = rphivalues[1];
    pphi = rphivalues[3];

    REAL8 x1,x2,x3;
    x1 = r*GET_COS(phi);
    x2 = r*GET_SIN(phi);
    x3 = 0;
    
    REAL8 v1,v2,v3;
    v1 = rdot*GET_COS(phi)-r*phidot*GET_SIN(phi);
    v2 = rdot*GET_SIN(phi)+r*phidot*GET_COS(phi);
    v3 = 0;
    
    REAL8 hr, hi;
    if(PNwaveformPRD544813rdotc_22mode(&hr, &hi,x1,x2,x3,v1,v2,v3,eta) == CEV_FAILURE)
    {
        print_warning("Error _%s: hECC is nan.\n", __func__);
        //print_err("\tphidot = %f, pphi = %f, rdot = %f\n", phidot, pphi, rdot);
        //print_err("\tv1 = %f, v2 = %f, v3 = %f\n", v1, v2, v3);
        return CEV_FAILURE;
    }
    //print_debug("hECC = %e + i%e\n", hr, hi);
    *hECC = hr + I*hi;
    return CEV_SUCCESS;
}

INT EOBCalculateEccCorrection(COMPLEX16 *hECC,
                              const REAL8Vector *values,
                              SpinEOBParams *seobParams)
{
    REAL8 r,phi,pr,pphi,rdot,phidot;
    
    r = values->data[0];
    phi = values->data[1];
    pr = values->data[2];
    pphi = values->data[3];
    
    REAL8 dydt[4],ttmp=0;
    if (XLALSpinAlignedHcapDerivative(ttmp, values->data, dydt, seobParams) == GSL_FAILURE)
        return CEV_FAILURE;
    
    rdot = dydt[0];
    phidot = dydt[1];
    
    REAL8 x1,x2,x3;
    x1 = r*GET_COS(phi);
    x2 = r*GET_SIN(phi);
    x3 = 0;
    
    REAL8 v1,v2,v3;
    v1 = rdot*GET_COS(phi)-r*phidot*GET_SIN(phi);
    v2 = rdot*GET_SIN(phi)+r*phidot*GET_COS(phi);
    v3 = 0;
    
    REAL8 hr, hi;
    if (PNwaveformPRD544813rdotc_22mode(&hr,&hi,x1,x2,x3,v1,v2,v3,seobParams->eobParams->eta) == CEV_FAILURE)
    {
        print_warning("Error _%s: hECC is nan.\n", __func__);
        //print_err("\tphidot = %f, pphi = %f, rdot = %f\n", phidot, pphi, rdot);
        return CEV_FAILURE;
    }
    *hECC = hr + I*hi;
    return CEV_SUCCESS;
}

