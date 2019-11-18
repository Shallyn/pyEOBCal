/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyHamiltonian.h"
#include <gsl/gsl_deriv.h>

#define LMAX 8
#define STEP_SIZE 1.0e-4


/* Wrapper for GSL to call the Hamiltonian function */
double
GSLSpinAlignedHamiltonianWrapper (double x, void *params)
{
  HcapDerivParams *dParams = (HcapDerivParams *) params;

  EOBParams *eobParams = dParams->params->eobParams;

  REAL8 tmpVec[6];

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
  REAL8Vector *s1Vec = dParams->params->s1VecOverMtMt;
  REAL8Vector *s2Vec = dParams->params->s2VecOverMtMt;
  REAL8Vector *sigmaKerr = dParams->params->sigmaKerr;
  REAL8Vector *sigmaStar = dParams->params->sigmaStar;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy (tmpVec, dParams->values, sizeof (tmpVec));

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = 3;
  r.data = tmpVec;
  p.data = tmpVec + 3;

  return SpinEOBHamiltonian(eobParams->eta, &r, &p, s1Vec, s2Vec,
				       sigmaKerr, sigmaStar, dParams->params->tortoise,
				       dParams->params->seobCoeffs) /
    eobParams->eta;
}



REAL8 SpinEOBHamiltonian(const REAL8 eta,
                         REAL8Vector *x,
                         REAL8Vector *p,
                         REAL8Vector *s1Vec,
                         REAL8Vector *s2Vec,
                         REAL8Vector *sigmaKerr,
                         REAL8Vector *sigmaStar,
                         int tortoise,
                         SpinEOBHCoeffs * coeffs)
{
    REAL8 r, r2, nx, ny, nz;
    REAL8 sKerr_x, sKerr_y, sKerr_z, a, a2;
    REAL8 sStar_x, sStar_y, sStar_z;
    REAL8 e3_x, e3_y, e3_z;
    REAL8 costheta;		/* Cosine of angle between Skerr and r */
    REAL8 xi2, xi_x, xi_y, xi_z;	/* Cross product of unit vectors in direction of Skerr and r */
    REAL8 vx, vy, vz, pxir, pvr, pn, prT, pr, pf, ptheta2;	/*prT is the tortoise pr */
    REAL8 w2, rho2;
    REAL8 u, u2, u3, u4, u5;
    REAL8 bulk, deltaT, deltaR, Lambda;
    REAL8 D, qq, ww, B, w, MU, nu, BR, wr, nur, mur;
    REAL8 wcos, nucos, mucos, ww_r, Lambda_r;
    REAL8 logTerms, deltaU, deltaU_u, Q, deltaT_r, pn2, pp;
    REAL8 deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z;
    REAL8 sx, sy, sz, sxi, sv, sn, s3;
    REAL8 H, Hns, Hs, Hss, Hreal, Hwcos, Hwr, HSOL, HSONL;
    REAL8 m1PlusetaKK;

    /* Terms which come into the 3.5PN mapping of the spins */
    //REAL8 aaa, bbb, a13P5, a23P5, a33P5, b13P5, b23P5, b33P5;
    REAL8 sMultiplier1, sMultiplier2;

    /*Temporary p vector which we will make non-tortoise */
    REAL8 tmpP[3];

    REAL8 csi;

    /* Spin gauge parameters. (YP) simplified, since both are zero. */
    // static const double aa=0., bb=0.;

    //printf( "In Hamiltonian:\n" );
    //printf( "x = %.16e\t%.16e\t%.16e\n", x->data[0], x->data[1], x->data[2] );
    //printf( "p = %.16e\t%.16e\t%.16e\n", p->data[0], p->data[1], p->data[2] );

    r2 =
        x->data[0] * x->data[0] + x->data[1] * x->data[1] +
        x->data[2] * x->data[2];
    r = sqrt (r2);
    nx = x->data[0] / r;
    ny = x->data[1] / r;
    nz = x->data[2] / r;

    sKerr_x = sigmaKerr->data[0];
    sKerr_y = sigmaKerr->data[1];
    sKerr_z = sigmaKerr->data[2];

    sStar_x = sigmaStar->data[0];
    sStar_y = sigmaStar->data[1];
    sStar_z = sigmaStar->data[2];

    a2 = sKerr_x * sKerr_x + sKerr_y * sKerr_y + sKerr_z * sKerr_z;
    a = sqrt (a2);

    if (a != 0.)
    {
        e3_x = sKerr_x / a;
        e3_y = sKerr_y / a;
        e3_z = sKerr_z / a;
    }
    else
    {
        e3_x = 0.;
        e3_y = 0.;
        e3_z = 1.;
    }

    costheta = e3_x * nx + e3_y * ny + e3_z * nz;

    xi2 = 1. - costheta * costheta;

    xi_x = -e3_z * ny + e3_y * nz;
    xi_y = e3_z * nx - e3_x * nz;
    xi_z = -e3_y * nx + e3_x * ny;

    vx = -nz * xi_y + ny * xi_z;
    vy = nz * xi_x - nx * xi_z;
    vz = -ny * xi_x + nx * xi_y;

    w2 = r2 + a2;
    rho2 = r2 + a2 * costheta * costheta;

    u = 1. / r;
    u2 = u * u;
    u3 = u2 * u;
    u4 = u2 * u2;
    u5 = u4 * u;
    //printf( "KK = %.16e\n", coeffs->KK );
    m1PlusetaKK = -1. + eta * coeffs->KK;
    /* Eq. 5.75 of BB1 */
    bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
    /* Eq. 5.73 of BB1 */
    logTerms =
        1. + eta * coeffs->k0 + eta * log (1. + coeffs->k1 * u + coeffs->k2 * u2 +
                        coeffs->k3 * u3 + coeffs->k4 * u4 +
                        coeffs->k5 * u5 +
                        coeffs->k5l * u5 * log (u));
    //printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
    /* Eq. 5.73 of BB1 */
    deltaU = bulk * logTerms;
    /* Eq. 5.71 of BB1 */
    deltaT = r2 * deltaU;
    /* ddeltaU/du */
    deltaU_u = 2. * (1. / m1PlusetaKK + a2 * u) * logTerms +
        bulk * (eta *
            (coeffs->k1 +
            u * (2. * coeffs->k2 +
            u * (3. * coeffs->k3 +
                u * (4. * coeffs->k4 +
                    5. * (coeffs->k5 +
                    coeffs->k5l * log (u)) * u))))) / (1. +
                                        coeffs->
                                        k1 * u +
                                        coeffs->
                                        k2 * u2 +
                                        coeffs->
                                        k3 * u3 +
                                        coeffs->
                                        k4 * u4 +
                                        (coeffs->
                                        k5 +
                                        coeffs->
                                        k5l *
                                        log (u))
                                        * u5);
    /* ddeltaT/dr */
    deltaT_r = 2. * r * deltaU - deltaU_u;
    /* Eq. 5.39 of BB1 */
    Lambda = w2 * w2 - a2 * deltaT * xi2;
    /* Eq. 5.83 of BB1, inverse */
    D = 1. + log (1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
    /* Eq. 5.38 of BB1 */
    deltaR = deltaT * D;
    /* See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
    qq = 2. * eta * (4. - 3. * eta);
    /* See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
    ww = 2. * a * r + coeffs->b3 * eta * a2 * a * u + coeffs->bb3 * eta * a * u;

    /* We need to transform the momentum to get the tortoise co-ord */
    if (tortoise)
    {
        csi = sqrt (deltaT * deltaR) / w2;	/* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    }
    else
    {
        csi = 1.0;
    }
    //printf( "csi(miami) = %.16e\n", csi );
    prT = p->data[0] * nx + p->data[1] * ny + p->data[2] * nz;
    /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    tmpP[0] = p->data[0] - nx * prT * (csi - 1.) / csi;
    tmpP[1] = p->data[1] - ny * prT * (csi - 1.) / csi;
    tmpP[2] = p->data[2] - nz * prT * (csi - 1.) / csi;

    pxir = (tmpP[0] * xi_x + tmpP[1] * xi_y + tmpP[2] * xi_z) * r;
    pvr = (tmpP[0] * vx + tmpP[1] * vy + tmpP[2] * vz) * r;
    pn = tmpP[0] * nx + tmpP[1] * ny + tmpP[2] * nz;

    pr = pn;
    pf = pxir;
    ptheta2 = pvr * pvr / xi2;
    //printf( "pr = %.16e, prT = %.16e\n", pr, prT );

    //printf( " a = %.16e, r = %.16e\n", a, r );
    //printf( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n",
    //D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );
    /* Eqs. 5.36 - 5.46 of BB1 */
    /* Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
    Hns =
        sqrt (1. + prT * prT * prT * prT * qq * u2 + ptheta2 / rho2 +
        pf * pf * rho2 / (Lambda * xi2) +
        pr * pr * deltaR / rho2) / sqrt (Lambda / (rho2 * deltaT)) +
        pf * ww / Lambda;

    //printf( "term 1 in Hns: %.16e\n",  prT*prT*prT*prT*qq*u2 );
    //printf( "term 2 in Hns: %.16e\n", ptheta2/rho2 );
    //printf( "term 3 in Hns = %.16e\n", pf*pf*rho2/(Lambda*xi2) );
    //printf( "term 4 in Hns = %.16e\n", pr*pr*deltaR/rho2 );
    //printf( "term 5 in Hns = %.16e\n", Lambda/(rho2*deltaT) );
    //printf( "term 6 in Hns = %.16e\n", pf*ww/Lambda );
    /* Eqs. 5.30 - 5.33 of BB1 */
    B = sqrt (deltaT);
    w = ww / Lambda;
    nu = 0.5 * log (deltaT * rho2 / Lambda);
    MU = 0.5 * log (rho2);
    /* dLambda/dr */
    Lambda_r = 4. * r * w2 - a2 * deltaT_r * xi2;
    ww_r =
        2. * a - (a2 * a * coeffs->b3 * eta) * u2 - coeffs->bb3 * eta * a * u2;
    /* Eqs. 5.47a - 5.47d of BB1 */
    BR =
        (-2. * deltaT + sqrt (deltaR) * deltaT_r) / (2. * sqrt (deltaR * deltaT));
    wr = (-Lambda_r * ww + Lambda * ww_r) / (Lambda * Lambda);
    nur =
        (r / rho2 +
        (w2 * (-4. * r * deltaT + w2 * deltaT_r)) / (2. * deltaT * Lambda));
    mur = (r / rho2 - 1. / sqrt (deltaR));
    /* Eqs. 5.47f - 5.47h of BB1 */
    wcos = -2. * a2 * costheta * deltaT * ww / (Lambda * Lambda);
    nucos = a2 * costheta * w2 * (w2 - deltaT) / (rho2 * Lambda);
    mucos = a2 * costheta / rho2;
    /* Eq. 5.52 of BB1, (YP) simplified */
    //Q = 1. + pvr*pvr/(exp(2.*MU)*xi2) + exp(2.*nu)*pxir*pxir/(B*B*xi2) + pn*pn*deltaR/exp(2.*MU);
    Q =
        1. + pvr * pvr / (rho2 * xi2) +
        deltaT * rho2 / Lambda * pxir * pxir / (B * B * xi2) +
        pn * pn * deltaR / rho2;

    pn2 = pr * pr * deltaR / rho2;
    pp = Q - 1.;
    //printf( "pn2 = %.16e, pp = %.16e\n", pn2, pp );
    //printf( "sigmaKerr = %.16e, sigmaStar = %.16e\n", sKerr_z, sStar_z );
    /* Eq. 5.68 of BB1, (YP) simplified for aa=bb=0. */
    /*
        deltaSigmaStar_x=(- 8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_x - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_x +
        eta*(-8.*sKerr_x - 36.*pn2*r*sKerr_x + 3.*pp*r*sKerr_x + 14.*sStar_x - 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x))/(12.*r);

        deltaSigmaStar_y=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_y - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_y +
        eta*(-8.*sKerr_y - 36.*pn2*r*sKerr_y + 3.*pp*r*sKerr_y + 14.*sStar_y - 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y))/(12.*r);

        deltaSigmaStar_z=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_z - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_z +
        eta*(-8.*sKerr_z - 36.*pn2*r*sKerr_z + 3.*pp*r*sKerr_z + 14.*sStar_z - 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z))/(12.*r);
    */
    deltaSigmaStar_x =
        eta * (-8. * sKerr_x - 36. * pn2 * r * sKerr_x + 3. * pp * r * sKerr_x +
        14. * sStar_x - 30. * pn2 * r * sStar_x +
        4. * pp * r * sStar_x) / (12. * r);

    deltaSigmaStar_y =
        eta * (-8. * sKerr_y - 36. * pn2 * r * sKerr_y + 3. * pp * r * sKerr_y +
        14. * sStar_y - 30. * pn2 * r * sStar_y +
        4. * pp * r * sStar_y) / (12. * r);

    deltaSigmaStar_z =
        eta * (-8. * sKerr_z - 36. * pn2 * r * sKerr_z + 3. * pp * r * sKerr_z +
        14. * sStar_z - 30. * pn2 * r * sStar_z +
        4. * pp * r * sStar_z) / (12. * r);


    /* Now compute the additional 3.5PN terms. */
    /* The following gauge parameters correspond to those given by
    * Eqs. (69) and (70) of BB2 (aaa -> a0, bbb -> b0).
    * In SEOBNRv1 model, we chose to set all of them to zero,
    * described between Eqs. (3) and (4).
    */
    /*
        aaa = -3./2.*eta;
        bbb = -5./4.*eta;
        a1 = eta*eta/2.;
        a2 = -(1./8.)*eta*(-7. + 8.*eta);
        a3 = -((9.*eta*eta)/16.);
        b1 = 1./16.*eta*(9. + 5.*eta);
        b2 = -(1./8.)*eta*(-17. + 5.*eta);
        b3 = -3./8.*eta*eta;
    */
    /*aaa = 0.;
        bbb = 0.;
        a13P5 = 0.;
        a23P5 = 0.;
        a33P5 = 0.;
        b13P5 = 0.;
        b23P5 = 0.;
        b33P5 = 0.;
    */
    /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
    /*
        sMultiplier1 =-(2.*(24.*b23P5 + eta*(-353. + 27.*eta) + bbb*(56. + 60.*eta)) +
        2.*(24.*b13P5 - 24.*b23P5 + bbb*(14. - 66.*eta) + 103.*eta - 60.*eta*eta)*pp*
        r + 120.*(2.*b33P5 - 3.*eta*(bbb + eta))*pn2*pn2*r*r +
        (-48.*b13P5 + 4.*bbb*(1. + 3.*eta) + eta*(23. + 3.*eta))*pp*pp*
        r*r + 6.*pn2*r*(16.*b13P5 + 32.*b23P5 + 24.*b33P5 - 47.*eta +
        54.*eta*eta + 24.*bbb*(1. + eta) +
        (24.*b13P5 - 24.*b33P5 - 16.*eta + 21.*eta*eta + bbb*(-2. + 30.*eta))*pp*
        r))/(72.*r*r);
    */
    sMultiplier1 =
        -(2. * eta * (-353. + 27. * eta) +
        2. * (103. * eta - 60. * eta * eta) * pp * r +
        120. * (-3. * eta * eta) * pn2 * pn2 * r * r +
        (eta * (23. + 3. * eta)) * pp * pp * r * r +
        6. * pn2 * r * (-47. * eta + 54. * eta * eta +
                (-16. * eta +
                21. * eta * eta) * pp * r)) / (72. * r * r);
    /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
    /*
        sMultiplier2 = (-16.*(6.*a23P5 + 7.*eta*(8. + 3.*eta) + aaa*(14. + 15.*eta)) +
        4.*(-24.*a13P5 + 24.*a23P5 - 109.*eta + 51.*eta*eta + 2.*aaa*(-7. + 33.*eta))*
        pp*r + 30.*(-16.*a33P5 + 3.*eta*(8.*aaa + 9.*eta))*pn2*pn2*r*r +
        (96.*a13P5 - 45.*eta - 8.*aaa*(1. + 3.*eta))*pp*pp*r*r -
        6.*pn2*r*(32.*a13P5 + 64.*a23P5 + 48.*a33P5 + 16.*eta + 147.*eta*eta +
        48.*aaa*(1. + eta) + (48.*a13P5 - 48.*a33P5 - 6.*eta + 39.*eta*eta +
        aaa*(-4. + 60.*eta))*pp*r))/(144.*r*r);
    */
    sMultiplier2 =
        (-16. * (7. * eta * (8. + 3. * eta)) +
        4. * (-109. * eta + 51. * eta * eta) * pp * r +
        810. * eta * eta * pn2 * pn2 * r * r - 45. * eta * pp * pp * r * r -
        6. * pn2 * r * (16. * eta + 147. * eta * eta +
                (-6. * eta +
                39. * eta * eta) * pp * r)) / (144. * r * r);
    /* Eq. 52 of BB2 */
    deltaSigmaStar_x +=
        sMultiplier1 * sigmaStar->data[0] + sMultiplier2 * sigmaKerr->data[0];
    deltaSigmaStar_y +=
        sMultiplier1 * sigmaStar->data[1] + sMultiplier2 * sigmaKerr->data[1];
    deltaSigmaStar_z +=
        sMultiplier1 * sigmaStar->data[2] + sMultiplier2 * sigmaKerr->data[2];

    /* And now the (calibrated) 4.5PN term */
    deltaSigmaStar_x += coeffs->d1 * eta * sigmaStar->data[0] / (r * r * r);
    deltaSigmaStar_y += coeffs->d1 * eta * sigmaStar->data[1] / (r * r * r);
    deltaSigmaStar_z += coeffs->d1 * eta * sigmaStar->data[2] / (r * r * r);
    //deltaSigmaStar_x += coeffs->d1v2 * eta * sigmaKerr->data[0] / (r * r * r);
    //deltaSigmaStar_y += coeffs->d1v2 * eta * sigmaKerr->data[1] / (r * r * r);
    //deltaSigmaStar_z += coeffs->d1v2 * eta * sigmaKerr->data[2] / (r * r * r);


    //printf( "deltaSigmaStar_x = %.16e, deltaSigmaStar_y = %.16e, deltaSigmaStar_z = %.16e\n",
    //   deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z );

    sx = sStar_x + deltaSigmaStar_x;
    sy = sStar_y + deltaSigmaStar_y;
    sz = sStar_z + deltaSigmaStar_z;


    sxi = sx * xi_x + sy * xi_y + sz * xi_z;
    sv = sx * vx + sy * vy + sz * vz;
    sn = sx * nx + sy * ny + sz * nz;

    s3 = sx * e3_x + sy * e3_y + sz * e3_z;
    /* Eq. 3.45 of BB1, second term */
    Hwr =
        (exp (-3. * MU - nu) * sqrt (deltaR) *
        (exp (2. * (MU + nu)) * pxir * pxir * sv -
        B * exp (MU + nu) * pvr * pxir * sxi +
        B * B * xi2 * (exp (2. * MU) * (sqrt (Q) + Q) * sv +
                pn * pvr * sn * sqrt (deltaR) -
                pn * pn * sv * deltaR))) / (2. * B * (1. +
                                sqrt (Q)) *
                            sqrt (Q) * xi2);
    /* Eq. 3.45 of BB1, third term */
    Hwcos =
        (exp (-3. * MU - nu) *
        (sn *
        (-(exp (2. * (MU + nu)) * pxir * pxir) +
        B * B * (pvr * pvr - exp (2. * MU) * (sqrt (Q) + Q) * xi2)) -
        B * pn * (B * pvr * sv -
            exp (MU +
                nu) * pxir * sxi) * sqrt (deltaR))) / (2. * B * (1. +
                                        sqrt
                                        (Q)) *
                                    sqrt (Q));
    /* Eq. 3.44 of BB1, leading term */
    HSOL =
        (exp (-MU + 2. * nu) * (-B + exp (MU + nu)) * pxir * s3) / (B * B *
                                    sqrt (Q) *
                                    xi2);
    /* Eq. 3.44 of BB1, next-to-leading term */
    HSONL =
        (exp (-2. * MU + nu) *
        (-(B * exp (MU + nu) * nucos * pxir * (1. + 2. * sqrt (Q)) * sn * xi2) +
        (-(BR * exp (MU + nu) * pxir * (1. + sqrt (Q)) * sv) +
        B * (exp (MU + nu) * nur * pxir * (1. + 2. * sqrt (Q)) * sv +
            B * mur * pvr * sxi + B * sxi * (-(mucos * pn * xi2) +
                            sqrt (Q) * (mur * pvr -
                                nur * pvr + (-mucos +
                                        nucos) *
                                pn * xi2)))) *
        sqrt (deltaR))) / (B * B * (sqrt (Q) + Q) * xi2);
    /* Eq. 3.43 and 3.45 of BB1 */
    Hs = w * s3 + Hwr * wr + Hwcos * wcos + HSOL + HSONL;
    /* Eq. 5.70 of BB1, last term */
    Hss = -0.5 * u3 * (sx * sx + sy * sy + sz * sz - 3. * sn * sn);
    /* Add correction for leading-order spin-induced quadrupole, relevant for BNS - this is zero when kappa_1,2=1 */
    /* Eq. 5.70 of BB1 */
    H = Hns + Hs + Hss;
    /* Add the additional calibrated term */
    /*H +=
        coeffs->dheffSS * eta * (sKerr_x * sStar_x + sKerr_y * sStar_y +
                    sKerr_z * sStar_z) / (r * r * r * r);*/
    /* One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric expressions of m1,m2 and S1,S2 */
    /*H += coeffs->dheffSSv2 * eta / (r*r*r*r) / (1.-4.*eta)
    * ( (sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z)*(1.-4.*eta+2.*eta*eta)
    +(sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z)*(-2.*eta+4.*eta*eta)
    +(sStar_x*sStar_x + sStar_y*sStar_y + sStar_z*sStar_z)*(2.*eta*eta) );*/
    H += coeffs->dheffSS * eta / (r * r * r * r)
        * (s1Vec->data[0] * s1Vec->data[0] + s1Vec->data[1] * s1Vec->data[1] +
        s1Vec->data[2] * s1Vec->data[2] + s2Vec->data[0] * s2Vec->data[0] +
        s2Vec->data[1] * s2Vec->data[1] + s2Vec->data[2] * s2Vec->data[2]);
    //printf( "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
    //printf( "H = %.16e\n", H );
    /* Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
    Hreal = sqrt (1. + 2. * eta * (H - 1.));
    return Hreal;
}

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  )
{

  REAL8 a2;
  REAL8 u, u2, u3, u4, u5;
  REAL8 m1PlusetaKK;

  REAL8 bulk;
  REAL8 logTerms;
  REAL8 deltaU;
  REAL8 deltaT;

  u = 1. / r;
  u2 = u * u;
  u3 = u2 * u;
  u4 = u2 * u2;
  u5 = u4 * u;

  a2 = a * a;

  m1PlusetaKK = -1. + eta * coeffs->KK;

  bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;

  logTerms =
    1. + eta * coeffs->k0 + eta * log (1. + coeffs->k1 * u + coeffs->k2 * u2 +
				       coeffs->k3 * u3 + coeffs->k4 * u4 +
				       coeffs->k5 * u5 +
				       coeffs->k5l * u5 * log (u));
  /*printf(" a = %.16e, u = %.16e\n",a,u);
     printf( "k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e , k4 = %.16e, k5 = %.16e, k5l = %.16e\n",coeffs->k0,
     coeffs->k1,coeffs->k2,coeffs->k3,coeffs->k4,coeffs->k5,coeffs->k5l);
     printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms ); */
  deltaU = bulk * logTerms;

  deltaT = r * r * deltaU;


  return deltaT;
}


/**
 * This function calculates the function \f$\Delta_r(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 10a and 10b
 */
REAL8
XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  )
{
  REAL8 u2, u3;
  REAL8 D;
  REAL8 deltaT;			/* The potential function, not a time interval... */
  REAL8 deltaR;

  u2 = 1. / (r * r);
  u3 = u2 / r;

  D = 1. + log (1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);

  deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT (coeffs, r, eta, a);

  deltaR = deltaT * D;
  return deltaR;
}


/**
 * Function to calculate the value of omega for the spin-aligned EOB waveform.
 * Can NOT be used in precessing cases. This omega is defined as \f$\dot{y}/r\f$ by setting \f$y=0\f$.
 * The function calculates omega = v/r, by first converting (r,phi,pr,pphi) to Cartesian coordinates
 * in which rVec={r,0,0} and pVec={0,pphi/r,0}, i.e. the effective-test-particle is positioned at x=r,
 * and its velocity along y-axis. Then it computes omega, which is now given by dydt/r = (dH/dp_y)/r.
 */
REAL8
XLALSimIMRSpinAlignedEOBCalcOmega (const REAL8 values[],/**<< Dynamical variables */
				   SpinEOBParams * funcParams,
							/**<< EOB parameters */
                    REAL8 STEP_SIZE_x /**<< Step size for numerical derivation of H */
  )
{

  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
  REAL8 cartValues[6];

  gsl_function F;
  INT gslStatus;

  REAL8 omega;
  REAL8 r;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */
  params.values = cartValues;
  params.params = funcParams;

  F.function = &GSLSpinAlignedHamiltonianWrapper;
  F.params = &params;

  /* Populate the Cartesian values vector */
  /* We can assume phi is zero wlog */
  memset (cartValues, 0, sizeof (cartValues));
  cartValues[0] = r = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate omega. In the chosen co-ordinate system, */
  /* we need dH/dpy to calculate this, i.e. varyParam = 4   */
  params.varyParam = 4;
  gslStatus = gsl_deriv_central (&F, cartValues[4],
					       STEP_SIZE_x, &omega, &absErr);

  if (gslStatus != GSL_SUCCESS)
    {
      print_warning ("XLAL Error - %s: Failure in GSL function\n", __func__);
      return REAL8_FAIL_NAN;
    }

  omega = omega / r;

  return omega;
}

/**
 * Function to calculate the non-Keplerian coefficient for the spin-aligned EOB model.
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$ defined in Eq. A2.
 * i.e. the function returns \f$(r_{\Omega} / r)^3\f$.
 */
REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff (const REAL8 values[],
					SpinEOBParams * funcParams
  )
{
  REAL8 omegaCirc;

  REAL8 tmpValues[4];

  REAL8 r3;

  /* We need to find the values of omega assuming pr = 0 */
  memcpy (tmpValues, values, sizeof (tmpValues));
  tmpValues[2] = 0.0;

  omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmega (tmpValues, funcParams, STEP_SIZE);
  if (IS_REAL8_FAIL_NAN (omegaCirc))
    {
      return REAL8_FAIL_NAN;
    }

  r3 = values[0] * values[0] * values[0];

  return 1.0 / (omegaCirc * omegaCirc * r3);
}

