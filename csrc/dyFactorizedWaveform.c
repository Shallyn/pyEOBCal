/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyFactorizedWaveform.h"
#include "dyHamiltonian.h"
#include "dyNewtonianMultipole.h"

INT CalculateSpinFactorizedWaveformCoefficients(FacWaveformCoeffs *const coeffs,
                                                SpinEOBParams *restrict params,
                                                const REAL8 m1,
                                                const REAL8 m2,
                                                const REAL8 eta,
                                                REAL8 a,
                                                const REAL8 chiS,
                                                const REAL8 chiA)
{
    REAL8 eta2 = eta * eta;
    REAL8 eta3 = eta2 * eta;

    REAL8 dM, dM2, chiA2, chiS2;		//dM3;
    REAL8 aDelta, a2, a3;

    /* Combination which appears a lot */
    REAL8 m1Plus3eta, m1Plus3eta2, m1Plus3eta3;

    dM2 = 1. - 4. * eta;
    chiA2 = chiA * chiA;
    chiS2 = chiS * chiS;
    REAL8 chiA3 = chiA2 * chiA;
    REAL8 chiS3 = chiS2 * chiS;

    if (dM2 < 0. && dM2 > -1e-4 ) 
    {
        dM2 = 0.;
    }
    if (dM < 0)
    {
        return CEV_FAILURE;
    }

    dM = sqrt(dM2);
    if (m1 < m2)
    {
        dM = -dM;
    }

    aDelta = 0.;			// a value in delta_lm is 0 in SEOBNRv1, SEOBNRv2, and SEOBNRv4
    a2 = a * a;
    a3 = a2 * a;

    m1Plus3eta = -1. + 3. * eta;
    m1Plus3eta2 = m1Plus3eta * m1Plus3eta;
    m1Plus3eta3 = m1Plus3eta * m1Plus3eta2;

    memset (coeffs, 0, sizeof (FacWaveformCoeffs));

    coeffs->delta22vh3 = 7. / 3.;
    coeffs->delta22vh6 = (-4. * aDelta) / 3. + (428. * CST_PI) / 105.;

    coeffs->delta22vh6 =
	    -4. / 3. * (dM * chiA + chiS * (1 - 2 * eta)) +
	    (428. * CST_PI) / 105.;

    coeffs->delta22v8 = (20. * aDelta) / 63.;
    coeffs->delta22vh9 = -2203. / 81. + (1712. * CST_PI * CST_PI) / 315.;
    coeffs->delta22v5 = -24. * eta;
    coeffs->delta22v6 = 0.0;

    coeffs->rho22v2 = -43. / 42. + (55. * eta) / 84.;
    coeffs->rho22v3 = (-2. * (chiS + chiA * dM - chiS * eta)) / 3.;

    coeffs->rho22v4 =
	    -20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM) -
	    (33025. * eta) / 21168. + (19583. * eta2) / 42336.;

    coeffs->rho22v5 = (-34. * a) / 21.;

    coeffs->rho22v5 =
	    (-34. / 21. + 49. * eta / 18. + 209. * eta2 / 126.) * chiS +
	    (-34. / 21. - 19. * eta / 42.) * dM * chiA;

    coeffs->rho22v6 =
        1556919113. / 122245200. + (89. * a2) / 252. -
        (48993925. * eta) / 9779616. - (6292061. * eta2) / 3259872. +
        (10620745. * eta3) / 39118464. + (41. * eta * CST_PI * CST_PI) / 192.;
    coeffs->rho22v6l = -428. / 105.;
    coeffs->rho22v7 = (18733. * a) / 15876. + a * a2 / 3.;
    /* See https://dcc.ligo.org/T1600383 */
    coeffs->rho22v7 =
	    a3 / 3. + chiA * dM * (18733. / 15876. + (50140. * eta) / 3969. +
		(97865. * eta2) / 63504.) +
	    chiS * (18733. / 15876. + (74749. * eta) / 5292. -
		(245717. * eta2) / 63504. + (50803. * eta3) / 63504.);
    coeffs->rho22v8 =
	    -387216563023. / 160190110080. + (18353. * a2) / 21168. -
	    a2 * a2 / 8.;
    coeffs->rho22v8l = 9202. / 2205.;
    coeffs->rho22v10 = -16094530514677. / 533967033600.;
    coeffs->rho22v10l = 439877. / 55566.;


     //RC: The delta coefficient before were put inside the if(dM2). This is wrong.
     //it didn't effect the models before because this function is used only to
     //calculate the 22 mode and the flux of the other modes, but for the flux
     //you only need |h_lm|. This error was present in all the modes, and now is fixed.
     //Is not a problem for the 22 mode because there is no if(dM2) for the 22
    coeffs->delta21vh3 = 2. / 3.;
    coeffs->delta21vh6 = (-17. * aDelta) / 35. + (107. * CST_PI) / 105.;
    coeffs->delta21vh7 = (3. * aDelta * aDelta) / 140.;
    coeffs->delta21vh9 = -272. / 81. + (214. * CST_PI * CST_PI) / 315.;
    coeffs->delta21v5 = -493. * eta / 42.;
    if (dM2)
    {

        //coeffs->rho21v1   = (-3.*(chiS+chiA/dM))/(4.);
        coeffs->rho21v1 = 0.0;
        //coeffs->rho21v2   = -59./56 - (9.*chiAPlusChiSdM*chiAPlusChiSdM)/(32.*dM2) + (23.*eta)/84.;
	    coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
	    coeffs->rho21v3 = 0.0;
        /*coeffs->rho21v3   = (-567.*chiA*chiA*chiA - 1701.*chiA*chiA*chiS*dM
         + chiA*(-4708. + 1701.*chiS*chiS - 2648.*eta)*(-1. + 4.*eta)
         + chiS* dM3 *(4708. - 567.*chiS*chiS
         + 1816.*eta))/(2688.*dM3); */
        coeffs->rho21v4 =
	        -47009. / 56448. - (865. * a2) / 1792. - (405. * a2 * a2) / 2048. -
	        (10993. * eta) / 14112. + (617. * eta2) / 4704.;
        coeffs->rho21v5 =
	        (-98635. * a) / 75264. + (2031. * a * a2) / 7168. -
	        (1701. * a2 * a3) / 8192.;
        coeffs->rho21v6 =
	        7613184941. / 2607897600. + (9032393. * a2) / 1806336. +
	        (3897. * a2 * a2) / 16384. - (15309. * a3 * a3) / 65536.;
        coeffs->rho21v6l = -107. / 105.;
        coeffs->rho21v7 =
	        (-3859374457. * a) / 1159065600. - (55169. * a3) / 16384. +
	        (18603. * a2 * a3) / 65536. - (72171. * a2 * a2 * a3) / 262144.;
        coeffs->rho21v7l = 107. * a / 140.;
        coeffs->rho21v8 = -1168617463883. / 911303737344.;
        coeffs->rho21v8l = 6313. / 5880.;
        coeffs->rho21v10 = -63735873771463. / 16569158860800.;
        coeffs->rho21v10l = 5029963. / 5927040.;

        coeffs->f21v1 = (-3. * (chiS + chiA / dM)) / (2.);
	    coeffs->f21v3 =
	        (chiS * dM * (427. + 79. * eta) +
	        chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
        /*RC: New terms for SEOBNRv4HM, they are put to zero if use_hm == 0 */
        coeffs->f21v4 = 0.0;
        coeffs->f21v5 = 0.0;
        coeffs->f21v6 = 0.0;
        coeffs->f21v7c = 0;
        /* End new terms for SEOBNRv4HM */
    }
    else
    {
        coeffs->f21v1 = -3. * chiA / 2.;
	    coeffs->f21v3 =
	        (chiS * dM * (427. + 79. * eta) +
	        chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
        /* New terms for SEOBNRv4HM, they are put to zero if use_hm == 0 */
        coeffs->f21v4 = 0.0;
        coeffs->f21v5 = 0.0;
        coeffs->f21v6 = 0.0;
        /* End new terms for SEOBNRv4HM */
    }


  /* l = 3, Eqs. A9a - A9c for rho, Eqs. A15b and A15c for f,
     Eqs. 22 - 24 of DIN and Eqs. 27c - 27e of PBFRT for delta */
        coeffs->delta33vh3 = 13. / 10.;
        coeffs->delta33vh6 = (-81. * aDelta) / 20. + (39. * CST_PI) / 7.;
        coeffs->delta33vh9 = -227827. / 3000. + (78. * CST_PI * CST_PI) / 7.;
        coeffs->delta33v5 = -80897. * eta / 2430.;
    if (dM2)
    {
        coeffs->rho33v2 = -7. / 6. + (2. * eta) / 3.;
        //coeffs->rho33v3 = (chiS*dM*(-4. + 5.*eta) + chiA*(-4. + 19.*eta))/(6.*dM);
        coeffs->rho33v3 = 0.0;
        coeffs->rho33v4 =
            -6719. / 3960. + a2 / 2. - (1861. * eta) / 990. +
            (149. * eta2) / 330.;
        coeffs->rho33v5 = (-4. * a) / 3.;
        coeffs->rho33v6 = 3203101567. / 227026800. + (5. * a2) / 36.;
        coeffs->rho33v6l = -26. / 7.;
        coeffs->rho33v7 = (5297. * a) / 2970. + a * a2 / 3.;
        coeffs->rho33v8 = -57566572157. / 8562153600.;
        coeffs->rho33v8l = 13. / 3.;
        coeffs->rho33v10 = 0;
        coeffs->rho33v10l = 0;
        coeffs->f33v3 =
	        (chiS * dM * (-4. + 5. * eta) + chiA * (-4. + 19. * eta)) / (2. * dM);
        coeffs->f33v4 = 0;
        coeffs->f33v5 = 0;
        coeffs->f33v6 = 0;
        coeffs->f33vh6 = 0;
    }
    else
    {
        coeffs->f33v3 = chiA * 3. / 8.;
        coeffs->f33v4 = 0;
        coeffs->f33v5 = 0;
        coeffs->f33v6 = 0;
        coeffs->f33vh6 = 0;
    }

    coeffs->delta32vh3 = (10. + 33. * eta) / (-15. * m1Plus3eta);
    coeffs->delta32vh4 = 4. * aDelta;
    coeffs->delta32vh6 = (-136. * aDelta) / 45. + (52. * CST_PI) / 21.;
    coeffs->delta32vh9 = -9112. / 405. + (208. * CST_PI * CST_PI) / 63.;    
    coeffs->rho32v = (4. * chiS * eta) / (-3. * m1Plus3eta);
    coeffs->rho32v2 = (-4. * a2 * eta2) / (9. * m1Plus3eta2) + (328. - 1115. * eta +
                       320. * eta2) / (270. *
                                       m1Plus3eta);
    coeffs->rho32v2 = (328. - 1115. * eta +
					      320. * eta2) / (270. *
							      m1Plus3eta);
    coeffs->rho32v3 =
        (2. *
        (45. * a * m1Plus3eta3 -
        a * eta * (328. - 2099. * eta + 5. * (733. + 20. * a2) * eta2 -
		960. * eta3))) / (405. * m1Plus3eta3);
    coeffs->rho32v3 = 2. / 9. * a;
    coeffs->rho32v4 = a2 / 3. + (-1444528.
			       + 8050045. * eta - 4725605. * eta2 -
			       20338960. * eta3 +
			       3085640. * eta2 * eta2) / (1603800. *
							  m1Plus3eta2);
    coeffs->rho32v5 = (-2788. * a) / 1215.;
    coeffs->rho32v6 = 5849948554. / 940355325. + (488. * a2) / 405.;
    coeffs->rho32v6l = -104. / 63.;
    coeffs->rho32v8 = -10607269449358. / 3072140846775.;
    coeffs->rho32v8l = 17056. / 8505.;

    coeffs->delta31vh3 = 13. / 30.;
    coeffs->delta31vh6 = (61. * aDelta) / 20. + (13. * CST_PI) / 21.;
    coeffs->delta31vh7 = (-24. * aDelta * aDelta) / 5.;
    coeffs->delta31vh9 = -227827. / 81000. + (26. * CST_PI * CST_PI) / 63.;
    coeffs->delta31v5 = -17. * eta / 10.;
    if (dM2)
    {

        coeffs->rho31v2 = -13. / 18. - (2. * eta) / 9.;
        //coeffs->rho31v3  = (chiA*(-4. + 11.*eta) + chiS*dM*(-4. + 13.*eta))/(6.*dM);
        coeffs->rho31v3 = 0.0;
        coeffs->rho31v4 = 101. / 7128.
	        - (5. * a2) / 6. - (1685. * eta) / 1782. - (829. * eta2) / 1782.;
        coeffs->rho31v5 = (4. * a) / 9.;
        coeffs->rho31v6 = 11706720301. / 6129723600. - (49. * a2) / 108.;
        coeffs->rho31v6l = -26. / 63.;
        coeffs->rho31v7 = (-2579. * a) / 5346. + a * a2 / 9.;
        coeffs->rho31v8 = 2606097992581. / 4854741091200.;
        coeffs->rho31v8l = 169. / 567.;

        coeffs->f31v3 =
	        (chiA * (-4. + 11. * eta) +
	        chiS * dM * (-4. + 13. * eta)) / (2. * dM);
    }
    else
    {
        coeffs->f31v3 = -chiA * 5. / 8.;
    }

    /* l = 4, Eqs. A10a - A10d for delta, Eq. A15d for f
     Eqs. 25 - 28 of DIN and Eqs. 27f - 27i of PBFRT for delta */

    coeffs->delta44vh3 = (112. + 219. * eta) / (-120. * m1Plus3eta);
    coeffs->delta44vh6 = (-464. * aDelta) / 75. + (25136. * CST_PI) / 3465.;
    coeffs->delta44vh9 = 0.;

    coeffs->rho44v2 =
        (1614. - 5870. * eta + 2625. * eta2) / (1320. * m1Plus3eta);
    coeffs->rho44v3 =
        (chiA * (10. - 39. * eta) * dM +
        chiS * (10. - 41. * eta + 42. * eta2)) / (15. * m1Plus3eta);
    coeffs->rho44v4 =
        a2 / 2. + (-511573572. + 2338945704. * eta - 313857376. * eta2 -
        6733146000. * eta3 +
        1252563795. * eta2 * eta2) / (317116800. * m1Plus3eta2);
    coeffs->rho44v5 = (-69. * a) / 55.;
    coeffs->rho44v8 = 0.;
    coeffs->rho44v8l = 0.;
    coeffs->rho44v10 = 0.;
    coeffs->rho44v10l = 0;
    coeffs->rho44v6 = 16600939332793. / 1098809712000. + (217. * a2) / 3960.;
    coeffs->rho44v6l = -12568. / 3465.;

    coeffs->delta43vh3 = (486. + 4961. * eta) / (810. * (1. - 2. * eta));
    coeffs->delta43vh4 = (11. * aDelta) / 4.;
    coeffs->delta43vh6 = 1571. * CST_PI / 385.;
  if (dM2)
    {

        //coeffs->rho43v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
        coeffs->rho43v = 0.0;
        coeffs->rho43v2 =
	        (222. - 547. * eta + 160. * eta2) / (176. * (-1. + 2. * eta));
        coeffs->rho43v4 = -6894273. / 7047040. + (3. * a2) / 8.;
        coeffs->rho43v5 = (-12113. * a) / 6160.;
        coeffs->rho43v6 = 1664224207351. / 195343948800.;
        coeffs->rho43v6l = -1571. / 770.;

        coeffs->f43v =
	        (5. * (chiA - chiS * dM) * eta) / (2. * dM * (-1. + 2. * eta));
    }
  else
    {
        coeffs->f43v = -5. * chiA / 4.;
    }

    coeffs->delta42vh3 = (7. * (1. + 6. * eta)) / (-15. * m1Plus3eta);
    coeffs->delta42vh6 = (212. * aDelta) / 75. + (6284. * CST_PI) / 3465.;

    coeffs->rho42v2 =
        (1146. - 3530. * eta + 285. * eta2) / (1320. * m1Plus3eta);
    coeffs->rho42v3 =
        (chiA * (10. - 21. * eta) * dM +
        chiS * (10. - 59. * eta + 78. * eta2)) / (15. * m1Plus3eta);
    coeffs->rho42v4 =
        a2 / 2. + (-114859044. + 295834536. * eta + 1204388696. * eta2 -
	       3047981160. * eta3 -
	       379526805. * eta2 * eta2) / (317116800. * m1Plus3eta2);
    coeffs->rho42v5 = (-7. * a) / 110.;
    coeffs->rho42v6 = 848238724511. / 219761942400. + (2323. * a2) / 3960.;
    coeffs->rho42v6l = -3142. / 3465.;

    coeffs->delta41vh3 = (2. + 507. * eta) / (10. * (1. - 2. * eta));
    coeffs->delta41vh4 = (11. * aDelta) / 12.;
    coeffs->delta41vh6 = 1571. * CST_PI / 3465.;

    if (dM2)
    {

        //coeffs->rho41v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
        coeffs->rho41v = 0.0;
        coeffs->rho41v2 =
	        (602. - 1385. * eta + 288. * eta2) / (528. * (-1. + 2. * eta));
        coeffs->rho41v4 = -7775491. / 21141120. + (3. * a2) / 8.;
        coeffs->rho41v5 = (-20033. * a) / 55440. - (5 * a * a2) / 6.;
        coeffs->rho41v6 = 1227423222031. / 1758095539200.;
        coeffs->rho41v6l = -1571. / 6930.;

        coeffs->f41v =
	        (5. * (chiA - chiS * dM) * eta) / (2. * dM * (-1. + 2. * eta));
    }
    else
    {
        coeffs->f41v = -5. * chiA / 4.;
    }

    /* l = 5, Eqs. A11a - A11e for rho,
     Eq. 29 of DIN and Eqs. E1a and E1b of PBFRT for delta */
        coeffs->delta55vh3 =
            (96875. + 857528. * eta) / (131250. * (1 - 2 * eta));
        coeffs->delta55vh6 = 0;
        coeffs->delta55vh9 = 0;
    if (dM2)
    {

        coeffs->rho55v2 =
	        (487. - 1298. * eta + 512. * eta2) / (390. * (-1. + 2. * eta));
        coeffs->rho55v3 = (-2. * a) / 3.;
        coeffs->rho55v4 = -3353747. / 2129400. + a2 / 2.;
        coeffs->rho55v5 = -241. * a / 195.;
        coeffs->rho55v6 = 0.;
        coeffs->rho55v6l = 0.;
        coeffs->rho55v8 = 0.;
        coeffs->rho55v8l = 0.;
        coeffs->rho55v10 = 0.;
        coeffs->rho55v10l = 0.;
        coeffs->f55v3 = 0.;
        coeffs->f55v4 = 0.;
        coeffs->f55v5c = 0;

    }
    else
    {
        coeffs->f55v3 = 0;
        coeffs->f55v4 = 0;
        coeffs->f55v5c = 0;
    }



    coeffs->delta54vh3 = 8. / 15.;
    coeffs->delta54vh4 = 12. * aDelta / 5.;

    coeffs->rho54v2 = (-17448. + 96019. * eta - 127610. * eta2
		     + 33320. * eta3) / (13650. * (1. - 5. * eta +
						   5. * eta2));
    coeffs->rho54v3 = (-2. * a) / 15.;
    coeffs->rho54v4 = -16213384. / 15526875. + (2. * a2) / 5.;

    coeffs->delta53vh3 = 31. / 70.;
    if (dM2)
    {

        coeffs->rho53v2 =
	        (375. - 850. * eta + 176. * eta2) / (390. * (-1. + 2. * eta));
        coeffs->rho53v3 = (-2. * a) / 3.;
        coeffs->rho53v4 = -410833. / 709800. + a2 / 2.;
        coeffs->rho53v5 = -103. * a / 325.;
    }

    coeffs->delta52vh3 = 4. / 15.;
    coeffs->delta52vh4 = 6. * aDelta / 5.;

    coeffs->rho52v2 = (-15828. + 84679. * eta - 104930. * eta2
		     + 21980. * eta3) / (13650. * (1. - 5. * eta +
						   5. * eta2));
    coeffs->rho52v3 = (-2. * a) / 15.;
    coeffs->rho52v4 = -7187914. / 15526875. + (2. * a2) / 5.;

    coeffs->delta51vh3 = 31. / 210.;
    if (dM2)
    {

        coeffs->rho51v2 =
	        (319. - 626. * eta + 8. * eta2) / (390. * (-1. + 2. * eta));
        coeffs->rho51v3 = (-2. * a) / 3.;
        coeffs->rho51v4 = -31877. / 304200. + a2 / 2.;
        coeffs->rho51v5 = 139. * a / 975.;
    }

    /* l = 6, Eqs. A12a - A12f for rho, Eqs. E1c and E1d of PBFRT for delta */

    coeffs->delta66vh3 = 43. / 70.;

    coeffs->rho66v2 = (-106. + 602. * eta - 861. * eta2
		     + 273. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
    coeffs->rho66v3 = (-2. * a) / 3.;
    coeffs->rho66v4 = -1025435. / 659736. + a2 / 2.;

    coeffs->delta65vh3 = 10. / 21.;
    if (dM2)
    {

        coeffs->rho65v2 = (-185. + 838. * eta - 910. * eta2
			 + 220. * eta3) / (144. * (dM2 + 3. * eta2));
        coeffs->rho65v3 = -2. * a / 9.;
    }

    coeffs->delta64vh3 = 43. / 105.;

    coeffs->rho64v2 = (-86. + 462. * eta - 581. * eta2
		     + 133. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
    coeffs->rho64v3 = (-2. * a) / 3.;
    coeffs->rho64v4 = -476887. / 659736. + a2 / 2.;

    coeffs->delta63vh3 = 2. / 7.;
    if (dM2)
    {

        coeffs->rho63v2 = (-169. + 742. * eta - 750. * eta2
			 + 156. * eta3) / (144. * (dM2 + 3. * eta2));
        coeffs->rho63v3 = -2. * a / 9.;
    }

    coeffs->delta62vh3 = 43. / 210.;

    coeffs->rho62v2 = (-74. + 378. * eta - 413. * eta2
		     + 49. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
    coeffs->rho62v3 = (-2. * a) / 3.;
    coeffs->rho62v4 = -817991. / 3298680. + a2 / 2.;

    coeffs->delta61vh3 = 2. / 21.;
    if (dM2)
    {

        coeffs->rho61v2 = (-161. + 694. * eta - 670. * eta2
			 + 124. * eta3) / (144. * (dM2 + 3. * eta2));
        coeffs->rho61v3 = -2. * a / 9.;
    }

    /* l = 7, Eqs. A13a - A13g for rho, Eqs. E1e and E1f of PBFRT for delta */
    coeffs->delta77vh3 = 19. / 36.;
    if (dM2)
    {

        coeffs->rho77v2 = (-906. + 4246. * eta - 4963. * eta2
			 + 1380. * eta3) / (714. * (dM2 + 3. * eta2));
        coeffs->rho77v3 = -2. * a / 3.;
    }

    coeffs->rho76v2 = (2144. - 16185. * eta + 37828. * eta2 - 29351. * eta3
		     + 6104. * eta2 * eta2) / (1666. * (-1 + 7 * eta -
							14 * eta2 +
							7 * eta3));

    coeffs->delta75vh3 = 95. / 252.;
    if (dM2)
    {

        coeffs->rho75v2 = (-762. + 3382. * eta - 3523. * eta2
			 + 804. * eta3) / (714. * (dM2 + 3. * eta2));
        coeffs->rho75v3 = -2. * a / 3.;
    }

    coeffs->rho74v2 = (17756. - 131805. * eta + 298872. * eta2 - 217959. * eta3
		     + 41076. * eta2 * eta2) / (14994. * (-1. + 7. * eta -
							  14. * eta2 +
							  7. * eta3));

    coeffs->delta73vh3 = 19. / 84.;
    if (dM2)
    {

        coeffs->rho73v2 = (-666. + 2806. * eta - 2563. * eta2
			 + 420. * eta3) / (714. * (dM2 + 3. * eta2));
        coeffs->rho73v3 = -2. * a / 3.;
    }

    coeffs->rho72v2 = (16832. - 123489. * eta + 273924. * eta2 - 190239. * eta3
		     + 32760. * eta2 * eta2) / (14994. * (-1. + 7. * eta -
							  14. * eta2 +
							  7. * eta3));

    coeffs->delta71vh3 = 19. / 252.;
    if (dM2)
    {

        coeffs->rho71v2 = (-618. + 2518. * eta - 2083. * eta2
			 + 228. * eta3) / (714. * (dM2 + 3. * eta2));
        coeffs->rho71v3 = -2. * a / 3.;
    }

    /* l = 8, Eqs. A14a - A14h */

    coeffs->rho88v2 = (3482. - 26778. * eta + 64659. * eta2 - 53445. * eta3
		     + 12243. * eta2 * eta2) / (2736. * (-1. + 7. * eta -
							 14. * eta2 +
							 7. * eta3));

    if (dM2)
    {
        coeffs->rho87v2 =
	        (23478. - 154099. * eta + 309498. * eta2 - 207550. * eta3 +
	        38920 * eta2 * eta2) / (18240. * (-1 + 6 * eta - 10 * eta2 +
					   4 * eta3));
    }

    coeffs->rho86v2 = (1002. - 7498. * eta + 17269. * eta2 - 13055. * eta3
		     + 2653. * eta2 * eta2) / (912. * (-1. + 7. * eta -
						       14. * eta2 +
						       7. * eta3));

    if (dM2)
    {
        coeffs->rho85v2 = (4350. - 28055. * eta + 54642. * eta2 - 34598. * eta3
			 + 6056. * eta2 * eta2) / (3648. * (-1. + 6. * eta -
							    10. * eta2 +
							    4. * eta3));
    }

    coeffs->rho84v2 = (2666. - 19434. * eta + 42627. * eta2 - 28965. * eta3
		     + 4899. * eta2 * eta2) / (2736. * (-1. + 7. * eta -
							14. * eta2 +
							7. * eta3));

    if (dM2)
    {
        coeffs->rho83v2 =
	        (20598. - 131059. * eta + 249018. * eta2 - 149950. * eta3 +
	        24520. * eta2 * eta2) / (18240. * (-1. + 6. * eta - 10. * eta2 +
					    4. * eta3));
    }

        coeffs->rho82v2 = (2462. - 17598. * eta + 37119. * eta2 - 22845. * eta3
		     + 3063. * eta2 * eta2) / (2736. * (-1. + 7. * eta -
							14. * eta2 +
							7. * eta3));

    if (dM2)
    {
        coeffs->rho81v2 =
	        (20022. - 126451. * eta + 236922. * eta2 - 138430. * eta3 +
	        21640. * eta2 * eta2) / (18240. * (-1. + 6. * eta - 10. * eta2 +
					    4. * eta3));
    }
    return CEV_SUCCESS;
}



/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables. This is optimized for flux calculation,
 * by ignoring complex arguments and keeping only absolute values.
 * Changes:
 * (i) Complex Argument of Tlm not exponentiated.
 * (ii) \f$exp(\ii deltalm)\f$ set to 1.
 * Eq. 17 and the entire Appendix of the paper https://journals.aps.org/prd/abstract/10.1103/PhysRevD.86.024011.
 */
INT
XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform (COMPLEX16 * restrict hlm,
						      /**< OUTPUT, hlm waveforms */
						REAL8Vector * restrict values,
						      /**< dyanmical variables */
						const REAL8 v,
						      /**< velocity */
						const REAL8 Hreal,
						      /**< real Hamiltonian */
						const INT l,
						      /**< l mode index */
						const INT m,
						      /**< m mode index */
						SpinEOBParams * restrict params)
{
    /* Status of function calls */
    INT status;
    INT i;

    REAL8 eta;
    REAL8 r, pp, Omega, v2, /*vh, vh3, */ k, hathatk, eulerlogxabs;	//pr
    REAL8 Slm, rholm, rholmPwrl;
    REAL8 auxflm = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 hNewton;
    gsl_sf_result z2;
    REAL8 Tlm;

    /* Non-Keplerian velocity */
    REAL8 vPhi, vPhi2;

    /* Pre-computed coefficients */
    FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

    eta = params->eobParams->eta;


    /* Check our eta was sensible */
    if (eta > 0.25 && eta < 0.25 +1e-4) {
            eta = 0.25;
    }
    /*else if ( eta == 0.25 && m % 2 )
        {
        // If m is odd and dM = 0, hLM will be zero
        memset( hlm, 0, sizeof( COMPLEX16 ) );
        return XLAL_SUCCESS;
        } */

    r = values->data[0];
    //pr    = values->data[2];
    pp = values->data[3];

    v2 = v * v;
    Omega = v2 * v;
    //vh3     = Hreal * Omega;
    //vh    = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log (2.0 * (REAL8) m * v);

    /* Calculate the non-Keplerian velocity */
    vPhi =
        XLALSimIMRSpinAlignedEOBNonKeplerCoeff (values->data, params);

    if (IS_REAL8_FAIL_NAN (vPhi))
    {
        return CEV_FAILURE;
    }

    vPhi = r * cbrt (vPhi);
    vPhi *= Omega;
    vPhi2 = vPhi * vPhi;

    /* Calculate the newtonian multipole, 1st term in Eq. 17, given by Eq. A1 */
    status = XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole (&hNewton,
                                    vPhi2, r,
                                    values->data[1],
                                    (UINT) l, m,
                                    params->
                                    eobParams);
    if (status == CEV_FAILURE)
    {
        return CEV_FAILURE;
    }

    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
    if (((l + m) % 2) == 0)
    {
        Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
    else
    {
        Slm = v * pp;
    }
    //printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

    /* Calculate the absolute value of the Tail term,
    * 3rd term in Eq. 17, given by Eq. A6, and Eq. (42) of
    * http://arxiv.org/pdf/1212.4357.pdf */
    k = m * Omega;
    hathatk = Hreal * k;
    hathatksq4 = 4. * hathatk * hathatk;
    hathatk4pi = 4. * CST_PI * hathatk;
    /*
        gsl_sf_result lnr1, arg1;
        XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
        if (status != GSL_SUCCESS)
        {
        XLALPrintError("XLAL Error - %s: Error in GSL function\n", __func__ );
        XLAL_ERROR( XLAL_EFUNC );
        }
    */
    status = gsl_sf_fact_e (l, &z2);
    if (status != GSL_SUCCESS)
    {
        print_warning("XLAL Error - %s: Error in GSL function: %s\n",
                        __func__, gsl_strerror (status));
        return CEV_FAILURE;
    }
    /*
        COMPLEX16 Tlmold;
        Tlmold = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * (
        arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
        Tlmold /= z2.val;
    */
    /* Calculating the prefactor of Tlm, outside the multiple product */
    Tlmprefac = sqrt (hathatk4pi / (1. - exp (-hathatk4pi))) / z2.val;

    /* Calculating the multiple product factor */
    for (Tlmprodfac = 1., i = 1; i <= l; i++)
        {
        Tlmprodfac *= (hathatksq4 + (REAL8) i * i);
        }

    Tlm = Tlmprefac * sqrt (Tlmprodfac);

    /* Calculate the residue phase and amplitude terms */
    /* deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15, others  */
    /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
    /* auxflm is a special part of the 5th term in Eq. 17, given by Eq. A15 */
    /* Actual values of the coefficients are defined in the next function of this file */
    switch (l)
        {
        case 2:
        switch (abs (m))
        {
        case 2:
        rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
                                + v * (hCoeffs->rho22v4
                                    +
                                    v *
                                    (hCoeffs->
                                    rho22v5 +
                                    v *
                                    (hCoeffs->
                                    rho22v6 +
                                    hCoeffs->
                                    rho22v6l *
                                    eulerlogxabs +
                                    v *
                                    (hCoeffs->
                                    rho22v7 +
                                    v *
                                    (hCoeffs->
                                    rho22v8 +
                                    hCoeffs->
                                    rho22v8l *
                                    eulerlogxabs +
                                    (hCoeffs->
                                    rho22v10 +
                                    hCoeffs->
                                    rho22v10l *
                                    eulerlogxabs)
                                    * v2)))))));
        break;
        case 1:
        {
            rholm = 1. + v * (hCoeffs->rho21v1
                    + v * (hCoeffs->rho21v2 +
                        v * (hCoeffs->rho21v3 +
                        v * (hCoeffs->rho21v4 +
                            v * (hCoeffs->rho21v5 +
                                v * (hCoeffs->rho21v6 +
                                hCoeffs->rho21v6l *
                                eulerlogxabs +
                                v *
                                (hCoeffs->rho21v7 +
                                hCoeffs->rho21v7l *
                                eulerlogxabs +
                                v *
                                (hCoeffs->rho21v8 +
                                hCoeffs->rho21v8l *
                                eulerlogxabs +
                                (hCoeffs->
                                    rho21v10 +
                                    hCoeffs->
                                    rho21v10l *
                                    eulerlogxabs) *
                                v2))))))));
            auxflm = v * hCoeffs->f21v1 + v2 * v * hCoeffs->f21v3;
        }
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        case 3:
        switch (m)
        {
        case 3:
        rholm =
            1. + v2 * (hCoeffs->rho33v2 +
                v * (hCoeffs->rho33v3 +
                    v * (hCoeffs->rho33v4 +
                    v * (hCoeffs->rho33v5 +
                        v * (hCoeffs->rho33v6 +
                        hCoeffs->rho33v6l * eulerlogxabs +
                        v * (hCoeffs->rho33v7 +
                            (hCoeffs->rho33v8 +
                            hCoeffs->rho33v8l *
                            eulerlogxabs) * v))))));
        auxflm = v * v2 * hCoeffs->f33v3;
        break;
        case 2:
        rholm = 1. + v * (hCoeffs->rho32v
                    + v * (hCoeffs->rho32v2 +
                    v * (hCoeffs->rho32v3 +
                        v * (hCoeffs->rho32v4 +
                            v * (hCoeffs->rho32v5 +
                            v * (hCoeffs->rho32v6 +
                                hCoeffs->rho32v6l *
                                eulerlogxabs +
                                (hCoeffs->rho32v8 +
                                hCoeffs->rho32v8l *
                                eulerlogxabs) *
                                v2))))));
        break;
        case 1:
        rholm =
            1. + v2 * (hCoeffs->rho31v2 +
                v * (hCoeffs->rho31v3 +
                    v * (hCoeffs->rho31v4 +
                    v * (hCoeffs->rho31v5 +
                        v * (hCoeffs->rho31v6 +
                        hCoeffs->rho31v6l * eulerlogxabs +
                        v * (hCoeffs->rho31v7 +
                            (hCoeffs->rho31v8 +
                            hCoeffs->rho31v8l *
                            eulerlogxabs) * v))))));
        auxflm = v * v2 * hCoeffs->f31v3;
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        case 4:
        switch (m)
        {
        case 4:

        rholm = 1. + v2 * (hCoeffs->rho44v2
                    + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                                    +
                                    v *
                                    (hCoeffs->
                                    rho44v5 +
                                    (hCoeffs->
                                    rho44v6 +
                                    hCoeffs->
                                    rho44v6l *
                                    eulerlogxabs) *
                                    v))));
        break;
        case 3:
        rholm = 1. + v * (hCoeffs->rho43v
                    + v * (hCoeffs->rho43v2
                    + v2 * (hCoeffs->rho43v4 +
                        v * (hCoeffs->rho43v5 +
                            (hCoeffs->rho43v6 +
                            hCoeffs->rho43v6l *
                            eulerlogxabs) * v))));
        auxflm = v * hCoeffs->f43v;
        break;
        case 2:
        rholm = 1. + v2 * (hCoeffs->rho42v2
                    + v * (hCoeffs->rho42v3 +
                        v * (hCoeffs->rho42v4 +
                        v * (hCoeffs->rho42v5 +
                            (hCoeffs->rho42v6 +
                            hCoeffs->rho42v6l *
                            eulerlogxabs) * v))));
        break;
        case 1:
        rholm = 1. + v * (hCoeffs->rho41v
                    + v * (hCoeffs->rho41v2
                    + v2 * (hCoeffs->rho41v4 +
                        v * (hCoeffs->rho41v5 +
                            (hCoeffs->rho41v6 +
                            hCoeffs->rho41v6l *
                            eulerlogxabs) * v))));
        auxflm = v * hCoeffs->f41v;
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        case 5:
        switch (m)
        {
        case 5:
        rholm = 1. + v2 * (hCoeffs->rho55v2
                    + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4
                                    +
                                    v *
                                    (hCoeffs->
                                    rho55v5 +
                                    hCoeffs->
                                    rho55v6 * v))));
        break;
        case 4:
        rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
                                + hCoeffs->rho54v4 * v));
        break;
        case 3:
        rholm = 1. + v2 * (hCoeffs->rho53v2
                    + v * (hCoeffs->rho53v3 +
                        v * (hCoeffs->rho53v4 +
                        hCoeffs->rho53v5 * v)));
        break;
        case 2:
        rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
                                + hCoeffs->rho52v4 * v));
        break;
        case 1:
        rholm = 1. + v2 * (hCoeffs->rho51v2
                    + v * (hCoeffs->rho51v3 +
                        v * (hCoeffs->rho51v4 +
                        hCoeffs->rho51v5 * v)));
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        case 6:
        switch (m)
        {
        case 6:
        rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
                                + hCoeffs->rho66v4 * v));
        break;
        case 5:
        rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
        break;
        case 4:
        rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
                                + hCoeffs->rho64v4 * v));
        break;
        case 3:
        rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
        break;
        case 2:
        rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
                                + hCoeffs->rho62v4 * v));
        break;
        case 1:
        rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        case 7:
        switch (m)
        {
        case 7:
        rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
        break;
        case 6:
        rholm = 1. + hCoeffs->rho76v2 * v2;
        break;
        case 5:
        rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
        break;
        case 4:
        rholm = 1. + hCoeffs->rho74v2 * v2;
        break;
        case 3:
        rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
        break;
        case 2:
        rholm = 1. + hCoeffs->rho72v2 * v2;
        break;
        case 1:
        rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        case 8:
        switch (m)
        {
        case 8:
        rholm = 1. + hCoeffs->rho88v2 * v2;
        break;
        case 7:
        rholm = 1. + hCoeffs->rho87v2 * v2;
        break;
        case 6:
        rholm = 1. + hCoeffs->rho86v2 * v2;
        break;
        case 5:
        rholm = 1. + hCoeffs->rho85v2 * v2;
        break;
        case 4:
        rholm = 1. + hCoeffs->rho84v2 * v2;
        break;
        case 3:
        rholm = 1. + hCoeffs->rho83v2 * v2;
        break;
        case 2:
        rholm = 1. + hCoeffs->rho82v2 * v2;
        break;
        case 1:
        rholm = 1. + hCoeffs->rho81v2 * v2;
        break;
        default:
        return CEV_FAILURE;
        break;
        }
        break;
        default:
        return CEV_FAILURE;
        break;
        }

    /* Raise rholm to the lth power */
    rholmPwrl = 1.0;
    i = l;
    while (i--)
        {
        rholmPwrl *= rholm;
        }
    /* In the equal-mass odd m case, there is no contribution from nonspin terms,
    * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
    * In this case, we must ignore the nonspin terms directly, since the leading term defined by
    * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
    */
    if (eta == 0.25 && m % 2)
        {
        rholmPwrl = auxflm;
        }
    else
        {
        rholmPwrl += auxflm;
        }

    /*if (r > 8.5)
        {
        printf("YP::dynamics variables in waveform: %i, %i, %e, %e\n",l,m,r,pp);
        printf( "rholm^l = %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e\n", rholmPwrl, Tlm, 0., Slm, creal(hNewton), cimag(hNewton) );} */
    /* Put all factors in Eq. 17 together */
    *hlm = Tlm * Slm * rholmPwrl;
    *hlm *= hNewton;
    /*if (r > 8.5)
        {
        printf("YP::FullWave: %.16e,%.16e, %.16e\n",hlm->re,hlm->im,sqrt(hlm->re*hlm->re+hlm->im*hlm->im));
        } */
    return CEV_SUCCESS;
}


/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of the paper.
 */
INT
XLALSimIMRSpinEOBGetSpinFactorizedWaveform (COMPLEX16 * restrict hlm,
						      /**< OUTPUT, hlm waveforms */
					    REAL8Vector * restrict values,
						      /**< dyanmical variables */
					    const REAL8 v,
						      /**< velocity */
					    const REAL8 Hreal,
						      /**< real Hamiltonian */
					    const INT l,
						      /**< l mode index */
					    const INT m,
						      /**< m mode index */
					    SpinEOBParams * restrict params
						       /**< Spin EOB parameters */
  )
{
  /* Status of function calls */
  INT status;
  INT i;


  REAL8 eta;
  REAL8 r, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;	//pr
  REAL8 Slm, deltalm, rholm;
  COMPLEX16 auxflm = 0.0;
  COMPLEX16 Tlm, rholmPwrl;
  COMPLEX16 hNewton;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 vPhi, vPhi2;

  /* Pre-computed coefficients */
  FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

  if (abs (m) > l)
    {
      return CEV_FAILURE;
    }
  if (m == 0)
    {
      return CEV_FAILURE;
    }
  eta = params->eobParams->eta;

  /* Check our eta was sensible */
  if (eta > 0.25 && eta < 0.25 +1e-4) {
      eta = 0.25;
  }
  if (eta > 0.25)
    {
      print_warning
	("XLAL Error - %s: Eta = %.16f seems to be > 0.25 - this isn't allowed!\n",
	 __func__,eta);
      return CEV_FAILURE;
    }
  /*else if ( eta == 0.25 && m % 2 )
     {
     // If m is odd and dM = 0, hLM will be zero
     memset( hlm, 0, sizeof( COMPLEX16 ) );
     return XLAL_SUCCESS;
     } */

  r = values->data[0];
  //pr    = values->data[2];
  pp = values->data[3];

  v2 = v * v;
  Omega = v2 * v;
  vh3 = Hreal * Omega;
  vh = cbrt (vh3);
  eulerlogxabs = CST_GAMMA + log (2.0 * (REAL8) m * v);

  /* Calculate the non-Keplerian velocity */
  //params->alignedSpins
  if (1)
    {
      // YP: !!!!! SEOBNRv3devel temporary change !!!!!

	  vPhi =
	    XLALSimIMRSpinAlignedEOBNonKeplerCoeff (values->data, params);
      // YP: !!!!! SEOBNRv3devel temporary change !!!!!

      if (IS_REAL8_FAIL_NAN (vPhi))
	{
	  return CEV_FAILURE;
	}

      vPhi = r * cbrt (vPhi);
      vPhi *= Omega;
      vPhi2 = vPhi * vPhi;
    }
  else
    {
      vPhi = v;
      vPhi2 = v2;
    }
  /* Calculate the newtonian multipole, 1st term in Eq. 17, given by Eq. A1 */
  // YP: !!!!! SEOBNRv3devel temporary change !!!!!
  status = XLALSimIMRSpinEOBCalculateNewtonianMultipole (&hNewton, vPhi2, r,
							 values->data[1],
							 (UINT) l, m,
							 params->eobParams);

  if (status == CEV_FAILURE)
    {
      return CEV_FAILURE;
    }

  /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
  if (((l + m) % 2) == 0)
    {
      Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
  else
    {
      Slm = v * pp;
    }
  //printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

  /* Calculate the Tail term, 3rd term in Eq. 17, given by Eq. A6 */
  k = m * Omega;
  hathatk = Hreal * k;
  status = gsl_sf_lngamma_complex_e (l + 1.0, -2.0 * hathatk, &lnr1,
					  &arg1);
  if (status != GSL_SUCCESS)
    {
      print_warning ("XLAL Error - %s: Error in GSL function: %s\n",
                      __func__, gsl_strerror (status));
      return CEV_FAILURE;
    }
  status = gsl_sf_fact_e (l, &z2);
  if (status != GSL_SUCCESS)
    {
      print_warning ("XLAL Error - %s: Error in GSL function: %s\n",
                      __func__, gsl_strerror (status));
      return CEV_FAILURE;
    }
  Tlm =
    cexp ((lnr1.val + CST_PI * hathatk) +
	  I * (arg1.val + 2.0 * hathatk * log (4.0 * k / sqrt (CST_E))));
  Tlm /= z2.val;


  /* Calculate the residue phase and amplitude terms */
  /* deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15, others  */
  /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
  /* auxflm is a special part of the 5th term in Eq. 17, given by Eq. A15 */
  /* Actual values of the coefficients are defined in the next function of this file */
  switch (l)
    {
    case 2:
      switch (abs (m))
	{
	case 2:
	  deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
							+
							vh * vh *
							(hCoeffs->delta22vh9 *
							 vh))) +
	    hCoeffs->delta22v5 * v * v2 * v2 +
	    hCoeffs->delta22v6 * v2 * v2 * v2 +
	    hCoeffs->delta22v8 * v2 * v2 * v2 * v2;
	  rholm =
	    1. + v2 * (hCoeffs->rho22v2 +
		       v * (hCoeffs->rho22v3 +
			    v * (hCoeffs->rho22v4 +
				 v * (hCoeffs->rho22v5 +
				      v * (hCoeffs->rho22v6 +
					   hCoeffs->rho22v6l * eulerlogxabs +
					   v * (hCoeffs->rho22v7 +
						v * (hCoeffs->rho22v8 +
						     hCoeffs->rho22v8l *
						     eulerlogxabs +
						     (hCoeffs->rho22v10 +
						      hCoeffs->rho22v10l *
						      eulerlogxabs) *
						     v2)))))));
	  break;
	case 1:
	  {
	    deltalm = vh3 * (hCoeffs->delta21vh3 + vh3 * (hCoeffs->delta21vh6
							  +
							  vh *
							  (hCoeffs->
							   delta21vh7 +
							   (hCoeffs->
							    delta21vh9) * vh *
							   vh))) +
	      hCoeffs->delta21v5 * v * v2 * v2 +
	      hCoeffs->delta21v7 * v2 * v2 * v2 * v;
	    rholm =
	      1. + v * (hCoeffs->rho21v1 +
			v * (hCoeffs->rho21v2 +
			     v * (hCoeffs->rho21v3 +
				  v * (hCoeffs->rho21v4 +
				       v * (hCoeffs->rho21v5 +
					    v * (hCoeffs->rho21v6 +
						 hCoeffs->rho21v6l *
						 eulerlogxabs +
						 v * (hCoeffs->rho21v7 +
						      hCoeffs->rho21v7l *
						      eulerlogxabs +
						      v * (hCoeffs->rho21v8 +
							   hCoeffs->rho21v8l *
							   eulerlogxabs +
							   (hCoeffs->
							    rho21v10 +
							    hCoeffs->
							    rho21v10l *
							    eulerlogxabs) *
							   v2))))))));
                   auxflm = v * hCoeffs->f21v1;
               }
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 3:
      switch (m)
	{
	case 3:
    deltalm =
        vh3 * (hCoeffs->delta33vh3 +
             vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
             hCoeffs->delta33v5 * v * v2 * v2;
            rholm =
                1. + v2 * (hCoeffs->rho33v2 +
                           v * (hCoeffs->rho33v3 +
                                v * (hCoeffs->rho33v4 +
                                     v * (hCoeffs->rho33v5 +
                                          v * (hCoeffs->rho33v6 +
                                               hCoeffs->rho33v6l * eulerlogxabs +
                                               v * (hCoeffs->rho33v7 +
                                                    v * (hCoeffs->rho33v8 +
                                                     hCoeffs->rho33v8l *
                                                     eulerlogxabs)))))));
            auxflm = v * v2 * hCoeffs->f33v3;
	  break;
	case 2:
	  deltalm =
	    vh3 * (hCoeffs->delta32vh3 +
		   vh * (hCoeffs->delta32vh4 +
			 vh * vh * (hCoeffs->delta32vh6 +
				    hCoeffs->delta32vh9 * vh3)));
	  rholm =
	    1. + v * (hCoeffs->rho32v +
		      v * (hCoeffs->rho32v2 +
			   v * (hCoeffs->rho32v3 +
				v * (hCoeffs->rho32v4 +
				     v * (hCoeffs->rho32v5 +
					  v * (hCoeffs->rho32v6 +
					       hCoeffs->rho32v6l *
					       eulerlogxabs +
					       (hCoeffs->rho32v8 +
						hCoeffs->rho32v8l *
						eulerlogxabs) * v2))))));
	  break;
	case 1:
	  deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
							+
							vh *
							(hCoeffs->delta31vh7 +
							 hCoeffs->delta31vh9 *
							 vh * vh))) +
	    hCoeffs->delta31v5 * v * v2 * v2;
	  rholm =
	    1. + v2 * (hCoeffs->rho31v2 +
		       v * (hCoeffs->rho31v3 +
			    v * (hCoeffs->rho31v4 +
				 v * (hCoeffs->rho31v5 +
				      v * (hCoeffs->rho31v6 +
					   hCoeffs->rho31v6l * eulerlogxabs +
					   v * (hCoeffs->rho31v7 +
						(hCoeffs->rho31v8 +
						 hCoeffs->rho31v8l *
						 eulerlogxabs) * v))))));
	  auxflm = v * v2 * hCoeffs->f31v3;
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 4:
      switch (m)
	{
	case 4:
	  deltalm = vh3 * (hCoeffs->delta44vh3 + hCoeffs->delta44vh6 * vh3)
	    + hCoeffs->delta44v5 * v2 * v2 * v;

    rholm = 1. + v2 * (hCoeffs->rho44v2
             + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                    +
                    v *
                    (hCoeffs->
                     rho44v5 +
                     (hCoeffs->
                      rho44v6 +
                      hCoeffs->
                      rho44v6l *
                      eulerlogxabs) *
                     v))));
	  break;
	case 3:
	  deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
						       +
						       hCoeffs->delta43vh6 *
						       vh * vh));
	  rholm =
	    1. + v * (hCoeffs->rho43v +
		      v * (hCoeffs->rho43v2 +
			   v2 * (hCoeffs->rho43v4 +
				 v * (hCoeffs->rho43v5 +
				      (hCoeffs->rho43v6 +
				       hCoeffs->rho43v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f43v;
	  break;
	case 2:
	  deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
	  rholm = 1. + v2 * (hCoeffs->rho42v2
			     + v * (hCoeffs->rho42v3 +
				    v * (hCoeffs->rho42v4 +
					 v * (hCoeffs->rho42v5 +
					      (hCoeffs->rho42v6 +
					       hCoeffs->rho42v6l *
					       eulerlogxabs) * v))));
	  break;
	case 1:
	  deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
						       +
						       hCoeffs->delta41vh6 *
						       vh * vh));
	  rholm =
	    1. + v * (hCoeffs->rho41v +
		      v * (hCoeffs->rho41v2 +
			   v2 * (hCoeffs->rho41v4 +
				 v * (hCoeffs->rho41v5 +
				      (hCoeffs->rho41v6 +
				       hCoeffs->rho41v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f41v;
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 5:
      switch (m)
	{
	case 5:
	  deltalm =
	    hCoeffs->delta55vh3 * vh3 + hCoeffs->delta55v5 * v2 * v2 * v;
	  rholm =
	    1. + v2 * (hCoeffs->rho55v2 +
		       v * (hCoeffs->rho55v3 +
			    v * (hCoeffs->rho55v4 +
				 v * (hCoeffs->rho55v5 +
				      hCoeffs->rho55v6 * v))));
	  break;
	case 4:
	  deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
	  rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
						     + hCoeffs->rho54v4 * v));
	  break;
	case 3:
	  deltalm = hCoeffs->delta53vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho53v2
			     + v * (hCoeffs->rho53v3 +
				    v * (hCoeffs->rho53v4 +
					 hCoeffs->rho53v5 * v)));
	  break;
	case 2:
	  deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
	  rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
						     + hCoeffs->rho52v4 * v));
	  break;
	case 1:
	  deltalm = hCoeffs->delta51vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho51v2
			     + v * (hCoeffs->rho51v3 +
				    v * (hCoeffs->rho51v4 +
					 hCoeffs->rho51v5 * v)));
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 6:
      switch (m)
	{
	case 6:
	  deltalm = hCoeffs->delta66vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
						     + hCoeffs->rho66v4 * v));
	  break;
	case 5:
	  deltalm = hCoeffs->delta65vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
	  break;
	case 4:
	  deltalm = hCoeffs->delta64vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
						     + hCoeffs->rho64v4 * v));
	  break;
	case 3:
	  deltalm = hCoeffs->delta63vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
	  break;
	case 2:
	  deltalm = hCoeffs->delta62vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
						     + hCoeffs->rho62v4 * v));
	  break;
	case 1:
	  deltalm = hCoeffs->delta61vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 7:
      switch (m)
	{
	case 7:
	  deltalm = hCoeffs->delta77vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
	  break;
	case 6:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho76v2 * v2;
	  break;
	case 5:
	  deltalm = hCoeffs->delta75vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
	  break;
	case 4:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho74v2 * v2;
	  break;
	case 3:
	  deltalm = hCoeffs->delta73vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
	  break;
	case 2:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho72v2 * v2;
	  break;
	case 1:
	  deltalm = hCoeffs->delta71vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    case 8:
      switch (m)
	{
	case 8:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho88v2 * v2;
	  break;
	case 7:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho87v2 * v2;
	  break;
	case 6:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho86v2 * v2;
	  break;
	case 5:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho85v2 * v2;
	  break;
	case 4:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho84v2 * v2;
	  break;
	case 3:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho83v2 * v2;
	  break;
	case 2:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho82v2 * v2;
	  break;
	case 1:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho81v2 * v2;
	  break;
	default:
	  return CEV_FAILURE;
	  break;
	}
      break;
    default:
      return CEV_FAILURE;
      break;
    }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  for(i = 0 ; i < l ; i++) rholmPwrl *= rholm;
  /* In the equal-mass odd m case, there is no contribution from nonspin terms,
   * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
   * In this case, we must ignore the nonspin terms directly, since the leading term defined by
   * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
   */
  if (eta == 0.25 && m % 2)
    {
      rholmPwrl = auxflm;
    }
  else
    {
      rholmPwrl += auxflm;
    }
  /* Put all factors in Eq. 17 together */
  *hlm = Tlm * cexp (I * deltalm) * Slm * rholmPwrl;
  *hlm *= hNewton;
  /*if (r > 8.5)
     {
     printf("YP::FullWave: Reh = %.16e, Imh = %.16e, hAmp = %.16e, hPhi = %.16e\n",creal(*hlm),cimag(*hlm),cabs(*hlm),carg(*hlm));
     } */
  return CEV_SUCCESS;
}

