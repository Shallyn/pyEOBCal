/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyHcapNumericalDerivative.h"
#include "dyFactorizedFlux.h"
#include "dyHamiltonian.h"
#include <gsl/gsl_deriv.h>

#define STEP_SIZE 1.0e-3

/**
 * Calculate the derivative of the Hamiltonian w.r.t. a specific parameter
 * Used by generic spin EOB model, including initial conditions solver.
 */
REAL8 XLALSpinHcapNumDerivWRTParam(
                 const INT paramIdx,      /**<< Index of the parameters */
                 const REAL8 values[],     /**<< Dynamical variables */
                 SpinEOBParams *funcParams /**<< EOB Parameters */
                 )
{
     HcapDerivParams params;

     REAL8 result;

     gsl_function F;
     INT         gslStatus;

     REAL8 mass1, mass2;

     /* The error in a derivative as measured by GSL */
     REAL8 absErr;

     /* Set up pointers for GSL */
     params.values  = values;
     params.params  = funcParams;

     F.function       = &GSLSpinHamiltonianWrapper;
     F.params         = &params;
     params.varyParam = paramIdx;

     mass1 = params.params->eobParams->m1;
     mass2 = params.params->eobParams->m2;

     /* Now calculate derivatives w.r.t. the required parameter */
     if ( paramIdx >=6 && paramIdx < 9 )
     {
          gslStatus = gsl_deriv_central( &F, values[paramIdx],
                         STEP_SIZE*mass1*mass1, &result, &absErr );
     }
     else if ( paramIdx >= 9 )
     {
          gslStatus = gsl_deriv_central( &F, values[paramIdx],
                         STEP_SIZE*mass2*mass2, &result, &absErr );
     }
     else
     {
          gslStatus = gsl_deriv_central( &F, values[paramIdx],
                         STEP_SIZE, &result, &absErr );
     }
     if ( gslStatus != GSL_SUCCESS )
     {
          return REAL8_FAIL_NAN;
     }

  //printf( "Abserr = %e\n", absErr );

  return result;
}

/**
 * Wrapper for GSL to call the Hamiltonian function
 */
REAL8 GSLSpinHamiltonianWrapper( double x, void *params )
{
     HcapDerivParams *dParams = (HcapDerivParams *)params;

     EOBParams *eobParams = dParams->params->eobParams;

     REAL8 tmpVec[12];
     REAL8 s1normData[3], s2normData[3], sKerrData[3], sStarData[3];

     /* These are the vectors which will be used in the call to the Hamiltonian */
     REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
     REAL8Vector sigmaKerr, sigmaStar;

     int i;
     REAL8 a;
     REAL8 m1 = eobParams->m1;
     REAL8 m2 = eobParams->m2;
     REAL8 mT2 = (m1+m2)*(m1+m2);

     /* Use a temporary vector to avoid corrupting the main function */
     memcpy( tmpVec, dParams->values, sizeof(tmpVec) );

     /* Set the relevant entry in the vector to the correct value */
     tmpVec[dParams->varyParam] = x;

     /* Set the LAL-style vectors to point to the appropriate things */
     r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
     sigmaKerr.length = sigmaStar.length = 3;
     r.data     = tmpVec;
     p.data     = tmpVec+3;
     spin1.data = tmpVec+6;
     spin2.data = tmpVec+9;
     spin1norm.data = s1normData;
     spin2norm.data = s2normData;
     sigmaKerr.data = sKerrData;
     sigmaStar.data = sStarData;

     memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) );
     memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) );

     for ( i = 0; i < 3; i++ )
     {
          s1normData[i] /= mT2;
          s2normData[i] /= mT2;
     }

     /* Calculate various spin parameters */
     CalculateSigmaKerr( &sigmaKerr, eobParams->m1,
                         eobParams->m2, &spin1, &spin2 );
     CalculateSigmaStar( &sigmaStar, eobParams->m1,
                         eobParams->m2, &spin1, &spin2 );
     a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0] + sigmaKerr.data[1]*sigmaKerr.data[1]
               + sigmaKerr.data[2]*sigmaKerr.data[2] );
     //printf( "a = %e\n", a );
     //printf( "aStar = %e\n", sqrt( sigmaStar.data[0]*sigmaStar.data[0] + sigmaStar.data[1]*sigmaStar.data[1] + sigmaStar.data[2]*sigmaStar.data[2]) );
     if ( isnan( a ) )
     {
          print_warning( "a is nan!!\n");
     }
     //XLALSimIMRCalculateSpinEOBHCoeffs( dParams->params->seobCoeffs, eobParams->eta, a );

     //printf( "Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs ) );
     return SpinEOBHamiltonian( eobParams->eta, &r, &p, 
          &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, 
          dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
}


#define lMax 8
int XLALSpinAlignedHcapDerivative(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
  HcapDerivParams params;

  /* Since we take numerical derivatives wrt dynamical variables */
  /* but we want them wrt time, we use this temporary vector in  */
  /* the conversion */
  REAL8           tmpDValues[6];

  /* Cartesian values for calculating the Hamiltonian */
  REAL8           cartValues[6];

  REAL8           H; //Hamiltonian
  REAL8           flux;

  gsl_function F;
  INT4         gslStatus;
  UINT i;

  REAL8Vector rVec, pVec;
  REAL8 rData[3], pData[3];

  /* We need r, phi, pr, pPhi to calculate the flux */
  REAL8       r;
  REAL8Vector polarDynamics;
  REAL8       polData[4];

  REAL8 Mtotal, eta;

  /* Spins */
  REAL8Vector *s1Vec = NULL;
  REAL8Vector *s2Vec = NULL;
  REAL8Vector *sKerr = NULL;
  REAL8Vector *sStar = NULL;

  REAL8 a;

  REAL8 omega;

  /* EOB potential functions */
  REAL8 DeltaT, DeltaR;
  REAL8 csi;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Declare NQC coefficients */
  EOBNonQCCoeffs *nqcCoeffs = NULL;

  /* Set up pointers for GSL */ 
  params.values  = cartValues;
  params.params  = (SpinEOBParams *)funcParams;
  nqcCoeffs = params.params->nqcCoeffs;

  s1Vec = params.params->s1VecOverMtMt;
  s2Vec = params.params->s2VecOverMtMt;
  sKerr = params.params->sigmaKerr;
  sStar = params.params->sigmaStar;

  F.function = &GSLSpinAlignedHamiltonianWrapper;
  F.params   = &params;

  Mtotal = params.params->eobParams->Mtotal;
  eta   = params.params->eobParams->eta;

  r = values[0];

  /* Since this is spin aligned, I make the assumption */
  /* that the spin vector is along the z-axis.         */
  a  = sKerr->data[2];

  /* Calculate the potential functions and the tortoise coordinate factor csi,
     given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
  DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );

  DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );

  csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
  //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
  //printf( "csi in derivatives function = %.16e\n", csi );

  /* Populate the Cartesian values vector, using polar coordinate values */
  /* We can assume phi is zero wlog */
  memset( cartValues, 0, sizeof( cartValues ) );
  cartValues[0] = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate derivatives w.r.t. each Cartesian variable */
  for ( i = 0; i < 6; i++ )
  {
    params.varyParam = i;
    gslStatus = gsl_deriv_central( &F, cartValues[i], 
                    STEP_SIZE, &tmpDValues[i], &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
      print_warning( "XLAL Error - %s: Failure in GSL function\n", __func__ );
      return CEV_FAILURE;
    }
  }

  /* Calculate the Cartesian vectors rVec and pVec */
  polarDynamics.length = 4;
  polarDynamics.data   = polData;

  memcpy( polData, values, sizeof( polData ) );

  rVec.length = pVec.length = 3;
  rVec.data   = rData;
  pVec.data   = pData;

  memset( rData, 0, sizeof(rData) );
  memset( pData, 0, sizeof(pData) );

  rData[0] = values[0];
  pData[0] = values[2];
  pData[1] = values[3] / values[0];
  /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
  H =  SpinEOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );

  //printf( "csi = %.16e, ham = %.16e ( tortoise = %d)\n", csi, H, params.params->tortoise );
  //exit(1);
  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
  H = H * Mtotal;

  /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
      tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/
  /* Now calculate omega, and hence the flux */

  omega = tmpDValues[4] / r;
     dvalues[0] = csi * tmpDValues[3];
    dvalues[1] = omega;
    
    dvalues[2] = -tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    dvalues[2] = dvalues[2] * csi;
    dvalues[3] = 0;

  //flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, nqcCoeffs, omega, params.params, H/Mtotal, lMax ,1);
    flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, values, dvalues, nqcCoeffs, omega, params.params, H/Mtotal, lMax);
    if (IS_REAL8_FAIL_NAN(flux) || isnan(flux) )
    {
        print_warning("Failed to calculate flux.\n");
        //print_err("\tdvalues = [%f, %f, %f]\n", dvalues[0], dvalues[1], dvalues[2]);
        return CEV_FAILURE;
    }

  /* Looking at the non-spinning model, I think we need to divide the flux by eta */
  flux = flux / eta;

  //printf( "Flux in derivatives function = %.16e\n", flux );

  /* Now we can calculate the final (spherical) derivatives */
  /* csi is needed because we use the tortoise co-ordinate */
  /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
  //dvalues[0] = csi * tmpDValues[3];
  //dvalues[1] = omega;
  /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
  dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
  dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux / omega;
  dvalues[3] = - flux / omega;

  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Values:\n%f %f %f %f\n", values[0], values[1], values[2], values[3] );

  //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Derivatives:\n%f %f %f %f\n", dvalues[0], r*dvalues[1], dvalues[2], dvalues[3] );

  if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
  {
    print_warning( "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
    return CEV_FAILURE;
  }
  return CEV_SUCCESS;
}
#undef lMax
