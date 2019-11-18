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

