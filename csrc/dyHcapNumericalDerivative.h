/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYHCAPNUMERICALDERIVATIVE__
#define __INCLUDE_DYHCAPNUMERICALDERIVATIVE__

#include "dyUtils.h"

REAL8 XLALSpinHcapNumDerivWRTParam(
                 const INT paramIdx,      /**<< Index of the parameters */
                 const REAL8 values[],     /**<< Dynamical variables */
                 SpinEOBParams *funcParams /**<< EOB Parameters */
                 );

REAL8 GSLSpinHamiltonianWrapper( double x, void *params );

int XLALSpinAlignedHcapDerivative(
                  double  t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

#endif

