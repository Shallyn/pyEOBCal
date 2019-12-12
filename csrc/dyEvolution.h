/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYEVOLUTION__
#define __INCLUDE_DYEVOLUTION__

#include "dyUtils.h"
#include "dyIntegrator.h"
#include "dyNewtonianMultipole.h"
#include "dyFactorizedWaveform.h"
#include "dyInitialCondition.h"
#include "dyHamiltonian.h"
#include "dyFactorizedFlux.h"


INT EvolutionCore(const REAL8 m1,
                  const REAL8 m2,
                  const REAL8 fMin,
                  const REAL8 ecc,
                  const REAL8 deltaT,
                  const REAL8 spin1z,
                  const REAL8 spin2z,
                  COMPLEX16TimeSeries **hout,
                  AdjParams    *adjParams,
                  CtrlParams   *ctrlParams);

INT applyDefaultAdjustableParameters(AdjParams *adjParams,
                                     const REAL8 m1,
                                     const REAL8 m2,
                                     const REAL8 spin1z,
                                     const REAL8 spin2z);

INT IterateNQCCorrectionCoeffs(const REAL8 m1,
                               const REAL8 m2,
                               const REAL8 fMin,
                               const REAL8 ecc,
                               const REAL8 deltaT,
                               const REAL8 spin1z,
                               const REAL8 spin2z,
                               REAL8Vector *tNR,
                               REAL8Vector *hrealNR,
                               REAL8Vector *himagNR,
                               EOBNonQCCoeffs    *nqcCoeffs_output,
                               COMPLEX16TimeSeries **hLM_output,
                               AdjParams    *adjParams,
                               CtrlParams   *ctrlParams,
                               REAL8 eps,
                               INT maxstep,
                               gboolean get_waveform);

#endif

