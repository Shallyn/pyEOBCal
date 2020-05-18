/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_AUXFUNCTIONS__
#define __INCLUDE_AUXFUNCTIONS__

#include "../csrc/dyUtils.h"

INT auxCalculateCompareFlux(const REAL8 values[],
                            SpinEOBParams *seobParams,
                            REAL8 *Eflux,
                            REAL8 *FphiEOB,
                            REAL8 *EfluxPN,
                            REAL8 *LfluxPN,
                            REAL8 *EfluxNagar,
                            REAL8 *LfluxNagar,
                            REAL8 *EfluxFac,
                            REAL8 *LfluxFac);


INT auxSpinEOBRadiationReactionForceAll(const REAL8 values[],
                                     SpinEOBParams *seobParams,
                                     REAL8 *Fr,
                                     REAL8 *Fphi,
                                     REAL8 *fluxout,
                                     REAL8 *fluxCircout,
                                     REAL8 *Fr2PN,
                                     REAL8 *Fphi2PN,
                                     REAL8 *FrFac,
                                     REAL8 *FphiFac,
                                     REAL8 *FrResum,
                                     REAL8 *FphiResum);

INT aux2PNEOBRadiationReactionForce(const REAL8 values[],
                                    SpinEOBParams *seobParams,
                                    REAL8 *Fr,
                                    REAL8 *Fphi);

INT aux2PNFactorizedRadiationReactionForce(const REAL8 values[],
                                        SpinEOBParams *seobParams,
                                        REAL8 *Fr,
                                        REAL8 *Fphi);

#endif

