/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYHYBRIDRINGDOWN__
#define __INCLUDE_DYHYBRIDRINGDOWN__

#include "dyUtils.h"

 INT XLALSimIMREOBAttachFitRingdown(
    REAL8Vector * signal1, /**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
    REAL8Vector * signal2, /**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
    const REAL8 dt,        /**<< Sample time step (in seconds) */
    const REAL8 mass1,     /**<< First component mass (in Solar masses) */
    const REAL8 mass2,     /**<< Second component mass (in Solar masses) */
    REAL8 spin1z,    /**<<The spin of the first object;  */
    REAL8 spin2z,    /**<<The spin of the first object;  */
    REAL8Vector *timeVec, /**<< Vector containing the time values */
    REAL8Vector *matchrange /**<< Time values chosen as points for performing comb matching */
    );
INT XLALSimIMREOBHybridAttachRingdown(
                                       REAL8Vector *signal1,    /**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
                                       REAL8Vector *signal2,    /**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
                                       const INT   l,          /**<< Current mode l */
                                       const INT   m,          /**<< Current mode m */
                                       const REAL8  dt,         /**<< Sample time step (in seconds) */
                                       const REAL8  mass1,      /**<< First component mass (in Solar masses) */
                                       const REAL8  mass2,      /**<< Second component mass (in Solar masses) */
                                       const REAL8  spin1z,     /**<<The spin of the first object; only needed for spin waveforms */
                                       const REAL8  spin2z,     /**<<The spin of the second object; only needed for spin waveforms */
                                       REAL8Vector *timeVec,    /**<< Vector containing the time values */
                                       REAL8Vector *matchrange /**<< Time values chosen as points for performing comb matching */
);

#endif

