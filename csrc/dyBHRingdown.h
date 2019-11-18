/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYBHRINGDOWN__
#define __INCLUDE_DYBHRINGDOWN__

#include "dyUtils.h"
INT XLALSimIMREOBGenerateQNMFreqV2(COMPLEX16Vector * modefreqs,
                              /**<< OUTPUT, complex freqs of overtones in unit of Hz */
    const REAL8 mass1,        /**<< The mass of the 1st component (in Solar masses) */
    const REAL8 mass2,        /**<< The mass of the 2nd component (in Solar masses) */
    const REAL8 spin1z,     /**<< The spin of the 1st object; only needed for spin waveforms */
    const REAL8 spin2z,     /**<< The spin of the 2nd object; only needed for spin waveforms */
    UINT l,                  /**<< The l value of the mode in question */
    INT m,                   /**<< The m value of the mode in question */
    UINT nmodes             /**<< The number of overtones that should be included (max 8) */
    );

#endif

