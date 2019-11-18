/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_DYINTEGRATOR__
#define __INCLUDE_DYINTEGRATOR__

#include "dyUtils.h"


typedef struct tagRK4GSLIntegrator
{
    gsl_odeiv_step    *step;
    gsl_odeiv_control *control;
    gsl_odeiv_evolve  *evolve;
    
    gsl_odeiv_system  *sys;
    
    int (* dydt) (double t, const double y[], double dydt[], void * params);
    int (* stop) (double t, const double y[], double dydt[], void * params);
    
    int retries;        /* retries with smaller step when derivatives encounter singularity */
    int stopontestonly;    /* stop only on test, use tend to size buffers only */
    
    int returncode;
} RK4GSLIntegrator;

RK4GSLIntegrator *InitRK4GSLIntegrator(int dim,
                                       int (* dydt) (double t, const double y[], double dydt[], void *params),
                                       int (* stop) (double t, const double y[], double dydt[], void *params),
                                       double eps_abs, double eps_rel );

int playRungeKutta4(RK4GSLIntegrator *integrator,
                    void *params,
                    REAL8 *yinit,
                    REAL8 tinit,
                    REAL8 tend,
                    REAL8 deltat,
                    REAL8Array **yout);


void DestroyIntegrator(RK4GSLIntegrator *integrator);


int playRungeKutta4(RK4GSLIntegrator    *integrator,
                    void                *params,
                    REAL8               *yinit,
                    REAL8               tinit,
                    REAL8               tend,
                    REAL8               deltat,
                    REAL8Array          **yout);


#endif

