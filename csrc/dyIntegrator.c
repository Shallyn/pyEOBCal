/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyIntegrator.h"

void DestroyIntegrator(RK4GSLIntegrator *integrator)
{
    if (!integrator) return;
    
    if (integrator->evolve)  gsl_odeiv_evolve_free(integrator->evolve);
    if (integrator->control) gsl_odeiv_control_free(integrator->control);
    if (integrator->step)    gsl_odeiv_step_free(integrator->step);
    
    free( integrator->sys );
    free( integrator );
    
    return;
}


RK4GSLIntegrator *InitRK4GSLIntegrator(int dim,
                                       int (* dydt) (double t, const double y[], double dydt[], void *params),
                                       int (* stop) (double t, const double y[], double dydt[], void *params),
                                       double eps_abs, double eps_rel )
{
    RK4GSLIntegrator *integrator;
    integrator = (RK4GSLIntegrator *) calloc(1, sizeof (RK4GSLIntegrator));
    integrator->step    = gsl_odeiv_step_alloc(gsl_odeiv_step_rk4,dim) ;
    integrator->control = gsl_odeiv_control_y_new(eps_abs,eps_rel);
    integrator->evolve  = gsl_odeiv_evolve_alloc(dim);

    /* allocate the GSL system (functions, etc.) */
    integrator->sys = (gsl_odeiv_system *)calloc(1, sizeof(gsl_odeiv_system));

    if ( !(integrator->step) || !(integrator->control) || !(integrator->evolve) || !(integrator->sys) )
    {
        DestroyIntegrator(integrator);
    }

    integrator->dydt = dydt;
    integrator->stop = stop;
    
    integrator->sys->function  = dydt;
    integrator->sys->jacobian  = NULL;
    integrator->sys->dimension = dim;
    integrator->sys->params    = NULL;
    
    integrator->retries = 6;
    integrator->stopontestonly = 0;


    return integrator;
}


int playRungeKutta4(RK4GSLIntegrator *integrator,
                    void *params,
                    REAL8 *yinit,
                    REAL8 tinit,
                    REAL8 tend,
                    REAL8 deltat,
                    REAL8Array **yout)
{
    int failed = 0;
    int status; /* used throughout */
    /* needed for the integration */
    size_t dim, bufferlength, cnt, retries;
    REAL8 t, tnew, h0;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */
    /* needed for the final interpolation */
    gsl_spline *interp = NULL;
    gsl_interp_accel *accel = NULL;
    int outputlen = 0;
    REAL8Array *output = NULL;
    REAL8 *times, *vector;      /* aliases */

    /* allocate the buffers!
     * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
     * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat) + 2;  /* allow for the initial value and possibly a final semi-step */
    buffers = CreateREAL8Array(2, 2*dim + 1, bufferlength);  /* 2-dimensional array, (dim+1) x bufferlength */
    temp = (REAL8*)calloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp)
    {
        print_warning("Error _%s: Failed to apply memory.\n");
        failed = 1;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;
    
    integrator->returncode = 0;
    
    cnt = 0;
    retries = integrator->retries;
    
    t = tinit;
    h0 = deltat;
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    UINT i;
    for (i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];
    

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    //print_debug("START: dydt_ini = [%.5e, %.5e, %.5e, %.5e]\n", dydt_in[0], dydt_in[1], dydt_in[2], dydt_in[3]);
    //print_debug("START: y_init = [%.5e, %.5e, %.5e, %.5e]\n", y[0], y[1], y[2], y[3]);

    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS)
    {
        print_warning("Initial derivative failed.\n");
        failed = 1;
        integrator->returncode = status;
        goto bail_out;
    }
    //print_debug("START: dydt_ini = [%.5e, %.5e, %.5e, %.5e]\n", dydt_in[0], dydt_in[1], dydt_in[2], dydt_in[3]);
    while (1) 
    {
        if (!integrator->stopontestonly && t >= tend)
        {
            break;
        }
        
        if (integrator->stop)
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS)
            {
                integrator->returncode = status;
                break;
            }
        }
        
        /* ready to try stepping! */
    try_step:
        /* if we would be stepping beyond the final time, stop there instead... */
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        //print_err("STEP: dydt_in = [%.5e, %.5e, %.5e, %.5e]\n", dydt_in[0], dydt_in[1], dydt_in[2], dydt_in[3]);
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));
        //print_err("STEP: dydt_in0 = [%.5e, %.5e, %.5e, %.5e]\n", dydt_in0[0], dydt_in0[1], dydt_in0[2], dydt_in0[3]);

        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
         * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
         * and the error code from the user-supplied function will be returned. */
        //print_err("STEP: dydt_out = [%.5e, %.5f, %.5e, %.5e]\n\n", dydt_out[0], dydt_out[1], dydt_out[2], dydt_out[3]);
        
        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS)
        {

            print_warning("GSL status is FAILURE, retry = %d.\n", retries);
            if (retries--)
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            }
            else
            {
                failed = 1;
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        }
        /*
        else
        {
            //print_err("DEBUG: exit.\n");
            //return CEV_FAILURE;
            retries = integrator->retries;      // we stepped successfully, reset the singularity retries 
        }*/
        
        tnew = t + h0;
        
        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);
        
        /* did the error-checker reduce the stepsize?
         * note: other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
         * GSL_ODEIV_HADJ_NIL if it was unchanged */
        if (status == GSL_ODEIV_HADJ_DEC)
        {
            //print_debug("GSL-HADJ-DEC. retry.\n");
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }
        
        /* update the current time and input derivatives */
        t = tnew;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        cnt++;
        
        /* check if interpolation buffers need to be extended */
        if (cnt >= bufferlength) 
        {
            REAL8Array *rebuffers;
            
            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
             * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, 2*dim + 1, 2 * bufferlength)))
            {
                print_warning("Cannot reallocate buffers.\n");
                failed = 1;
                goto bail_out;
            } else {
                for (i = 0; i <= 2*dim; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], cnt * sizeof(REAL8));
                DestroyREAL8Array(buffers);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }
        
        /* copy time and state into interpolation buffers */
        buffers->data[cnt] = t;
        for (i = 1; i <= dim; i++)
        {
            buffers->data[i * bufferlength + cnt] = y[i - 1];   /* y does not have time */
            buffers->data[(i + dim)*bufferlength + cnt] = dydt_out[i-1];
        }
    }
    /* copy the final state into yinit */
    
    memcpy(yinit, y, dim * sizeof(REAL8));
    
    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (cnt == 0)
    {
        goto bail_out;
    }
    
    interp = gsl_spline_alloc(gsl_interp_cspline, cnt + 1);
    accel = gsl_interp_accel_alloc();
    
    outputlen = (int)(t / deltat) + 1;
    output = CreateREAL8Array(2, 2*dim + 1, outputlen);
    
    if (!interp || !accel || !output)
    {
        print_warning("Failed to allocate GSL interpolator.\n");
        failed = 1;
        if (output)
            DestroyREAL8Array(output);
        outputlen = 0;
        goto bail_out;
    }
    
    /* make an array of times */
    times = output->data;
    INT j;
    for (j = 0; j < outputlen; j++)
        times[j] = tinit + deltat * j;
    
    /* interpolate! */
    for (i = 1; i <= 2*dim; i++) 
    {
        gsl_spline_init(interp, &buffers->data[0], &buffers->data[bufferlength * i], cnt + 1);
        
        vector = output->data + outputlen * i;
        for (j = 0; j < outputlen; j++) 
        {
            gsl_spline_eval_e(interp, times[j], accel, &(vector[j]));
        }
    }
    
    /* deallocate stuff and return */
bail_out:
    
    if (buffers)
        DestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    if (temp)
        free(temp);
    
    if (interp)
        gsl_spline_free(interp);
    if (accel)
        gsl_interp_accel_free(accel);
    
    if (failed)
    {
        return CEV_FAILURE;
    }
    *yout = output;
    return outputlen;
}

