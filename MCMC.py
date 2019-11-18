#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 00:38:27 2019

@author: drizl
"""

import numpy as np
from scipy import integrate
import emcee
import time
from .SXS import SXSAdjustor
from optparse import OptionParser
from pathlib import Path

def parseargs(argv):
    from WTestLib.SXS import DEFAULT_TABLE, DEFAULT_SRCLOC

    parser = OptionParser(description='Waveform Calibrator With SXS')
    parser.add_option('--f-ini', type = 'float', default = 0.002, help = 'Initial orbital frequency[M]')
    parser.add_option('--f-min', type = 'float', default = 20, help = 'Initial orbital frequency[Hz]')
    parser.add_option('--srate', type = 'float', default = 16384, help = 'sample rate')
    parser.add_option('--SXS', type = 'str', default = '0001', help = 'SXS template for calibration')
    parser.add_option('--table', type = 'str', default = str(DEFAULT_TABLE), help = 'Path of SXS table.')
    parser.add_option('--srcloc', type = 'str', default = str(DEFAULT_SRCLOC), help = 'Path of SXS waveform data.')
    parser.add_option('--ecc', type = 'float', default = 0, help = 'Eccentricity of waveform.')
    parser.add_option('--jobtag', type = 'str', help = 'JobID for the code run')
    parser.add_option('--prefix', type = 'str', default = '.', help = 'dir for results saving.')
    parser.add_option('--max-step', type = 'int', default = 1000, help = 'Max MCMC step.')
    args, _ = parser.parse_args(argv)
    return args



def main(argv = None):
    args = parseargs(argv)

    fini = args.f_ini
    fmin = args.f_min
    srate = args.srate
    SXSnum = args.SXS
    SRCLOC = args.srcloc
    SXSTable = args.table
    ecc = args.ecc
    prefix = Path(args.prefix)
    Max_Steps = args.max_step

    if not prefix.exists():
        prefix.mkdir(parents = True)
    
    fchain = prefix / 'chain.txt'

    ADJ = SXSAdjustor(SXSnum, f_min_dimless = fini, f_min = fmin, srate = srate, srcloc = SRCLOC, table = SXSTable)
    def lnprob(pms):
        return ADJ.get_lnprob(pms, ecc = ecc)
    # Initial val [KK, dSS, dSO, dtPeak] 
    p0 = [ADJ.adjParamsV4.KK, ADJ.adjParamsV4.dSS, ADJ.adjParamsV4.dSO, ADJ.adjParamsV4.dtPeak]
    ndim, nwalkers = 4, 8
    print("starting p0")
    p0 =[p0 + 0.05*np.abs(p0)*np.random.randn(ndim) for i in range(nwalkers)]
    print("ending p0")
    #p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)) * 1e-8 + result
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a=2, threads=4)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a=6, threads=4)
    """
    The easiest and simplest indi-
    cator that things are going well is the acceptance fraction; it should be in the 0.2 to 0.5 range
    (there are theorems about this for specific problems; for example Gelman, Roberts, & Gilks
    1996). In principle, if the acceptance fraction is too low, you can raise it by decreasing the "a"
    parameter; and if it is too high, you can reduce it by increasing the "a" parameter. However,
    in practice, we find that "a = 2" is good in essentially all situations.

    If you underestimate the errors in your parameters, acceptance rate will be 100%,
    which means that you never reach regions of lower likelihood. On the other hand, if
    you overestimate the errors, acceptance rate will be 0%, which means that the code
    will take a long time to get a useful number of points.
    """


    pos, prob, state = sampler.run_mcmc(p0, 100)
    print ("Before burn-in, when nsteps is 100,Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction))) 
    sampler.reset()

    f = open(fchain, "w")
    f.close()

    start_time = time.time()
    i_count=0      # used to record the sampling steps
    R_c = np.linspace(0, 1, ndim)      # set initial value for the convergence criterion
    converg_state = False     # set initial value for the convergence state
    for result in sampler.sample(pos, iterations=Max_Steps, storechain=False):
        #print(" the value of result is : %s" %(result))
        position = result[0]
        f = open(fchain, "a")
        for k in range(position.shape[0]):   # k = 0,1,2,..., nwalkers-1
            for j in range(ndim):
                f.write("%.8f " %position[k,j])
            f.write("\n")   
        f.close()
        i_count = i_count+1  # used to record the sampling steps
        
        if (i_count>100 and i_count % 100 == 0): #check the convergence state of the parameters once every 1000 steps, after N_steps > 1000.        
            chains = np.loadtxt(fchain)
            i_c = 0    # used to record the number of parameters that have satisfied the convergence criterion 
            for i in range(ndim):
                para = chains[:,i]
                para_reshape = para.reshape((i_count,nwalkers))
                # the second half of each chain is used to check the convergence state
                para_reshape = para_reshape[i_count/2:i_count,:]   
                # mean of each walker
                walker_mean = np.mean(para_reshape, axis=0, keepdims=True) 
                # variance between each walker
                var_mean = np.var(walker_mean)                            
                # variance of each walker
                walker_var = np.var(para_reshape, axis=0, keepdims=True)   
                mean_var = np.mean(walker_var)
                
            # sample from one walker ==> one chain
            # For multiple (nwalkers) chains  the code computes the Gelman and Rubin "R statistic"
            # Please See Page 38 of "eprint arXiv:0712.3028" for the definitions of "R statistic"
            
                R_c[i] = (mean_var*(1.0-2.0/i_count)+var_mean*(1.0+1.0/nwalkers))/mean_var 

                if abs(R_c[i]-1) < 0.01:
                    i_c= i_c+1
            
            if i_c == ndim:    # When every parameter has satisfied the convergence criterion,  we can stop sampling.
                converg_state = True
        if converg_state:
            print ("The chains have converged very well!\n The chain steps= %d" %i_count)
            break 

    if not converg_state:
        print ("The chains have not converged well enough. The values of R-1 for the model parameters are:")  
        for i in range(ndim):
            print ("%.3f" %(R_c[i]-1))  
    
    print ("Time elapsed: %.3f" %((time.time()-start_time)/60), "mins")
