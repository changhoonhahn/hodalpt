'''

module for different priors usd in the hodalpt project 

'''
import numpy as np 


def sample_NLB(seed, width=0.5):
    ''' sample non-local bias parameters centered on values based on best-fit
    to Quijote fiducial 


    returns
    -------
    dict with alpha, beta, nmean, and rsd parameters
    function to write pm 50 percent ALPT nlb priors centered on quijote fiducial best fit.
    alpha, beta, nmean are arrays of cenral best fit values (16,), width is desired prior width (percentile)
    returns dictionaries for alpha, beta, nmean 
    '''
    # read in best-fit NBL parameters
    #parfile columns: nmean,normalization,alpha,beta
    fpars = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'sims', 'dat',
                'parslist_quijote_bestfit.txt') 
    pars = np.loadtxt(fpars, delimiter=',', skiprows=1).T
    alpha_arr = pars[2] 
    beta_arr = pars[3]
    nmean_arr = pars[0]

    theta_rsd = {
    'bv': 0.7289, 
    'bb': 1.1652,
    'betarsd': 1.3136, 
    'gamma': 0.4944}

    rng = np.random.default_rng(seed)
    sample_alpha = np.zeros_like(alpha_arr)
    sample_beta = np.zeros_like(beta_arr)
    sample_nmean = np.zeros_like(nmean_arr)
    for i in range(len(alpha_arr)): 
        alpha = alpha_arr[i]
        beta = beta_arr[i]
        nmean = nmean_arr[i]

        alow = alpha - alpha*width
        ahigh = alpha + alpha*width

        blow = beta - beta*width
        bhigh = beta + beta*width

        nlow = nmean - nmean*width
        nhigh = nmean + nmean*width

        sample_alpha[i] = rng.uniform(alow, ahigh)
        sample_beta[i] = rng.uniform(blow, bhigh)
        sample_nmean[i] = rng.uniform(nlow, nhigh)

    # sample rsd params
    theta = {}
    theta['alpha'] = sample_alpha 
    theta['beta'] = sample_beta
    theta['nmean'] = sample_beta
    for key, val in theta_rsd.items():
        val_low = val - val*width
        val_high = val + val*width
        theta[key] = rng.uniform(val_low, val_high)
    
    return sample_alpha, sample_beta, sample_nmean, dict_rsd


def sample_HOD(seed): 
    ''' sample HOD parameters from Gaussian priors set around SIMBIG CMASS
    constraints  
    '''
    rng = np.random.default_rng(seed)

    hod = {
        'logMmin': rng.normal(12.97, 0.11),
        'sigma_logM': max(rng.normal(0.40, 0.1), 1e-3),
        'logM0': rng.normal(13.67, 0.3),
        'logM1': rng.normal(13.68, 0.31),
        'alpha': max(rng.normal(0.79, 0.26), 1e-3),
        'Abias': rng.normal(0.01, 0.16),
        'eta_conc': rng.normal(1.11,0.40),
        'eta_cen':rng.normal(0.31, 0.13),
        'eta_sat': rng.normal(0.85, 0.27) 
        }
    return hod

