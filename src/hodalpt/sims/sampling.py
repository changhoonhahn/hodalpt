import numpy as np
import h5py
from pathlib import Path
_DAT_DIR = Path(__file__).parent / 'dat'

############################################################################################
# REPRODUCIBLE SAMPLING F'Ns
############################################################################################
def sample_NLB(seed, width=0.5):
    # function to write pm 50 percent ALPT nlb priors centered on quijote fiducial best fit.
    # alpha, beta, nmean are arrays of cenral best fit values (16,), width is desired prior width (percentile)
    # returns dictionaries for alpha, beta, nmean 
    pars = np.loadtxt(_DAT_DIR / 'parslist_quijote_bestfit.txt', delimiter=',', skiprows=1).T
    #parfile columns: nmean,normalization,alpha,beta
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
    dict_rsd = {}
    for key, val in theta_rsd.items():
        val_low = val - val*width
        val_high = val + val*width
        dict_rsd[key] = rng.uniform(val_low, val_high)

    return sample_alpha, sample_beta, sample_nmean, dict_rsd

def index_to_seed(i, master_seed=42):
    """Map sample index → unique deterministic seed, independent of n_samples."""
    return (master_seed * 2654435761 + i) & 0xFFFFFFFF

def generate_and_save_NLB(outfile, n_samples, width=0.5, master_seed=42):   
    """
    Reproducible by construction: sample i always has the same seed
    regardless of n_samples. Running with n=10 then n=100 gives identical
    first-10 dictionaries.
    """
    with h5py.File(outfile, 'w') as f:
        f.attrs['master_seed'] = master_seed
        f.attrs['width'] = width
        samples_grp = f.create_group('samples')
        for i in range(n_samples):
            s = index_to_seed(i, master_seed)
            sa, sb, sn, s_rsd= sample_NLB(seed=s, width=width)
            grp = samples_grp.create_group(str(i))
            grp.create_dataset('alpha',  data=sa.reshape(4, 4))
            grp.create_dataset('beta',   data=sb.reshape(4, 4))
            grp.create_dataset('nmean',  data=sn.reshape(4, 4))
            rsd_grp = grp.create_group('theta_rsd')
            for key, val in s_rsd.items():
                rsd_grp.create_dataset(key, data=val)
    
    print(f"Wrote {n_samples} samples → {outfile}")

def sample_HOD(seed): 
    ''' sample HOD parameters from gaussian around SIMBIG CMASS priors
    '''
    rng = np.random.default_rng(seed)
    hod = {
        'logMmin': rng.normal(12.97, 0.11),
        'sigma_logM': max(rng.normal(0.40, 0.1), 1e-3),
        'logM0': rng.normal(13.67, 0.3),
        'logM1': rng.normal(13.68, 0.31),
        'alpha': max(rng.normal(0.79, 0.26), 1e-3),
        'Abias': rng.normal(0.01, 0.16), 
        'eta_conc': max(rng.normal(1.11,0.40),1e-3),
        'eta_cen':max(rng.normal(0.31, 0.13),1e-3),
        'eta_sat': max(rng.normal(0.85, 0.27),1e-3)
        }
    return hod
def generate_and_save_HOD(outfile, n_samples, master_seed=42):  
    """
    Reproducible by construction: sample i always has the same seed
    regardless of n_samples. Running with n=10 then n=100 gives identical
    first-10 dictionaries.
    """
    with h5py.File(outfile, 'w') as f:

        f.attrs['master_seed'] = master_seed

        samples_grp = f.create_group('samples')
        for i in range(n_samples):
            s = index_to_seed(i, master_seed)
            hod = sample_HOD(s)
            grp = samples_grp.create_group(str(i))
            for key, val in hod.items():
                grp.create_dataset(key, data=val)
    
    print(f"Wrote {n_samples} samples → {outfile}")

#############################################################
# THETA VECTORIZATION  (defines canonical parameter ordering)
#############################################################
_HOD_KEYS = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha',
             'Abias', 'eta_conc', 'eta_cen', 'eta_sat']

_RSD_KEYS = ['bv', 'bb', 'betarsd', 'gamma']

def hod_to_vec(hod):
    """Flatten HOD dict to 1-D array in canonical order (_HOD_KEYS)."""
    return np.array([hod[k] for k in _HOD_KEYS])

def nlb_to_vec(theta_gal, theta_rsd):
    """Flatten NLB dicts to 1-D array: alpha(16), beta(16), nmean(16), rsd(4)."""
    return np.concatenate([
        theta_gal['alpha'].ravel(),
        theta_gal['beta'].ravel(),
        theta_gal['nmean'].ravel(),
        [theta_rsd[k] for k in _RSD_KEYS],
    ])


#############################################################
# READ-IN
#############################################################
def read_samps_NLB(fn, idx):
    with h5py.File(fn, 'r') as f:
        grp = f[f'samples/{idx}']
        theta_gal = {
            'alpha': grp['alpha'][()],   # shape (4,4)
            'beta':  grp['beta'][()],
            'nmean': grp['nmean'][()],
        }
        theta_rsd = {k: grp['theta_rsd'][k][()] for k in grp['theta_rsd']}
    return theta_gal, theta_rsd

def read_samps_HOD(fn, idx):
    with h5py.File(fn, 'r') as f:
        grp = f[f'samples/{idx}']
        hod = {k: grp[k][()] for k in grp}
    return hod