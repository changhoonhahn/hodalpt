'''

script to generate the sobol sequence for cosmological parameters. 

'''
import numpy as np
from scipy.stats import qmc


def sobol_sample_5d(n_samples, seed=42):
    ''' Generate Sobol samples for LCDM cosmological parameters

    Parameters
    ----------
    n_samples : int
    '''
    # quijote LHC bounds for consistency 
    bounds = [
            [0.1, 0.5], # Omega_m
            [0.03, 0.07], # Omega_b
            [0.5, 0.9], # h
            [0.8, 1.2], # n_s
            [0.6, 1.0]] # sigma_8

    sampler = qmc.Sobol(d=len(bounds), scramble=True, seed=seed)

    # Sobol performs best with powers of 2
    m = int(np.ceil(np.log2(n_samples)))
    u = sampler.random_base2(m=m)[:n_samples]

    lower = np.array([b[0] for b in bounds])
    upper = np.array([b[1] for b in bounds])

    samples = qmc.scale(u, lower, upper)

    return samples


if __name__ == "__main__":
    n_samples = 2048 # this *should* be plenty 
    samples = sobol_sample_5d(n_samples) 
    
    np.savetxt('cosmology_alpt_sobol%i.dat' % n_samples, samples, 
               header='Omega_m, Omega_b, h, n_s, sigma_8') 
        
