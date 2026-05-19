'''

script for constructing the Quijote Rockstar + HOD galaxy samples for test sets


'''
import os, sys
import time
import warnings
import h5py
import numpy as np
from hodalpt.sims import quijote as Q
from nbodykit.lab import ArrayCatalog, FFTPower
import h5py

warnings.filterwarnings('ignore', message='You have selected 18 bins')

i_lhc = int(sys.argv[1]) # LHC index
_t0 = time.time()

np.random.seed(i_lhc) 

path_quij = os.path.join(sys.argv[2], str(i_lhc))


def sample_HOD(): 
    ''' sample HOD parameters
    '''
    # logMmin, sigma_logM, logM0, logM1, alpha
    # hod_lower_bound = np.array([12.0, 0.1, 13.0, 13.0, 0.])
    # hod_upper_bound = np.array([14.0, 0.6, 15.0, 15.0, 1.5])    
    hod = {
        'logMmin': np.random.normal(12.97, 0.11, size=1)[0],
        'sigma_logM': max(np.random.normal(0.40, 0.1, size=1)[0], 1e-3),
        'logM0': np.random.normal(13.67, 0.3, size=1)[0],
        'logM1': np.random.normal(13.68, 0.31, size=1)[0],
        'alpha': np.random.normal(0.79, 0.26, size=1)[0],
        'Abias': np.random.normal(0.01, 0.16, size=1)[0],
        'eta_conc': max(np.random.normal(1.11,0.40, size=1)[0], 1e-3),
        'eta_cen': max(np.random.normal(0.31, 0.13, size=1)[0], 1e-3),
        'eta_sat': max(np.random.normal(0.85, 0.27, size=1)[0], 1e-3)
        }
    return hod

# sample HOD parameters from the original SimBIG prior 
hod = sample_HOD()

# populate Rockstar halos with galaxies using HOD 
gals = Q.HODgalaxies(hod, path_quij, z=0.5)

# impose RSD
xyz = Q.Box_RSD(gals, LOS=[0,0,1], Lbox=1000.)


cat = ArrayCatalog({'Position': xyz}, BoxSize=1000.) 
r   = FFTPower(cat, mode='2d', Nmesh=256, dk=0.005, kmin=0.008,
               Nmu=10, los=[0,0,1], poles=[0, 2])

fname = os.path.join(path_quij, 'spec.hdf5')

poles = r.poles
k      = poles['k']
p0     = poles['power_0'].real - poles.attrs['shotnoise']
p2     = poles['power_2'].real
nmodes = poles['modes']
theta_hod = np.array(list(hod.values()))
theta_c = np.loadtxt(path_quij+'/Cosmo_params.dat')

with h5py.File(fname, 'w') as f:
    f['xyz']      = xyz
    f['k']        = k
    f['p0k']      = p0
    f['p2k']      = p2
    f['nmodes']   = nmodes
    f['hod']      = theta_hod
    f['cosmo']    = theta_c
    f['shotnoise'] = poles.attrs['shotnoise']
    # save useful metadata
    f.attrs['N']       = poles.attrs['N1']
    f.attrs['BoxSize'] = 1000.
    f.attrs['Nmesh']   = 256
    f.attrs['kmin']    = 0.008
    f.attrs['dk']      = 0.005
print(f'Wrote sample {i_lhc} to {fname}. ({time.time() - _t0:.1f}s)')
