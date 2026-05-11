'''

script for constructing the Quijote Rockstar + HOD galaxy samples for test sets


'''
import os, sys
import numpy as np 
from hodalpt.sims import quijote as Q
from nbodykit.lab import ArrayCatalog, FFTPower


i_lhc = int(sys.argv[1]) # LHC index 

np.random.seed(i_lhc) 

path_quij = os.path.join(sys.argv[2], str(i_lhc))


def sample_HOD(): 
    ''' sample HOD parameters
    '''
    # logMmin, sigma_logM, logM0, logM1, alpha
    hod_lower_bound = np.array([12.0, 0.1, 13.0, 13.0, 0.])
    hod_upper_bound = np.array([14.0, 0.6, 15.0, 15.0, 1.5])    

    hod = {
        'logMmin': np.random.uniform(12.0, 14.0, size=1)[0],
        'sigma_logM': np.random.uniform(0.1, 0.6, size=1)[0],
        'logM0': np.random.uniform(13.0, 15.0, size=1)[0],
        'logM1': np.random.uniform(13.0, 15.0, size=1)[0],
        'alpha': np.random.uniform(0.0, 1.5, size=1)[0],
        'Abias': np.clip(0.2 * np.random.normal(size=1), -1., 1.)[0], 
        'eta_conc': np.random.uniform(0.2, 2.0, size=1)[0], 
        'eta_cen': np.random.uniform(0., 0.7, size=1)[0], 
        'eta_sat': np.random.uniform(0.2, 2.0, size=1)[0]
    }
    return hod

# sampel HOD parameters from the original SimBIG prior 
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

with h5py.File(fname, 'w') as f:
    f['xyz']      = xyz
    f['k']        = k
    f['p0k']      = p0
    f['p2k']      = p2
    f['nmodes']   = nmodes
    f['shotnoise'] = poles.attrs['shotnoise']
    # save useful metadata
    f.attrs['N']       = poles.attrs['N1']
    f.attrs['BoxSize'] = 1000.
    f.attrs['Nmesh']   = 256
    f.attrs['kmin']    = 0.008
    f.attrs['dk']      = 0.005
