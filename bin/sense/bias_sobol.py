#!/bin/python
''' python script to run nonlocal bias on ALPT  mocks from the sobol sequence
    saves two spectra to bias/0 bias/1
'''
# directory where test sobol is located: /corral/utexas/AST25023/simbig/alpt/sobol/300
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from hodalpt import stats
from hodalpt.sims import alpt as CS
import h5py
from nbodykit.lab import ArrayCatalog, FFTPower
import time

t0 = time.time()

dir_sobol = sys.argv[1] #/corral/utexas/AST25023/simbig/alpt/sobol
i_sobol = int(sys.argv[2]) 
path_priors = sys.argv[3] #'/work/11053/mcasas/ls6/hodalpt/bin/sense/alpt_nlb_priors_50p.hdf5' the full path to priors 
silent = (sys.argv[4] == 'True') # silent = True

idx_priors = [(2*i_sobol), (2*i_sobol) + 1]
# bias directories 
dm_dir = os.path.join(dir_sobol, str(i_sobol), 'alpt')+'/'
outdir = os.path.join(dir_sobol, str(i_sobol), 'bias')
os.makedirs(outdir, exist_ok=True)

def read_samps(fn, idx):
    with h5py.File(fn, 'r') as f:
        i = idx  # sample index
        grp = f[f'samples/{i}']
        
        theta_gal = {
            'alpha': grp['alpha'][()],   # shape (4,4)
            'beta':  grp['beta'][()],
            'nmean': grp['nmean'][()],
        }
        
        theta_rsd = {k: grp['theta_rsd'][k][()] for k in grp['theta_rsd']}
    return theta_gal, theta_rsd

theta_gal_0, theta_rsd_0 = read_samps(path_priors, idx_priors[0])
theta_gal_1, theta_rsd_1 = read_samps(path_priors, idx_priors[1])

print('galaxy assignment running for sobol%i' % i_sobol)
outdir_bias0 = os.path.join(outdir, '0')
outdir_bias1 = os.path.join(outdir, '1')
os.makedirs(outdir_bias0, exist_ok=True)
os.makedirs(outdir_bias1, exist_ok=True)

xyz_g0 = CS.CSbox_galaxy(theta_gal_0, theta_rsd_0, dm_dir, bias_model='nonlocal0', subgrid=True, silent=silent)
xyz_g1 = CS.CSbox_galaxy(theta_gal_1, theta_rsd_1, dm_dir, bias_model='nonlocal0', subgrid=True, silent=silent)

# # write out xyz_g1 to bias/0 positions (remove later) 
# outfn_xyz = os.path.join(outdir_bias1, 'xyz_g.npy')
# np.save(outfn_xyz, xyz_g1)


print('calculating spectra for sobol%i' % i_sobol)
# spec_0 = stats.Pk_periodic(xyz_g0.T, Lbox=1000, Ngrid=256, Nmubin=20, fft='pyfftw', silent=silent,rsd=2)
# spec_1 = stats.Pk_periodic(xyz_g1.T, Lbox=1000, Ngrid=256, Nmubin=20, fft='pyfftw', silent=silent,rsd=2)

# nbodykit
cat0 = ArrayCatalog({'Position': xyz_g0}, BoxSize=1000.)
r0   = FFTPower(cat0, mode='2d', Nmesh=256, dk=0.005, kmin=0.008,
               Nmu=10, los=[0,0,1], poles=[0,2])

cat1 = ArrayCatalog({'Position': xyz_g1}, BoxSize=1000.)
r1   = FFTPower(cat1, mode='2d', Nmesh=256, dk=0.005, kmin=0.008,
               Nmu=10, los=[0,0,1], poles=[0,2])

def save_spectrum(fname, r, overwrite=True):
    if not overwrite and os.path.isfile(fname):
        raise FileExistsError(f"{fname} already exists.")
        
    """Save FFTPower multipoles to HDF5."""
    poles = r.poles
    k      = poles['k']
    p0     = poles['power_0'].real - poles.attrs['shotnoise']
    p2     = poles['power_2'].real
    nmodes = poles['modes']
    
    with h5py.File(fname, 'w') as f:
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

outfn0 = os.path.join(outdir_bias0, 'spec.hdf5')
outfn1 = os.path.join(outdir_bias1, 'spec.hdf5')

save_spectrum(outfn0, r0, overwrite=True)
save_spectrum(outfn1, r1, overwrite=True)

print('bias computed & spectra written for sobol%i completed in %.1f min' % (i_sobol, (time.time() - t0) / 60.))
# print('bias computed & spectra written for sobol%i' % i_sobol)
