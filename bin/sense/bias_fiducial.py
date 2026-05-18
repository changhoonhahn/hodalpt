import os, sys
import numpy as np
import h5py
from hodalpt.sims.sampling import read_samps_NLB, read_samps_HOD, nlb_to_vec, hod_to_vec
# from hodalpt import stats
from hodalpt.sims import alpt as CS
from hodalpt.sims import quijote as Q
from nbodykit.lab import ArrayCatalog, FFTPower
import time

t0 = time.time()
'''
for one fiducial_hr realization:
    read in samples
    generate galaxy catalog from theta_bias
        save positions xyz_g
    compute power spectra
        save k p0k p2k
'''
i0 = int(sys.argv[1])
i1 = int(sys.argv[2])

bias_fn = 'NLB_samples_10000.hdf5'
hod_fn = 'HOD_samples_1000.hdf5'

dm_dir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/alpt/'
outdir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias'
path_quij = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0'
outdir_NLB = os.path.join(outdir,'NLB')
outdir_HOD =  os.path.join(outdir,'HOD')
os.makedirs(outdir_NLB, exist_ok=True)
os.makedirs(outdir_HOD, exist_ok=True)

def save_spectrum(fname, xyz, theta, overwrite=True):
    """Save FFTPower multipoles to HDF5."""
    if not overwrite and os.path.isfile(fname):
        raise FileExistsError(f"{fname} already exists.")
    
    cat = ArrayCatalog({'Position': xyz}, BoxSize=1000.) 
    r   = FFTPower(cat, mode='2d', Nmesh=256, dk=0.005, kmin=0.008,
               Nmu=10, los=[0,0,1], poles=[0, 2])

    poles = r.poles
    k      = poles['k']
    p0     = poles['power_0'].real - poles.attrs['shotnoise']
    p2     = poles['power_2'].real
    nmodes = poles['modes']
    
    with h5py.File(fname, 'w') as f:
        f['theta']    = theta
        f['xyz']      = xyz
        f['k']        = k
        f['p0']       = p0
        f['p2']       = p2
        f['nmodes']   = nmodes
        f['shotnoise'] = poles.attrs['shotnoise']
        # save useful metadata
        f.attrs['N']       = poles.attrs['N1']
        f.attrs['BoxSize'] = 1000.
        f.attrs['Nmesh']   = 256
        f.attrs['kmin']    = 0.008
        f.attrs['dk']      = 0.005



with h5py.File(hod_fn, 'r') as f:
    n_hod = len(f['samples'])




for i in range(i0, i1):
    print('computing bias for sample%i' % i)

    fname_NLB = outdir_NLB+'/spec.%i.h5' % i
    theta_gal, theta_rsd = read_samps_NLB(bias_fn, i)
    xyz_nlb = CS.CSbox_galaxy(theta_gal, theta_rsd, dm_dir, bias_model='nonlocal0', subgrid=True, silent=True)
    save_spectrum(fname_NLB, xyz_nlb, nlb_to_vec(theta_gal, theta_rsd))

    if i < n_hod:
        fname_HOD = outdir_HOD+'/spec.%i.h5' % i
        hod = read_samps_HOD(hod_fn, i)
        gals =  Q.HODgalaxies(hod, path_quij, z=0.5)
        xyz_hod = Q.Box_RSD(gals, LOS=[0,0,1], Lbox=1000.)
        save_spectrum(fname_HOD, xyz_hod, hod_to_vec(hod))
    print('sample%i done in %.1f min' % (i, (time.time() - t0) / 60.))

    

        