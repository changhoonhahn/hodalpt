import numpy as np
import h5py
from hodalpt import stats
import os
import sys
import time
import pyfftw
pyfftw.config.NUM_THREADS = 1
# import pyspectrum as pyS

i0 = int(sys.argv[1])
i1 = int(sys.argv[2])

NLB_dir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/NLB'
HOD_dir =  '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/HOD'

def load_pos(specfn):
    with h5py.File(specfn, 'r') as f:
        xyz = f['xyz'][:]
    return xyz

for i in range(i0,i1):
    t0 = time.time()
    specfn = os.path.join(NLB_dir, 'spec.%i.h5' %i)
    if not os.path.isfile(specfn):
        print('MISSING %i' % i)
        continue
    xyz = load_pos(specfn)
    bispec = stats.B0_periodic(xyz.T, w=None, Lbox=1000., fft='pyfftw', silent=True)
    with h5py.File(specfn, 'a') as f:
        if 'b123' not in f:
            f['i_k1'] = bispec['i_k1']
            f['i_k2'] = bispec['i_k2']
            f['i_k3'] = bispec['i_k3']
            f['b123'] = bispec['b123']
            f['q123'] = bispec['q123']
    print('bispectrum added for NLB sample %i' % i)
    if i < 1000:
        specfn = os.path.join(HOD_dir, 'spec.%i.h5' %i)
        if not os.path.isfile(specfn):
            print('MISSING HOD %i' % i)
        else:
            xyz = load_pos(specfn)
            bispec = stats.B0_periodic(xyz.T, w=None, Lbox=1000., fft='pyfftw', silent=True)
            with h5py.File(specfn, 'a') as f:
                if 'b123' not in f:
                    f['i_k1'] = bispec['i_k1']
                    f['i_k2'] = bispec['i_k2']
                    f['i_k3'] = bispec['i_k3']
                    f['b123'] = bispec['b123']
                    f['q123'] = bispec['q123']
            print('bispectrum added for HOD sample %i' % i)
    print('bias sample %i finished in %.1f s' % (i, time.time() - t0))
