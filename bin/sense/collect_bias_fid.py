import h5py
import numpy as np
import os
import sys

NLB_bias_dir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/NLB'
HOD_bias_dir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/HOD'
kmax = 1.0

out_fn_NLB = os.path.expandvars('$WORK/hodalpt/bin/npe/bias_fid_NLB_data.hdf5')
out_fn_HOD = os.path.expandvars('$WORK/hodalpt/bin/npe/bias_fid_HOD_data.hdf5')


def load_data(outfn, kmax):
    with h5py.File(outfn, 'r') as f:
        k     = f['k'][:]
        p0    = f['p0'][:]
        p2    = f['p2'][:]
        theta = f['theta'][:]
        mask  = k <= kmax
    return k[mask], p0[mask], p2[mask], theta


def collect(bias_dir, kmax, out_fn):
    available = sorted(
        [f for f in os.listdir(bias_dir) if f.startswith('spec.') and f.endswith('.h5')],
        key=lambda f: int(f.split('.')[1])
    )
    print(f'Found {len(available)} spectra in {bias_dir}')

    all_theta, all_p0, all_p2 = [], [], []
    k = None

    for fname in available:
        fn = os.path.join(bias_dir, fname)
        k, p0, p2, theta = load_data(fn, kmax)
        all_theta.append(theta)
        all_p0.append(p0)
        all_p2.append(p2)

    with h5py.File(out_fn, 'w') as f:
        f.create_dataset('theta', data=np.array(all_theta))
        f.create_dataset('p0',    data=np.array(all_p0))
        f.create_dataset('p2',    data=np.array(all_p2))
        f.create_dataset('k',     data=k)

    print(f'Wrote {out_fn}: theta {np.array(all_theta).shape}, p0 {np.array(all_p0).shape}')


collect(NLB_bias_dir, kmax, out_fn_NLB)
collect(HOD_bias_dir, kmax, out_fn_HOD)
