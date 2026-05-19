import h5py
import numpy as np
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

NLB_bias_dir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/NLB'
HOD_bias_dir = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/HOD'
kmax = 1.0
N_WORKERS = 16  # tune to your allocation's core count

out_fn_NLB = os.path.expandvars('$WORK/hodalpt/bin/npe/bias_fid_NLB_data.hdf5')
out_fn_HOD = os.path.expandvars('$WORK/hodalpt/bin/npe/bias_fid_HOD_data.hdf5')


def load_data(args):
    fn, kmax = args
    with h5py.File(fn, 'r') as f:
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
    n = len(available)
    print(f'Found {n} spectra in {bias_dir}')

    paths = [(os.path.join(bias_dir, fname), kmax) for fname in available]

    t0 = time.time()
    results = [None] * n

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(load_data, args): i for i, args in enumerate(paths)}
        done = 0
        for fut in as_completed(futures):
            results[futures[fut]] = fut.result()
            done += 1
            if done % 500 == 0 or done == n:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (n - done) / rate if rate > 0 else float('inf')
                print(f'  {done}/{n}  ({rate:.0f} files/s, ETA {eta:.0f}s)')

    k = results[0][0]
    all_p0    = np.array([r[1] for r in results])
    all_p2    = np.array([r[2] for r in results])
    all_theta = np.array([r[3] for r in results])

    with h5py.File(out_fn, 'w') as f:
        f.create_dataset('theta', data=all_theta)
        f.create_dataset('p0',    data=all_p0)
        f.create_dataset('p2',    data=all_p2)
        f.create_dataset('k',     data=k)

    print(f'Wrote {out_fn}: theta {all_theta.shape}, p0 {all_p0.shape}')


collect(NLB_bias_dir, kmax, out_fn_NLB)
collect(HOD_bias_dir, kmax, out_fn_HOD)
