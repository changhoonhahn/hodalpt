# utils.py
from sbi import utils as Ut
import numpy as np
import torch
import h5py

def training_data(kmax, infn, sumstat):
    with h5py.File(infn, 'r') as f:
        k     = f['k'][:]
        p0    = f['p0'][:]
        p2    = f['p2'][:]
        theta = f['theta'][:]
        mask  = k <= kmax
        if sumstat == 'p0':
            x = p0[:,mask]
        elif sumstat == 'p2':
            x = p2[:,mask]
        elif sumstat == 'pk_all':
            x = np.concatenate([p0[:,mask], p2[:,mask]], axis=1)
    return theta, x

def train_val_split(theta, x, seed, n_test=100, save_path=None):
    rng = np.random.default_rng(seed)
    shuffled_idx = rng.permutation(len(x))

    x_all = x[shuffled_idx]
    y_all = theta[shuffled_idx]

    test_idx = shuffled_idx[:n_test]
    if save_path is not None:
        np.savez(save_path, 
                 test_idx=test_idx,
                 x_test=x_all[:n_test],
                 y_test=y_all[:n_test])

    x_train = torch.tensor(x_all[n_test:],       dtype=torch.float32)
    y_train = torch.tensor(y_all[n_test:],       dtype=torch.float32)

    return x_train, y_train



def get_prior_bounds():
    # from sobol
    cosmo_lower = torch.tensor([0.1, 0.03, 0.5, 0.8, 0.6])
    cosmo_upper = torch.tensor([0.5, 0.07, 0.9, 1.2, 1.0])
    # central from parslist_quijote_bestfit, 50% bounds
    alpt_lower = torch.tensor([3.8557e-01, 6.8758e-01, 5.0618e-01, 1.3292e+00, 5.5732e-01, 9.3604e-01,
            1.3567e+00, 1.4422e+00, 4.8722e-01, 1.6048e+00, 2.1101e+00, 1.4812e+00,
            4.2917e-01, 4.0281e-01, 4.0621e-01, 2.3636e-01, 1.2361e+00, 2.2217e+00,
            2.9492e-01, 2.4092e+01, 4.4946e+01, 4.6910e+00, 3.3169e+01, 3.5682e+01,
            8.2116e+00, 3.2208e+00, 4.9426e+01, 2.3394e+01, 5.2163e+01, 4.0580e+00,
            5.0913e+01, 2.1794e+01, 7.5953e-05, 3.2730e-05, 2.5717e-06, 7.5000e-09,
            4.4805e-05, 1.2004e-04, 2.9497e-05, 3.7650e-07, 3.1065e-06, 1.3229e-05,
            1.0483e-05, 2.9475e-07, 1.2750e-08, 5.0250e-08, 5.4000e-08, 3.0000e-09,
            5.4668e-01, 8.7390e-01, 9.8520e-01, 3.7080e-01])
    alpt_upper = torch.tensor([6.4262e-01, 1.1460e+00, 8.4363e-01, 2.2153e+00, 9.2887e-01, 1.5601e+00,
            2.2612e+00, 2.4036e+00, 8.1204e-01, 2.6747e+00, 3.5169e+00, 2.4687e+00,
            7.1528e-01, 6.7135e-01, 6.7702e-01, 3.9393e-01, 2.0601e+00, 3.7028e+00,
            4.9153e-01, 4.0153e+01, 7.4911e+01, 7.8184e+00, 5.5282e+01, 5.9470e+01,
            1.3686e+01, 5.3680e+00, 8.2376e+01, 3.8990e+01, 8.6938e+01, 6.7634e+00,
            8.4856e+01, 3.6324e+01, 1.2659e-04, 5.4550e-05, 4.2862e-06, 1.2500e-08,
            7.4675e-05, 2.0007e-04, 4.9162e-05, 6.2750e-07, 5.1775e-06, 2.2049e-05,
            1.7471e-05, 4.9125e-07, 2.1250e-08, 8.3750e-08, 9.0000e-08, 5.0000e-09,
            9.1113e-01, 1.4565e+00, 1.6420e+00, 6.1800e-01])

    lower_bounds = torch.cat([cosmo_lower, alpt_lower])  # (57,)
    upper_bounds = torch.cat([cosmo_upper, alpt_upper])  # (57,)
    return lower_bounds, upper_bounds

def get_prior():
    lower_bounds, upper_bounds = get_prior_bounds()
    return Ut.BoxUniform(low=lower_bounds, high=upper_bounds)