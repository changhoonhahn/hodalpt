import os
import numpy as np
import glob

def check_bias_runs(i_sobol, dir_sobol, overwrite=False):
    if overwrite:
        return False
    """Check that both bias spectra exist for realization i_sobol."""
    expected = [
        os.path.join(dir_sobol, str(i_sobol), 'bias', '0', 'spec.hdf5'),
        os.path.join(dir_sobol, str(i_sobol), 'bias', '1', 'spec.hdf5'),
    ]
    for f in expected:
        if not os.path.isfile(f):
            print('%i not complete (missing %s)' % (i_sobol, f))
            return False
    return True

def run_bias_sobol(i0, i1, path_priors, time=5, queue='normal', silent=True):
    """Run bias_sobol.py for a batch of sobol realizations [i0, i1)."""
    dir_sobol = '/corral/utexas/AST25023/simbig/alpt/sobol'
    scriptdir = os.path.dirname(os.path.abspath(__file__))

    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    a = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J bias.sobol.%i_%i' % (i0, i1),
        '#SBATCH -o o/bias.sobol.%i_%i' % (i0, i1),
        '#SBATCH -p %s' % queue,
        '#SBATCH -N 1',
        '#SBATCH -n 1',
        '#SBATCH --time=%s:%s:00' % (str(hr).zfill(2), str(mn).zfill(2)),
        '#SBATCH -A AST25023',
        '',
        'module purge',
        'module load intel',
        'module load impi',
        'module load fftw3/3.3.10',
        'module load gsl',
        '',
        'unset PYTHONPATH',
        'source ~/.bashrc',
        '',
        'conda activate simbig',
        '',
        ''
    ])

    any_to_run = False
    for i_sobol in range(i0, i1):
        if not check_bias_runs(i_sobol, dir_sobol):
            a += "python %s/bias_sobol.py %s %i %s %s\n" % (
                scriptdir, dir_sobol, i_sobol, path_priors, str(silent))
            any_to_run = True 
        else:
            print('%i already complete' % i_sobol)
    
    if not any_to_run:
        print('All realizations %i-%i already complete, skipping submission.' % (i0, i1))
        return None

    slurm_path = os.path.join(os.environ['WORK'], 'script.slurm')
    with open(slurm_path, 'w') as f:
        f.write(a)
    os.system('sbatch %s' % slurm_path)
    return None


def submit_all_batches(path_priors, n_total=2048, batch_size=128,
                       time_per_run=2.5, queue='normal', silent=True):
    """
    Submit all batches covering [0, n_total) in chunks of batch_size.
    time_per_run: minutes per bias_sobol.py call (used to estimate wall time).
    """
    time_hrs = (batch_size * time_per_run) / 60.0 * 1.2  # 20% buffer
    time_hrs = max(time_hrs, 0.5)

    for i0 in range(0, n_total, batch_size):
        i1 = min(i0 + batch_size, n_total)
        print(f'Submitting batch [{i0}, {i1})  wall time = {time_hrs:.2f}h')
        run_bias_sobol(i0, i1, path_priors, time=time_hrs, queue=queue, silent=silent)


if __name__ == '__main__':
    PATH_PRIORS = '/work/11053/mcasas/ls6/hodalpt/bin/sense/alpt_nlb_priors_50p.hdf5'
    dir_sobol = '/corral/utexas/AST25023/simbig/alpt/sobol'
    # debug
    for i_sobol in range(0, 4):
        result = check_bias_runs(i_sobol, dir_sobol)
        print('i_sobol=%i check_bias_runs=%s' % (i_sobol, result))
    # --- Single batch test ---
    # run_bias_sobol(0, 4, PATH_PRIORS, time=0.5, queue='normal', silent=False)

    # --- Submit all 2048 in batches of 128 ---
    submit_all_batches(PATH_PRIORS, n_total=2048, batch_size=128,
                       time_per_run=1.5, queue='normal', silent=True)