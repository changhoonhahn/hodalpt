import os
import numpy as np

def run_npe_optuna(sumstat, kmax, queue='gpu-a100', time=24, n_trials=100):
    '''
    Submit npe_optuna.py as a SLURM job on Lonestar6

    Parameters
    ----------
    sumstat : str
        Summary statistic, e.g. 'p0', 'p2', 'pk_all'
    kmax : float
        Maximum k value
    queue : str
        SLURM partition. Use 'gpu-a100' for GPU nodes
    time : float
        Wall time in hours
    n_trials : int
        Number of Optuna trials
    '''
    scriptdir = os.path.dirname(os.path.abspath(__file__))

    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    # format kmax string to match study_name convention
    if int(kmax * 100) % 10 == 0:
        kmax_str = 'kmax%.1f' % kmax
    else:
        kmax_str = 'kmax%.2f' % kmax

    job_name = 'npe.maf.%s.%s' % (sumstat, kmax_str)

    a = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J %s' % job_name,
        '#SBATCH -o o/%s' % job_name,
        '#SBATCH -p %s' % queue,
        '#SBATCH -N 1',
        '#SBATCH -n 1',
        #'#SBATCH --gpus=1',                          # request one GPU
        '#SBATCH --time=%s:%s:00' % (str(hr).zfill(2), str(mn).zfill(2)),
        '#SBATCH -A AST25023',
        '',
        'module purge',
        'module load intel',
        'module load impi',
        "module load fftw3/3.3.10",
        "module load gsl", 
        '',
        'unset PYTHONPATH',
        'source ~/.bashrc',
        '',
        'conda activate simbig',
        '',
        'python %s/npe_optuna.py %s %s %s' % (scriptdir, sumstat, kmax, n_trials),
        ''])

    f = open(os.path.join(os.environ['WORK'], 'script.slurm'), 'w')
    f.write(a)
    f.close()
    os.system('sbatch %s' % os.path.join(os.environ['WORK'], 'script.slurm'))
    return None


if __name__ == '__main__':
    run_npe_optuna('p0', 0.2, queue='development', time=0.5, n_trials=1)
