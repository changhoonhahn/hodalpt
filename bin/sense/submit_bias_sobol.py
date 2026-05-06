import os
import numpy as np

def run_bias_sobol_test(i_sobol, path_priors, time=1, queue='development', silent=True):
    ''' run bias_sobol.py for a single sobol realization (test run)
    '''
    dir_sobol = '/corral/utexas/AST25023/simbig/alpt/sobol'
    scriptdir = os.path.dirname(__file__)
    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    a = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J bias.sobol.%i' % i_sobol,
        '#SBATCH -o o/bias.sobol.%i' % i_sobol,
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
        ''])

    a += "python %s/bias_sobol.py %s %i %s %s\n" % (scriptdir, dir_sobol, i_sobol, path_priors, str(silent))

    f = open(os.path.join(os.environ['WORK'], 'script.slurm'), 'w')
    f.write(a)
    f.close()
    os.system('sbatch %s' % os.path.join(os.environ['WORK'], 'script.slurm'))
    return None

if __name__ == '__main__':
    run_bias_sobol_test(
        i_sobol=0,
        path_priors='/work/11053/mcasas/ls6/hodalpt/bin/sense/alpt_nlb_priors_50p.hdf5',
        time=1,
        queue='normal',
        silent=True
    )
