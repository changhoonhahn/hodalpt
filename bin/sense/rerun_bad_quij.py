'''
Rerun quij_test.py for the bad-eta/sigma_logM indices identified in inspect_quij_eta.ipynb.
Submits a single SLURM job on Lonestar that reruns only those indices.
'''
import os
import numpy as np

BAD_IDX = [39, 66, 70, 151, 155, 257, 318, 340, 361]

def run_quijote_bad(indices=BAD_IDX, time=1, queue='development', silent=True):
    _dir = '/corral/utexas/AST25023/simbig/quijote/latinhypercube_hr/'
    scriptdir = os.path.dirname(os.path.abspath(__file__))

    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    label = 'quij.hod.rerun'
    a = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J %s' % label,
        '#SBATCH -o o/%s' % label,
        '#SBATCH -p %s' % queue,
        '#SBATCH -N 1',
        '#SBATCH -n 1',
        '#SBATCH --time=%s:%s:00' % (str(hr).zfill(2), str(mn).zfill(2)),
        '#SBATCH -A AST25022',
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

    for i_lhc in indices:
        a += "python %s/quij_test.py %i %s\n" % (scriptdir, i_lhc, _dir)

    slurm_path = os.path.join(os.environ['WORK'], 'script_rerun_bad.slurm')
    with open(slurm_path, 'w') as f:
        f.write(a)
    os.system('sbatch %s' % slurm_path)
    return None


if __name__ == '__main__':
    run_quijote_bad(BAD_IDX, time=1, queue='development')
