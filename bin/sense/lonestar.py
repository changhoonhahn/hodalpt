'''

generate and run slurm scripts for TACC Lonestar

'''
import os, sys
import numpy as np


def run_alpt_sobol(i0, i1, time=1, queue='development', silent=True):
    ''' run ALPT for specific LHC realization 
    '''
    _dir= '/corral/utexas/AST25023/simbig/alpt/sobol/'
    scriptdir = os.path.dirname(__file__)
    
    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    # write slurm file for submitting the job
    a = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J alpt.sobol.%i_%i' % (i0, i1),
        '#SBATCH -o o/alpt.sobol.%i_%i' % (i0, i1),
        '#SBATCH -p %s' % queue, 
        '#SBATCH -N 1',               
        '#SBATCH -n 1',               
        '#SBATCH --time=%s:%s:00' % (str(hr).zfill(2), str(mn).zfill(2)),
        '#SBATCH -A AST25023', 
        '',
        "module purge ",
        "module load intel",  
        "module load impi", 
        "module load fftw3/3.3.10",
        "module load gsl", 
        "", 
        "unset PYTHONPATH", 
        "source ~/.bashrc", 
        "", 
        "conda activate simbig",
        '',
        ''])
    
    for i_lhc in range(i0, i1): 
        if not check_alpt_runs(i_lhc): 
            a += "python %s/alpt.py %s %i %s\n" % (scriptdir, _dir, i_lhc, str(silent))
        else: 
            print('%i already complete' % i_lhc) 

    # create the script.sh file, execute it and remove it
    f = open(os.path.join(os.environ['WORK'], 'script.slurm'),'w')
    f.write(a)
    f.close()
    os.system('sbatch %s' % os.path.join(os.environ['WORK'], 'script.slurm'))
    #os.system('rm script.slurm')
    return None


def check_alpt_runs(i_lhc): 
    ''' check that 
    '''
    import glob
    dir_lhc = '/corral/utexas/AST25023/simbig/alpt/sobol/'

    _has_file = True 
    if len(glob.glob(dir_lhc+str(i_lhc)+'/alpt/super_BOXpos*')) != 3: 
        print('%i did not finish' % i_lhc) 
        _has_file = False
    return _has_file 


if __name__=="__main__": 
    #for i in range(100): check_alpt_runs(i) 
    #run_alpt_sobol(0, 50, queue='development', time=1)

    run_alpt_sobol(0, 1, queue='development', time=0.5, silent=False)
