'''

generate and run slurm scripts for TACC Lonestar

'''
import os, sys
import numpy as np


def run_alpt_lhc(i0, i1, time=1, queue='development'):
    ''' run ALPT for specific LHC realization 
    '''
    dir_lhc = '/corral/utexas/AST25023/simbig/quijote/latinhypercube_hr/'

    # write slurm file for submitting the job
    a = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J alpt.lhc.%i_%i' % (i0, i1),
        '#SBATCH -o o/alpt.lhc.%i_%i' % (i0, i1),
        '#SBATCH -p %s' % queue, 
        '#SBATCH -N 1',               
        '#SBATCH -n 1',               
        '#SBATCH --time=%s:00:00' % (str(time).zfill(2)),
        '#SBATCH -A AST25023', 
        '',
        "module purge ",
        "module load fftw3/3.3.10 gsl", 
        "module load intel",  
        "module load impi", 
        "", 
        "unset PYTHONPATH", 
        "source ~/.bashrc", 
        "", 
        "conda activate simbig",
        '',
        ''])
    
    for i_lhc in range(i0, i1): 
        if not check_alpt_runs(i_lhc): 
            a += "python /home1/11004/chahah/projects/hodalpt/bin/lhc_alpt/alpt.py %s %i True\n" % (dir_lhc, i_lhc)
        else: 
            print('%i already complete' % i_lhc) 

    # create the script.sh file, execute it and remove it
    f = open(os.path.join(os.environ['WORK'], 'script.slurm'),'w')
    f.write(a)
    f.close()
    os.system('sbatch %s' % os.path.join(os.environ['WORK'], 'script.slurm'))
    #os.system('rm script.slurm')
    return None


def check_ics(i_lhc): 
    ''' check that the ICs are available 
    '''
    import glob
    dir_lhc = '/corral/utexas/AST25023/simbig/quijote/latinhypercube_hr/'
    
    _has_file = True 
    if len(glob.glob(dir_lhc+str(i_lhc)+'/CAMB.params')) != 1: 
        print('LHC %i does not have CAMB.params' % i_lhc) 
        _has_file = False 
    if len(glob.glob(dir_lhc+str(i_lhc)+'/Cosmo_params.dat')) != 1: 
        print('LHC %i does not have Cosmo_params.dat' % i_lhc) 
        _has_file = False 
    if len(glob.glob(dir_lhc+str(i_lhc)+'/ICs/ics*hdf5')) != 8: 
        print('LHC %i does not have ICs' % i_lhc) 
        _has_file = False 
    return _has_file 


def check_alpt_runs(i_lhc): 
    ''' check that 
    '''
    import glob
    dir_lhc = '/corral/utexas/AST25023/simbig/quijote/latinhypercube_hr/'

    _has_file = True 
    if len(glob.glob(dir_lhc+str(i_lhc)+'/alpt/super_BOXpos*')) != 3: 
        print('%i did not finish' % i_lhc) 
        _has_file = False
    return _has_file 


def clean_ics(i_lhc): 
    ''' remove heavy ics
    '''
    dir_lhc = '/corral/utexas/AST25023/simbig/quijote/latinhypercube_hr/'

    if check_alpt_runs(i_lhc): 
        print('removing %i/ICs/ics*hdf5' % i_lhc)
        os.system('rm %s' % dir_lhc+str(i_lhc)+'/ICs/ics*hdf5')

    return None 


if __name__=="__main__": 
    #for i in range(2, 9):
    #    run_alpt_lhc(10+10*i, 20+10*i, queue='normal')
    #run_alpt_lhc(80, 100, time=3, queue='normal')
    #for i in range(100): clean_ics(i)
    #for i in range(100): check_alpt_runs(i) 
    run_alpt_lhc(0, 50, queue='development', time=1)

