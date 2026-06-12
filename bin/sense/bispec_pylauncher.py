'''
bispec_pylauncher.py

Parallelizes bispectrum computation for fiducial_HR NLB and HOD samples
across TACC nodes using pylauncher. Each commandline computes bispectrum
for a single sample index, so pylauncher can saturate all available cores.

Scale (measured):
    i < 1000  (NLB + HOD): ~623 s per sample
    i >= 1000 (NLB only):  ~312 s per sample (estimated ~half)
    8 nodes on Lonestar6 (128 cores/node = 1024 slots):
        wall time ≈ (1000×623 + 9000×312) / 1024 ≈ 57 min
    6h wall time default is intentionally conservative.

Usage (run locally to submit):
    python bispec_pylauncher.py 0 10000
    python bispec_pylauncher.py 0 10000 --nodes 16 --time 6 --queue normal
    python bispec_pylauncher.py 0 10000 --skip-done   # skip already-done samples
    python bispec_pylauncher.py resume 0 10000         # restart interrupted job

Do NOT call the 'launch' / 'resume-run' subcommands directly;
they are invoked by the SLURM script.
'''
import os
import sys
import numpy as np


NLB_DIR = '/corral/utexas/AST25023/simbig/quijote/fiducial_HR/0/bias/NLB'


def _pending_indices(i0, i1, skip_done=False):
    '''Return list of sample indices that still need bispectrum computed.'''
    indices = list(range(i0, i1))
    if not skip_done:
        return indices

    import h5py
    pending = []
    for i in indices:
        fn = os.path.join(NLB_DIR, 'spec.%i.h5' % i)
        try:
            with h5py.File(fn, 'r') as f:
                if 'b123' not in f:
                    pending.append(i)
        except Exception:
            pending.append(i)
    n_done = len(indices) - len(pending)
    if n_done:
        print('skipping %i already-complete samples' % n_done)
    return pending


def run_bispec_fiducial_pylauncher(i0, i1, nodes=8, time=6, queue='normal',
                                   skip_done=False):
    '''Submit bispectrum computation for fiducial_HR NLB/HOD samples via pylauncher.

    Generates one task per sample index. pylauncher distributes tasks across
    all cores on the requested nodes, dispatching new tasks as cores free up.
    HOD bispectrum is computed alongside NLB for i < 1000 (handled inside
    bispec_fid.py).

    Parameters
    ----------
    i0, i1 : int
        Sample index range [i0, i1).
    nodes : int
        Number of compute nodes. Lonestar6 has 128 cores/node;
        8 nodes = 1024 parallel slots, finishing ~10k samples in ~1 h wall time.
    time : float
        Requested wall-clock hours.
    queue : str
        SLURM partition ('normal', 'development').
    skip_done : bool
        Check h5 files and omit samples with b123 already written.
    '''
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    workdir   = os.environ.get('WORK', os.getcwd())
    scratch   = os.environ.get('SCRATCH', workdir)

    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    indices = _pending_indices(i0, i1, skip_done=skip_done)
    if not indices:
        print('all samples already complete, nothing to submit')
        return None
    print('%i samples to process' % len(indices))

    # per-sample log directory on scratch (pylauncher also captures output in
    # its workdir, but per-sample files make debugging easier)
    logdir = os.path.join(scratch, 'bispec_logs_%i_%i' % (i0, i1))
    os.makedirs(logdir, exist_ok=True)

    # commandlines file: one task per index
    # pylauncher dispatches these across cores as they become free
    cmdfile = os.path.join(scratch, 'bispec_cmds_%i_%i.txt' % (i0, i1))
    with open(cmdfile, 'w') as f:
        for i in indices:
            logfile = os.path.join(logdir, 'sample_%i.log' % i)
            f.write('python %s/bispec_fid.py %i %i > %s 2>&1\n'
                    % (scriptdir, i, i + 1, logfile))
    print('wrote %i commandlines to %s' % (len(indices), cmdfile))

    # pylauncher workdir persists queuestate for --resume
    pyl_workdir = os.path.join(scratch, 'pylauncher_bispec_%i_%i' % (i0, i1))

    slurm = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J bispec.pyl.%i_%i' % (i0, i1),
        '#SBATCH -o %s/bispec.pyl.%i_%i.%%j.out' % (workdir, i0, i1),
        '#SBATCH -e %s/bispec.pyl.%i_%i.%%j.err' % (workdir, i0, i1),
        '#SBATCH -p %s' % queue,
        '#SBATCH -N %i' % nodes,          # pylauncher uses -N; ignores -n
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
        'module load pylauncher',
        '',
        # call back into this script with 'launch' subcommand
        'python %s launch %s %s' % (os.path.abspath(__file__), cmdfile, pyl_workdir),
        '',
    ])

    slurmfile = os.path.join(workdir, 'bispec_pylauncher_%i_%i.slurm' % (i0, i1))
    with open(slurmfile, 'w') as f:
        f.write(slurm)
    print('submitting %s' % slurmfile)
    os.system('sbatch %s' % slurmfile)
    return None


def resume_bispec_fiducial_pylauncher(i0, i1, nodes=8, time=4, queue='normal'):
    '''Resume an interrupted pylauncher bispectrum job.

    Reads the queuestate file from the original pylauncher workdir to skip
    tasks that already completed and re-run the rest.
    '''
    workdir = os.environ.get('WORK', os.getcwd())
    scratch = os.environ.get('SCRATCH', workdir)

    hr = int(np.floor(time))
    mn = int((time * 60) % 60)

    pyl_workdir = os.path.join(scratch, 'pylauncher_bispec_%i_%i' % (i0, i1))
    queuestate  = os.path.join(pyl_workdir, 'queuestate')

    if not os.path.exists(queuestate):
        raise FileNotFoundError(
            'no queuestate at %s — has this job been run before?' % queuestate)

    slurm = '\n'.join([
        '#!/bin/bash',
        '#SBATCH -J bispec.pyl.resume.%i_%i' % (i0, i1),
        '#SBATCH -o %s/bispec.pyl.resume.%i_%i.%%j.out' % (workdir, i0, i1),
        '#SBATCH -e %s/bispec.pyl.resume.%i_%i.%%j.err' % (workdir, i0, i1),
        '#SBATCH -p %s' % queue,
        '#SBATCH -N %i' % nodes,
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
        'module load pylauncher',
        '',
        'python %s resume-run %s' % (os.path.abspath(__file__), queuestate),
        '',
    ])

    slurmfile = os.path.join(workdir, 'bispec_pylauncher_resume_%i_%i.slurm' % (i0, i1))
    with open(slurmfile, 'w') as f:
        f.write(slurm)
    print('resuming from %s' % queuestate)
    os.system('sbatch %s' % slurmfile)
    return None


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    subcmd = sys.argv[1]

    if subcmd == 'launch':
        # called from SLURM — runs pylauncher to distribute commandlines
        import pylauncher
        cmdfile     = sys.argv[2]
        pyl_workdir = sys.argv[3]
        pylauncher.ClassicLauncher(cmdfile, workdir=pyl_workdir, debug='job', delay=0.01, cores=9)

    elif subcmd == 'resume-run':
        # called from SLURM — resumes from queuestate
        import pylauncher
        queuestate = sys.argv[2]
        pylauncher.ResumeClassicLauncher(queuestate, debug='job', delay=0.01, cores=9)

    elif subcmd == 'resume':
        # local submission of a resume SLURM job
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument('subcmd')
        p.add_argument('i0', type=int)
        p.add_argument('i1', type=int)
        p.add_argument('--nodes', type=int,   default=8)
        p.add_argument('--time',  type=float, default=4.0)
        p.add_argument('--queue', type=str,   default='normal')
        args = p.parse_args()
        resume_bispec_fiducial_pylauncher(
            args.i0, args.i1, nodes=args.nodes, time=args.time, queue=args.queue)

    else:
        # local submission of a new job: python bispec_pylauncher.py i0 i1 [opts]
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument('i0',         type=int)
        p.add_argument('i1',         type=int)
        p.add_argument('--nodes',     type=int,   default=8)
        p.add_argument('--time',      type=float, default=6.0)
        p.add_argument('--queue',     type=str,   default='normal')
        p.add_argument('--skip-done', action='store_true')
        args = p.parse_args()
        run_bispec_fiducial_pylauncher(
            args.i0, args.i1,
            nodes=args.nodes, time=args.time, queue=args.queue,
            skip_done=args.skip_done)
