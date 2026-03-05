#!/bin/python
''' python script to run 
'''
import os, sys
import numpy as np 
import matplotlib.pyplot as plt
from hodalpt.sims import alpt as CS

# directory where Quijote LHC is located /corral/utexas/AST25023/simbig/quijote/latinhypercube_hr
dir_lhc = sys.argv[1] 
i_lhc = int(sys.argv[2]) # LHC realization
silent = (sys.argv[3] == 'True') 

# directories 
ic_path = os.path.join(dir_lhc, str(i_lhc), 'ICs')
outdir = os.path.join(dir_lhc, str(i_lhc), 'alpt')

print('ALPT running for LHC%i' % i_lhc) 
if not os.path.isfile(os.path.join(outdir, 'Quijote_ICs_delta_z127_n256_CIC.DAT')): 
    CS.CSbox_alpt(ic_path, outdir, seed=0, dgrowth_short=5., make_ics=True,
            subgrid=True, return_pos=False, silent=silent)
else: 
    # run subgrid on 
    CS.CSbox_alpt(ic_path, outdir, seed=0, dgrowth_short=5., make_ics=False,
            subgrid=True, return_pos=False, silent=silent)
