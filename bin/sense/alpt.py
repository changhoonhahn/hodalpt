#!/bin/python
''' python script to run ALPT mocks from the sobol sequence 
'''
import os, sys
import numpy as np 
import matplotlib.pyplot as plt
from hodalpt.sims import alpt as CS

# directory where ALPT mocks will be located 
dir_alpt = sys.argv[1] 
i_sobol = int(sys.argv[2]) # i_sobol realization
silent = (sys.argv[3] == 'True') 

# read cosmology 
cosmos = np.loadtxt('cosmology_alpt_sobol2048.dat', skiprows=1) 

cosmo = {}
cosmo['Omega_m']    = cosmos[i_sobol][0] 
cosmo['Omega_b']    = cosmos[i_sobol][1] 
cosmo['h']          = cosmos[i_sobol][2]
cosmo['n_s']        = cosmos[i_sobol][3]
cosmo['sigma_8']    = cosmos[i_sobol][4]

# directories 
outdir = os.path.join(dir_alpt, str(i_sobol), 'alpt')

print('ALPT running for %i' % i_sobol) 
CS.CSbox_alpt(cosmo, outdir, seed=0, dgrowth_short=5., 
              Nmesh_ic=512, Nsample_ic=256, 
              make_ics=True, 
              subgrid=True,
              return_pos=False, 
              silent=silent)
