import numpy as np
import h5py
from hodalpt.sims import sampling as S


outfile_NLB = 'NLB_samples_10000.hdf5'
outfile_HOD = 'HOD_samples_1000.hdf5'

S.generate_and_save_NLB(outfile_NLB, 10000)
S.generate_and_save_HOD(outfile_HOD, 1000)
