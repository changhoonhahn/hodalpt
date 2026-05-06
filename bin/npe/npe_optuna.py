import os, sys
import numpy as np
import h5py
import torch
from torch import nn
# from torch.utils.tensorboard.writer import SummaryWriter
import optuna
# from utils.py (in hodalpt/bin/npe)
from utils import training_data, train_val_split, get_prior, get_prior_bounds
from sbi import utils as Ut
from sbi.neural_nets import posterior_nn
from sbi import inference as Inference

import matplotlib.pyplot as plt

# change for TACC
#output_dir = '/Users/mcc3842/CosmicSim2025/hodalpt/bin/npe' # local
#output_dir = '/corral/utexas/AST25023/simbig/npe'
output_dir = '/work/11053/mcasas/ls6/hodalpt/bin/npe' # for now, writing to work
#infn = '/Users/mcc3842/CosmicSim2025/hodalpt/bin/sense/alpt_dataset.hdf5'
infn = '/work/11053/mcasas/ls6/hodalpt/bin/sense/alpt_dataset.hdf5'

sumstat     = sys.argv[1]
kmax        = float(sys.argv[2]) 
nf_model    = 'maf'
cosmo_only  = False
n_trials    = int(sys.argv[3])

##################################################################################
cuda = torch.cuda.is_available()
device = ("cuda:0" if cuda else "cpu")

seed = 12387
torch.manual_seed(seed)
if cuda:
    torch.cuda.manual_seed(seed)
##################################################################################
# load/split data 
##################################################################################

if int(kmax * 100) % 10 == 0: 
    kmax_str = 'kmax%.1f' % kmax
else:  
    kmax_str = 'kmax%.2f' % kmax


test_name = 'test_data.exp1.%s.%s.%s.npz' % (
    sumstat,                                    # p0, p2, pk_all
    kmax_str,                                   # kmax0.2, kmax0.15, etc
    ['full', 'cosmo_only'][cosmo_only]          # full (all 57 params) or cosmo only
)
theta, x = training_data(kmax, infn, sumstat)

test_data_path = os.path.join(output_dir, test_name)
x_train, y_train = train_val_split(theta, x, seed, n_test=100, save_path=test_data_path)
#x_train = torch.log(x_train)
print('reserving 100 sims for testing at '+str(test_data_path))
##################################################################################
# set prior 
##################################################################################
# om, ob, h, ns, s8, alpha, beta, nmean, rsd
prior = get_prior()
##################################################################################
# OPTUNA
##################################################################################
# Optuna Parameters

cosmo_only = False
study_name = 'qphi.exp1.%s.%s.%s.%s' % (
    nf_model,                                   # maf
    sumstat,                                    # p0, p2, pk_all
    kmax_str,                                   # kmax0.2, kmax0.15, etc
    ['full', 'cosmo_only'][cosmo_only]          # full (all 57 params) or cosmo only
)
n_jobs     = 1
if not os.path.isdir(os.path.join(output_dir, study_name)): 
    os.system('mkdir %s' % os.path.join(output_dir, study_name))
storage    = 'sqlite:///%s/%s/%s.db' % (output_dir, study_name, study_name)
n_startup_trials = 20

n_blocks_min, n_blocks_max = 2, 10
n_transf_min, n_transf_max = 2, 10
n_hidden_min, n_hidden_max = 64, 512
n_lr_min,     n_lr_max     = 5e-6, 1e-3
# remove p_drop pars for now ?
# p_drop_min, p_drop_max = 0., 1. 



def Objective(trial):
    ''' objective function for optuna 
    '''
    # Generate the model                                         
    n_blocks = trial.suggest_int("n_blocks", n_blocks_min, n_blocks_max)
    n_transf = trial.suggest_int("n_transf", n_transf_min,  n_transf_max)
    n_hidden = trial.suggest_int("n_hidden", n_hidden_min, n_hidden_max, log=True)
    lr = trial.suggest_float("lr", n_lr_min, n_lr_max, log=True) 
    p_drop = trial.suggest_float("p_drop", 0., 0.3)

    mlp = nn.Identity()
    
    neural_posterior = posterior_nn(nf_model, 
            embedding_net=mlp, 
            hidden_features=n_hidden, 
            num_transforms=n_transf, 
            num_blocks=n_blocks, 
            dropout_probability=p_drop, 
            use_batch_norm=False
            )

    anpe = Inference.SNPE(prior=prior,
            density_estimator=neural_posterior,
            device=device)#, 
            #summary_writer=SummaryWriter('%s/%s/%s.%i' % 
               # (output_dir, study_name, study_name, trial.number)))

    anpe.append_simulations(y_train.to(device), x_train.to(device))

    p_theta_x_est = anpe.train(
            validation_fraction=0.08,
            training_batch_size=200,
            learning_rate=lr, #clip_max_norm=clip_max,
            stop_after_epochs=20, 
            show_train_summary=True)

    
    # save trained NPE  
    qphi    = anpe.build_posterior(p_theta_x_est)
    fqphi   = os.path.join(output_dir, study_name, '%s.%i.pt' % (study_name, trial.number))
    torch.save(qphi, fqphi)
        
    best_validation_loss = anpe._summary['best_validation_loss'][0]

    # anpe._summary_writer.add_hparams(
    #         {'n_blocks': n_blocks, 'n_transf': n_transf, 'n_hidden': n_hidden, 'lr': lr, 'p_drop': p_drop},
    #         {'best_validation_loss': best_validation_loss}
    #         )
  

    return best_validation_loss


sampler     = optuna.samplers.TPESampler(n_startup_trials=n_startup_trials) # multivariate=True)
study       = optuna.create_study(study_name=study_name, sampler=sampler, storage=storage, directions=["minimize"], load_if_exists=True) # , "minimize", "minimize"

study.optimize(Objective, n_trials=n_trials, n_jobs=n_jobs)
print("  Number of finished trials: %i" % len(study.trials))
