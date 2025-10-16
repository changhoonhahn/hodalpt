'''

module to construct CosmicSignal ALPT galaxy mocks 


author(s):
    * Francesco Sinigaglia 
    * ChangHoon Hahn: minor modifications 

'''
import os, glob 
import numpy as np 
import configparser
from . import cwc as C


def CSbox_galaxy(theta_gal, theta_rsd, dm_dir, Ngrid=256, Lbox=1000.,
                 zsnap=0.5, lambdath_tweb=0.0, lambdath_twebdelta = 0.0,
                 seed=123456, silent=True): 
    ''' construct CosmicSignal galaxy mock given DM box. Applies the bias model
    in hodalpt.sims.cwc to specified ALPT DM output 


    paramters
    ---------
    dm_dir 


    return 
    ------


    '''
    np.random.seed(seed)

    assert os.path.isdir(dm_dir), "specify correct directory for the DM files"

    # find file suffix and check for consistency 
    suffix = 'OM'+(glob.glob(os.path.join(dm_dir, 'deltaBOXOM*'))[0].split('deltaBOXOM')[-1]).split('.gz')[0]
    omega_m = float(suffix.split('OM')[1].split('OL')[0])
    if not silent: print(suffix)
    if not silent: print('Omega_m %f' % omega_m) 
    def Fname(prefix): return os.path.join(dm_dir, '%s%s' % (prefix, suffix))

    dm_filename         = Fname('deltaBOX')
    tweb_filename       = Fname('Tweb_')
    twebdelta_filename  = Fname('TwebDelta_')

    vx_filename = dm_dir + 'VExEULz%3.3f.dat' % zsnap
    vy_filename = dm_dir + 'VExEULz%3.3f.dat' % zsnap
    vz_filename = dm_dir + 'VExEULz%3.3f.dat' % zsnap

    posx_filename = Fname('BOXposx')
    posy_filename =	Fname('BOXposy')
    posz_filename = Fname('BOXposz')
    
    assert os.path.isfile(dm_filename), 'missing %s' % dm_filename
    assert os.path.isfile(tweb_filename), 'missing %s' % tweb_filename
    assert os.path.isfile(twebdelta_filename), 'missing %s' % twebdelta_filename
    assert os.path.isfile(vx_filename), 'missing %s' % vx_filename
    assert os.path.isfile(vy_filename), 'missing %s' % vy_filename 
    assert os.path.isfile(vz_filename), 'missing %s' % vz_filename 
    assert os.path.isfile(posx_filename), 'missing %s' % posx_filename
    assert os.path.isfile(posy_filename), 'missing %s' % posy_filename
    assert os.path.isfile(posz_filename), 'missing %s' % posz_filename

    # parse galaxy bias parameters 
    alpha   = theta_gal['alpha']  
    beta    = theta_gal['beta']
    dth     = theta_gal['dth'] 
    rhoeps  = theta_gal['rhoeps']
    eps     = theta_gal['eps']     
    rhoepsprime = 0.
    epsprime    = 0.
    nmean   = theta_gal['nmean'] # yes, this is weird but lets not overthink it for now 


    # parse rsd parameters
    bv      = theta_rsd['bv'] 
    bb      = theta_rsd['bb']
    betarsd = theta_rsd['betarsd']
    gamma   = theta_rsd['gamma'] 
    
    # Observer positions            
    obspos = [Lbox/2., Lbox/2., Lbox/2.]

    lcell = Lbox/Ngrid

    xobs = obspos[0]
    yobs = obspos[1]
    zobs = obspos[2]
    
    if not silent: print('Reading input ...')
    # read inputs 
    delta = np.fromfile(dm_filename, dtype=np.float32)  # In real space
    delta = np.reshape(delta, (Ngrid,Ngrid,Ngrid))

    tweb = np.fromfile(tweb_filename, dtype=np.float32)  # In real space
    twebdelta = np.fromfile(twebdelta_filename, dtype=np.float32)  # In real space 

    # Positions
    posx = np.fromfile(posx_filename, dtype=np.float32)   
    posy = np.fromfile(posy_filename, dtype=np.float32) 
    posz = np.fromfile(posz_filename, dtype=np.float32)

    # Now they are velocity vectors
    vx = np.fromfile(vx_filename, dtype=np.float32)  
    vy = np.fromfile(vy_filename, dtype=np.float32) 
    vz = np.fromfile(vz_filename, dtype=np.float32) 

    # Reshape arrays from 1D to 3D --> reshape only arrays which have mesh structure, e.g. NOT positions
    #delta = np.reshape(delta, (Ngrid,Ngrid,Ngrid))
    tweb = np.reshape(tweb, (Ngrid,Ngrid,Ngrid))
    twebdelta = np.reshape(twebdelta, (Ngrid,Ngrid,Ngrid))

    vx = np.reshape(vx, (Ngrid,Ngrid,Ngrid))
    vy = np.reshape(vy, (Ngrid,Ngrid,Ngrid))
    vz = np.reshape(vz, (Ngrid,Ngrid,Ngrid))

    # Apply the bias and get halo/galaxy number counts
    if not silent: print('Getting number counts via parametric bias ...')
    ncounts = C.biasmodel_local_box(Ngrid, Lbox, delta,  nmean, alpha, beta, dth, rhoeps, eps, rhoepsprime, epsprime, xobs, yobs, zobs)
    ncountstot = np.sum(ncounts) # total number of objects
    if not silent: print('Number counts diagnostics (min, max, mean): ', np.amin(ncounts), np.amax(ncounts), np.mean(ncounts))

    ncounts = np.reshape(ncounts, (Ngrid, Ngrid, Ngrid))

    # Now assign positions
    if not silent: print('Preparing galaxy positions ...')
    posxarr_prep, posyarr_prep, poszarr_prep = C.prepare_indices_array(posx, posy, posz, Ngrid, Lbox)
    if not silent: print('Sampling galaxy positions ...')
    posx, posy, posz = C.sample_galaxies(Lbox, Ngrid, posxarr_prep, posyarr_prep, poszarr_prep, ncounts)


    if not silent: print('apply RSD ...')
    posx, posy, posz = C.real_to_redshift_space_local_box(delta, tweb, posx,
                                                          posy, posz, vx, vy,
                                                          vz, Ngrid, Lbox,
                                                          xobs, yobs, zobs, bv,
                                                          bb, betarsd, gamma,
                                                          zsnap, omega_m) 
    
    return np.vstack([posx, posy, posz]).T
