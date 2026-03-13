'''

module to construct CosmicSignal ALPT galaxy mocks 


author(s):
    * Francesco Sinigaglia 
    * ChangHoon Hahn: minor modifications 

'''
import numpy as np 
import subprocess
import os, sys, glob 

import camb

from . import cwc as C
from . import util as U 
from . import quijote as Q


def CSbox_galaxy(theta_gal, theta_rsd, dm_dir, Ngrid=256, Lbox=1000.,
                 zsnap=0.5, lambdath_tweb=0.0, lambdath_twebdelta=0.0,
                 bias_model='local', subgrid=False, seed=123456, silent=True): 
    ''' construct CosmicSignal galaxy mock given DM box. Applies the bias model
    in hodalpt.sims.cwc to specified ALPT DM output 

    .. code-block:: python
        
       from hodalpt.sims import alpt as Alpt
       theta_gal = {'alpha': 1.9230, 'beta': 2.0253, 'dth': -0.7889, 'rhoeps':
                    14.6874, 'eps': 0.5616, 'nmean': 3.3e-4}
       theta_rsd = {'bv': 0.7289, 'bb': 1.1652, 'betarsd': 1.3136, 'gamma':
                    0.4944}

       xyz = Alpt.CSbox_galaxy(theta_gal, theta_rsd, '/Users/hahnchanghoon/data/simbig/quijote/fiducial/0/alpt/', silent=False)


    paramters
    ---------
    theta_gal : dict
        Dictionary specify the galaxy bias parameters: 'alpha', 'beta', 'dth',
        'rhoeps', 'eps'. For some default values use: `theta_gal = { 'alpha':
        1.9230, 'beta': 2.0253, 'dth': -0.7889, 'rhoeps': 14.6874, 'eps':
        0.5616, 'nmean': 3.3e-4}`
        
    theta_rsd : dict
        Dictionary specifying the RSD parameters: bv, bb, betarsd, gamma. For
        some default values use `theta_rsd = {'bv': 0.7289, 'bb': 1.1652,
        'betarsd': 1.3136, 'gamma': 0.4944}`

    dm_dir : str
        Directory with the ALPT DM files 

    bias_model : str
        specify which bias model to use 


    return 
    ------
    xyz : array 
        Ngal x 3 array specifying the x, y, z position of galaxies.
    '''
    if bias_model not in ['local', 'nonlocal0']: 
        raise NotImplementedError('%s bias model not implemented yet' % bias_model) 

    np.random.seed(seed)

    assert os.path.isdir(dm_dir), "specify correct directory for the DM files"

    def Fname(prefix): return os.path.join(dm_dir, '%s.dat' % prefix)
    
    prefix_subgrid = ''
    if subgrid: prefix_subgrid = 'super_'
    dm_filename         = Fname(prefix_subgrid+'deltaBOX')
    tweb_filename       = Fname('Tweb_')
    twebdelta_filename  = Fname('TwebDelta_')

    vx_filename = dm_dir + 'VExEULz%3.3f.dat' % zsnap
    vy_filename = dm_dir + 'VEyEULz%3.3f.dat' % zsnap
    vz_filename = dm_dir + 'VEzEULz%3.3f.dat' % zsnap
   
    if not silent: print('reading %s' % Fname(prefix_subgrid+'BOXpos*'))
    posx_filename = Fname(prefix_subgrid+'BOXposx')
    posy_filename =	Fname(prefix_subgrid+'BOXposy')
    posz_filename = Fname(prefix_subgrid+'BOXposz')
    
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
    nmean   = theta_gal['nmean'] # yes, this is weird but lets not overthink it for now 
    if bias_model  == 'local': 
        dth     = theta_gal['dth'] 
        rhoeps  = theta_gal['rhoeps']
        eps     = theta_gal['eps']     
        rhoepsprime = theta_gal['rhoepsprime'] 
        epsprime    = theta_gal['epsprime']

    # parse rsd parameters
    bv      = theta_rsd['bv'] 
    bb      = theta_rsd['bb']
    betarsd = theta_rsd['betarsd']
    gamma   = theta_rsd['gamma'] 
    
    lcell = Lbox/Ngrid
    
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

    # impose boundary conditions 
    posx = (posx + Lbox) % Lbox
    posy = (posy + Lbox) % Lbox
    posz = (posz + Lbox) % Lbox

    # Now they are velocity vectors
    vx = np.fromfile(vx_filename, dtype=np.float32)  
    vy = np.fromfile(vy_filename, dtype=np.float32) 
    vz = np.fromfile(vz_filename, dtype=np.float32) 

    # Reshape arrays from 1D to 3D --> reshape only arrays which have mesh structure, e.g. NOT positions
    #delta = np.reshape(delta, (Ngrid,Ngrid,Ngrid))
    tweb        = np.reshape(tweb, (Ngrid,Ngrid,Ngrid))
    twebdelta   = np.reshape(twebdelta, (Ngrid,Ngrid,Ngrid))

    vx = np.reshape(vx, (Ngrid,Ngrid,Ngrid))
    vy = np.reshape(vy, (Ngrid,Ngrid,Ngrid))
    vz = np.reshape(vz, (Ngrid,Ngrid,Ngrid))

    # Apply the bias and get halo/galaxy number counts
    if not silent: print('Getting number counts via parametric bias ...')
    if bias_model == 'local': 
        ncounts = C.biasmodel_local_box(Ngrid, Lbox, delta,  nmean, alpha, beta, dth, rhoeps, eps, 
                                        rhoepsprime, epsprime)
    elif bias_model == 'nonlocal0': 
        ncounts = C.biasmodel_nonlocal_box(Ngrid, Lbox, delta, tweb, twebdelta, 
                                  nmean, alpha, beta)#, dth, rhoeps, eps)
    else: 
        raise NotImplementedError('%s bias model not implemented yet' % bias_model) 
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
                                                          bv, bb, betarsd, gamma,
                                                          zsnap, omega_m) 
    
    return np.vstack([posx, posy, posz]).T


def CSbox_alpt(cosmo, outdir, seed=0, dgrowth_short=5., ngrid=256,
               lbox=1000., zsnap=0.5, lambdath_tweb=0., lambdath_twebdelta=0., 
               Nmesh_ic=512, Nsample_ic=256, subgrid=True, return_pos=True,
               silent=True): 
    ''' run CosmicSignal ALPT code for given LCDM cosmological parameters
    '''
    zsnap = [zsnap] 
    assert len(zsnap) == 1
    zmin = zsnap[0] # hardcoded
    zmax = zsnap[0] # hardcoded

    # cosmological parameters
    Omega_m = cosmo['Omega_m'] 
    Omega_b = cosmo['Omega_b']
    n_s     = cosmo['n_s']
    sigma8  = cosmo['sigma_8']
    h       = cosmo['h']
    w0      = -1 # hardcoded to LCDM for now 
    wa      = 0 # hardcoded to LCDM for now 

    ic_paramfile = '2LPT.param' # hardcoded IC filename in Quijote 

    # webonx executable should be in same directory as this script
    scriptdir = os.path.dirname(__file__)
    webonx = os.path.join(scriptdir, 'webonx')
    assert os.path.isfile(webonx), "webonx executable not found"

    os.makedirs(outdir, exist_ok=True) # make output directory in case it doesn't exist 

    # Write input redshift snapshots file
    _write_z_input_file(zsnap, outdir)

    # Write cosmological parameters file 
    _write_cosmology_par_input_file(Omega_m, Omega_b, w0, n_s, wa, sigma8, h, outdir)

    # write input parameter file 
    _write_input_par_file(ngrid, lbox, seed, 1, lambdath_tweb,
                          lambdath_twebdelta, Omega_m, dgrowth_short, zmin, zmax, outdir)
    
    os.chdir(outdir)

    if not silent: print(f'Generating ICs and writing out delta IC')
    # grid ICs and get delta
    delta = generate_IC(Omega_m, Omega_b, h, n_s, sigma8, lbox, ngrid, seed,
                        ic_dir=outdir, Nmesh=Nmesh_ic, Nsample=Nsample_ic,
                        silent=silent)

    # write delta to outdir 
    delta.astype('float32').tofile(os.path.join(outdir, 'Quijote_ICs_delta_z127_n256_CIC.DAT'))


    # Compute displacement fields at different redshifts
    if not silent: print(f'Computing displacement fields at z=%s' % zsnap[0])
    sys.stdout.flush()
    if not silent: subprocess.run([webonx,])
    else: subprocess.run([webonx,], stdout=subprocess.DEVNULL)
    sys.stdout.flush()

    prefix_subgrid = ''
    if subgrid: 
        prefix_subgrid = 'super_'

        # write input parameter file for subgrid model 
        _write_input_par_file(ngrid, lbox, seed, -70, lambdath_tweb,
                              lambdath_twebdelta, Omega_m, dgrowth_short, zmin, zmax, outdir)
        # run subgrid model 
        if not silent: print(f'Running subgrid model')
        sys.stdout.flush()
        if not silent: subprocess.run([webonx,])
        else: subprocess.run([webonx,], stdout=subprocess.DEVNULL)
        sys.stdout.flush()

    if not return_pos: 
        return None 
    else: 
        posx = np.fromfile(os.path.join(outdir, prefix_subgrid+'BOXposx.dat'), dtype=np.float32)
        posy = np.fromfile(os.path.join(outdir, prefix_subgrid+'BOXposy.dat'), dtype=np.float32)
        posz = np.fromfile(os.path.join(outdir, prefix_subgrid+'BOXposz.dat'), dtype=np.float32)

        # impose boundary conditions 
        posx = (posx + lbox) % lbox
        posy = (posy + lbox) % lbox
        posz = (posz + lbox) % lbox
        
        return np.vstack([posx, posy, posz]).T


def CSbox_alpt_Q(ic_path, outdir, seed=0, dgrowth_short=5., ngrid=256,
               lbox=1000., zsnap=0.5, lambdath_tweb=0., lambdath_twebdelta=0., 
               make_ics=True, subgrid=True, return_pos=True, silent=True): 
    ''' run CosmicSignal ALPT for Quijote initial condiitons 
    '''
    zsnap = [zsnap] 
    assert len(zsnap) == 1
    zmin = zsnap[0] # hardcoded
    zmax = zsnap[0] # hardcoded

    ic_paramfile = '2LPT.param' # hardcoded IC filename in Quijote 

    # webonx executable should be in same directory as this script
    scriptdir = os.path.dirname(__file__)
    webonx = os.path.join(scriptdir, 'webonx')
    assert os.path.isfile(webonx), "webonx executable not found"

    os.makedirs(outdir, exist_ok=True) # make output directory in case it doesn't exist 

    # Set cosmological parameters for this run from Quijote IC  
    omega_m, omega_b, w0, n_s, wa, sigma8, hh = _read_cosmo_pars_from_config(ic_path, ic_paramfile)

    # Write input redshift snapshots file
    _write_z_input_file(zsnap, outdir)

    # Write cosmological parameters file 
    _write_cosmology_par_input_file(omega_m, omega_b, w0, n_s, wa, sigma8, hh, outdir)

    # write input parameter file 
    _write_input_par_file(ngrid, lbox, seed, 1, lambdath_tweb,
                          lambdath_twebdelta, omega_m, dgrowth_short, zmin, zmax, outdir)
    
    if make_ics: 
        if not silent: print(f'Computing and writing out delta IC')
        # grid ICs and get delta
        delta = _make_ics_quijote(ic_path, float(lbox), ngrid)
        # write delta to outdir 
        delta.astype('float32').tofile(os.path.join(outdir, 'Quijote_ICs_delta_z127_n256_CIC.DAT'))

    os.chdir(outdir)

    # Compute displacement fields at different redshifts
    if not silent: print(f'Computing displacement fields at z=%s' % zsnap[0])
    sys.stdout.flush()
    if not silent: subprocess.run([webonx,])
    else: subprocess.run([webonx,], stdout=subprocess.DEVNULL)
    sys.stdout.flush()

    prefix_subgrid = ''
    if subgrid: 
        prefix_subgrid = 'super_'

        # write input parameter file for subgrid model 
        _write_input_par_file(ngrid, lbox, seed, -70, lambdath_tweb,
                              lambdath_twebdelta, omega_m, dgrowth_short, zmin, zmax, outdir)
        # run subgrid model 
        if not silent: print(f'Running subgrid model')
        sys.stdout.flush()
        if not silent: subprocess.run([webonx,])
        else: subprocess.run([webonx,], stdout=subprocess.DEVNULL)
        sys.stdout.flush()

    if not return_pos: 
        return None 
    else: 
        posx = np.fromfile(os.path.join(outdir, prefix_subgrid+'BOXposx.dat'), dtype=np.float32)
        posy = np.fromfile(os.path.join(outdir, prefix_subgrid+'BOXposy.dat'), dtype=np.float32)
        posz = np.fromfile(os.path.join(outdir, prefix_subgrid+'BOXposz.dat'), dtype=np.float32)

        # impose boundary conditions 
        posx = (posx + lbox) % lbox
        posy = (posy + lbox) % lbox
        posz = (posz + lbox) % lbox
        
        return np.vstack([posx, posy, posz]).T


def _write_input_par_file(ngrid, lbox, seed, sfmodel, lambdath_tweb, lambdath_twebdelta, omegam, dgrowth_short, zmin, zmax, outdir):
    ff = open(os.path.join(outdir, 'input.par'), 'w')
    
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('# BOX: parameter file\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('seed = %d\n' %seed)
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('Nx = %d # Number of pixels x\n' %ngrid)
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('Lx = %s # Box length in x-direction [Mpc/h]\n' %str(lbox))
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('fnameIC = Quijote_ICs_delta_z127_n256_CIC.DAT# Attention no space after = unless you give a name\n')
    ff.write('fnameDM = deltaBOX.dat\n')
    ff.write('sfmodel = %s\n' %str(sfmodel))
    ff.write('filter = 1\n')
    ff.write('dgrowth_short = %s\n' % str(dgrowth_short))
    ff.write('rsml = 1.0\n')
    ff.write('dtol = 0.005\n')
    ff.write('curlfrac = 0\n')
    ff.write('write_box = true\n')
    ff.write('check_calibration = false\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('z = %s\n' % zmin)
    ff.write('zref = 127.\n')
    ff.write('slength = 5\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('mpro = 1\n')
    ff.write('nminperc = 20\n')
    ff.write('tcw = 4\n')
    ff.write('tcwD = 1234\n')
    ff.write('deltathL = -.5\n')
    ff.write('deltathH = -.5\n')
    ff.write('lth = %s # lambdath phi-web\n' %lambdath_tweb)
    ff.write('lthD = %s # lambdath delta-web\n' %lambdath_twebdelta)
    ff.write('fridge = 0.1\n')
    ff.write('boostk = 1.7\n')
    ff.write('boostf = 1.7\n')
    ff.write('boosts = 1.7\n')
    ff.write('alpha = 0.5\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('NsnapOMP = 1\n')
    ff.write('zsnapsmin = %s\n' % zmin)
    ff.write('zsnapsmax = %s\n' % zmax)
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('xllc = %s\n' %str(-np.float32(lbox)/2))
    ff.write('yllc = %s\n' %str(-np.float32(lbox)/2))
    ff.write('zllc = %s\n' %str(-np.float32(lbox)/2))

    ff.close()
    return None 


def _write_z_input_file(zsnap, outdir):
    ''' write z input file used by webonx 
    '''
    ff = open(os.path.join(outdir, 'z_input.par'), 'w')

    for ii in range(len(zsnap)):
        ff.write(str(zsnap[ii]) + '\n')
    ff.close()
    return None 


def _write_cosmology_par_input_file(omega_m, omega_b, w0, n_s, wa, sigma8, hpar, outdir):
    ''' write cosmology input file for webonx
    '''
    ff = open(os.path.join(outdir, 'cosmology.par'), 'w')

    ff.write('#---------------------------------------------------------------------\n')
    ff.write('# CosmicSignal: cosmology parameter file\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('omega_m = %s\n' %str(omega_m))
    ff.write('omega_b = %s\n' %str(omega_b))
    ff.write('wpar    = %s\n' %str(w0))
    ff.write('n_s     = %s\n' %str(n_s))
    ff.write('wprime  = %s\n' %str(wa))
    ff.write('sigma8  = %s\n' %str(sigma8))
    ff.write('rsmooth = 8.0\n')
    ff.write('hpar    = %s\n' %str(hpar))
    ff.write('betapar = 1.5\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('readPS = true\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('fnamePS = linear_power_z0.DAT\n')
    ff.write('#---------------------------------------------------------------------\n')
    ff.write('#---------------------------------------------------------------------\n')

    ff.close()
    return None 


def _read_cosmo_pars_from_config(ic_path, ic_paramfile):
    ''' read cosmological parameters from config file 
    '''
    fn = os.path.join(ic_path, ic_paramfile)

    raw = np.genfromtxt(fn, comments='%')

    omega_m = raw[7,1]
    omega_b = raw[9,1]
    hh = raw[11,1]
    #zz = raw[12,1]
    sigma8 = raw[13,1]
    n_s = raw[19,1]
    w0 = -1
    wa = 0.

    return omega_m, omega_b, w0, n_s, wa, sigma8, hh


def _make_ics_quijote(ic_path, lbox, ngrid):
    ''' grid quijote initial conditions 
    '''
    ics = Q.IC(ic_path.split('ICs')[0])
    posx = ics.pos[:,0] 
    posy = ics.pos[:,1] 
    posz = ics.pos[:,2] 
    
    weight = np.ones(len(posx))                                                                                                         

    delta = U.get_cic(posx, posy, posz, weight, float(lbox), ngrid)
    delta = delta.flatten()
    delta = delta/np.mean(delta)-1.

    return delta 


def generate_IC(Omega_m, Omega_b, h, ns, s8, lbox, ngrid, seed, ic_dir=None, Nmesh=512,
                Nsample=256, silent=True): 
    ''' generate ICs using 2LPTicQ in a way consistent with Quijote simulations
    '''
    # get CAMB matter power spectrum 
    k, zs, Pkmm = camb_Pkmm0(Omega_m, Omega_b, h, ns, s8)

    # write power spectrum to file 
    for j,z in enumerate(zs):
        np.savetxt(os.path.join(ic_dir, 'Pk_mm_z=%.3f.txt' % z), np.transpose([k,Pkmm[j,:]]))

    # write input file for 2LPTic
    _write_2lptic_param(Omega_m, Omega_b, h, ns, s8, seed, 
                        Lbox=lbox, 
                        Nmesh=Nmesh,
                        Nsample=Nsample, 
                        ic_dir=ic_dir)
    f2lptparam = os.path.join(ic_dir, '2LPT.param') 
    
    # run 2LPT 
    scriptdir = os.path.dirname(__file__)
    twolpt = os.path.join(scriptdir, '2LPTic')
    assert os.path.isfile(twolpt), "2LPTic executable not found"
    
    if not silent: print(f'Computing 2LPT IC')
    sys.stdout.flush()
    if not silent: subprocess.run([twolpt, f2lptparam])
    else: subprocess.run([twolpt, f2lptparam], stdout=subprocess.DEVNULL)
    sys.stdout.flush()

    # read ICs
    if not silent: print(f'reading in IC')
    _, pos, vel, _ = U.read_gadget_ic(os.path.join(ic_dir, 'ics'))
    pos = pos/1000.

    weight = np.ones(pos.shape[0])                                                                                                         

    delta = U.get_cic(pos[:,0], pos[:,1], pos[:,2], weight, float(lbox), ngrid)
    delta = delta.flatten()
    delta = delta/np.mean(delta)-1.
    return delta


def _write_2lptic_param(Omega_m, Omega_b, h, ns, s8, seed, Lbox=1000.,
                        Nmesh=512, Nsample=256, ic_dir=None): 
    ''' write 2LPT.param file for 2LPTicQ 
    '''
    if (Nsample % 64) != 0: raise ValueError('should be divisible by 64') 
    fglass = os.path.join(os.path.dirname(__file__), 'GLASS',
                          'dummy_glass_dmonly_64.dat')
    fpk = os.path.join(ic_dir, 'Pk_mm_z=0.000.txt')
    
    a = """
    Nmesh            %i      
    Nsample          %i       
    Box              %f
    FileBase         ics         
    OutputDir        %s          
    GlassFile        %s
    GlassTileFac     %i         
    Omega            %.4f    
    OmegaLambda      %.4f    
    OmegaBaryon      0.0000    
    OmegaDM_2ndSpecies  0.0    
    HubbleParam      %.4f    
    Redshift         127       
    Sigma8           %.4f       
    SphereMode       0         
    WhichSpectrum    2         
    FileWithInputSpectrum  %s
    InputSpectrum_UnitLength_in_cm  3.085678e24 
    ShapeGamma       0.201     
    PrimordialIndex  1.0       
    
    Phase_flip          0      
    RayleighSampling    1      
    Seed                %d      
    
    NumFilesWrittenInParallel 1  
    UnitLength_in_cm          3.085678e21  
    UnitMass_in_g             1.989e43     
    UnitVelocity_in_cm_per_s  1e5          
    
    WDM_On               0      
    WDM_Vtherm_On        0      
    WDM_PartMass_in_kev  10.0   
    """ % (Nmesh, Nsample, Lbox*1000., ic_dir+'/', 
           fglass, Nsample/64, 
           Omega_m, 1.0-Omega_m, h, s8, fpk, seed) 

    f = open('%s/2LPT.param' % ic_dir, 'w')
    f.write(a)
    f.close() 
    return None 


def camb_Pkmm0(Omega_m, Omega_b, h, ns, s8):
    ''' get z=0 matter power spectrum given LCDM cosmological parameters using
    CAMB. This is based on Paco's script:
    https://github.com/franciscovillaescusa/Quijote-simulations/blob/master/latin_hypercube/parameters_file.py 
    '''
    # hardcoded CAMB parameters
    hierarchy    = 'degenerate'
    Mnu          = 0.0 #eV
    Nnu          = 0   #number of massive neutrinos
    Neff         = 3.046
    As           = 2.13e-9
    tau          = None
    Omega_k      = 0.0
    pivot_scalar = 0.05
    pivot_tensor = 0.05
    kmax         = 10.0
    k_per_logint = 20
    redshifts    = [0]
    
    # run CAMB
    Omega_c  = Omega_m - Omega_b
    pars     = camb.CAMBparams()

    # set accuracy of the calculation
    # HighAccuracyDefault is not defunct
    pars.set_accuracy(AccuracyBoost=5.0, lSampleBoost=5.0, lAccuracyBoost=5.0, #HighAccuracyDefault=True,
                      DoLateRadTruncation=True)

    # set value of the cosmological parameters
    pars.set_cosmology(H0=h*100.0, ombh2=Omega_b*h**2, omch2=Omega_c*h**2, 
                       mnu=Mnu, omk=Omega_k, neutrino_hierarchy=hierarchy, 
                       num_massive_neutrinos=Nnu, nnu=Neff, tau=tau)
                   
    # set the value of the primordial power spectrum parameters
    pars.InitPower.set_params(As=As, ns=ns, 
                              pivot_scalar=pivot_scalar, pivot_tensor=pivot_tensor)

    # set redshifts, k-range and k-sampling
    pars.set_matter_power(redshifts=redshifts, kmax=kmax, k_per_logint=k_per_logint)

    # compute results
    results = camb.get_results(pars)

    # interpolate to get Pmm, Pcc...etc
    k, zs, Pkmm = results.get_matter_power_spectrum(minkh=2e-5, maxkh=kmax, 
                                                    npoints=400, var1=7, var2=7, 
                                                    have_power_spectra=True, 
                                                    params=None)
    return k, zs, Pkmm 
