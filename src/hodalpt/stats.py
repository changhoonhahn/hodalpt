'''

module for calculating summary statistics given DM particles, halos, or
galaxies in a box 


'''
from pyspectrum import pyspecturm as PyS


def Pk_periodic(xyz, w=None, Lbox=2600, Ngrid=360, rsd=False, Nmubin=10, fft='pyfftw', silent=True):  
    ''' calculate the powerspectrum given xyz positions 

    Parameters
    ----------
    xyz : 2d array 
        3xN array of object positions. 
    
    w : 1d array, optional 
        N-dim array of weights for objects (default: None) 

    Lbox : float, optional 
        box size in Mpc/h (default: 2600) 

    Ngrid : int, optional 
        FFT grid size (default:360)  

    rsd : bool or int, optional 
        If rsd == False, assumes there are no redshift-space distoritons and
        calculates the power spectrum monopole. 
        If rsd == 0, 1, 2, assumes there are redshift-space distortions and
        calculates power specturm multipoles along `rsd` axes. The value
        specifies the redshift-space LOS direction {0: 'x', 1: 'y', 2: 'z'}.
        (default: False) 

    Nmubin : int, optional
        number of Mu bins for (k, mu) binning. (default: 10)

    fft : string optional 
        fftw version to use. Options are 'pyfftw' and 'fortran'. (default: pyfftw) 

    silent : boolean, optional 
        if True nothing is printed. 


    Return
    ------
    spec : dictionary containing 
        * k : average k in k bin
        * p0k : monopole
        * p2k : quadrupole
        * p4k : hexadecapole
        * counts : number of modes in k bin
        * k_kmu : average k of (k,mu) bin
        * mu_kmu : average mu of (k,mu) bin
        * p_kmu : P(k,mu)
        * counts_kmu : number of modes in (k, mu) bin
    '''
    if isinstance(rsd, bool): 
        if rsd: 
            raise ValueError("If you want RSD, specify the LOS direction {0: 'x', 1: 'y', 2: 'z'}.")
        else: 
            # calculate power spectrum monopole using pyspectrum
            spec = PyS.Pk_periodic(xyz, w=w, Lbox=Lbox, Ngrid=Ngrid, fft=fft, silent=silent)
    else: 
        spec = PyS.Pk_periodic_rsd(xyz, w=w, Lbox=Lbox, Ngrid=Ngrid, rsd=rsd,
                            Nmubin=Nmubin, fft=fft, 
                            code='fortran', # by default uses the wrapped fortran code. 
                            silent=silent)
    return spec 


def B0_periodic(xyz, w=None, Lbox=2600, fft='pyfftw', silent=True):
    ''' calculate the bispectrum monopole for periodic box. 

    This function wraps `pyspectrum.Bk_periodic` with hardcoded choices 
    (Ngrid=360, step=3, Ncut=3, Nmax=40). Do not modify these choices as it
    will require recomputing the triangle counts which take forever. 

    :param xyz: 
        3xN dimensional array of the positions of objects (e.g. DM, halos,
        galaxies) 

    :param w: 
        N dimensional array of weights

    :param Lbox: (default: 2600.)
        box size 

    :param Ngrid: (default: Ngrid) 
        grid size 
    
    :return bispec: 
        dictionary containing the bispectrum. 
        bispec.keys() = 'meta', 'i_k1','i_k2','i_k3','p0k1','p0k2','p0k3','b123','q123','counts'
        bispec['meta'] is a dictionary with meta data regarding how the bispectrum is calculated
        and such.  i_k1, i_k2, i_k3 are the triangle side lenghts in units of fundmental mode
        p0k1, p0k2, p0k3 are corresponding powerspectrum values shot noise corrected. b123 is 
        the bispectrum also shot noise corrected and q123 is the reduced bispectrum (where both
        b123 and pk1pk2pk3 are shot noise corrected. 'counts' are the number of modes. 
    '''
    return PyS.Bk_periodic(xyz, w=w, Lbox=Lbox, Ngrid=360, step=3, Ncut=3, Nmax=40, fft=fft)
