import numpy as np 


class Snap(object): 
    ''' class object for particle snapshots 
    '''
    def __init__(self): 
        self.BoxSize  = None #Mpc/h
        self.Nall     = None #Total number of particles
        self.Masses   = None #Masses of the particles in Msun/h
        self.Omega_m  = None #value of Omega_m
        self.Omega_l  = None #value of Omega_l
        self.h        = None #value of h
        self.redshift = None #redshift of the snapshot
        self.Hubble   = None #Value of H(z) in km/s/(Mpc/h)

        self.pos      = None #positions in Mpc/h
        self.vel      = None #peculiar velocities in km/s

    def _read_quijote_header(self, header): 
        ''' given Quijote simulation header, 
        '''
        self.BoxSize  = header['boxsize']/1e3  #Mpc/h
        self.Nall     = header['nall']         #Total number of particles
        self.Masses   = header['massarr']*1e10 #Masses of the particles in Msun/h
        self.Omega_m  = header['omega_m']      #value of Omega_m
        self.Omega_l  = header['omega_l']      #value of Omega_l
        self.h        = header['hubble']       #value of h
        self.redshift = header['redshift']     #redshift of the snapshot
        self.Hubble   = 100.0*np.sqrt(self.Omega_m*(1.0+self.redshift)**3+self.Omega_l)#Value of H(z) in km/s/(Mpc/h)
        return None 

    def save_cs(self):
        ''' save snapshot to binary format that can be ready by CosmicSignals 
        '''
        return None 
