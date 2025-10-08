'''


module for reading in Quijote data 


'''
import os
import h5py 
import hdf5plugin
import numpy as np 
# --- emanu --- 
from . import sims 


quijote_zsnap_dict = {0.: 4, 0.5: 3, 1.:2, 2.: 1, 3.: 0}


def Nbody(_dir, z=0.5): 
    ''' read CDM particles from Quijote n-body simulation output snapshot 

    Parameters
    ----------
    _dir : string
        directory that contains all the snapshots 

    z : float
        redshift of the snapshot you want read 
    '''

    # redshift snapshot 
    assert z in quijote_zsnap_dict.keys(), 'snapshots are available at z=0, 0.5, 1, 2, 3'
    snapnum = quijote_zsnap_dict[z]

    snapshot = os.path.join(_dir, 'snapdir_%s' % str(snapnum).zfill(3), 'snap_%s' % str(snapnum).zfill(3))
    return _read_snap(snapshot)


def IC(_dir): 
    ''' read initial conditions of CDM particles from Quijote. The IC is
    generated at z=127. 

    Parameters
    ----------
    _dir : string
        directory that contains the snapshots and IC files 
    '''
    return _read_snap(os.path.join(_dir, 'ICs', 'ICs'))


def _read_snap(snapshot):
    ''' read snapshot in hdf5 format. This is streamlined and modified version
    of readgadget from pylians. 
    '''
    ptype    = [1] #[1](CDM), [2](neutrinos) or [1,2](CDM+neutrinos)

    # read header
    header = _read_header(snapshot)
    snap = sims.Snap()
    snap._read_quijote_header(header)

    # read positions, velocities and IDs of the particles
    pos = _read_block(snapshot, "POS ", ptype)/1e3 #positions in Mpc/h
    vel = _read_block(snapshot, "VEL ", ptype)     #peculiar velocities in km/s
    #ids = readgadget.read_block(snapshot, "ID  ", ptype)-1   #IDs starting from 0

    snap.pos = pos 
    snap.vel = vel
    return snap 


def _read_header(snapshot): 
    ''' read header of snapshot file 
    '''
    filename, fformat = _fname_format(snapshot)

    hdr = {} 

    f = h5py.File(filename, 'r')
    hdr['time']     = f['Header'].attrs[u'Time']
    hdr['redshift'] = f['Header'].attrs[u'Redshift']
    hdr['boxsize']  = f['Header'].attrs[u'BoxSize']
    hdr['filenum']  = f['Header'].attrs[u'NumFilesPerSnapshot']
    hdr['omega_m']  = f['Header'].attrs[u'Omega0']
    hdr['omega_l']  = f['Header'].attrs[u'OmegaLambda']
    hdr['hubble']   = f['Header'].attrs[u'HubbleParam']
    hdr['massarr']  = f['Header'].attrs[u'MassTable']
    hdr['npart']    = f['Header'].attrs[u'NumPart_ThisFile']
    hdr['nall']     = f['Header'].attrs[u'NumPart_Total']
    hdr['cooling']  = f['Header'].attrs[u'Flag_Cooling']
    hdr['format']   = 'hdf5'
    f.close()

    return hdr 


# This function reads a block from an entire gadget snapshot (all files)
# it can read several particle types at the same time. 
# ptype has to be a list. E.g. ptype=[1], ptype=[1,2], ptype=[0,1,2,3,4,5]
def _read_block(snapshot, block, ptype, verbose=False):

    # find the format of the file and read header
    filename, fformat = _fname_format(snapshot)
    head    = _read_header(filename)    
    Nall    = head['nall']
    filenum = head['filenum']

    # find the total number of particles to read
    Ntotal = 0
    for i in ptype:
        Ntotal += Nall[i]

    # find the dtype of the block
    if   block=="POS ":  dtype=np.dtype((np.float32,3))
    elif block=="VEL ":  dtype=np.dtype((np.float32,3))
    elif block=="MASS":  dtype=np.float32
    elif block=="ID  ":  dtype=_read_field(filename, block, ptype[0]).dtype
    else: raise Exception('block not implemented in readgadget!')

    # define the array containing the data
    array = np.zeros(Ntotal, dtype=dtype)


    # do a loop over the different particle types
    offset = 0
    for pt in ptype:

        if filenum==1:
            array[offset:offset+Nall[pt]] = _read_field(snapshot, block, pt)
            offset += Nall[pt]

        # multi-file hdf5 snapshot
        else:

            # do a loop over the different files
            for i in range(filenum):
                
                # find the name of the file to read
                filename = '%s.%d.hdf5'%(snapshot,i)

                # read number of particles in the file and read the data
                npart = _read_header(filename)['npart'][pt]
                array[offset:offset+npart] = _read_field(filename, block, pt)
                offset += npart   

    if offset!=Ntotal:  raise Exception('not all particles read!!!!')
            
    return array


# This function reads a block of an individual file of a gadget snapshot
def _read_field(snapshot, block, ptype):

    filename, fformat = _fname_format(snapshot)
    head              = _read_header(filename)

    prefix = 'PartType%d/'%ptype
    f = h5py.File(filename, 'r')
    if   block=="POS ":  suffix = "Coordinates"
    elif block=="MASS":  suffix = "Masses"
    elif block=="ID  ":  suffix = "ParticleIDs"
    elif block=="VEL ":  suffix = "Velocities"
    else: raise Exception('block not implemented in readgadget!')
    array = f[prefix+suffix][:]
    f.close()

    if block=="VEL ":  array *= np.sqrt(head['time'])
    if block=="POS " and array.dtype==np.float64:
        array = array.astype(np.float32)
    return array


def _fname_format(snapshot):
    # find snapshot name and format
    if os.path.exists(snapshot):
        if snapshot[-4:]=='hdf5':  filename, fformat = snapshot, 'hdf5'
        else:                      filename, fformat = snapshot, 'binary'
    elif os.path.exists(snapshot+'.0'):
        filename, fformat = snapshot+'.0', 'binary'
    elif os.path.exists(snapshot+'.hdf5'):
        filename, fformat = snapshot+'.hdf5', 'hdf5'
    elif os.path.exists(snapshot+'.0.hdf5'):
        filename, fformat = snapshot+'.0.hdf5', 'hdf5'
    else:  
        raise Exception('File (%s) not found!' % snapshot)
    return filename,fformat

