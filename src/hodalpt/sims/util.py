import struct
import numpy as np

from numba import njit, prange


@njit(parallel=False, cache=True, fastmath=True)
def get_cic(posx, posy, posz, weight, lbox, ngrid):
    ''' cloud in cell gridding of x, y, z positions 
    '''
    lcell = lbox/ngrid

    delta = np.zeros((ngrid,ngrid,ngrid))

    for ii in prange(len(posx)):
        xx = posx[ii]
        yy = posy[ii]
        zz = posz[ii]

        if xx<0:
            xx += lbox
        if xx>=lbox:
            xx -= lbox

        if yy<0:
            yy += lbox
        if yy>=lbox:
            yy -= lbox

        if zz<0:
            zz += lbox
        if zz>=lbox:
            zz -= lbox


        indxc = int(xx/lcell)
        indyc = int(yy/lcell)
        indzc = int(zz/lcell)

        wxc = xx/lcell - indxc
        wyc = yy/lcell - indyc
        wzc = zz/lcell - indzc

        if wxc <=0.5:
            indxl = indxc - 1
            if indxl<0:
                indxl += ngrid
            wxc += 0.5
            wxl = 1 - wxc
        elif wxc >0.5:
            indxl = indxc + 1
            if indxl>=ngrid:
                indxl -= ngrid
            wxl = wxc - 0.5
            wxc = 1 - wxl

        if wyc <=0.5:
            indyl = indyc - 1
            if indyl<0:
                indyl += ngrid
            wyc += 0.5
            wyl = 1 - wyc
        elif wyc >0.5:
            indyl = indyc + 1
            if indyl>=ngrid:
                indyl -= ngrid
            wyl = wyc - 0.5
            wyc = 1 - wyl

        if wzc <=0.5:
            indzl = indzc - 1
            if indzl<0:
                indzl += ngrid
            wzc += 0.5
            wzl = 1 - wzc
        elif wzc >0.5:
            indzl = indzc + 1
            if indzl>=0:
                indzl -= ngrid
            wzl = wzc - 0.5
            wzc = 1 - wzl                                                                                                 

        ww = weight[ii]

        delta[indxc,indyc,indzc] += ww * wxc*wyc*wzc
        delta[indxl,indyc,indzc] += ww * wxl*wyc*wzc
        delta[indxc,indyl,indzc] += ww * wxc*wyl*wzc
        delta[indxc,indyc,indzl] += ww * wxc*wyc*wzl
        delta[indxl,indyl,indzc] += ww * wxl*wyl*wzc
        delta[indxc,indyl,indzl] += ww * wxc*wyl*wzl
        delta[indxl,indyc,indzl] += ww * wxl*wyc*wzl
        delta[indxl,indyl,indzl] += ww * wxl*wyl*wzl

    return delta


def read_gadget_ic(filename):
    """
    Read a Gadget format-1 IC/snapshot file with blocks:
      header, positions, velocities, IDs

    Returns
    -------
    header : dict
    pos : (N, 3) float32 array
    vel : (N, 3) float32 array
    ids : (N,) int32 array
    """
    with open(filename, "rb") as f:
        # -----------------
        # Header block
        # -----------------
        header_bytes = _read_block(f)
        if len(header_bytes) != 256:
            raise ValueError(f"Expected 256-byte header, got {len(header_bytes)} bytes")

        # Gadget header layout
        # 6i, 6d, d, d, i, i, 6i, i, i, d, d, d, d
        fmt = "<6i6dddii6iiidddd"
        header_size = struct.calcsize(fmt)
        raw = struct.unpack(fmt, header_bytes[:header_size])

        i = 0
        npart = np.array(raw[i:i+6], dtype=np.int32); i += 6
        mass = np.array(raw[i:i+6], dtype=np.float64); i += 6
        time = raw[i]; i += 1
        redshift = raw[i]; i += 1
        flag_sfr = raw[i]; i += 1
        flag_feedback = raw[i]; i += 1
        npartTotal = np.array(raw[i:i+6], dtype=np.uint32); i += 6
        flag_cooling = raw[i]; i += 1
        num_files = raw[i]; i += 1
        boxsize = raw[i]; i += 1
        omega0 = raw[i]; i += 1
        omegaLambda = raw[i]; i += 1
        hubble = raw[i]; i += 1

        header = {
            "npart": npart,
            "mass": mass,
            "time": time,
            "redshift": redshift,
            "flag_sfr": flag_sfr,
            "flag_feedback": flag_feedback,
            "npartTotal": npartTotal,
            "flag_cooling": flag_cooling,
            "num_files": num_files,
            "BoxSize": boxsize,
            "Omega0": omega0,
            "OmegaLambda": omegaLambda,
            "HubbleParam": hubble,
        }

        N = int(npart.sum())

        # -----------------
        # Positions block
        # -----------------
        pos_bytes = _read_block(f)
        pos = np.frombuffer(pos_bytes, dtype=np.float32).reshape(N, 3)

        # -----------------
        # Velocities block
        # -----------------
        vel_bytes = _read_block(f)
        vel = np.frombuffer(vel_bytes, dtype=np.float32).reshape(N, 3)

        # -----------------
        # IDs block
        # -----------------
        id_bytes = _read_block(f)
        ids = np.frombuffer(id_bytes, dtype=np.int32)

        if len(ids) != N:
            raise ValueError(f"Expected {N} particle IDs, found {len(ids)}")

    return header, pos, vel, ids


def _read_block(f):
    """Read one Fortran-style unformatted block."""
    block_size_start = struct.unpack("<I", f.read(4))[0]
    data = f.read(block_size_start)
    block_size_end = struct.unpack("<I", f.read(4))[0]

    if block_size_start != block_size_end:
        raise IOError(
            f"Block size mismatch: start={block_size_start}, end={block_size_end}"
        )
    return data

