"""
Uses CRTM to compute surface transmittance.
"""

import crtmmodis_ as crtm

from numpy import array, zeros

def getSfcTrans(sample):
    """
    Uses CRTM to compute surface transmittance.

       tau_21, tau_31 = getSfcTrans(sample)
    
    """
    
    N = len(sample.tsh)
    
    u = sample.u.T
    v = sample.v.T 
    t = sample.t.T
    q = sample.qv.T
    o3 = sample.o3.T
    delp = sample.delp.T
    ts = sample.tsh
    u10m = sample.u[:,-1] # fix this, should be u10m from file
    v10m = sample.v[:,-1] # fix this, should be v10m from file

    source_zenith_angle = zeros(N)
    sensor_zenith_angle = zeros(N)
    sensor_scan_angle   = zeros(N)
    
    channels = array([2, 11])
    
    sfctrans,rc = crtm.getsfctrans(u,v,t,q,o3,delp,ts,u10m,v10m,
                               sensor_zenith_angle,source_zenith_angle,sensor_scan_angle,
                               channels,
                               coefpath='/users/adasilva/data/crtm/rel-2.1.3/crtm_coeffs/little_endian/')

    if rc:
        raise ValueError('on return from getsfctrans, rc=<%d>'%rc)

    tau_21, tau_31 = sfctrans

    return (tau_21, tau_31)
    
#-------------------------------
        
if __name__ == "__main__":

    from pyobs import NPZ
    
    print("ok")

    sample = NPZ('/Users/adasilva/workspace/Data_Analysis/AGU-2014/seac4rs_01.npz')

    tau_21, tau_31 = getSfcTrans(sample)
