"""
   Implements Python interface to Mie Calculator. Based on f2py extension Mie_.

"""

from numpy import array, isfortran, ones, float32, size, zeros
from pylab import plot, axis, title
from MAPL  import config

import scat_MieObs_ 

#...........................................................................

def getFresnelLUT(filename):
    """
    Given a filename, return the ocean Fresnel reflectance LUT
    """
#   Call Fortran for the actual read
#   ---------------------------------------
    oceanler = scat_MieObs_.readfresnel(filename)

    return(oceanler)

#...........................................................................

def getMieDims(rcfile='Aod_EOS.rc'):
    """
        Return dimensions of Mie-table like nPol and nMom.
    """
    cf = config.Config(rcfile)
    dutable = cf('filename_optical_properties_DU')
    nCh, nRh, nBin, nMom, nPol, rc = scat_MieObs_.getmiedims(dutable) 
    if rc != 0:
       raise ValueError("on return from getMieDims, rc = %d"%rc)

    return (nMom, nPol)

#---
def _toArray(channels):
    """Return numpy array, if *channels* is a scalar."""

    if type(channels) is int:
        channels = (channels, )
    elif type(channels) is float:
        channels = (channels, )
    channels = array(channels)
    return channels
#---
def aerToUpper(aer):
    """
    Create upper case aliases for aer variables to cope with changes
    in filename.
    """
    vnames = list(aer.__dict__.keys())
    for v in vnames:
        V = v.upper()
        if v != V:
            aer.__dict__[V] = aer.__dict__[v] # make alias
#---
def getEdgeVars(aer,ptop=1.):
    """
    Given aer object with (airdens,delp) attributes
    returns
           pe --- layer edge pressure [Pa]
           ze --- layer edge height above sfc [m]
           te --- temperature [K] at layer edge

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    """
    if needs_transpose(aer):
        pe, ze, te = scat_MieObs_.getedgevars(aer.AIRDENS.T,aer.DELP.T,ptop)
    else: 
        pe, ze, te = scat_MieObs_.getedgevars(aer.AIRDENS,aer.DELP,ptop)

    return (pe,ze,te)

#--
def getAOPscalar(aer,channels,Verbose=False,rcfile='Aod_EOS.rc'):
    """
    Compute (tau,ssa,g) given aer object.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    """

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)

    # Pack inputs
    # -----------
    vnames = [ 'du001', 'du002', 'du003', 'du004', 'du005',
               'ss001', 'ss002', 'ss003', 'ss004', 'ss005',
               'BCphobic', 'BCphilic',
               'OCphobic', 'OCphilic',
               'SO4' ]

    nq, nch = len(vnames), len(channels)
    nobs = size(aer.PS)
    if needs_transpose(aer): # aer is (nobs,km)
        km = aer.DELP.shape[1]
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH.T
#        rh = zeros((km,nobs))
        for n, v in zip(list(range(nq)),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V].T * aer.DELP.T / 9.81 
    else:                    # aer is (km,nobs)
        km = aer.DELP.shape[0]
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH
#        rh = zeros((km,nobs))
        for n, v in zip(list(range(nq)),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V] * aer.DELP / 9.81 

    # Do the Mie calculation
    # ----------------------
    tau,ssa,g,rc = scat_MieObs_.getaopscalar(rcfile,channels,pad(vnames),Verbose,qm,rh)

    if rc!=0:
        print("<<<ERROR>>> on return from scat_MieObs_.getaopscalar, rc = ", rc)
        raise ValueError('cannot get Aerosol Optical Properties (scalar version)')

    return (tau,ssa,g)
#--
def getAOPvector(aer,channels,I=None,Verbose=False,rcfile='Aod_EOS.rc',nMom=301):
    """
    Compute (tau,ssa,g,pmom) given aer object.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    
    J --- index of subset of observations to process
    """

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)

    # Pack inputs
    # -----------
    vnames = [ 'du001', 'du002', 'du003', 'du004', 'du005',
               'ss001', 'ss002', 'ss003', 'ss004', 'ss005',
               'BCphobic', 'BCphilic',
               'OCphobic', 'OCphilic',
               'SO4' ]

    nq, nch = len(vnames), len(channels) 
    if I is None:
        nobs = size(aer.PS)
        I = list(range(0,nobs))
    else :
        nobs = len(I)
    
    if needs_transpose(aer): # aer is (nobs,km)
        km = aer.DELP.shape[1]
        
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[I].T
        for n, v in zip(list(range(nq)),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V][I].T * aer.DELP[I].T / 9.81 
    else:                    # aer is (km,nobs)
        km = aer.DELP.shape[0]
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[:,I]
        for n, v in zip(list(range(nq)),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V][:,I] * aer.DELP[:,I] / 9.81 
    

    # Do the Mie calculation
    # ----------------------
    nMom_mieTable,nPol_mieTable = getMieDims(rcfile=rcfile)  # return nMom & nPol of the Mie tables
    nPol = 6                                                 # for dust non spherical
    tau,ssa,g,pmom,rc = scat_MieObs_.getaopvector(rcfile,channels,pad(vnames),Verbose,qm,rh,nMom,nPol)

    if rc!=0:
        print("<<<ERROR>>> on return from scat_MieObs_.getaopvector, rc = ", rc)
        raise ValueError('cannot get Aerosol Optical Properties (vector version)')

    return (tau,ssa,g,pmom)
#---
def getO3tau(aer,channels,Verbose=False,I=None,OMPS=False,):
    """
    Compute (tau) due to ozone absorption given aer object including temperature.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    
    J --- index of subset of observations to process
    """

    # Can Li provides O3 cross section coefficients at OMI wavelengths convolved with
    # OMI slit function.  Given ozone in VMR and the pressure thickness of each layer
    # we need to calculate the # of molecules cm-2 and multiply by the cross section
    # derived from the following coefficients as:
    #  xsec = 1.e-20 * a0 * (1. + a1*T + a2*T*T)
    # where T is the temperature [deg C] and xsec is [cm2]

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)


    nch = len(channels) 
    if I is None:
        nobs = size(aer.PS)
        I = list(range(0,nobs))
    else :
        nobs = len(I)
    
    if needs_transpose(aer): # aer is (nobs,km)
        km    = aer.DELP.shape[1]
        delp  = aer.DELP[I].T
        o3    = aer.O3[I].T    #O3 is vmr [mol mol-1]
        t     = aer.T[I].T     # [K]
    else:                    # aer is (km,nobs)
        km    = aer.DELP.shape[0]
        delp  = aer.DELP[:,I]
        o3    = aer.O3[:,I]
        t     = aer.T[:,I]

    if OMPS:
        tau,rc = scat_MieObs_.geto3tauomps(channels,Verbose,o3,t,delp)
    else:
        tau,rc = scat_MieObs_.geto3tauomi(channels,Verbose,o3,t,delp)

    return (tau)
#---
def getSO2tau(aer,channels,Verbose=False,I=None,OMPS=False,):
    """
    Compute (tau) due to SO2 absorption given aer object including temperature.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    
    J --- index of subset of observations to process
    """

    # Can Li provides SO2 cross section coefficients at OMI wavelengths convolved with
    # OMI slit function.  Given SO2 in MMR and the pressure thickness of each layer
    # we need to calculate the # of molecules cm-2 and multiply by the cross section
    # derived from the following coefficients as:
    #  xsec = 1.e-20 * a0 * (1. + a1*T + a2*T*T)
    # where T is the temperature [deg C] and xsec is [cm2]

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)


    nch = len(channels) 
    if I is None:
        nobs = size(aer.PS)
        I = list(range(0,nobs))
    else :
        nobs = len(I)
    
    if needs_transpose(aer): # aer is (nobs,km)
        km    = aer.DELP.shape[1]
        delp  = aer.DELP[I].T
        so2   = aer.SO2S[I].T + aer.SO2V[I].T   #O3 is vmr [mol mol-1]
        t     = aer.T[I].T     # [K]
    else:                    # aer is (km,nobs)
        km    = aer.DELP.shape[0]
        delp  = aer.DELP[:,I]
        so2   = aer.SO2S[:,I] + aer.SO2V[:,I]
        t     = aer.T[:,I]

    if OMPS:
        tau,rc = scat_MieObs_.getso2tauomps(channels,Verbose,so2,t,delp)
    else:
        tau,rc = scat_MieObs_.getso2tauomi(channels,Verbose,so2,t,delp)

    return (tau)
#---
def getAOPext(aer,channels,I=None,vnames=None,Verbose=False):
    """
    Compute (ext,backscat,aback_sfc,aback_toa) given aer object.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    
    I --- index of subset of observations to process
    """

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)

    # Pack inputs
    # -----------
    if vnames is None:
       vnames = [ 'du001', 'du002', 'du003', 'du004', 'du005',
               'ss001', 'ss002', 'ss003', 'ss004', 'ss005',
               'BCphobic', 'BCphilic',
               'OCphobic', 'OCphilic',
               'SO4' ]

    nq, nch = len(vnames), len(channels) 

    if I is None:
        nobs = size(aer.PS)
        I = list(range(0,nobs))
    else :
        nobs = len(I)
    
    if needs_transpose(aer): # aer is (nobs,km)
        km = aer.DELP.shape[1]
        
        qc = ones((km,nq,nobs),dtype=float32)
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[I].T
#        rh = zeros((km,nobs))
        for n, v in zip(list(range(nq)),vnames):
            V = v.upper()
            qc[:,n,:] = aer.__dict__[V][I].T * aer.AIRDENS[I].T 
            qm[:,n,:] = aer.__dict__[V][I].T * aer.DELP[I].T / 9.81 
    else:                    # aer is (km,nobs)
        km = aer.DELP.shape[0]
        
        qc = ones((km,nq,nobs),dtype=float32)
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[I]
#        rh = zeros((km,nobs))
        for n, v in zip(list(range(nq)),vnames):
            V = v.upper()
            qc[:,n,:] = aer.__dict__[V][:,I] * aer.AIRDENS[:,I] 
            qm[:,n,:] = aer.__dict__[V][:,I] * aer.DELP[:,I] / 9.81 

    # Do the Mie calculation
    # ----------------------
    
    ext,sca,backscat,aback_sfc,aback_toa,rc = scat_MieObs_.getext(channels,pad(vnames),Verbose,qc,qm,rh)

    if rc!=0:
        print("<<<ERROR>>> on return from scat_MieObs_.getaopvector, rc = ", rc)
        raise ValueError('cannot get Aerosol Optical Properties (vector version)')

    return (ext,sca,backscat,aback_sfc,aback_toa)

#........................................................................
def pad(names):
    """
    Make all strings in list *names* the same size for f2py's benefit.
    """
    return [ "%-16s"%v for v in names ]

#........................................................................
def needs_transpose(aer):
    """
    Returns
         True  if aer object has shapes (nobs,km)
         False if aer object has shapes (km,nobs)
    This is needed to cope with the way the fortran expects
    the data arrays.
    """
    if aer.PS.shape[0] == aer.DELP.shape[0]: # (nobs,km)
        return True
    else:
        return False
    
#........................................................................

if __name__ == "__main__":

    from pyobs import LIDAR_L2, NPZ
    from mie   import getTopo
    
    channels = array([532.,])

    c_dn = '/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4/' # Pete's Calipso interp
    a_dn = '/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4/'         # aer_Nv + interp

    c_fn = 'dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc' # Pete's interp
    a_fn = 'dR_Fortuna-2-4-b4.calipso_aer.20090715.npz' # aer collocation

    # Read relevant data
    # ------------------
    c = LIDAR_L2(c_dn+c_fn)  # Pete's interp with CALPSO coordinates
    a = NPZ(a_dn+a_fn)        # aer_v interpolated to obs location
    aerToUpper(a)
#   ---------------------------------------------------------------------------    
#   NOTE:
#
#   For creating the npz file above with the aer_v data interpolated to obs ,
#   location you do this:
#      c.sampleFile(aer_v,npzFile=npzFile)   
#   where aer_v is the full, gridded, netcdf file.
#   In principle, you do not need to write the npzFile as the result of the
#   sampling is also available as attribute *sample*. So, instead of reading
#   the npzFile you could do
#     a = c.sample
#   But sampling takes time, so saving a "sampling  file" is convenient.
#   BTW, npzFile is too python specific. Soon, I'll implement a NetCDF option,
#   so you will be able to say
#     c.sample(aer_v,ncFile=NetCDF_filename)
#   See bottom of GMAO_pyobs/pyobs.lidar_l2.py file for an example.
#   ---------------------------------------------------------------------------    

    pe, ze, te = getEdgeVars(a)
    
    tau, ssa, g = getAOPscalar(a,channels,rcfile='Aod_EOS.rc')
    ext,sca,backscat,aback_sfc,aback_toa = getAOPext(a,channels)
