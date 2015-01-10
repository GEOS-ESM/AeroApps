"""
   Implements Python interface to Mie Calculator. Based on f2py extension Mie_.

"""

from numpy import array, isfortran
from pylab import plot, axis, title
from MAPL  import config

from Mie_ import *

def getDims(filename):
    """Return dimensions of GFIO-compatible file. The output tuple is
    (tm,km,jm,im)."""
    im, jm, km, tm, rc = getdims(filename)
    if rc != 0:
        raise ValueError, "on return from getDIms, rc = %d"%rc

    return (tm, km, jm, im)

def getMieDims(rcfilename='Aod_EOS.rc'):
    """Return dimensions of Mie-table like nPol and nMom."""
    cf = config.Config(rcfilename)
    dutable = cf('filename_optical_properties_DU')
    nCh, nRh, nBin, nMom_mieTab, nPol_mieTab, rc = getmiedims(dutable) 
    if rc != 0:
       raise ValueError, "on return from getMieDIms, rc = %d"%rc

    return (nMom_mieTab, nPol_mieTab)

#def setnMom(rcfilename='Aod_EOS.rc'):
#    nMom = 6      
#    cf = config.Config(rcfilename)
#    cf.set('n_moments', nMom)
#    return (nMom)

#...............................................................

def getMieGriddedScal(filename,channels,nymd=-1,nhms=-1,verbose=0):
    """
    Given a file name with a Chem Bundle containing aerosol concentrations,
    along with a tuple (or single value) of wavenumbers (channels, in nm), 
    this function returns 3D gridded fields with the following parameters:
    
    tau --- aerosol optical  
    ssa --- single scattering albedo
    g   --- asymmetry factor


    of shape (km,jm,im,nch), C ordering. In addition, it also returns the
    following met fields of shape(km+1,jm,im), C ordering:

    pe  --- pressure at layer edges [Pa]    
    ze  --- height above surface at layer edges [Pa]    
    te  --- temperature at layer edges [K]    

    On input, nymd, nhms are the date and time to be read from file; the
    default is the last time stamp on file.
    
    """

#   Make sure channels is a numpy array
#   -------------------------------------
    channels = _toArray(channels)

#   Call Fortran for the actual calculation
#   ---------------------------------------

    tm, km, jm, im = getDims(filename)
    nMom, nPol = getMieDims(rcfilename='Aod_EOS.rc')
    
    tau, ssa, g, pmom, pe, ze, te, rc = \
          getmiegridded(filename,im,jm,km,nMom,nPol,nymd,nhms,channels,verbose)

    if rc != 0:
        raise ValueError, "on return from getMie extension, rc = "+str(rc)

#   Transpose the array, converting it to C order in the process;
#   remove dimensions of length 1, most like the channel dimension
#   --------------------------------------------------------------
    tau = tau.transpose().squeeze()
    ssa = ssa.transpose().squeeze()
    g   =   g.transpose().squeeze()
    pe  =  pe.transpose().squeeze()
    ze  =  ze.transpose().squeeze()
    te  =  te.transpose().squeeze()

    return (tau, ssa, g, pe, ze, te)
#........................................................................

def getMieGriddedVect(filename,channels,nymd=-1,nhms=-1,verbose=0):
    """
    Given a file name with a Chem Bundle containing aerosol concentrations,
    along with a tuple (or single value) of wavenumbers (channels, in nm), 
    this function returns 3D gridded fields with the following parameters:
    
    tau --- aerosol optical  
    ssa --- single scattering albedo
    g   --- asymmetry factor

    pmom --- components of the phase function matrix
    nMom --- nombre de moments que l'on souhaite pour Vlidort

    of shape (km,jm,im,nch), C ordering. In addition, it also returns the
    following met fields of shape(km+1,jm,im), C ordering:

    pe  --- pressure at layer edges [Pa]    
    ze  --- height above surface at layer edges [Pa]    
    te  --- temperature at layer edges [K]    

    On input, nymd, nhms are the date and time to be read from file; the
    default is the last time stamp on file.
    
    """

#   Make sure channels is a numpy array
#   -------------------------------------
    channels = _toArray(channels)

#   Call Fortran for the actual calculation
#   ---------------------------------------

    tm, km, jm, im = getDims(filename)
    nMom_mieTab, nPol_mieTab = getMieDims(rcfilename='Aod_EOS.rc')
    
    tau, ssa, g, pmom, pe, ze, te, rc = \
          getmiegridded(filename,im,jm,km,nMom,nPol_mieTab,nymd,nhms,channels,verbose)

    if rc != 0:
        raise ValueError, "on return from getMie extension, rc = "+str(rc)

#   Transpose the array, converting it to C order in the process;
#   remove dimensions of length 1, most like the channel dimension
#   --------------------------------------------------------------
    tau = tau.transpose().squeeze()
    ssa = ssa.transpose().squeeze()
    g   =   g.transpose().squeeze()
    pmom =  pmom.transpose().squeeze()
    pe  =  pe.transpose().squeeze()
    ze  =  ze.transpose().squeeze()
    te  =  te.transpose().squeeze()

    return (tau, ssa, g, pmom, pe, ze, te)

#.......................................................................

def getMieScal(filename,channels,lon,lat,nymd=-1,nhms=-1,verbose=0):
    """
    Given a file name with a Chem Bundle containing aerosol concentrations,
    along with a tuple (or single value) of wavenumbers (channhels, in nm), 
    and a list of coordinates:

    lon --- longitudes in degrees, [-180,180]
    lat --- latitudes in degrees, [-90,90]

    this function returns
    
    tau --- aerosol optical  
    ssa --- single scattering albedo
    g   --- asymmetry factor
  
   

    of shape = (km,nch,nobs), in Fortran ordering. In addition, it also returns the
    following met fields of shape(km+1,nobs), Fortran ordering:

    pe  --- pressure at layer edges [Pa]    
    ze  --- height above surface at layer edges [Pa]    
    te  --- temperature at layer edges [K]    

    On input, nymd, nhms are the date and time to be read from file; the
    default is the last time stamp on file.
    
    """

#   Make sure channels is a numpy array
#   -------------------------------------
    channels = _toArray(channels)

#   Call Fortran for the actual calculation
#   ---------------------------------------
    tm, km, jm, im = getDims(filename)
    nMom_mieTab, nPol_mieTab = getMieDims(rcfilename='Aod_EOS.rc')
    nMom = 4
#    print nMom
    tau, ssa, pe, ze, te, g, rc =scalar(filename,km,nMom,nPol_mieTab,nymd,nhms,channels,lon,lat,verbose)
    
    if rc != 0:
        raise ValueError, "on return from getMie extension, rc = "+str(rc)

    return (tau, ssa, pe, ze, te, g)

#........................................................................

def getMieVect(filename,channels,lon,lat,nymd=-1,nhms=-1,verbose=0):
    """
    Given a file name with a Chem Bundle containing aerosol concentrations,
    along with a tuple (or single value) of wavenumbers (channhels, in nm), 
    and a list of coordinates:

    lon --- longitudes in degrees, [-180,180]
    lat --- latitudes in degrees, [-90,90]

    this function returns
    
    tau --- aerosol optical  
    ssa --- single scattering albedo
    g   --- asymmetry factor
  
    pmom --- components of the phase function matrix
   

    of shape = (km,nch,nobs), in Fortran ordering. In addition, it also returns the
    following met fields of shape(km+1,nobs), Fortran ordering:

    pe  --- pressure at layer edges [Pa]    
    ze  --- height above surface at layer edges [Pa]    
    te  --- temperature at layer edges [K]    

    On input, nymd, nhms are the date and time to be read from file; the
    default is the last time stamp on file.
    
    """

#   Make sure channels is a numpy array
#   -------------------------------------
    channels = _toArray(channels)

#   Call Fortran for the actual calculation
#   ---------------------------------------
    tm, km, jm, im = getDims(filename)
    nMom_mieTab, nPol_mieTab = getMieDims(rcfilename='Aod_EOS.rc')
    nMom = 500
    
#    print nMom_mieTab, nPol_mieTab, tm, km, jm, im
    nPol = nPol_mieTab
    tau, ssa, pe, ze, te, g, pmom,  rc = \
         vector(filename,km,nMom,nPol_mieTab,nymd,nhms,channels,lon,lat,verbose)
    
    if rc != 0:
        raise ValueError, "on return from getMie extension, rc = "+str(rc)

    return (nMom, nPol, tau, ssa, pe, ze, te, g, pmom)

#...........................................................................

def getTopo(filename,lon,lat):
    """
    Given a file containng topography, returns the topographic height
    interpolated to obs locations.
    """
#   Call Fortran for the actual calculation
#   ---------------------------------------
    tm, km, jm, im = getDims(filename)
    zs, rc = \
         gettopobs(filename,lon,lat)

    if rc != 0:
        raise ValueError, "on return from getTopObs extension, rc = "+str(rc)

    return zs

#...........................................................................

def getAlbedo(filename,time,lon,lat):
    from grads import GrADS
    from grads.gacore import dt2gat
    from datetime import datetime 
    time_ = datetime(2000,time.month,time.day,12)
    ga = GrADS(Window=False,Echo=False)
    fh = ga.open(filename)
    ga("set time %s"%dt2gat(time_))

    albedo, levs = ga.interp('albedo',lon,lat)

    return albedo
    
#........................................................................

def _toArray(channels):
    """Return numpy array, if *channels* is a scalar."""

    if type(channels) is int:
        channels = (channels, )
    elif type(channels) is float:
        channels = (channels, )
    channels = array(channels)
    return channels

def plot_tau(tau,pe,n):
    """Simple plot of tau for location n."""
    pm = (pe[1:,n] + pe[0:-1,n] ) / 200.
    plot(tau[:,0,n],pm)
    axis([0, 0.025, 1000, 100])
    title('Aerosol Optical Thickness')

def plot_ssa(ssa,pe,n):
    """Simple plot of tau for location n."""
    pm = (pe[1:,n] + pe[0:-1,n] ) / 200.
    plot(ssa[:,0,n],pm)
    axis([0.75, 1, 1000, 100])
    title('Single Scattering Albedo')

#........................................................................

if __name__ == "__main__":
   
    import pyods

    nymd = 20080629
    nhms = 0
    ods = pyods.ODS('ods/omi.aero_tc8.obs.20080629.ods',nymd,nhms).select(kt=45)

#    filename = 'data/a5_arctas_02.inst3d_aer_v.20080629_0000z.nc'
#    channels = 550.
#    tau, ssa, g, pm = getMieGridded(filename,channels,nymd,nhms)
#    tau, ssa, g, pe, ze, te = \
#           getMieScal(filename,channels,ods.lon,ods.lat,nymd,nhms,verbose=1)
#    zs = getTopo('data/topography.1152x721.nc',ods.lon,ods.lat)

#    he = ze + zs.reshape((1,)+zs.shape)

