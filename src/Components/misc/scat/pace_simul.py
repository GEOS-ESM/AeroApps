"""

  Quick and dirty PACE simulator

"""

import os
from numpy         import savez, load, tile, zeros, pi, sin, cos, transpose
from scipy.special import legendre
from datetime      import datetime
from time          import clock 

#import VLIDORT_

from gfio   import GFIO, GFIOctl
from pyobs  import NPZ
from mieobs import *  

MISSING = -1.2676506E+030

PACE_channels = [  350, 360, 385, 412, 425, 443, 460, 475, 490, 510, 532, 555, 583,
                   617, 640, 655, 665, 678, 710, 748, 765, 820, 865, 1245, 1640, 2135,]

PACE_scat = array([180.,178.29,176.07,173.84,171.61,169.37,167.14,164.90,162.67,160.43,
                   158.20,155.96,153.72,151.49,149.25,147.02,144.78,142.55,140.31,138.07,
                   135.84,133.60,131.37,129.13,126.89,124.66,122.42,120.19,117.95,115.71,
                   113.48,111.24,109.01,106.77,104.53,102.30,100.06,97.83,95.59,93.35,
                   91.12,90.00,88.88,86.65,84.41,82.17,79.94,77.70,75.47,73.23,70.99,
                   68.76,66.52,64.29,62.05,59.81,57.58,55.34,53.11,50.87,48.63,46.40,
                   44.16,41.93,39.69,37.45,35.22,32.98,30.75,28.51,26.28,24.04,21.80,
                   19.57,17.33,15.10,12.86,10.63,8.39,6.16,3.93,1.71,0.])

#....................................................

def szangle(lon,lat,tyme):

    # find the sza
    # ------------
    a0 = 0.006918
    a1 = 0.399912
    a2 = 0.006758
    a3 = 0.002697
    b1 = 0.070257
    b2 = 0.000907
    b3 = 0.000148

    doy, xhour = [], [] 
    for t in tyme:
        t0 = datetime(t.year,1,1)
        doy.append((t - t0).days + 1)
        xhour.append(t.hour)
    doy, xhour = aray(doy), array(xhour)
    
#   solar declination in radians
#   ----------------------------
    r  = 2.*pi*(doy-1) / 365.
    dec = a0 - a1*cos(r) + b1*sin(r) - a2*cos(2.*r) + b2*sin(2.*r)- a3*cos(3.*r) + b3*sin(3.*r)  

    # local time in hours
    # -------------------
    timloc = xhour + lon[:]/15.              
    I = (timloc<0)
    timloc[I] = timloc[I] = timloc[I]+24.
    I = (timloc>24)
    timloc[I] = timloc[I] = timloc[I] - 24

    ahr = abs(timloc - 12.)*15.*pi/180.      # hour angle in radians
   
    latRad = lat[:]*2*pi/360.                # latitude in radians

    cossza = sin(latRad)*sin(dec)+cos(latRad)*cos(dec)*cos(ahr)
    cossza[cossza<0] = 0.0
    sza = acos(cossza) * 180. / pi
    sza[sza>=90] = 89.999

    return sza   


#----
def vlidort_vector (channels, solar_zenith, relat_azymuth, sensor_zenith,
                    albedo, tau, ssa, g, pmom, pe, he, te, verbose=0 ):
        """
        
        Uses VLIDORT to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo   --- surface albedo
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pmom     --- scatter phase matrix
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           radiance       ---  radiance (normalized)
           reflectance    ---  reflectance
        """

        radiance_, reflectance_, rc = VLIDORT_.vector(channels,
                                                      tau, ssa, g, pmom,
                                                      pe, he, te, 
                                                      albedo, solar_zenith,
                                                      relat_azymuth,
                                                      sensor_zenith, MISSING, verbose)
        if rc != 0:
            raise ValueError("on return from VLIDORT_.vector/vector, rc = "+str(rc))

        return (radiance_, reflectance_)

#---
def synthLegendre(pmom,scat_angles=None,N=None):
    """
    Reconstruct phase function in physical space as a function of
    scattering angle.
    """
    if scat_angles==None:
        if N==None:
            raise ValueError("scattering angles not specified")
        else:
            scat_angles = linspace(-180.,180.,N)

    mu = cos(pi*scat_angles/180.)
    mu[mu<-1.] = -1.
    mu[mu>+1.] = +1.
    
    nScat = len(mu)

    pmom = pmom.T
    nPol, nMom, nobs, nch, km = pmom.shape

    # Pre-compute Legendre polynomials
    # --------------------------------
    L = zeros((nMom,nScat))
    for j in range(nMom):
        L[j,:] = legendre(j)(mu)

    L = L.T
    
    # Accumulate
    # ----------
    P = zeros((nPol,nScat,nobs,nch,km))
    for i in range(nPol):
        for k in range(nScat):
#            print 'Polarization, mu: ', i, mu[k]
            P_ik = zeros((nobs,nch,km))
            for j in range(nMom):
                P_ik += pmom[i,j,:,:,:] * L[k,j] # can be optimized!
            P[i,k,:,:,:] += P_ik

    pmom = pmom.T # leave it as we found it

    return P.T # (km,nch,nobs,nScat,nPol)
        
#-------------------------------------------------------------------------------
def toFortran():
    """
    Read NPZ file and write out a simple Fortran file with AOD, SSA, P, etc.
    """
    import g5pace_ as g

    npz = NPZ('pace_g5aero.npz')

    d = g.g5pace

    # Transpose phase function
    # ------------------------
    km, nch, nobs, nscat, npol = list(range(5)) # axis ordering for P
    naxes = (nobs,npol,km,nscat,nch)      # new axis ordering
    P = transpose(npz.P,axes=naxes)

    for i in range(len(npz.lon)):

        t = npz.tyme[i]
        d.year, d.month, d.day, d.hour, d.min, d.sec = \
        t.year, t.month, t.day, t.hour, t.minute, t.second
        
        d.lon = npz.lon[i]
        d.lat = npz.lat[i]

        d.u10m = npz.u10m[i]
        d.v10m = npz.v10m[i]

        d.pe = npz.pe[:,i]
        d.he = npz.he[:,i]
        d.te = npz.te[:,i]

        d.qv = npz.qv[:,i]
        d.o3 = npz.o3[:,i]

        d.channels = npz.channels[:]
        d.scat_angles = npz.scat_angles[:]
        
        d.tau[:,:] = npz.tau[:,:,i]
        d.ssa[:,:] = npz.ssa[:,:,i]
        d.g[:,:]   = npz.g[:,:,i]

        d.p11 = P[i,0]
        d.p12 = P[i,1]
        d.p33 = P[i,2]
        d.p34 = P[i,3]

        g.write_('g5aero_p%02d.bin'%(i+1))

if __name__ == "__main__":

    # Site dependency
    # ---------------
    site = os.uname()[1]
    if ('discover' in site) or ('borg' in site):
        dirn = '/discover/nobackup/projects/gmao/yotc/pub/e572p1_tst_01/opendap/assim/'
        aer_fn = dirn + 'inst3_3d_aer_Nv'
        asm_fn = dirn + 'inst3_3d_asm_Nv'
        sfc_fn = dirn + 'inst3_2d_asm_Nx'
        topo_fn  = '/home/adasilva/iesa/aerosol/data/VLIDORT/topo/topography.1152x721.nc'
    elif 'calculon' in site:
        aer_fn = '/nobackup/MERRAero/inst3d_aer_v.ddf'
        topo_fn  = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    else:
        raise ValueError('unknown site <%s>'%site)

    # Channels, moments
    # -----------------
    nMom = 255
    channels = array(PACE_channels[:])

    # Coordinates
    # -----------
    lon = array([-103.,
                 133.,
                   0,
                 -20.,
                 -40.,
                  ])
    lat = array([-48.,
                 41.5,
                 0,
                 22.5,
                 15.,
                  ])
    tyme = array([datetime(2011,0o7,12,0),
                  datetime(2011,0o7,12,6),
                  datetime(2011,0o7,17,12),
                  datetime(2011,0o7,17,12),
                  datetime(2011,0o7,17,12),
                  ])
    
    # Get aer_Nv variables interpolated to obs location
    # -------------------------------------------------
    a = GFIOctl(aer_fn).sampleVars(lon,lat,tyme,Verbose=True)
    A = GFIOctl(asm_fn).sampleVars(lon,lat,tyme,onlyVars=('O3', 'QV','DELP'), Verbose=True)
    s = GFIOctl(sfc_fn).sampleVars(lon,lat,tyme,onlyVars=('U10M', 'V10M'), Verbose=True)
    a.DELP, a.O3, a.QV, a.U10M, a.V10M = A.DELP, A.O3, A.QV, s.U10M, s.V10M  # consolidate
    zs = GFIO(topo_fn).interp('zs',lon,lat)

    # For robustness, create upper case aliases for variables
    # -------------------------------------------------------
    aerToUpper(a)

    # Get the VLIDORT inputs
    # ----------------------
    pe, ze, te = getEdgeVars(a)
    km = a.DELP.shape[0]
    he = ze + tile(zs,(km+1,1))
    tau, ssa, g, pmom = getAOPvector(a,channels,nMom=nMom,rcfile='Aod_PACE.rc')

    P = synthLegendre(pmom,scat_angles=PACE_scat)

    # Save data
    # ---------
    savez('pace_g5aero.npz',
           lon=lon, lat=lat, tyme=tyme,
           u10m=a.U10M, v10m=a.V10M,
           o3 = a.O3, qv = a.QV,
           channels=channels, scat_angles=PACE_scat,
           pe=pe, he=he, te=te,
           tau=tau, ssa=ssa, g=g, P=P)

    # Save Fortran binaries
    # ---------------------
    toFortran()
    
def hold():

    # get albedo
    # ----------
    albedo = TBD

    # find solar zenith angle
    # -----------------------
    sza = szangle(lat,tyme)

    # Make up viewwing geometry
    # -------------------------
     
    # call vlidort
    # ------------
    radiance, reflectance = vlidort_vector(c, channels, sza, albedo, tau_, ssa_, g_,pmom_, pe, he, te)

