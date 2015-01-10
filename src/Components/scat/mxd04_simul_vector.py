
from numpy import savez, load, tile, zeros, size, transpose, pi, interp, sqrt
from time  import time
from datetime import datetime

import VLIDORT_BRDF_
from pyobs  import NPZ, mxd04
from gfio   import GFIO, GFIOctl
from mieobs import *    

MISSING = mxd04.MISSING

def vlidort_vector (m, tau, ssa, g, pmom, pe, he, te, albedo,U10m, V10m,mr, I=None, verbose=0):
        """
        
        Uses VLIDORT to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]  

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo         --- surface albedo
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pmom     --- scatter phase matrix
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           reflectance
	   
        """
        relat_azimuth = abs(m.SensorAzimuth[m.qa_flag==3][I]-m.SolarAzimuth[m.qa_flag==3][I])
	radiance, reflectance,refl_cx, BRDF, rc = VLIDORT_BRDF_.vector(m.channels,
						    tau, ssa, g,pmom,
						    pe, he, te,
						    albedo[I],U10m[I],V10m[I],mr,
						    m.SolarZenith[m.qa_flag==3][I],
						    relat_azimuth, # check!
						    m.SensorZenith[m.qa_flag==3][I],
						    MISSING, verbose)
        
        if rc != 0:
            raise ValueError, "on return from VLIDORT_BRDF_.scalar/vector, rc = "+str(rc)

        return reflectance,refl_cx, BRDF

#---
def getAlbedo(albe_fn,channels,u10m,v10m):
    """
    Returns ocean albedo.
    """

    # Get precomputed albedo LUT
    # --------------------------
    npz = NPZ(albe_fn)

    # Trimmed wind speed
    # ------------------
    w10m = sqrt(u10m*u10m+v10m*v10m)
    w10m[w10m<0] = 0
    w10m[w10m>50.] = 50.

    # Interpolate albedo
    # ------------------
    albedo = zeros((len(w10m),len(channels)))
    for ch,i in zip(channels,range(len(channels))):
        j = list(npz.channels).index(ch)
        albedo[:,i] = interp(w10m,npz.speed,npz.albedo[:,j])

    return albedo, w10m

#---
def gridBox(m,dLon=0.3125,dLat=0.25,Verbose=False):
    """
    Group observations by gridbox.
    """
    nLon = int(360. / dLon)
    nLat = int(180. / dLat) + 1
    I = ((m.Longitude+180.)/dLon).round().astype(int) # I = 0,...,nLon
    J = ((m.Latitude+90.)/dLat).round().astype(int)   # J = 0,...,nLat
    N = I*nLon + J

    G = dict()
    for n in N:
	G[n] = True
    for n in G.keys():
	if Verbose:
	   N_ = (N==n)
	   i = I[N_][0]
	   j = J[N_][0]
	   lon = round(-180 + i * dLon,2)
	   lat = round( -90 + j * dLat,2)
	   print "Grid box at (%7.2f,%6.2f) has %2d observations"%(lon,lat,size(I[N_]))

    return dict(G=G.keys(),N=N,I=I,J=J) 

#-------------------------------------------------------------------------------
if __name__ == "__main__":

    # External files
    # --------------
    aero_fn = '/nobackup/MERRAero/inst3d_aer_v.ddf'
    topo_fn = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    wind_fn = '/nobackup/MERRA/opendap/slv_Nx.ctl'
    albe_fn = 'ocean_albedo.npz'
    
    sat = 'MYD04'
    Algo = 'OCEAN'
    path = '/nobackup/MODIS/Level2'
    syn_time = datetime(2007,07,10,0)

    # Read MODIS L2 data for this synoptic time
    # -----------------------------------------
    Files = mxd04.granules(path,sat,syn_time,nsyn=24) 
    print 'files', Files   
    m = mxd04.MxD04_L2(Files,Algo,nsyn=24)

#    G = gridBox(m,Verbose=True)
#def HOLD():
		
    # Get Aerosol/Met variables interpolated to obs location
    # ------------------------------------------------------
    a = GFIOctl(aero_fn).sampleVars(m.Longitude[m.qa_flag==3][4000:6000],m.Latitude[m.qa_flag==3][4000:6000],m.Time[m.qa_flag==3][4000:6000],Verbose=True)
    w = GFIOctl(wind_fn).sampleVars(m.Longitude[m.qa_flag==3][4000:6000],m.Latitude[m.qa_flag==3][4000:6000],m.Time[m.qa_flag==3][4000:6000],
                                    onlyVars=['U10M','V10M'],Verbose=True)
    zs = GFIO(topo_fn).interp('zs',m.Longitude[m.qa_flag==3][4000:6000],m.Latitude[m.qa_flag==3][4000:6000])

    aerToUpper(a) # For robustness, create upper case aliases for variables
#def hold():
    a.SS001[:,:] = a.SS002[:,:]
    a.SS002[:,:] = a.SS003[:,:]
    a.SS003[:,:] = a.SS004[:,:]
    a.SS004[:,:] = a.SS005[:,:]
    a.SS005[:,:] = 0.0

    # Estimate ocean albedo
    # ---------------------
    albedo,wind = getAlbedo(albe_fn,m.channels,w.U10M,w.V10M)


    # Get the VLIDORT inputs
    # ----------------------
    pe, ze, te = getEdgeVars(a)
    km = a.DELP.shape[0]
    he = ze + tile(zs,(km+1,1))
    mr = [1.336,1.333,1.331,1.328,1.324,1.317,1.306]  # refractive index fct m.channels
    nobs=m.Longitude[m.qa_flag==3][4000:6000].shape[0]
    nch = len(m.channels)
    nMom = 300
    mobs = 100
    reflectance_ = ones((nobs,nch))
    refl_cx_ = ones((nobs,nch))
    BRDF_ = ones((nobs,nch))
    for i in range(0,nobs,mobs):
        I = range(i,min(i+mobs,nobs))
        print 'I', i, i+mobs, 'nobs=', nobs
        tau_, ssa_, g_, pmom_ = getAOPvector(a,m.channels,\
        I=I,nMom=nMom,rcfile='Aod_MODIS.rc')

        pe_=pe[:,I].copy()
        he_=he[:,I].copy()
        te_=te[:,I].copy()
    # Vector reflectance
    # ------------------
        t0 = time()
        reflectance,refl_cx,BRDF = vlidort_vector (m, tau_, ssa_, g_,pmom_, pe_, he_, te_,\
        albedo,w.U10M,w.V10M,mr,I=I,verbose=1)
        print 'VLIDORT (Vector): %4.2f minutes/1000 obs'%((time()-t0)/(60.*m.nobs/1000.))
        reflectance_[I,:] = reflectance.astype('float32') 
        refl_cx_[I,:] = refl_cx.astype('float32')
        BRDF_[I,:] = BRDF.astype('float32')
    savez('TEST_2_6/reflectance_vector_cx_MODIS_July_10_2007_NSTREAM_6_4000_6000_old_2.5_with0.5sigma2.npz',\
    refl_MODIS=m.reflectance[m.qa_flag==3][4000:6000],refl_VL=reflectance_,refl_cx=refl_cx_,\
    albedo=albedo, BRDF=BRDF_, Wind=wind,sza=m.SolarZenith[m.qa_flag==3],\
    view=m.SensorZenith[m.qa_flag==3][4000:6000],lat=m.Latitude[m.qa_flag==3][4000:6000],\
    lon=m.Longitude[m.qa_flag==3][4000:6000],\
    relat=abs(m.SensorAzimuth[m.qa_flag==3][4000:6000]-m.SolarAzimuth[m.qa_flag==3][4000:6000]))
