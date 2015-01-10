
from numpy import savez, load, tile, zeros, size, transpose, pi, interp, sqrt
from time  import time
from datetime import datetime

#import VLIDORT_OMI_
import VLIDORT_BRDF_
from pyobs  import NPZ, mxd04
from gfio   import GFIO, GFIOctl
from mieobs import *    

MISSING = mxd04.MISSING

def vlidort_scalar (m, tau, ssa, g, pe, he, te, albedo,U10m, V10m,mr,verbose=0):
        """
        
        Uses VLIDORT to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]  

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo         --- surface albedo
           U10m           --- Wind speed components
           V10m
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           reflectance
	   
        """
        relat_azimuth = abs(m.SensorAzimuth[m.qa_flag==3]-m.SolarAzimuth[m.qa_flag==3])
	radiance, reflectance,refl_cx, BRDF,rc = VLIDORT_BRDF_.scalar(m.channels,
						    tau, ssa, g,
						    pe, he, te,
						    albedo,U10m,V10m,mr,
						    m.SolarZenith[m.qa_flag==3],
						    relat_azimuth, # check!
						    m.SensorZenith[m.qa_flag==3],
						    MISSING, verbose)
        
        if rc != 0:
            raise ValueError, "on return from VLIDORT_.scalar/vector, rc = "+str(rc)

        return reflectance, refl_cx, BRDF

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

    return albedo

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
    syn_time = datetime(2007,06,28,0)

    # Read MODIS L2 data for this synoptic time
    # -----------------------------------------
    Files = mxd04.granules(path,sat,syn_time,nsyn=24)    
    m = mxd04.MxD04_L2(Files,Algo,nsyn=24)

#    G = gridBox(m,Verbose=True)
#def HOLD():
		
    # Get Aerosol/Met variables interpolated to obs location
    # ------------------------------------------------------
    a = GFIOctl(aero_fn).sampleVars(m.Longitude[m.qa_flag==3],m.Latitude[m.qa_flag==3],m.Time[m.qa_flag==3],Verbose=True)
    w = GFIOctl(wind_fn).sampleVars(m.Longitude[m.qa_flag==3],m.Latitude[m.qa_flag==3],m.Time[m.qa_flag==3],
                                    onlyVars=['U10M','V10M'],Verbose=True)
    zs = GFIO(topo_fn).interp('zs',m.Longitude[m.qa_flag==3],m.Latitude[m.qa_flag==3])

    aerToUpper(a) # For robustness, create upper case aliases for variables

    a.SS001[:,:] = a.SS002[:,:]
    a.SS002[:,:] = a.SS003[:,:]
    a.SS003[:,:] = a.SS004[:,:]
    a.SS004[:,:] = a.SS005[:,:]
    a.SS005[:,:] = 0.0

    # Estimate ocean albedo
    # ---------------------
    albedo = getAlbedo(albe_fn,m.channels,w.U10M,w.V10M)
    albedo = albedo
    # Get the VLIDORT inputs
    # ----------------------
    pe, ze, te = getEdgeVars(a)
    km = a.DELP.shape[0]
    he = ze + tile(zs,(km+1,1))
    tau, ssa, g = getAOPscalar(a,m.channels,rcfile='Aod_MODIS.rc')
    mr = [1.336,1.333,1.331,1.328,1.324,1.317,1.306]  # refractive index fct m.channels

    # Scalar reflectance
    # ------------------
    t0 = time()
#    reflectance = vlidort_scalar (m, tau, ssa, g, pe, he, te, albedo, verbose=1)
    # Using my Cox Munk
    reflectance,refl_cx,BRDF = vlidort_scalar (m, tau, ssa, g, pe, he, te, albedo,w.U10M,w.V10M,mr,verbose=1)
    print 'VLIDORT (Scalar): %4.2f minutes/1000 obs'%((time()-t0)/(60.*m.nobs/1000.))
    path = '/nobackup/2/vbuchard/VLIDORT_TEST/MODIS_refl'
    savez(path+'reflectance_scalar_MODIS_Giss_Cox_Munk_VL_f90_all.npz',\
    refl_MODIS=m.reflectance[m.qa_flag==3],refl_VL=reflectance,refl_cx=refl_cx,\
    albedo=albedo, BRDF=BRDF, sza=m.SolarZenith[m.qa_flag==3],\
    view=m.SensorZenith[m.qa_flag==3],\
    relat=abs(m.SensorAzimuth[m.qa_flag==3]-m.SolarAzimuth[m.qa_flag==3]))
