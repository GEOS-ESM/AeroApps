from pyobs import LIDAR_L2, NPZ, OMAERUV_L2
from numpy import savez, load, tile,zeros,size,transpose

from time import time
import VLIDORT_OMI_
#import VLIDORT90_
from math import *
from mie    import getTopo, getAlbedo
from mieobs import *    
from time import clock 

MISSING = -1.2676506E+030
pi = 3.1415926


def vlidort_scalar (c,tau, ssa, g, pe, he, te,iGood=None, verbose=0 ):
        """
       
        Uses VLIDORT f90 to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]  

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo         --- surface albedo
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           radiance       ---  radiance (normalized)
           reflectance    ---  reflectance 
        """

        nch = len(c.channels)               # wavelength number
  
        radiance_, reflectance_, rc = VLIDORT_OMI_.scalar(c.channels,
                                                     tau, ssa, g, pe, he, te, 
                                                     c.albedo[iGood],
                                                     c.solar_zenith[iGood], 
                                                     c.relat_azymuth[iGood],
                                                     c.sensor_zenith[iGood], 
                                                     MISSING,verbose)
        
        if rc != 0:
            raise ValueError, "on return from VLIDORT_OMI_.scalar/vector, rc = "+str(rc)

        return (radiance_, reflectance_)
#---------------------------------------------------------------------
def vlidort_ai_vector (c,tau, ssa, g, pmom, pe, he, te, I=None,iGood=None,verbose=0 ):
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
           radiance       ---  radiance (normalized)
           ai             ---  aerosol index
           reflectivity   ---  reflectivity R*
        """
       

        nch = len(c.channels)               # wavelength number
  
        if I is None:
            I = range(c.ps[iGood].shape[0])
        
        radiance,AI, reflectivity,rc = VLIDORT_OMI_.ai_vector(c.channels,
                                                     tau, ssa, g,pmom, pe, he, te, 
                                                     c.albedo[iGood][I],
                                                     c.solar_zenith[iGood][I], 
                                                     c.relat_azymuth[iGood][I],
                                                     c.sensor_zenith[iGood][I], 
                                                     MISSING,c.radiance[iGood][I],
                                                     c.ai[iGood][I],verbose)
        if rc != 0:
            raise ValueError, "on return from VLIDORT_OMI_.scalar/vector, rc = "+str(rc)
        

        return (radiance, reflectivity,AI)




#-------------------------------------------------------------------------------
if __name__ == "__main__":


#   read OMI files level 2 
    topo_fn  = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    omi = OMAERUV_L2('/nobackup/OMI/Level2/OMAERUV/2007/06/11/OMI-Aura_L2-OMAERUV_2007m0611t0158-o15448_v003-2009m0305t013016.he5')
    omi.sampleLoadz('/nobackup/2/vbuchard/dR_MERRA-AA-r1_OMI/2007/06/11/dR_MERRA-AA-r1_aer_L2-OMAERUV_2007m0611t0158.npz') # AERv file interpolated to obs location 

    aerToUpper(omi.sample)

    # Get the VLIDORT inputs
    # ----------------------
    rank = len(omi.sample.delp.shape)
    if rank > 2 :
       nt, nr, nz = omi.sample.delp.shape 
       raise IOError, "rank > 2 to do - > code not adapted for this case"
    else :
       nobs, nz = omi.sample.delp.shape 

    km = nz
              
    lons = omi.lon.ravel()
    lats = omi.lat.ravel()
    zs = getTopo(topo_fn,lons,lats)


    pe, ze, te = getEdgeVars(omi.sample)
    
    Igood_ = omi.sample.iGood.ravel()    # omi.sample.iGood -> only True values after \
                                         # omi_sampleReduce -> shape (nt,nr) -> 1D array 

    zs_=zs[Igood_]                       # surface height only for good values in npz file         
    he = ze + tile(zs_,(km+1,1))         # height above sea level

    # call vlidort scalar       all in one time
    # -----------------
    print 'rh', omi.sample.RH[0], omi.sample.RH.shape
#    print 'qm', omi.sample.qm[0], omi.sample.qm.shape
    tau, ssa, g= getAOPscalar(omi.sample,omi.channels)
#    radiance_,reflectance_=vlidort_scalar(omi,tau, ssa, g, pe, he, te,omi.sample.iGood)
"""
    # call vlidort vector       all in one time
    # -----------------
             
    nch = len(omi.channels)
    nMom = 300 
    
    mobs=200 
    radiance_ = ones((nobs,nch))
    AI_ = ones((nobs))
    refl_=  ones((nobs))  
    AOT_=ones((nobs,nch))   
    SSA_=ones((nobs,nch))      
    for i in range(0,nobs,mobs):
        
        I = range(i,min(i+mobs,nobs))
        
        print 'channels', omi.channels, omi.sample.RH
        tau_, ssa_, g_, pmom_ = getAOPvector(omi.sample,omi.channels,\
        I= I,nMom=nMom,rcfile='Aod_EOS.rc')

        pe_=pe[:,I].copy()
        he_=he[:,I].copy()
        te_=te[:,I].copy()
        t0 = time()
        print 'time beginning', t0
        radiance,reflectivity,AI=vlidort_ai_vector(omi,tau_, ssa_, g_,pmom_,\
                                                 pe_, he_, te_,I,omi.sample.iGood)
        print 'timmmmeeeee end', (time()-t0)
        radiance_[I,:] = radiance.astype('float32')
        AI_[I] = AI.astype('float32')
        refl_[I] = reflectivity.astype('float32')
#        AOT_[I,:]= AOT
#        SSA_[I,:]= Single
"""    
#    print 'save in a file', radiance_
#    path = '/nobackup/2/vbuchard/VLIDORT_TEST/'
#    savez('TEST_2_6/test_AI_2007_06_11_f90_old2.5.npz',lon=lons[Igood_],\
#    lat=lats[Igood_],\
#    rad_OMI=omi.radiance[omi.sample.iGood],rad_VL=radiance_,\
#    AI_VL=AI_,AI_OMI=omi.ai[omi.sample.iGood], refl_VL=refl_,\
#    refl_OMI=omi.Reflectivity[omi.sample.iGood])

#    savez(path+'test_AI_2007_06_11_f90_n16.npz',lon=lons, lat=lats,\
#    rad_OMI=omi.radiance[omi.sample.iGood],rad_VL=radiance_,\
#    )
