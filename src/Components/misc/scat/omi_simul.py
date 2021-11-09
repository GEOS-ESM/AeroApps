#!/usr/bin/env python

"""
  A Python script to create VLIDORT/OMI Level 2a files.
"""


from pyobs import LIDAR_L2, NPZ, OMAERUV_L2, aura
from numpy import savez, load, tile,zeros,size,transpose,where

from time import time
from optparse   import OptionParser   # Command-line args  
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
def vlidort_ai_vector (c,tau, ssa, g, pmom, pe, he, te, oceanler, I=None,iGood=None,verbose=0 ):
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
           oceanler --- ocean Fresnel reflectance lookup table
          
        Outputs :
           radiance       ---  radiance (normalized)
           ai             ---  aerosol index
           reflectivity   ---  reflectivity R*
        """
       

        nch = len(c.channels)               # wavelength number
        if I is None:
            I = range(c.ps[iGood].shape[0])
        
        radiance,AI,Single,AOT, reflectivity,rc,i388calc, trans, spher, residue, BRDF = VLIDORT_OMI_.ai_vector(c.channels,
                                                     tau, ssa, g,pmom, pe, he, te, oceanler, 
                                                     c.albedo[iGood][I],
                                                     c.solar_zenith[iGood][I], 
                                                     c.relat_azymuth[iGood][I],
                                                     c.sensor_zenith[iGood][I], 
                                                     c.qa_grd[iGood][I], 
                                                     MISSING,verbose)
        if rc != 0:
            raise ValueError, "on return from VLIDORT_OMI_.scalar/vector, rc = "+str(rc)
        

        return (radiance, reflectivity,Single,AOT,AI, i388calc, trans, spher, residue, BRDF)




#-------------------------------------------------------------------------------
if __name__ == "__main__":

#-------------------------------------------------------
#   Parse and create filenames, open and read files
#-------------------------------------------------------

#   CHECK INPUT ARGUMENTS
    parser = OptionParser(usage="Usage: %prog [options] OMI file",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()
 
#  GET OMI FILE FROM INPUT ARGUMENT LIST
    if len(args) == 1:
        omi_file = args[0]
    else:
        parser.error("must have 1 argument: OMI file")

# PARSE INPUT OMI FILENAME
    instr, sat, lev, prod, vtime, orb,ver, ptime, t = aura.parseFilename(omi_file)
    model_filen = 'c180R_v202_aura_schill_aer_L2-OMAERUV_'+vtime+'Full.npz'
    print model_filen

# INITIALIZE OMI OBJECT
    omi = OMAERUV_L2(omi_file)

# INPUT MODEL DATAFILE
    omi.sampleLoadz('inst3d_aer_v/'+model_filen)

# READ THE OCEAN LER LUT 
    oceanlut = '/discover/nobackup/pcolarco/aura/Fresnel_OceanLUT_v1.bin'
    oceanler = getFresnelLUT(oceanlut)

#   Work on the QA flag to reduce to values 0 - 7
    dd = omi.ibits(omi.qa_grd,0,4)
    dd

    aerToUpper(omi.sample)

    # Get the VLIDORT inputs
    # ----------------------
    rank = len(omi.sample.DELP.shape)
    if rank > 2 :
       nt, nr, nz = omi.sample.DELP.shape 
       raise IOError, "rank > 2 to do - > code not adapted for this case"
    else :
       nobs, nz = omi.sample.DELP.shape 

    km = nz
              
    lons = omi.lon.ravel()
    lats = omi.lat.ravel()
    topo_fn  = '/discover/nobackup/pcolarco/aura/topography.1152x721.nc'
    zs = getTopo(topo_fn,lons,lats)

#  FROM AIRDENS,DELP IN SAMPLE FILE, GET LAYER EDGE INFO: pe,te,he (see above comments)
    pe, ze, te = getEdgeVars(omi.sample)
    
    Igood_ = omi.sample.iGood.ravel()    # omi.sample.iGood -> only True values after \
                                         # omi_sampleReduce -> shape (nt,nr) -> 1D array 

    zs_=zs[Igood_]                       # surface height only for good values in npz file         
    he = ze + tile(zs_,(km+1,1))         # height above sea level


# ---------------------------------------------
#  Call vlidort scalar       all in one time
# ---------------------------------------------

#    print 'rh', omi.sample.RH[0], omi.sample.RH.shape
#    print 'qm', omi.sample.qm[0], omi.sample.qm.shape
    tau, ssa, g= getAOPscalar(omi.sample,omi.channels)
#    radiance_,reflectance_=vlidort_scalar(omi,tau, ssa, g, pe, he, te,omi.sample.iGood)


# -------------------------------------------
# Call vlidort vector       all in one time
# -------------------------------------------
             
    nch  = len(omi.channels)
    nMom = 300 
    mobs = 200

    radiance_ = ones((nobs,nch))
    BRDF_     = ones((nobs,nch))
    AI_       = ones((nobs))
    i388calc_ = ones((nobs))
    trans_    = ones((nobs))
    spher_    = ones((nobs))
    refl_     = ones((nobs))  
    residue_  = ones((nobs))
    AOT_      = ones((nobs,nch))   
    SSA_      = ones((nobs,nch))
#   nobs      = 70200

#   For testing only
    print nobs
#    for i in range(0,nobs,mobs):
    for i in range(25000,50000,mobs):
        
        I = range(i,min(i+mobs,nobs))
        tau_, ssa_, g_, pmom_ = getAOPvector(omi.sample,omi.channels,\
        I= I,nMom=nMom,rcfile='Aod_EOS.rc')

#        print(tau_[71,0,:])
#        print(ssa_[71,0,:])

        pe_=pe[:,I].copy()
        he_=he[:,I].copy()
        te_=te[:,I].copy()
        t0 = time()

        print 'time beginning', t0, i
        radiance,reflectivity,Single,AOT,AI,i388calc,trans,spher,residue,BRDF= \
                                                 vlidort_ai_vector(omi,tau_, ssa_, g_,pmom_,\
                                                 pe_, he_, te_,oceanler,I,omi.sample.iGood)

        print 'timmmmeeeee end', (time()-t0), i

        radiance_[I,:] = radiance.astype('float32')
        BRDF_[I,:]     = BRDF.astype('float32')
        AI_[I]         = AI.astype('float32')
        i388calc_[I]   = i388calc.astype('float32')
        trans_[I]      = trans.astype('float32')
        spher_[I]      = spher.astype('float32')
        residue_[I]    = residue.astype('float32')
        refl_[I]       = reflectivity.astype('float32')
        AOT_[I,:]      = AOT
        SSA_[I,:]      = Single
        print(AOT_[I,0])
        print(SSA_[I,0])
        print(AI_[I])
        print(lons[I])
        print(lats[I])

    print 'save in a file', radiance_
    savez('./out/ai.'+model_filen,lon=lons[Igood_], \
                              lat=lats[Igood_], \
                              rad_OMI=omi.radiance[omi.sample.iGood], \
                              rad_VL=radiance_,\
                              AI_VL=AI_, \
                              AI_OMI=omi.ai[omi.sample.iGood], \
                              refl_VL=refl_,\
                              refl_OMI=omi.Reflectivity[omi.sample.iGood],\
                              SSA = SSA_, \
                              AOT = AOT_, \
                              i388calc = i388calc_, \
                              trans=trans_, \
                              spher=spher_, \
                              residue=residue_, \
                              brdf=BRDF_)

    savez('./out/radiance.'+model_filen,lon=lons, \
                                    lat=lats, \
                                    rad_OMI=omi.radiance[omi.sample.iGood], \
                                    rad_VL=radiance_ )
