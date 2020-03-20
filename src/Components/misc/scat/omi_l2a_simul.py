#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
  A Python script to create VLIDORT/OMI Level 2a files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
import sys
import VLIDORT_OMI_

from time       import clock
from optparse   import OptionParser   # Command-line args  
from pyobs      import OMAERUV_L2,aura
from mieobs     import getEdgeVars, getAOPvector, getAOPscalar,aerToUpper,aerToLower 
from mie        import getTopo
from MAPL       import strTemplate
from numpy      import tile, ones, savez


MISSING = -1.2676506E+030

#---------------------------------------------------------------------
def makethis_dir(filename):
    """Creates the relevant directory if necessary."""
    path, filen = os.path.split(filename)
    if path != '':
        rc = os.system('mkdir -p '+path)
        if rc:
            raise IOError, "could not create directory "+path
#---------------------------------------------------------------------
def vlidort_scalar (c,tau, ssa, g, pe, he, te,iGood=None, verbose=0 ):
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
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           radiance       ---  radiance (normalized)
           ai             ---  aerosol index
           reflectivity   ---  reflectivity
        """

        nch = len(c.channels)               # wavelength number
  
        radiance, ai, reflectivity, rc = VLIDORT_OMI_.ai_scalar(c.channels,
                                                     tau, ssa, g, pe, he, te, 
                                                     c.albedo[iGood],
                                                     c.solar_zenith[iGood], 
                                                     c.relat_azymuth[iGood],
                                                     c.sensor_zenith[iGood], 
                                                     MISSING,c.radiance[iGood], 
                                                     c.ai[iGood], verbose)
        
        if rc != 0:
            raise ValueError, "on return from OMI_.scalar/vector, rc = "+str(rc)

        return (radiance, ai, reflectivity)
        
#---------------------------------------------------------------------
def vlidort_vector (c,tau, ssa, g, pmom, pe, he, te, I=None,iGood=None,verbose=0 ):
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
           reflectivity   ---  reflectivity
        """
       

        nch = len(c.channels)               # wavelength number
  
        if I is None:
            I = range(c.ps[iGood].shape[0])
        
        radiance, ai, reflectivity, rc = VLIDORT_OMI_.ai_vector(c.channels,
                                                     tau, ssa, g,pmom, pe, he, te, 
                                                     c.albedo[iGood][I],
                                                     c.solar_zenith[iGood][I], 
                                                     c.relat_azymuth[iGood][I],
                                                     c.sensor_zenith[iGood][I], 
                                                     MISSING,c.radiance[iGood][I], 
                                                     c.ai[iGood][I], verbose)
        if rc != 0:
            raise ValueError, "on return from OMI_.scalar/vector, rc = "+str(rc)
        

        return (radiance, ai, reflectivity)



#------------------------------------------------------------------------------------
if __name__ == "__main__":

    expid        = 'dR_MERRA-AA-r2'

#   Defaults may be platform dependent
#   ----------------------------------
    if os.path.exists('/nobackup/2/vbuchard/dR_MERRA-AA-r1_OMI/'): # New calculon
        aer_dir = '/nobackup/2/vbuchard/dR_MERRA-AA-r1_OMI/%y4/%m2/%d2'     
        topography   = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    elif os.path.exists('/home/vbuchard/workspace/OMI_AI/dR_MERRA-AA-r2_2007/'): # Discover
        aer_dir = '/discover/nobackup/vbuchard/VLIDORT/AI_exp_2007/OMI_1.4.7/SAMPLE/seac4rs_01/%y4/%m2/%d2' 
        topography   = '/discover/nobackup/projects/gmao/iesa/aerosol/data/VLIDORT/topo/topography.1152x721.nc'
    else: # Must be somewhere else, no good defaults
        aer_dir = '/'  # to complete
        topography   = 'topography.1152x721.nc'

    out_dir      = '/discover/nobackup/vbuchard/VLIDORT/AI_exp_2007/OMI_1.4.7/Results_seac4rs_01/SSA_dec_0.05/.'
#    out_fname    = '%s.omi_l3a.%y4%m2%d2%s.nc4'
    out_fname    = 'SEAC4RS_01.OMAERUV1.4.7_l2a.SSA_OC_d_0.05_%y4%m2%d2_%h2%n2z'
    scalar       = False
    res          = 'c'
    
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] OMI file",
                          version='omi_level2a-1.0.0' )

    parser.add_option("-A", "--aer", dest="aer_dir", default=aer_dir,
                      help="aerosol sample directory (default=%s)"\
                           %aer_dir )

    parser.add_option("-t", "--topo", dest="topography", default=topography,
                      help="File name for topography data (default=%s)"\
                           %topography )

    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="experiment id (default=%s)"\
                           %expid )

    parser.add_option("-d", "--dir", dest="out_dir", default=out_dir,
                      help="output directory (default=%s)"\
                           %out_dir )

    parser.add_option("-f", "--fname", dest="out_fname", default=out_fname,
                      help="output file tempalte  (default=%s)"\
                           %out_fname )

    parser.add_option("-r", "--res", dest="res", default=res,
                      help="grid resolution of VLIDORT simulations   (default=%s)"\
                           %out_fname )

    parser.add_option("-u", "--uncompressed",
                      action="store_true", dest="uncompressed",
                      help="do not use n4zip to compress gridded output file (default=False)")

    parser.add_option("-S", "--skipVLIDORT",
                      action="store_true", dest="skipVLIDORT",
                      help="do not call VLIDORT (default=False)")

    parser.add_option("-F", "--force",
                      action="store_true", dest="force",
                      help="overwrites output file")

    parser.add_option("-s", "--scalar",
                      action="store_true", dest="scalar",
                      help="runs VLIDORT in scalar mode")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="turn on verbosity.")

    (options, args) = parser.parse_args()
    
    if len(args) == 1:
       omi_file = args[0]
    else:
       parser.error("must have 1 argument: OMI file")

#   Parse the OMI file name (ex OMI-Aura_L2-OMAERUV_2007m0611t0019-o15447_v003-2009m0305t013119.he5)
#   ------------------------      
    instr, sat, lev, prod, vtime, orb,ver, ptime, t = aura.parseFilename(omi_file)


#   GEOS-5 Aerosol sample file name (npz file)
#   -------------------------
    aer_dir = strTemplate(options.aer_dir,dtime=t)
    aer_file = aer_dir+'/'+options.expid+'_aer_'\
               + lev + '-'+ prod + '_' + vtime + '.npz'


    if options.verbose:
        print ""
        print "                          VLIDORT/OMI Level 2A Processing"
        print "                          -------------------------------"
        print ""
        t0 = clock()
            
    if options.skipVLIDORT:
        print "IMPORTANT: Skipping VLIDORT"

#   Output ungridded file
#   -------------------
    out_file = strTemplate(options.out_dir+'/'+options.out_fname,
                           expid=options.expid, dtime=t)
    if os.path.exists(out_file) and (options.force is not True):
        print "omi_l2a: Output file <%s> exists --- cannot proceed."%out_file
        raise IOError, "Specify --force to overwrite existing output file."
    else:
        makethis_dir(out_file) # make sure directory exists

#   Read OMI and GEOS-5 aerosol data at OMI obs
#   -------------
    if options.verbose:
        t = clock() - t0
        print "- Reading OMI measurements at t=%f"%t

    omi = OMAERUV_L2(omi_file)
    omi.sampleLoadz(aer_file)
    aerToUpper(omi.sample)
    aerToLower(omi.sample)
#   Compute Mie parameters at the OMI locations
#   -------------------------------------------
#    aerosols = strTemplate(options.aer_template,nymd=nymd,nhms=nhms)
    if not options.skipVLIDORT:
        rank = len(omi.sample.delp.shape)
        if rank > 2 :
           nt, nr, nz = omi.sample.delp.shape 
           raise IOError, "rank > 2 to do - > code not adapted for this case"
        else :
           nobs, nz = omi.sample.delp.shape 

        nch = len(omi.channels)
        km = nz              
        lons = omi.lon.ravel()
        lats = omi.lat.ravel()

        zs = getTopo(options.topography,lons,lats)

        pe, ze, te = getEdgeVars(omi.sample)
    
        Igood_ = omi.sample.iGood.ravel()    # omi.sample.iGood -> only True values after \
                                             # omi_sampleReduce -> shape (nt,nr) -> 1D array 

        zs_=zs[Igood_]                       # surface height only for good values in npz file         
        he = ze + tile(zs_,(km+1,1))         # height above sea level


#   Loop over observations in batches of mobs observations if vector
#   ------------------------------------------------------
        radiance_ = ones((nobs,nch))
        ai_ = ones(nobs)
        g5_aaod=ones(nobs)
        if not options.scalar : # vector version
           nMom = 300  
           mobs = 200
           for i in range(0, nobs, mobs):

              if options.verbose:
                t = clock() - t0
                print "- Processing batch number %d"%i
                print "  [] Reading and interpolating IOPs at t=%f"%t

              I = range(i,min(i+mobs,nobs))

              tau_, ssa_, g_, pmom_=  getAOPvector(omi.sample,omi.channels,\
                                              I=I,nMom=nMom,rcfile='Aod_EOS.rc')
              aaod = tau_*(1-ssa_)
	      aaod_ = aaod.sum(axis=0)[1] # 388 nm

              # decrease SSA at only 354 nm
              # --------------------------
#              print 'ssa avant ', ssa_.shape, ssa_[45:72,0,0:2],ssa_[45:72,1,0:2] 
#              ssa_[:,0,:] = ssa_[:,0,:]*0.95
#              print 'ssa apres ', ssa_.shape, ssa_[45:72,0,0:2],ssa_[45:72,1,0:2]

              pe_=pe[:,I].copy()
              he_=he[:,I].copy()
              te_=te[:,I].copy()
            
              if options.verbose:
                t = clock() - t0
                print "  [] Doing RT calculation with VLIDORT at t=%f"%t

              radiance, ai,reflectivity=vlidort_vector(omi,tau_, ssa_, g_,pmom_,\
                                                 pe_, he_, te_,I,omi.sample.iGood)
              
              radiance_[I,:] = radiance.astype('float32')
              ai_[I]=ai.astype('float32')             
              g5_aaod[I] = aaod_.astype('float32') 

        else : # scalar version
             
           tau, ssa, g= getAOPscalar(omi.sample,omi.channels)
           radiance, ai,reflectivity_=vlidort_scalar(omi,tau, ssa, g, pe, he, te,omi.sample.iGood)
           radiance_[:,:] = radiance.astype('float32')
           ai_[:]=ai.astype('float32')

    print 'ai', ai_.shape, 'ssa', ssa_.shape

#   Results :
#   ---------
    iGood= omi.sample.iGood
    
#    results=dict(GEOS5_rad=radiance_, GEOS5_ai=ai_, OMI_rad=omi.radiance[iGood], 
#                 OMI_ai=omi.ai[iGood],OMI_aaod=omi.aaod[iGood],lat=omi.lat[iGood],
#                 lon=omi.lon[iGood],iGood=omi.sample.iGood)

#   Write output file : ungridded data 
#   -------------------------
#    savez(out_file,**results)
    savez(out_file,channels=omi.channels,GEOS5_rad=radiance_, GEOS5_ai=ai_,\
          aaod_g5388 = g5_aaod,\
          OMI_rad=omi.radiance[iGood],\
          OMI_ai=omi.ai[iGood],aaod_OMI388=omi.aaod[iGood][:,1],OMI_aod=omi.aod[iGood],\
          lat=omi.lat[iGood],\
          lon=omi.lon[iGood],iGood=omi.sample.iGood)  

   
    os.system('killme.csh %d &' %os.getpid())

#   Write ungridded data
#   --------------------
#    name, ext = os.path.splitext(out_file)
#    npz_file = name + '.npz'
#    omi.writeu(npz_file)
    
