#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
  A Python script to create VLIDORT/OMI Level 3a files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
import sys

from time       import clock
from optparse   import OptionParser   # Command-line args  
from omi        import OMI
from mie        import getMieScal, getMieVect, getTopo
from MAPL       import strTemplate

#---------------------------------------------------------------------
def makethis_dir(filename):
    """Creates the relevant directory if necessary."""
    path, filen = os.path.split(filename)
    if path != '':
        rc = os.system('mkdir -p '+path)
        if rc:
            raise IOError, "could not create directory "+path
        
#---------------------------------------------------------------------

if __name__ == "__main__":

    expid        = 'vlidort'

#   Defaults may be platform dependent
#   ----------------------------------
    if os.path.exists('/nobackup/OMI/Level2/ODS/'): # New calculon
        ods_template = '/nobackup/OMI/Level2/ODS/Y2008/M06/omi.aero_tc8.obs.%y4%m2%d2.ods'
        aer_template = '/nobackup/1/ARCTAS/Y%y4/M%m2/d5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc'
        topography   = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    elif os.path.exists('/discover/nobackup/projects/gmao/iesa/'): # Discover
        ods_template = '/discover/nobackup/projects/gmao/iesa/aerosol/data/OMI/Level2/ODS/Y%y4/M%m2/omi.aero_tc8.obs.%y4%m2%d2.ods'
        aer_template = '/discover/nobackup/projects/gmao/iesa/aerosol/data/ARCTAS/Y%y4/M%m2/D%d2/d5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc4'
        topography   = '/discover/nobackup/projects/gmao/iesa/aerosol/data/VLIDORT/topo/topography.1152x721.nc'
    else: # Must be somewhere else, no good defaults
        ods_template = 'omi.aero_tc8.obs.%y4%m2%d2.ods'
        aer_template = 'd5_arctas_02.inst3d_aer_v.%y4%m2%d2_%h2%n2z.nc'
        topography   = 'topography.1152x721.nc'

    out_dir      = './'
    out_fname    = '%s.omi_l3a.%y4%m2%d2_%h2%n2z.nc4'

    scalar       = False
    res          = 'c'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] nymd nhms",
                          version='omi_level3a-1.0.0' )

    parser.add_option("-A", "--aero", dest="aer_template", default=aer_template,
                      help="File name template for input Aerosol concentrations (default=%s)"\
                           %aer_template )

    parser.add_option("-O", "--ods", dest="ods_template", default=ods_template,
                      help="File name template for input ODS (default=%s)"\
                           %ods_template )

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
    
    if len(args) == 2:
        nymd, nhms = args
    else:
        parser.error("must have 2 arguments: date and time")
        
    if options.verbose:
        print ""
        print "                          VLIDORT/OMI Level 3A Processing"
        print "                          -------------------------------"
        print ""
        t0 = clock()
            
    if options.skipVLIDORT:
        print "IMPORTANT: Skipping VLIDORT"

#   Output gridded file
#   -------------------
    out_file = strTemplate(options.out_dir+'/'+options.out_fname,
                           expid=options.expid,nymd=nymd,nhms=nhms)
    if os.path.exists(out_file) and (options.force is not True):
        print "omi_l3a: Output file <%s> exists --- cannot proceed."%out_file
        raise IOError, "Specify --force to overwrite existing output file."
    else:
        makethis_dir(out_file) # make sure directory exists

#   Read OMI data
#   -------------
    if options.verbose:
        t = clock() - t0
        print "- Reading OMI measurements at t=%f"%t
    omi_file = strTemplate(options.ods_template,nymd=nymd,nhms=nhms)
    omi = OMI(omi_file,nymd,nhms)

#   Compute Mie parameters at the OMI locations
#   -------------------------------------------
    aerosols = strTemplate(options.aer_template,nymd=nymd,nhms=nhms)

#   Loop over observations in batches of mobs observations
#   ------------------------------------------------------
    if not options.skipVLIDORT:
        mobs = 100
        for i in range(0, omi.nobs, mobs):

            if options.verbose:
                t = clock() - t0
                print "- Processing batch number %d"%i
                print "  [] Reading and interpolating IOPs at t=%f"%t

            I = range(i,min(i+mobs,omi.nobs))

            nMom, nPol, tau, ssa, \
            pe, ze, te, g, pmom =  getMieVect(aerosols,omi.channels,\
                                              omi.lon[I],omi.lat[I],nymd,nhms)

            zs = getTopo(options.topography,omi.lon[I],omi.lat[I])
            he = ze + zs.reshape((1,)+zs.shape)

            if options.verbose:
                t = clock() - t0
                print "  [] Doing RT calculation with VLIDORT at t=%f"%t

            omi.vlidort(nMom,nPol,tau,ssa,g,pmom,pe,he,te,verbose=2,
                        scalar=options.scalar, I=I)

#   Write gridded output file
#   -------------------------
    omi.writeg(filename=out_file,res=options.res)

#   Compress it unless the user disabled it
#   ---------------------------------------
    if not options.uncompressed:
        rc = os.system("n4zip "+out_file)
        if rc:
            warnings.warn('cannot compress output file <%s>'%out_file)

#   Write ungridded data
#   --------------------
    name, ext = os.path.splitext(out_file)
    npz_file = name + '.npz'
    omi.writeu(npz_file)
    





