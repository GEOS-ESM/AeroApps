#!/usr/bin/env python3
#
# Sample GEOS-5 collections at OMI Level 2 orbits.
#

import os
from optparse   import OptionParser 
from pyobs      import OMAERUV_L2, aura
from numpy      import loadtxt
from mieobs     import aerToUpper 
#---------------------------------------------------------------------
def makethis_dir(filename):
    """Creates the relevant directory if necessary."""
    path, filen = os.path.split(filename)
    if path != '':
        rc = os.system('mkdir -p '+path)
        if rc:
            raise IOError("could not create directory "+path)

if __name__ == "__main__":

    expid   = 'ops'
    out_dir = '/discover/nobackup/vbuchard/VLIDORT/ENSEMBLE_exp/sample_dR_MERRA-AA-r2_ens10deg/%y4/%m2/%d2'
    file_dx_dy = 'delta_xy_list_10deg.dat'

    # Parse command line options
    # --------------------------
    parser = OptionParser(usage="Usage: %prog [options] omaeruv_filename g5_collection",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="npz_fn", default=None,
                      help="output NPZ file name (default=None)")
    
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="turn on verbosity.")

    parser.add_option("-x", "--expid", dest="expid", default=expid,
                       help="experiment id (default=%s)"\
                           %expid )

    parser.add_option("-F", "--Full", 
                      action="store_true",dest="Full_interp", 
                      help="Interpolation at all OMI obs location(default=%s)")

    parser.add_option("-d", "--dir", dest="out_dir", default=out_dir,
                      help="output directory (default=%s)"\
                           %out_dir )

    parser.add_option("-e", "--ensemble",action="store_true", dest="ensemble", default=None,
                      help="create an ensemble of 128 members (default=None)") 
         
    (options, args) = parser.parse_args()

    if len(args) == 2:
        omi_fn, g5_fn = args
    else:
        parser.error("must have 2 arguments: omaeruv_filename g5_collection")


    # Read OMI orbit
    # --------------
    omi = OMAERUV_L2(omi_fn)

    onlyVars = ['ps', 'delp', 'LWI', 'AIRDENS', 'RH',
                'du001', 'du002', 'du003', 'du004', 'du005',
                'ss001', 'ss002', 'ss003', 'ss004', 'ss005',
                'BCphobic', 'BCphilic', 'OCphobic', 'OCphilic', 'SO4']
#    onlyVars = ['PS', 'DELP', 'LWI', 'AIRDENS', 'RH',
#                'DU001', 'DU002', 'DU003', 'DU004', 'DU005',
#                'SS001', 'SS002', 'SS003', 'SS004', 'SS005',
#                'BCPHOBIC', 'BCPHILIC', 'OCPHOBIC', 'OCPHILIC', 'SO4']
    # Sample g5 collection at OMI orbit
    # ---------------------------------
    if options.npz_fn is None:
        instr, sat, lev, prod, vtime, orb,ver, ptime, t =aura.parseFilename(omi_fn)

    if options.Full_interp:
       options.npz_fn = options.out_dir+'/'+options.expid+'_aer_'+lev+ '-'+\
               prod + '_' + vtime +'Full'+'.npz'
       omi.sampleFile(g5_fn, onlyVars=onlyVars, npzFile=options.npz_fn, Verbose=options.verbose)
    else:
       if options.ensemble:   
           options.npz_fn = options.out_dir+'/'+options.expid+'_aer_'+lev+ '-'+\
           prod + '_' + vtime +'_obs.npz' # center point (OMI obs) 
           omi.sampleFile(g5_fn, onlyVars=onlyVars,npzFile=None,dx=None,dy=None, Verbose=options.verbose)
           Igood = omi.qa_flag <1         # OMI doc : FinalAlgorithmFlags - 0 - Most reliable( AAOD, SSA, AOD)
           omi.sampleReduce(I=Igood,npzFile=options.npz_fn)
def hold():
           data = loadtxt(file_dx_dy)         
           for i in range(len(data)) :         
              options.npz_fn = options.out_dir+'/'+options.expid+'_aer_'+lev+ '-'+\
              prod + '_' + vtime +'e'+('%03d')%i+'.npz'
              omi.sampleFile(g5_fn, onlyVars=onlyVars,npzFile=None,dx=data[:,0],dy=data[:,1], Verbose=options.verbose)
              Igood = omi.qa_flag <1         # OMI doc : FinalAlgorithmFlags - 0 - Most reliable( AAOD, SSA, AOD)
              omi.sampleReduce(I=Igood,npzFile=options.npz_fn)


