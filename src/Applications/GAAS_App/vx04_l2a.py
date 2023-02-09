#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
  A Python script to create NNR retrievals.
  It uses class VX04 to directly read VIIRS Aerosol Level 2 
  Files (AERDB and AERDT).

  This utility reads VIIRS Level2 files and creates an ODS file with
  NNR retrievals, as well as a *gritas* type gridded output.
  
  2023 Based on mxd04_l2a.py
  patricia.castellanos@nasa.gov
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
import sys
import subprocess

from time            import clock
from optparse        import OptionParser   # Command-line args  
from dateutil.parser import parse as isoparse
from vx04_nnr       import Vx04_NNR
from MAPL            import strTemplate

Ident = dict( vsnppdto = ('SNPP','dt_ocean'),
              vsnppdtl = ('SNPP','dt_land'),
              vsnppdbo = ('SNPP','db_ocean'),
              vsnppdbl = ('SNPP','db_land')
            )

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

    expid = 'nnr_001'
    ident = 'vsnppdbl'
    
#   Defaults may be platform dependent
#   ----------------------------------
    if os.path.exists('/nobackup/VIIRS/Level2/'): # New calculon
        l2_path = '/nobackup/VIIRS/'
        out_dir = '/nobackup/NNR/VIIRS/%coll/Level%lev/%prod/Y%y4/M%m2'
        nn_file = '/nobackup/NNR/Net/VIIRS/nnr_003.%ident_Tau.net'
        blank_ods = '/nobackup/NNR/Misc/blank.ods'
        aer_x   = '/nobackup/NNR/Misc/tavg1_2d_aer_Nx'
    else: # Must be somewhere else, no good defaults
        out_dir      = './'
        l2_path = './'
        nn_file = '%ident_Tau.net'
        blank_ods = 'blank.ods'
        aer_x   = 'tavg1_2d_aer_Nx'        

    out_tmpl = '%s.%prod_L%leva.%algo.%y4%m2%d2_%h2%n2z.%ext'
    coll = '002'
    res = 'c'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] ident isotime",
                          version='vx04_l2a-1.0.0' )


    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="Experiment id (default=%s)"\
                           %expid )

    parser.add_option("-d", "--dir", dest="out_dir", default=out_dir,
                      help="Output directory (default=%s)"\
                           %out_dir )

    parser.add_option("-A", "--aer_x", dest="aer_x", default=aer_x,
                      help="GrADS ctl for speciated AOD file (default=%s)"\
                           %aer_x )

    parser.add_option("-B", "--blank_ods", dest="blank_ods", default=blank_ods,
                      help="Blank ODS file name for fillers  (default=%s)"\
                           %blank_ods )

    parser.add_option("-C", "--collection", dest="coll", default=coll,
                      help="VIIRS collection (default=%s)"\
                           %coll )

    parser.add_option("-o", "--fname", dest="out_tmpl", default=out_tmpl,
                      help="Output file name template (default=%s); ODS file name will be derived from it by changing extension to '.ods' and replacing 'Level3' with 'Level2'."\
                           %out_tmpl )

    parser.add_option("-L", "--l2_dir", dest="l2_path", default=l2_path,
                      help="Top directory for VIIRS Level 2 files (default=%s)"\
                           %l2_path )

    parser.add_option("-N", "--net", dest="nn_file", default=nn_file,
                      help="Neural net file template  (default=%s)"\
                           %nn_file )

    parser.add_option("-r", "--res", dest="res", default=res,
                      help="Resolution for gridded output (default=%s)"\
                           %out_tmpl )

    parser.add_option("-u", "--uncompressed",
                      action="store_true", dest="uncompressed",default=False,
                      help="Do not use n4zip to compress gridded/ODS output file (default=False)")

    parser.add_option("-F", "--force",
                      action="store_true", dest="force",default=False,
                      help="Overwrites output file")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",default=False,
                      help="Turn on verbosity.")

    parser.add_option("--writenpz", dest="writenpz", default=False,
                      help="Write an ungridded npz file in addition to ODS and gridded files  (default=False)")    

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        ident, isotime = args
        sat, algo = Ident[ident]
        prod = sat + '04'
    else:
        parser.error("must have 3 arguments: ident, date and time")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
    if options.verbose:
        print ""
        print "                          VIIRS Level 2A Processing"
        print "                          -------------------------"
        print ""
        t0 = clock()

#   Time variables
#   --------------
    syn_time = isoparse(isotime)
    nymd     = str(syn_time.date()).replace('-','')
    nhms     = str(syn_time.time()).replace(':','')
            
#   Form output gridded file name
#   -----------------------------
    out_tmpl = options.out_dir+'/'+options.out_tmpl
    out_tmpl = out_tmpl.replace('%coll',options.coll).replace('%prod',prod).replace('%algo',algo).replace('%lev','3').replace('%ext','nc4')
    out_file = strTemplate(out_tmpl,expid=options.expid,nymd=nymd,nhms=nhms)
    name, ext = os.path.splitext(out_file)
    if os.path.exists(out_file) and (options.force is not True):
        print "vx04_l2a: Output Gridded file <%s> exists --- cannot proceed."%out_file
        raise IOError, "Specify --force to overwrite existing output file."    
    if os.path.exists(out_file) and options.force:
        os.remove(out_file)    

#   Form ODS file name
#   ------------------
    ods_tmpl = options.out_dir+'/'+options.out_tmpl
    ods_tmpl = ods_tmpl.replace('%coll',options.coll).replace('%prod',prod).replace('%algo',algo).replace('%lev','2').replace('%ext','ods')
    ods_file = strTemplate(ods_tmpl,expid=options.expid,nymd=nymd,nhms=nhms)
    if os.path.exists(ods_file) and (options.force is not True):
        print "vxd04_l2a: Output ODS file <%s> exists --- cannot proceed."%ods_file
        raise IOError, "Specify --force to overwrite existing output file."
    if os.path.exists(ods_file) and options.force:
        os.remove(ods_file)

#   Aerosol composition file name
#   -----------------------------
    if options.aer_x[-3:] == 'nc4':
      aer_x = strTemplate(options.aer_x,expid=options.expid,nymd=nymd,nhms=nhms)
    else:
      aer_x = options.aer_x

        
#   VIIRS Level 2 NNR Aerosol Retrievals
#   ------------------------------------
    if options.verbose:
        print "NNR Retrieving %s %s on "%(sat,algo.upper()),syn_time

    viirs = Vx04_NNR(options.l2_path,sat,algo.upper(),syn_time,aer_x,
                      coll=options.coll,
                      cloud_thresh=0.7,
                      cloudFree = 0.0,
                      aodmax = 1.0,
                      verbose=options.verbose)
    if viirs.nobs < 1:
        if options.verbose:
            print 'WARNING: no observation for this time in file <%s>'%ods_file
    
    elif any(viirs.iGood) == False:
        if options.verbose:
            print 'WARNING: no GOOD observation for this time in file <%s>'%ods_file
        viirs.nobs = 0

    nn_file = options.nn_file.replace('%ident',ident)
    viirs.apply(nn_file)

#   Write ODS
#   ---------
    makethis_dir(ods_file)
    if viirs.nobs>0:
        viirs.writeODS(ods_file,revised=True)
    else:
        if os.system('ods_blank.x %s %s %s %s'%(options.blank_ods,nymd,nhms,ods_file)):
            warnings.warn('cannot create empty output file <%s>'%ods_file)
        else:
            if options.verbose:
                print "[w] Wrote empty ODS file "+ods_file

#   Write gridded output file (revised channels only)
#   -------------------------------------------------
    makethis_dir(out_file)
    if viirs.nobs>0:
      if str.isdigit(options.res):
        viirs.writeg(filename=out_file,refine=int(options.res),channels=viirs.channels_)
      else:
        viirs.writeg(filename=out_file,res=options.res,channels=viirs.channels_)

#   Write ungridded data
#   --------------------
    if options.writenpz:
        name, ext = os.path.splitext(out_file)
        npz_file = name.replace('Level3','Level2') + '.npz'
        makethis_dir(npz_file)
        viirs.write(npz_file)
    
#   Compress nc output unless the user disabled it
#   ----------------------------------------------
    if viirs.nobs>0:
        if not options.uncompressed:
            if subprocess.call("n4zip " + out_file,shell=True):
                warnings.warn('cannot compress output file <%s>'%out_file)
            if subprocess.call("n4zip " + ods_file,shell=True):
                warnings.warn('cannot compress output ods file <%s>'%ods_file)
