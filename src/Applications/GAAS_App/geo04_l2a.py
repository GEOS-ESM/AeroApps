#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
  A Python script to create ODS file retrievals for geostationnary.
  It now uses class GEO04 to directly read GOES and HIMAWARI Aerosol Level 2 
  Files.
  Adapted from MODIS
  March 2021, virginie buchard
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
from MAPL            import strTemplate
from pyobs.geo04     import GEO04_L2, granules

Ident = dict( g16o = ('G16','ocean'),
              g16l = ('G16','land'),
              g17o = ('G17','ocean'),
              g17l = ('G17','land'),
              himo = ('HIM','ocean'),
              himl = ('HIM','land')
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

    expid = 'aerdt'
    out_dir = './'
    sat = 'G16'
    l2_path = '/css/geostationary'
    blank_ods = 'blank_syn8.ods'

    out_tmpl = '%s.%sat_L%leva.%algo.%y4%m2%d2_%h2%n2z.%ext'
    res = 'c'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] ident isotime",
                          version='geo04_l2a-1.0.0' )


    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="Experiment id (default=%s)"\
                           %expid )

    parser.add_option("-d", "--dir", dest="out_dir", default=out_dir,
                      help="Output directory (default=%s)"\
                           %out_dir )

#    parser.add_option("-A", "--aer_x", dest="aer_x", default=aer_x,
#                      help="GrADS ctl for speciated AOD file (default=%s)"\
#                           %aer_x )

    parser.add_option("-B", "--blank_ods", dest="blank_ods", default=blank_ods,
                      help="Blank ODS file name for fillers  (default=%s)"\
                           %blank_ods )

    parser.add_option("-S", "--sat", dest="sat", default=sat,
                      help="Satellite name: G16, G17 or HIM  (default=%s)"\
                           %sat )
    
    parser.add_option("-o", "--fname", dest="out_tmpl", default=out_tmpl,
                      help="Output file name template (default=%s); ODS file name will be derived from it by changing extension to '.ods' and replacing 'Level3' with 'Level2'."\
                           %out_tmpl )

    parser.add_option("-L", "--l2_dir", dest="l2_path", default=l2_path,
                      help="Top directory for ABI/AHI Level 2 files (default=%s)"\
                           %l2_path )

#    parser.add_option("-N", "--net", dest="nn_file", default=nn_file,
#                      help="Neural net file template  (default=%s)"\
#                           %nn_file )

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

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        ident, isotime = args
        sat, algo = Ident[ident]
    else:
        parser.error("must have 2 arguments: ident, and isotime")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
    if options.verbose:
        print ""
        print "                          GEO ABI/AHI Level 2A Processing"
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
    out_tmpl = out_tmpl.replace('%sat',options.sat).replace('%algo',algo).replace('%lev','3').replace('%ext','nc4')
    out_file = strTemplate(out_tmpl,expid=options.expid,nymd=nymd,nhms=nhms)
    name, ext = os.path.splitext(out_file)
    if os.path.exists(out_file) and (options.force is not True):
        print "geo04_l2a: Output Gridded file <%s> exists --- cannot proceed."%out_file
        raise IOError, "Specify --force to overwrite existing output file."    
    if os.path.exists(out_file) and options.force:
        os.remove(out_file)    

#   Form ODS file name
#   ------------------
    ods_tmpl = options.out_dir+'/'+options.out_tmpl
    ods_tmpl = ods_tmpl.replace('%sat',options.sat).replace('%algo',algo).replace('%lev','2').replace('%ext','ods')
    ods_file = strTemplate(ods_tmpl,expid=options.expid,nymd=nymd,nhms=nhms)
    if os.path.exists(ods_file) and (options.force is not True):
        print "geo04_l2a: Output ODS file <%s> exists --- cannot proceed."%ods_file
        raise IOError, "Specify --force to overwrite existing output file."
    if os.path.exists(ods_file) and options.force:
        os.remove(ods_file)

#   Aerosol composition file name (vb: no nedd for now)
#   -----------------------------
#    if options.aer_x[-3:] == 'nc4':
#      aer_x = strTemplate(options.aer_x,expid=options.expid,nymd=nymd,nhms=nhms)
#    else:
#      aer_x = options.aer_x

        
#   GOES/HIM Level 2 NNR Aerosol Retrievals skip for now...........
#   ------------------------------------
#    if options.verbose:
#        print "NNR Retrieving %s %s on "%(prod,algo.upper()),syn_time

#    geo = GEO04_NNR(options.l2_path,prod,algo.upper(),syn_time,aer_x,
#                      coll=options.coll,
#                      cloud_thresh=0.7,
#                      cloudFree = 0.0,
#                      aodmax = 1.0,
#                      verbose=options.verbose)
#    if modis.nobs < 1:
#        if options.verbose:
#            print 'WARNING: no observation for this time in file <%s>'%ods_file
    
#    elif any(modis.iGood) == False:
#        if options.verbose:
#            print 'WARNING: no GOOD observation for this time in file <%s>'%ods_file
#        geo.nobs = 0

#    nn_file = options.nn_file.replace('%ident',ident)
#    geo.apply(nn_file)



#   READ GOES/HIM Level 2 files using geo04.py under pyobs
#   --------------------------------------------
    Files  =  granules(options.l2_path, sat, syn_time, nsyn = 8)
    geo = GEO04_L2(Files, options.sat, algo.upper(), syn_time, nsyn =8,\
                Verb=options.verbose)
    # Filter obs (see in geo04.py for more details)
    #-------
    geo.filter()
    
#   Write ODS
#   ---------
    makethis_dir(ods_file)
    if geo.nobs>0:
        geo.writeODS(ods_file,revised=False,channels=[550,])
    else:
        if os.system('ods_blank.x %s %s %s %s'%(options.blank_ods,nymd,nhms,ods_file)):
            warnings.warn('cannot create empty output file <%s>'%ods_file)
        else:
            if options.verbose:
                print "[w] Wrote empty ODS file "+ods_file

#   Write gridded output file (revised channels only)
#   -------------------------------------------------
    makethis_dir(out_file)
    if geo.nobs>0:
      if str.isdigit(options.res):
        geo.writeg(filename=out_file,refine=int(options.res),channels=geo.aChannels)
      else:
        geo.writeg(filename=out_file,res=options.res,channels=geo.aChannels)

#   Write ungridded data
#   --------------------
#    name, ext = os.path.splitext(out_file)
#    npz_file = name.replace('Level3','Level2') + '.npz'
#    makethis_dir(npz_file)
#    modis.write(npz_file)
    
#   Compress nc output unless the user disabled it
#   ----------------------------------------------
    if geo.nobs>0:
        if not options.uncompressed:
            if subprocess.call("n4zip " + out_file,shell=True):
                warnings.warn('cannot compress output file <%s>'%out_file)
            if subprocess.call("n4zip_ods " + ods_file,shell=True):
                warnings.warn('cannot compress output ods file <%s>'%ods_file)
