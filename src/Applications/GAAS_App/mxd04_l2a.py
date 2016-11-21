#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
  A Python script to create NNR retrievals.
  It now uses class MXD04 to directly read MODIS Aerosol Level 2 
  Files (MOD04/MYD04).

  This utility reads MODIS Level2 files and creates an ODS file with
  NNR retrievals, as well as a *gritas* type gridded output.
  
  February 2011, revised Novembre 2016 for MODIS Collection 6.
  arlindo.dasilva@nasa.gov
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
import sys

from time         import clock
from optparse     import OptionParser   # Command-line args  
from datetime     import datetime

from mxd04_nnr    import MxD04_NNR
from MAPL         import strTemplate
from grads        import GrADS

Ident = dict( modo = ('MOD04','ocean'),
              modl = ('MOD04','land'),
              modd = ('MOD04','deep'),
              mydo = ('MYD04','ocean'),
              mydl = ('MYD04','land'),
              mydd = ('MYD04','deep')
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

    expid = 'nnr3'
    ident = 'modo'
    
#   Defaults may be platform dependent
#   ----------------------------------
    if os.path.exists('/nobackup/MODIS/Level2/'): # New calculon
        l2_path = '/nobackup/MODIS/Level2/'
        out_dir = '/nobackup/NNR/Level%lev/%prod/Y%y4/M%m2'
        nn_file = '/nobackup/NNR/Net/nnr_002.%ident_Tau.net'
        blank_ods = '/nobackup/NNR/Misc/blank.ods'
        coxmunk_lut = '/nobackup/NNR/Misc/coxmunk_lut.npz'
#    elif os.path.exists('/discover/nobackup/projects/gmao/iesa/'): # Discover
#        raise ValueError, 'not setup yet'
    else: # Must be somewhere else, no good defaults
        out_dir      = './'
        l2_path = './'
        nn_file = '%ident_Tau.net'
        blank_ods = 'blank.ods'
        coxmunk_lut = 'cox-munk_lut.npz'

    out_tmpl = '%s.%prod_l%leva.%algo.%y4%m2%d2_%h2%n2z.%ext'
    wind_file = 'merra_slv-hourly.ddf'
    albedo_file = 'albedo_clim.ctl'
    coll = '006'
    res = 'e'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] ident nymd nhms",
                          version='mxd04_l2a-1.0.0' )


    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="Experiment id (default=%s)"\
                           %expid )

    parser.add_option("-d", "--dir", dest="out_dir", default=out_dir,
                      help="Output directory (default=%s)"\
                           %out_dir )

    parser.add_option("-A", "--albedo", dest="albedo_file", default=albedo_file,
                      help="GrADS ctl for Albedo file (default=%s)"\
                           %albedo_file )

    parser.add_option("-B", "--blank_ods", dest="blank_ods", default=blank_ods,
                      help="Blank ODS file name for fillers  (default=%s)"\
                           %blank_ods )

    parser.add_option("-C", "--collection", dest="coll", default=coll,
                      help="MODIS collection (default=%s)"\
                           %coll )

    parser.add_option("-o", "--fname", dest="out_tmpl", default=out_tmpl,
                      help="Output file name template (default=%s); ODS file name will be derived from it by changing extension to '.ods' and replacing 'Level3' with 'Level2'."\
                           %out_tmpl )

    parser.add_option("-L", "--l2_dir", dest="l2_path", default=l2_path,
                      help="Top directory for MODIS Level 2 files (default=%s)"\
                           %l2_path )
    parser.add_option("-M", "--coxmunk", dest="coxmunk_lut", default=coxmunk_lut,
                      help="Blank ODS file name for fillers  (default=%s)"\
                           %coxmunk_lut )

    parser.add_option("-N", "--net", dest="nn_file", default=nn_file,
                      help="Neural net file template  (default=%s)"\
                           %nn_file )

    parser.add_option("-W", "--wind", dest="wind_file", default=wind_file,
                      help="GrADS ctl for Wind file (default=%s)"\
                           %wind_file )

    parser.add_option("-r", "--res", dest="res", default=res,
                      help="Resolution for gridded output (default=%s)"\
                           %out_tmpl )

    parser.add_option("-u", "--uncompressed",
                      action="store_true", dest="uncompressed",
                      help="Do not use n4zip to compress gridded/ODS output file (default=False)")

    parser.add_option("-F", "--force",
                      action="store_true", dest="force",
                      help="Overwrites output file")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Turn on verbosity.")

    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        ident, nymd, nhms = args
        prod, algo = Ident[ident]
    else:
        parser.error("must have 3 arguments: ident, date and time")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
    if options.verbose:
        print ""
        print "                          MODIS Level 2A Processing"
        print "                          -------------------------"
        print ""
        t0 = clock()

#   Time variables
#   --------------
    nymd_, nhms_ = ( int(nymd), int(nhms) )
    year, month, day = (nymd_/10000, (nymd_%10000)/100, nymd_%100)
    hour = nhms_/10000
    syn_time = datetime(year,month,day,hour,0,0) # no munte, second
            
#   Form output gridded file name
#   -----------------------------
    out_tmpl = options.out_dir+'/'+options.out_tmpl
    out_tmpl = out_tmpl.replace('%prod',prod).replace('%algo',algo).replace('%lev','3').replace('%ext','nc4')
    out_file = strTemplate(out_tmpl,expid=options.expid,nymd=nymd,nhms=nhms)
    name, ext = os.path.splitext(out_file)

#   Form ODS file name
#   ------------------
    ods_tmpl = options.out_dir+'/'+options.out_tmpl
    ods_tmpl = ods_tmpl.replace('%prod',prod).replace('%algo',algo).replace('%lev','2').replace('%ext','ods')
    ods_file = strTemplate(ods_tmpl,expid=options.expid,nymd=nymd,nhms=nhms)
    if os.path.exists(ods_file) and (options.force is not True):
        print "mxd04_l2a: Output ODS file <%s> exists --- cannot proceed."%ods_file
        raise IOError, "Specify --force to overwrite existing output file."

#   Gather Auxiliary data
#   ---------------------
    ga = GrADS(Echo=False,Window=False)
    if algo == 'ocean':
        if options.wind_file[-3:] == 'nc4':
            wind_file = strTemplate(options.wind_file,expid=options.expid,nymd=nymd,nhms=nhms)
            ga('sdfopen %s'%wind_file)
        else:
            ga('open %s'%options.wind_file)
        expr='mag(u10m,v10m)'
        vname = 'wind'
    elif (algo == 'land') or (algo == 'deep'):
        ga('open %s'%options.albedo_file)
        expr='albedo'
        vname = 'albedo'
    else:
        raise ValueError, 'unknown algo <%s>'%algo
        
#   MODIS Level 2 NNR Aerosol Retrievals
#   ------------------------------------
    if options.verbose:
        print "NNR Retrieving %s %s on "%(prod,algo.upper()),syn_time

    modis = MxD04_NNR(options.l2_path,prod,algo.upper(),syn_time,
                      ga,expr=expr,vname=vname,coll=options.coll,
                      cloud_thresh=0.7,coxmunk_lut=options.coxmunk_lut,
                      verbose=options.verbose)
    if modis.nobs < 1:
        if options.verbose:
            print 'WARNING: no observation for this time in file <%s>'%ods_file
    
    elif any(modis.iGood) == False:
        if options.verbose:
            print 'WARNING: no GOOD observation for this time in file <%s>'%ods_file
        modis.nobs = 0
    
    nn_file = options.nn_file.replace('%ident',ident)
    modis.apply(nn_file)

#   Write ODS
#   ---------
    makethis_dir(ods_file)
    if modis.nobs>0:
        modis.writeODS(ods_file,revised=True)
    else:
        if os.system('ods_blank.x %s %s %s %s'%(options.blank_ods,nymd,nhms,ods_file)):
            warnings.warn('cannot create empty output file <%s>'%ods_file)
        else:
            if options.verbose:
                print "[w] Wrote empty ODS file "+ods_file

#   Write gridded output file (revised channels only)
#   -------------------------------------------------
    makethis_dir(out_file)
    if modis.nobs>0:
        modis.writeg(filename=out_file,res=options.res,channels=modis.channels_)

#   Write ungridded data
#   --------------------
#    name, ext = os.path.splitext(out_file)
#    npz_file = name.replace('Level3','Level2') + '.npz'
#    makethis_dir(npz_file)
#    modis.write(npz_file)

#   Compress nc output unless the user disabled it
#   ----------------------------------------------
    if modis.nobs>0:
        if not options.uncompressed:
            if os.system("n4zip "+out_file):
                warnings.warn('cannot compress output file <%s>'%out_file)
            if os.system("n4zip "+ods_file):
                warnings.warn('cannot compress output ods file <%s>'%ods_file)
