#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to parse Prep config file and create ODS with unbiases AOD data.
  
  February 2011.
  arlindo.dasilva@nasa.gov
"""

from os       import system
from optparse import OptionParser
from MAPL     import Config

if __name__ == "__main__":
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog prep_config_file nymd nhms",
                          version='aod_l2a-1.0.0' )
    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry run.")

    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        prep_config, nymd, nhms = args
    else:
        parser.error("must have 3 arguments: prep_config_filename nymd nhms")

    # Parse prep config
    # -----------------
    cf = Config(prep_config,delim=' = ')
 
    Options =     " --expid=" + cf('AOD_L2A_EXPID')        + \
                 " --l2_dir=" + cf('AOD_L2A_L2_DIR')       + \
                     " --res=" + cf('AOD_L2A_RESOLUTION')   + \
                   "  --dir=" + cf('AOD_L2A_OUT_DIR')      + \
                  " --fname=" + cf('AOD_L2A_OUT_TEMPLATE') + \
                    " --net=" + cf('AOD_L2A_NN_FILE')      + \
                 " --albedo=" + cf('AOD_L2A_ALBEDO_FILE')  + \
                   " --wind=" + cf('AOD_L2A_WIND_FILE')    + \
              " --blank_ods=" + cf('AOD_L2A_BLANK_ODS')

    if cf('AOD_L2A_OVERWRITE').upper() == 'YES': Options += " --force"
    if   cf('AOD_L2A_VERBOSE').upper() == 'YES': Options += " -v"
              
    # Generate products
    # -----------------
    i = 0
    Coll = cf('AOD_L2A_COLLECTION').split(',')
    for ident in cf('AOD_L2A_IDENTS').split(','):
        coll = Coll[i] 
        cmd = "mxd04_l2a.py %s --collection=%s %s %s %s"%(Options,Coll[i],ident,nymd,nhms)
        print cmd
        if not options.dryrun:
            if system(cmd):
                raise ValueError, "mxd04.py failed for %s on %s %s"%(ident,nymd,nhms)

        i += 1
    

