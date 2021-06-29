#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to parse Prep config file and create ODS with GEOSTAIONARY aerosol retrievals.
  adapted from MODIS
  March 2021.
  virginie.buchard@nasa.gov
"""

from os       import system
from optparse import OptionParser
from MAPL     import Config

if __name__ == "__main__":
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog prep_config_file isotime",
                          version='aerdt_l2a-1.0.0' )
    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry run.")

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        prep_config, isotime = args
    else:
        parser.error("must have 2 arguments: prep_config_filename isotime")

    # Parse prep config
    # -----------------
    cf = Config(prep_config,delim=' = ')
 
    Options =     " --expid=" + cf('GEO_L2A_EXPID')        + \
                 " --l2_dir=" + cf('GEO_L2A_L2_DIR')       + \
                    " --res=" + cf('GEO_L2A_RESOLUTION')   + \
                   "  --dir=" + cf('GEO_L2A_OUT_DIR')      + \
                  " --fname=" + cf('GEO_L2A_OUT_TEMPLATE') + \
              " --blank_ods=" + cf('GEO_L2A_BLANK_ODS')   

    if cf('GEO_L2A_OVERWRITE').upper() == 'YES': Options += " --force"
    if cf('GEO_L2A_VERBOSE').upper() == 'YES': Options += " -v"
    if cf('GEO_L2A_UNCOMPRESSED').upper() == 'YES': Options += " -u "
              
    # Generate products
    # -----------------
    i = 0
    SAT = cf('GEO_L2A_SAT').split(',')
    for ident in cf('GEO_L2A_IDENTS').split(','):
        
        cmd = "geo04_l2a.py %s --sat=%s %s %s "%(Options,SAT[i],ident,isotime)
        if not options.dryrun:
            if system(cmd):
                raise ValueError, "geo04_l2a.py failed for %s on %s "%(ident,isotime)

        i += 1
    
 
