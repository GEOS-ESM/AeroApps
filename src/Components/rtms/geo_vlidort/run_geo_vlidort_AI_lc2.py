#!/usr/bin/env python3
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to parse Prep config file and 
  call geo_vlidort_lc2.py
  
  May 2017
  patricia.castellanos@nasa.gov
"""

from os       import system
from optparse import OptionParser
from MAPL     import Config

if __name__ == "__main__":
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog prep_config_file start_isotime end_isotime")
    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry run.")

    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        prep_config, start_isotime, end_isotime = args
    else:
        parser.error("must have 3 arguments: prep_config_filename start_isotime end_isotime")


    # Parse prep config
    # -----------------
    cf = Config(prep_config,delim=' = ')
 
    Options =           " --instname=" + cf('LC2_INSTNAME')      + \
                  " --version_string=" + cf('LC2_VERSION')       + \
                        " --channels=" + cf('LC2_CHANNELS')     + \
                         " --nodemax=" + cf('LC2_NODEMAX')      + \
                          " --layout=" + cf('LC2_LAYOUT')       + \
                          " --cldmax=" + cf('LC2_CLDMAX')       + \
                             " --dir=" + cf('LC2_DIR')          + \
                         " --runmode=" + cf('LC2_RUNMODE')       + \
                         " --runfile=" + cf('LC2_RUNFILE')      + \
                          " --execfile=" + cf('LC2_EXECFILE') 

    # Optional
    try:
        Options += " --c_band=" + cf('LC2_SURF_C_BAND')                      
    except:
        pass

    if cf('LC2_PROFILE').upper()     == 'YES': Options += " --profile"
    if cf('LC2_VERBOSE').upper()     == 'YES': Options += " --verbose"
    if cf('LC2_KEEP_LB').upper()     == 'YES': Options += " --keep_lb"
    if cf('LC2_ARCHIVE_LB').upper()  == 'YES': Options += " --archive_lb"

    # Generate LC2 products
    # ---------------------
    cmd = "geo_vlidort_AI_lc2.py %s %s %s "%(Options,start_isotime,end_isotime)
    print(cmd)
    if not options.dryrun:
        if system(cmd):
            raise ValueError("geo_vlidort_AI_lc2.py failed for %s to %s "%(start_isotime, end_isotime))

    




