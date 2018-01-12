#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to parse Prep config file and 
  call geo_vlidort_cloud_lc2.py
  
  Aug 2017
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
 
    Options =           " --instname=" + cf('LC2_INSTNAME')     + \
                         " --angname=" + cf('LC2_ANGNAME')     + \
                  " --version_string=" + cf('LC2_VERSION')      + \
                        " --channels=" + cf('LC2_CHANNELS')     + \
                        " --icldtable="+ cf('LC2_ICLDTABLE')    + \
                        " --lcldtable="+ cf('LC2_LCLDTABLE')    + \
                         " --surface=" + cf('LC2_SURFACE')      + \
                     " --surfversion=" + cf('LC2_SURF_VERSION') + \
                          " --interp=" + cf('LC2_SURF_INTERP')  + \
                         " --nodemax=" + cf('LC2_NODEMAX')      + \
                          " --layout=" + cf('LC2_LAYOUT')       + \
                             " --dir=" + cf('LC2_DIR')          + \
                         " --runmode=" + cf('LC2_RUNMODE')      + \
                         " --runfile=" + cf('LC2_RUNFILE')      + \
                          " --execfile=" + cf('LC2_EXECFILE') 

    # Optional
    try:
        Options += " --c_band=" + cf('LC2_SURF_C_BAND')                      
    except:
        pass

    try:
        Options += " --cld_band=" + cf('LC2_CLD_C_BAND')                      
    except:
        pass        

    if cf('LC2_PROFILE').upper()   == 'YES': Options += " --profile"
    if cf('LC2_VERBOSE').upper()   == 'YES': Options += " --verbose"
    if cf('LC2_ADDOUTPUT').upper() == 'YES': Options += " --additional"
    if cf('LC2_KEEP').upper()      == 'YES': Options += " --keep"



    # Generate LC2 products
    # ---------------------
    cmd = "geo_vlidort_cloud_lc2.py %s %s %s "%(Options,start_isotime,end_isotime)
    print cmd
    if not options.dryrun:
        if system(cmd):
            raise ValueError, "geo_vlidort_cloud_lc2.py failed for %s to %s "%(start_isotime, end_isotime)




    




