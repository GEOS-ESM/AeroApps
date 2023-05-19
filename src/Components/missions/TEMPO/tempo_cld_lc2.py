#!/usr/bin/env python3
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to run different episodes of TEMPO swath using
  run_geo_vlidort_lc2.py
  
  May 2017.
  patricia.castellanos@nasa.gov
"""
from os       import system
from optparse        import OptionParser

if __name__ == "__main__":

    # pcfFile   = 'geo_vlidort_cloud_lc2.pcf'
    isoStartTime = '2005-12-31T18:00:00Z'
    isoEndTime   = '2005-12-31T18:00:00Z'
    # episode   = None
    # dryrun    = False
    bin = 'run_geo_vlidort_cloud_lc2.py'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] pcfFile",
                          version='1.0.0' )

    parser.add_option("-b", "--bin", dest="bin",default=bin, 
              help="executable for run (default=%s)" %bin ) 

    parser.add_option("-s", "--isoStartTime", dest="isoStartTime", default=isoStartTime,
              help="start isotime" ) 

    parser.add_option("-e", "--isoEndTime", dest="isoEndTime", default=isoEndTime,
              help="end isotime" )                  

    parser.add_option("-E", "--episode", dest="episode",
              help="preset episode number" )

    parser.add_option("-d", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="do a dry run")        

    # Parse arguments
    # ---------------
    (options, args) = parser.parse_args()
    if len(args) == 1:
        pcfFile = args[0]
    else:
        parser.error("must have 1 arguments: pcfFile")  

    isoStartTime = options.isoStartTime
    isoEndTime   = options.isoEndTime
    episode      = options.episode
    dryrun       = options.dryrun
    bin          = options.bin        

    if episode is not None:
        if int(episode) == 1:
            isoStartTime  = '2005-12-31T00:00:00Z'
            isoEndTime    = '2006-01-01T23:00:00Z'
        elif int(episode) == 2:
            isoStartTime  = '2006-07-27T00:00:00Z'
            isoEndTime    = '2006-08-09T23:00:00Z'
        elif int(episode) == 3:
            isoStartTime  = '2007-03-28T00:00:00Z'
            isoEndTime    = '2007-03-29T23:00:00Z'
        elif int(episode) == 4:
            isoStartTime  = '2007-04-10T00:00:00Z'
            isoEndTime    = '2007-04-11T23:00:00Z'

    if dryrun:
        command = bin + ' --dryrun'
    else:
        command = bin

    command += ' ' + pcfFile + ' ' + isoStartTime + ' ' + isoEndTime

    print(command)
    if system(command):
        raise ValueError("%s failed for %s to %s "%(bin,isoStartTime,isoEndTime))
