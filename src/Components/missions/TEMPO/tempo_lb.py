#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to run different LevelB Sampling for TEMPO
  
  
  patricia.castellanos@nasa.gov
"""
from os              import system
from optparse        import OptionParser
from datetime        import timedelta
from dateutil.parser import parse         as isoparser
if __name__ == "__main__":
    # Defaults
    bin = 'geo_samplery.py'
    levelBdir = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/DATA/LevelB/'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog isoStartTime isoEndTime  rcFile",
                          version='1.0.0' )

    parser.add_option("-b", "--bin", dest="bin",default=bin, 
              help="executable for run (default=%s)" %bin ) 

    parser.add_option("--levelBdir", dest="levelBdir",default=levelBdir,
              help="LevelB Data Root Directory" )

    parser.add_option("-d", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="do a dry run")        

    # Parse arguments
    # ---------------
    (options, args) = parser.parse_args()
    if len(args) == 3:
        isoStartTime = args[0]
        isoEndTime   = args[1]
        rcFile       = args[2]
    else:
        parser.error("must have 1 arguments: pcfFile")    

    levelBdir    = options.levelBdir
    dryrun       = options.dryrun
    bin          = options.bin

    if dryrun:
        command = bin + ' --dryrun'
    else:
        command = bin



    command += '-v -C -r ' + rcFile 

    startdate = isoparser(isoStartTime)
    enddate   = isoparser(isoEndTime)

    outTemplate = '{}/Y{}/M{}/D{}/tempo-g5nr.lb2.{}.{}_{}z.nc4'
    while startdate <= enddate:
        Y = startdate.year
        M = str(startdate.month).zfill(2)
        D = str(startdate.day).zfill(2)
        H = str(startdate.hour).zfill(2)
        col = rcFile[:-3]
        outfile = outTemplate.format(levelBdir,Y,M,D,col,Y+M+D,H)


        run_command = command + ' -o' + outfile + ' ' + startdate.isoformat()

        print run_command
        if system(run_command):
            raise ValueError, "%s failed for %s to %s "%(bin,startdate.isoformat(),rcFile)

        startdate += timedelta(hours=1)
