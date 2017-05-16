#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""

  Simple wrapper script to run different episodes of TEMPO swath using
  run_geo_vlidort_lc2.py
  
  May 2017.
  patricia.castellanos@nasa.gov
"""
from os       import system

if __name__ == "__main__":

    pcfFile   = 'geo_vlidort_lc2.pcf'
    start_isotime = '2005-12-31T15:00:00Z'
    end_isotime   = '2005-12-31T15:00:00Z'
    episode   = None
    dryrun    = True

    if episode is not None:
        if episode == 1:
            start_isotime  = '2005-12-31T00:00:00Z'
            end_isotime    = '2006-01-01T23:00:00Z'
        elif episode == 2:
            start_isotime  = '2006-07-27T00:00:00Z'
            end_isotime    = '2006-08-09T23:00:00Z'
        elif episode == 3:
            start_isotime  = '2007-03-28T00:00:00Z'
            end_isotime    = '2007-03-29T23:00:00Z'
        elif episode == 4:
            start_isotime  = '2007-04-10T00:00:00Z'
            end_isotime    = '2007-04-11T23:00:00Z'


    bin = 'run_geo_vlidort_lc2.py'

    if dryrun:
        command = bin + ' --dryrun'
    else:
        command = bin

    command += ' ' + pcfFile + ' ' + start_isotime + ' ' + end_isotime

    print command
    if system(command):
        raise ValueError, "run_geo_vlidort_lc2.py failed for %s to %s "%(start_isotime,end_isotime)