
#!/usr/bin/env python
"""
  Implements direct emission estimates based on MODIS Deep Blue C6
  retrievals.

"""

import os
import sys

from numpy    import linspace, array, savez
from glob     import glob
from datetime import datetime, timedelta
from string   import Template

Force = False


def monthly_means(path,coll,outFile,year,month):
    """
    Given a time range in gatime, use lats4d to compute sampled mothly mean.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(outFile):
            print '<> File exists, skipping <%s>'%filename
            return

    if year<2010:
        stream = 'MERRA2_300'
    else:
        stream = 'MERRA2_400'
        
    here = os.getcwd()
    ldir = path+'/'+stream+'/Y%d/M%02d'%(year,month)
    os.chdir(ldir)
    patt = '*.%s.*.nc4'%coll
    inFiles = sorted(glob(patt))
    if len(inFiles)<1:
        raise ValueError, 'no files for %d-%02d'%(year,month)

    cmd = 'GFIO_mean_r8.x -o '+outFile+' '+' '.join(inFiles)

    print 30*'-'
    print cmd
    print 30*'-'
    if os.system(cmd):
        raise RuntimeError, 'error on return from %s'%cmd

    os.chdir(here)

#....................................................................

if __name__ == "__main__":

    inPath = '/home/adasilva/opendap/m2'

    # Parse command line
    # ------------------
    if len(sys.argv) == 3:
        year1 = sys.argv[1]
        year2 = sys.argv[2]
    elif len(sys.argv) == 2 :
        year1 = sys.argv[1]
        year2 = year1
    else:
        print "   Usage:   %s  year1 [year2]"%sys.argv[0]
        print "Examples:   %s   2003  2005"%sys.argv[0]
        print "            %s   2003"%sys.argv[0]
        raise RuntimeError, 'not enough parameters'

    year1, year2 = int(year1), int(year2)    

    # Loop over time and products
    # ---------------------------
    for year in range(year1,year2+1):
        for month in range(1,13):
            for coll in ( 'inst3_3d_gas_Nv',):
                
                mcoll = coll.replace('inst1','instM')\
                            .replace('inst3','instM')\
                            .replace('inst6','instM')\
                            .replace('tavg1','tavgM')\
                            .replace('tavg3','tavgM')\
                            .replace('tavg6','tavgM')

                dirn = os.getcwd()+'/%s/Y%d/M%02d'%(mcoll,year,month)
                os.system('/bin/mkdir -p '+dirn)

                outFile = dirn+'/MERRA2.%s.%d%02d.nc4'%(mcoll,year,month)

                monthly_means(inPath,coll,outFile,year,month)


#---------------





