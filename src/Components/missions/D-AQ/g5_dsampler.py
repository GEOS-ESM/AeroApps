#!/usr/bin/env python

from datetime import datetime, timedelta
from pyobs.dragon import DRAGON
from glob import glob
from numpy import sort

if __name__ == "__main__":

#    anet_dir = '/nobackup/AERONET/Level1.5/Dragon'
#    outTopDir = '/nobackup/DISCOVER-AQ/Archive'
    anet_dir = '/home/adasilva/iesa/aerosol/data/AERONET/Dragon/Level2'
    outTopDir = './Archive'

#    template = 'discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict'
#    top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim'

    template = 'discoveraq-geos5das-PRODUCT_SITE_DATE_RC.1001.ict'
#    top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/candidate_fp/assim'
    top_url = '/home/adasilva/opendap/daq/assim'

    coll_names = [ 'inst1_2d_hwl_Nx',
                   'tavg1_2d_slv_Nx',
                     ]

#    for fname in sort(glob("%s/dragon_aod_*-JUL*.txt"%anet_dir)):
    for fname in sort(glob("%s/dragon_aod_*[123]?-JUL*.txt"%anet_dir)):

        
        print "[] Processing ", fname
        d = DRAGON(fname)
        d.sample_N_writeICARTT(coll_names = coll_names,
                               template=template, top_url=top_url,
                               outTopDir=outTopDir, Verbose=True)
  
    
