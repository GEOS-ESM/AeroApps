#!/usr/bin/env python

from datetime import datetime, timedelta
from pyobs    import G5_ICARTT

if __name__ == "__main__":

#    template = 'discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict'
#    top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim'

    outTopDir = './Archive'

    template = 'discoveraq-geos5das-PRODUCT_SITE_DATE_RC.1001.ict'
#    top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/candidate_fp/assim'

    top_url = '/home/adasilva/opendap/daq/assim'

    coll_names = [ 'inst1_2d_hwl_Nx',
                   'inst3_2d_asm_Nx',
                   'tavg1_2d_flx_Nx',
                   'tavg1_2d_lnd_Nx',
                   'tavg1_2d_ocn_Nx',
                   'tavg1_2d_rad_Nx',
                   'tavg1_2d_slv_Nx',
                   'tavg3_2d_adg_Nx',
                   'tavg3_2d_aer_Nx',
                   'tavg3_2d_chm_Nx',
#                  'tavg3_2d_tag_Nx',
                 ]

    dt    = timedelta(seconds=86400)
    t_beg = datetime(2011,7,1)
    t_end = datetime(2011,7,31)
   
    g = G5_ICARTT('./d-aq_sites.rc', coll_names=coll_names,
                  template=template, top_url=top_url,
                  outTopDir = outTopDir,
                  Verbose=True)

    t = t_beg
    while t <= t_end:
        g.Sample_N_Write(t)
        t = t + dt
    
    
