#!/usr/bin/env python

from datetime import datetime, timedelta
from pyobs    import G5_AIRCRAFT

if __name__ == "__main__":

#    template = 'discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict'
#    top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim'

    template = 'discoveraq-geos5das-PRODUCT_SITE_DATE_RB.1001.ict'
    top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/candidate_fp/assim'

    coll_names = [ 'inst3_3d_asm_Nv',    # collection with altitude H
                   'inst3_3d_aer_Nv',
                   'inst3_3d_chm_Nv',
                   'inst3_3d_tag_Nv',
                 ]
    AltColl_name = coll_names[0]

    daq_dirn  = '/nobackup/DISCOVER-AQ/Archive/UMD-AIRCRAFT/HAO.HE'
    fname = daq_dirn + '/discoveraq-UMDAircraft_UMD-AIRCRAFT_20110710_RA_L1.ict' 

    a = G5_AIRCRAFT(fname,
                    coll_names = coll_names,
                    coll_dir = 'Collections/Aircraft',
                    top_url=top_url,
                    template=template,
                    outTopDir='./Archive',
                    Verbose=True,
                    )

    for c in coll_names:
        a.sample(c,AltVar='h',AltColl_name=AltColl_name)
        a.write()
    
    
