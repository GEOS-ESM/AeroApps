#!/usr/bin/env python
"""
Plot flight track on a background map
"""

import sys
import os
from numpy import array
from pyobs import ICARTT
from grads import GrADS

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Usage: "
        print "       plot_track track.ict"
        print "Example:" 
        print "       plot_track ../Tracks/routine/ORACLES-Flt-plan_p3_20160901_RA.ict "
        raise SystemExit, "Error: not enough arguments"
    else:
        track = sys.argv[1]

    # Tracks
    # ------
    m = ICARTT(track)

    # Start pygrads
    # -------------
    ga = GrADS(Window=False,Echo=False)

    # Domain
    # ------
    ga("""
       @ open $GADSET/model
       set lon -15 20
       set lat -30 0
       """)

    # Basemap
    # -------
    ga.basemap()
    ga.blue_marble(Show=True)
    ga.map.plot(m.Longitude,m.Latitude,'m',linewidth=2)
    
    # Hours
    # -----
    I = array([t.minute for t in m.tyme])==0
    ga.map.plot(m.Longitude[I],m.Latitude[I],'yo',ms=15)

    
    
    
    
