
#!/usr/bin/env python3
"""
Plot flight track on a background map
"""

import sys
import os
from numpy import array
from pyobs import ICARTT
from grads import GrADS
from matplotlib.pyplot import title, savefig

if __name__ == "__main__":

    nav = '/home/adasilva/iesa/aerosol/data/WE-CAN/Nav'
    if len(sys.argv) < 2:
        print("Usage: ")
        print("       plot_track flight_number")
        print("Example:") 
        print("       plot_track 01")
        raise SystemExit("Error: not enough arguments")
    else:
        rf = sys.argv[1]
        track = nav + '/WECANrf'+rf+'.ict.gz'
        figFile = './WECANrf'+rf+'.png'


    # Tracks
    # ------
    m = ICARTT(track)

    # Start pygrads
    # -------------
    ga = GrADS(Window=False,Echo=False)

    # Domain
    # -----
    ga("""
       @ open $GADSET/model
       set lon -126.5 -106.5'
       set lat 36 50'
       """)

    # Basemap
    # -------
    ga.basemap()
    ga.blue_marble(Show=True)
    ga.map.plot(m.Nav['Longitude'],m.Nav['Latitude'],'m',linewidth=2)
    ga.map.drawcountries(color='white')
    ga.map.drawstates(color='white')
    ga.map.drawcoastlines(color='white')


    # Hours
    # -----
    I = array([t.minute for t in m.tyme])==0
    ga.map.plot(m.Nav['Longitude'][I],m.Nav['Latitude'][I],'yo',ms=8)
    title('WE-CAN Research Flight %s on %s '%(rf,m.tyme[0].isoformat().split('T')[0]))

    savefig(figFile,bbox_inches='tight')
    
