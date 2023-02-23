#!/usr/bin/env python3
"""
Plots AOD analysis from opendap.
"""

import sys, os
from datetime import datetime, timedelta

from grads import GrADS, gacore, gacm

from matplotlib.pyplot import contourf, xlabel, ylabel, title, grid, plot, \
                              figure, gca, clf, cm, savefig, axes, \
                              colorbar, legend, show, subplot
from matplotlib.colors import LogNorm, Normalize

#-------------------------
def plot_aod(ga,tit,vname,fname,tyme,sub=None,
             Log=False,vmin=0.1,vmax=0.8,figFile=None):

    print("[] Plotting <%s>"%tit)

    # clf()
    ga('reinit')
    ga.open(fname)
    ga('set lon -130 -70')
    ga('set lat 5 55')
    # ga('set lon -120 -80')
    # ga('set lat 25 50')
    ga('set time %s'%gacore.dt2gat(tyme))

    ga.blue_marble('on')
    if Log:
        ga.imshow(vname,norm=LogNorm(vmin=vmin,vmax=vmax),
                  cmap=gacm.jet)
    else:
        v = ga.exp(vname)
        v[v>vmax] = vmax
        ga.imshow(v,vmin=vmin,vmax=vmax,dlon=360,dlat=180,
                  cmap=gacm.jet_l,sub=sub)
        
    ga.map.drawcountries()
    ga.map.drawstates()
    ga.map.drawcoastlines()

    # ga.map.plot(lons,lats,'w',linewidth=1.5)
    
    title(tit)
    
    #if figFile is not None:
        # savefig(figFile,dpi=180)

def get_dtag(t):
    dtag = '%4d%02d%02d_%02d'%(t.year,t.month,t.day,t.hour)
    return dtag
        
def get_tstr(t):
    dtag = '%02dZ %4d-%02d-%02d'%(t.hour,t.year,t.month,t.day)
    return dtag
        
#-------------------------
if __name__ == "__main__":

    # Parse command line
    # ------------------
    if len(sys.argv) < 2:
        print("Usage: ")
        print("       plot_verif yyyymmdd_hh")
        print("Example:") 
        print("       plot_verif 20130814_12")
        raise SystemExit("Error: not enough arguments")
    else:
        dtag = sys.argv[1]
        year, month, day, hour = dtag[0:4], dtag[4:6], dtag[6:8], dtag[9:11] 
        vtime = datetime(year=int(year), month=int(month), day=int(day), hour=int(hour))

    # Verification times
    # ------------------
    day = timedelta(seconds=24*60*60)
    fstart = vtime-day, vtime-2*day, vtime-3*day

    # GrADS
    # -----
    ga = GrADS(Window=False,Echo=False)
    gas_Nx = 'http://opendap.nccs.nasa.gov:80/dods/GEOS-5/fp/0.25_deg/assim/inst3_2d_gas_Nx'

    # Plot analysis
    # -------------
    clf()
    tv = get_tstr(vtime)
    td = get_dtag(vtime)
    plot_aod(ga,'AOD at %s'%tv,'aodana',gas_Nx,vtime,
             figFile='aod.ana.%s.png'%td,sub=221)

    # Forecasts
    # ---------
    root = 'http://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast'
    vdate = datetime(vtime.year,vtime.month,vtime.day)
    i = 0
    for dt in (0,1,2):
        i += 1
        t0 = vdate - dt * day
        hwl_Nx = root+'/inst1_2d_hwl_Nx/inst1_2d_hwl_Nx.'+get_dtag(t0)
        #tit = '%d-Day AOD Forecast Valid at %s'%(dt,tv)
        dh = vtime - t0
        nh = int(dh.total_seconds()/3600.)
        tit = '%d-Hour Forecast'%nh
        plot_aod(ga,tit,'totexttau',hwl_Nx,vtime,
                 sub=221+i,figFile='aod.f%02d.%s.png'%(dt,td))
   
    savefig('aod.verif.%s.png'%td,dpi=180)
