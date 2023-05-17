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
             Log=False,vmin=0.1,vmax=0.8,figFile=None,reverse=False):

    print("[] Plotting <%s>"%tit)

    clf()
    ga('reinit')
    ga.open(fname)
    ga('set lon -130 -70')
    ga('set lat 5 55')
    # ga('set lon -120 -80')
    # ga('set lat 25 50')
    ga('set time %s'%gacore.dt2gat(tyme))

    ga.blue_marble('on')
    if Log:
	if reverse == True:
		ga.imshow(vname,norm=LogNorm(vmin=vmin,vmax=vmax),
			  cmap=gacm.jet_r)
	else:
		ga.imshow(vname,norm=LogNorm(vmin=vmin,vmax=vmax),
			  cmap=gacm.jet)
    else:
	if reverse == True:
		v = ga.exp(vname)
		v[v>vmax] = vmax
		ga.imshow(v,vmin=vmin,vmax=vmax,dlon=360,dlat=180,
			  cmap=gacm.jet_r)
	else:
		v = ga.exp(vname)
		v[v>vmax] = vmax
		ga.imshow(v,vmin=vmin,vmax=vmax,dlon=360,dlat=180,
			  cmap=gacm.jet_l)
        
    ga.map.drawcountries()
    ga.map.drawstates()
    ga.map.drawcoastlines()

    # ga.map.plot(lons,lats,'w',linewidth=1.5)
    
    title(tit)
    
    if figFile is not None:
        savefig(figFile,dpi=180)

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
    if len(sys.argv) < 3:
        print("Usage: ")
        print("       anim_aod ending_yyyymmdd_hh ndays_back")
        print("Example:") 
        print("       anim_aod 20130819_21 7")
        raise SystemExit("Error: not enough arguments")
    else:
        dtag = sys.argv[1]
        ndays_back = int(sys.argv[2])
        year, month, day, hour = dtag[0:4], dtag[4:6], dtag[6:8], dtag[9:11] 
        etime = datetime(year=int(year), month=int(month), day=int(day), hour=int(hour))

    # Verification times
    # ------------------
    day = timedelta(seconds=24*60*60)
    hour = timedelta(seconds=60*60)
    btime = etime - ndays_back * day

    # GrADS
    # -----
    ga = GrADS(Window=False,Echo=False)
    hwl_Nx = 'http://opendap.nccs.nasa.gov:80/dods/GEOS-5/fp/0.25_deg/assim/inst1_2d_hwl_Nx'

    # Plot AOD
    # --------
    dt = etime - btime
    nhours = 1+int(dt.total_seconds()/3600.)

    # loop over time
    # --------------
    vnames = ['bcexttau', 'ocexttau', 'duexttau', 'suexttau', 'ssexttau']
    VN = ['BC', 'OC', 'DU', 'SU', 'SS']
    imin = [0.05, 0.3, 0.1, 0.3, 0.3]
    imax = [0.15, 0.8, 0.3, 0.8, 0.8]
    for iv, v in enumerate(VN):
	for h in range(hours):

		t = btime + h * hour
        
		clf()
		tv = get_tstr(t)
		td = get_dtag(t)
		name = 'aod.'+v+'.'
		ender = '%03d.png' % (h+1)
		figFile = name+ender
		print(figFile)
        
		plot_aod(ga,'GEOS-5 '+v+' fraction of AOD at %s'%tv,'('+vnames[iv]+'/totexttau)',hwl_Nx,t,figFile=figFile,vmin=imin[iv],vmax=imax[iv])

		if iv == 0:
			figFile = 'ssa.'+ender
			print(figFile)
			plot_aod(ga,'GEOS-5 SSA at %s'%tv,'(totscatau/totexttau)',hwl_Nx,t,figFile=figFile,vmin=0.84,vmax=0.93,reverse=True)

			figFile = 'ang.'+ender
			print(figFile)
			plot_aod(ga,'GEOS-5 Angstrom parameter [470-870 nm] at %s'%tv,'totangstr',hwl_Nx,t,figFile=figFile,vmin=0.2,vmax=1.8,reverse=True)

        
