#!/usr/bin/env python3

"""
    Plot the output of imager_swath
    Make sure it looks sane

    P. Castellanos 2024

"""

import os
import sys

import argparse
from   dateutil.parser import parse         as isoparser
from   datetime        import datetime, timedelta
import numpy           as np
from   MAPL.config     import Config
import xarray          as xr
from   glob            import glob
import cartopy
import matplotlib.pyplot as plt
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib as mpl

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    DT_mins = 60
    dt_mins = 5

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("start_isotime",
                        help="starting iso time")
    parser.add_argument("end_isotime",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with track input file names")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("var",
                        help="variable name: time, vza, vaa, sza, saa, or scatAngle")

    parser.add_argument("-D","--DT_mins", default=DT_mins, type=int,
                        help="Timestep in minutes for each plot (default=%i)"%DT_mins)

    parser.add_argument("-d","--dt_mins", default=dt_mins, type=int,
                        help="Timestep in minutes for each input granule (default=%i)"%dt_mins)

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')

    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()

    cf             = Config(args.inst_pcf,delim=' = ')
    instname       = cf('instname')

    args = parser.parse_args()


    # loop through start to end time, plotting
    # ------------------------------------------
    sdate = isoparser(args.start_isotime)
    edate = isoparser(args.end_isotime)
    DT    = timedelta(minutes=args.DT_mins)
    dt    = timedelta(minutes=args.dt_mins)

    transform=cartopy.crs.PlateCarree()
    cmap = mpl.cm.jet

    while sdate < edate:
        ed    = sdate + DT 
        sd    = sdate
        nymd  = str(sdate.date()).replace('-','')      
        year = sdate.strftime('%Y')
        month = sdate.strftime('%m')
        day   = sdate.strftime('%d')

        inFile = []
        while sd < ed:
            hour  = str(sd.hour).zfill(2)
            minute = str(sd.minute).zfill(2)

            # read 5-minute granules input files
            # ------------------------------------
            inFile.append(inTemplate.replace('%col',instname).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME))

            sd += dt

        ds = xr.open_mfdataset(inFile,combine='by_coords',coords=['time'])
        ds = ds.squeeze()

        lon0 = ds.trjLon.where((ds.trjLat <5) & (ds.trjLat >-5)).mean().values
        if np.isnan(lon0):
            lon0 = ds.trjLon.mean().values
        
#        projection=cartopy.crs.Orthographic(central_longitude=lon0,central_latitude=0)
        projection = cartopy.crs.PlateCarree(central_longitude=lon0)
        axes_class = (GeoAxes,
                  dict(map_projection=projection))
        # plot swath on sphere
        # ------------------
        fig = plt.figure(figsize=[8,10])
        axgr = AxesGrid(fig,  (0.1, 0.5, 0.8, 0.4) , axes_class=axes_class,
                nrows_ncols=(1, 1),
                axes_pad=0.6,
                cbar_location='right',
                cbar_mode='single',
                cbar_pad=0.2,
                cbar_size='3%',
                label_mode='')  # note the empty label_mode

        ax = axgr[0]
        sdata = ds.variables[args.var]
        bounds = np.linspace(sdata.min(),sdata.max(),100)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
        label = args.var
        ax.pcolormesh(ds.longitude,ds.latitude,sdata,transform=transform,norm=norm,cmap=cmap)
        ax.set_title('swath path')
        ax.coastlines()
        ax.set_global()

        # colorbar
        cb = axgr.cbar_axes[0].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label=label,ticks=bounds[::10])
#                            cb.ax.set_yticklabels(ticks.astype(str))


        plt.show()
        sys.exit()

        plt.suptitle("{} {} ialong={}".format(pdate.isoformat(),args.var,ialong))
        plt.savefig("{}_{}_{}.png".format(pdate.strftime("%Y%m%d_%H"),args.var,ialong))
                    
        sdate += DT
