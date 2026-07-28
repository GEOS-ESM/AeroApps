#!/usr/bin/env python3

from optparse   import OptionParser   # Command-line args
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import os
from traj_make_plot import make_plot, kdeplot
from datetime import datetime, timedelta

if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] filename",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        parcel_ = args[0]
    else:
        parser.error("must have 1 argument: filename")

    x = parcel_.split(".")
    label = x[1]
    print(label)
    timestr = label[-19:-3]
    firestr = label[0:-20]

    alts = ["10km","12km","14km"]
    nt = 73
    colors=["red","blue","green"]
    mcolor=["maroon","cornflowerblue","lime"]

#   Loop over valid times for plots
#   Haven't quite worked out how I want to do "now" -- probably should be
#   relative to forecast start time?
    now  = datetime.strptime(timestr,"%Y-%m-%dT%H:%M")
    now0 = datetime(now.year,now.month,now.day)
    print(now0)
    dp1  = now0+timedelta(days=0)
    dp2  = now0+timedelta(days=1)
    dp3  = now0+timedelta(days=2)

    for dp in [dp1,dp2,dp3]:
        vtime = dp.strftime("%Y-%m-%dT21:00:00")
        dstart = datetime(dp.year,dp.month,dp.day,19,59)
        dend   = datetime(dp.year,dp.month,dp.day,22,59)
        ax, ax2 = make_plot(label,vtime)

        for j,h in enumerate(alts):
            if j == 0:
                cmap = "YlOrRd"
                zorder = 97
            elif j == 1:
                cmap = "Blues"
                zorder = 98
            else:
                cmap = "YlGn"
                zorder = 99
            parcel = f"{parcel_}_{h}.nc"
            print(parcel)
            ds = xr.open_mfdataset(parcel)
            lat = ds["lat"].values
            lon = ds["lon"].values
            alt = ds["PAlt"].values
#           This converts np.datetime64 -> datetime.datetime; note [us] is microseconds
            time = ds["time"].values.astype('datetime64[us]').astype('O')
            ds.close()

            ind  = []
            ind_   = [i for i,d in enumerate(time.tolist()) if (d>=dstart) & (d<=dend)]
            if ind_ != []:
                ind.append(ind_)

            for i,it in enumerate(ind):
                kdeplot(ax,it,lon,lat,cmap,zorder=zorder)
            if j == 0:
                ax.plot(lon[0],lat[0],marker="p",c="k",transform=ccrs.PlateCarree(),zorder=100)
            for i in range(0,2500,25):
                cs = ax.plot(lon[:,i],lat[:,i],c=mcolor[j],transform=ccrs.PlateCarree(),zorder=50,alpha=.5)
            for i in range(0,2500,25):
                cs = ax2.plot(time,alt[:,i],color=mcolor[j])

            plt.savefig(f"{timestr}:00+{vtime}_{firestr}.density.png")
#            plt.close()
