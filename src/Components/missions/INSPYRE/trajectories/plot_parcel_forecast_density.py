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
    ax, ax2 = make_plot(label)
    nt = 73
    colors=["blue","red","green"]
    mcolor=["cornflowerblue","maroon","lime"]

    for j,h in enumerate(alts):
        if j == 0:
            cmap = "Blues"
            zorder = 97
        elif j == 1:
            cmap = "YlOrRd"
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
#       This converts np.datetime64 -> datetime.datetime; note [us] is microseconds
        time = ds["time"].values.astype('datetime64[us]').astype('O')
        now  = time[0]
        now0 = datetime(now.year,now.month,now.day)
        dp1  = now0+timedelta(days=0)
        dp2  = now0+timedelta(days=1)
        dp3  = now0+timedelta(days=2)
        ind  = []
        for dp in [dp1,dp2,dp3]:
            dstart = datetime(dp.year,dp.month,dp.day,19,59)
            dend   = datetime(dp.year,dp.month,dp.day,22,59)
            ind_   = [i for i,d in enumerate(time.tolist()) if (d>=dstart) & (d<=dend)]
            if ind_ != []:
                ind.append(ind_)
        ds.close()

        for i,it in enumerate(ind):
            kdeplot(ax,it,lon,lat,cmap,zorder=zorder)
        if j == 0:
            ax.plot(lon[0],lat[0],marker="p",c="k",transform=ccrs.PlateCarree(),zorder=100)
        for i in range(0,2500,25):
            cs = ax2.plot(time,alt[:,i],color=colors[j])

    plt.savefig(f"{timestr}:00_{firestr}.density.png")
