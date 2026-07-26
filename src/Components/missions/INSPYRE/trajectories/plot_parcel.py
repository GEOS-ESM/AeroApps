#!/usr/bin/env python3

from optparse   import OptionParser   # Command-line args
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import os
from traj_make_plot import make_plot



if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] filename",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        parcel = args[0]
    else:
        parser.error("must have 1 argument: filename")

    x = parcel.split(".")
    label = x[1]
    print(label)
    timestr = label[-19:-3]
    firestr = label[0:-20]

    ax, ax2 = make_plot(label)
    nt = 73
    print(parcel)
    ds = xr.open_mfdataset(parcel)
    lat = ds["lat"].values
    lon = ds["lon"].values
    alt = ds["PAlt"].values
    time = ds["time"].values
    ds.close()

    for i in range(0,2500,25):
        cs = ax.plot(lon[:,i],lat[:,i],c="blue",transform=ccrs.PlateCarree())

    for i in range(0,2500,25):
        cs2 = ax.scatter(x=lon[24,i],y=lat[24,i],c="k",marker="o",transform=ccrs.PlateCarree(),zorder=100)
        cs2 = ax.scatter(x=lon[48,i],y=lat[48,i],c="k",marker="^",transform=ccrs.PlateCarree(),zorder=100)
        cs2 = ax.scatter(x=lon[72,i],y=lat[72,i],c="k",marker="P",transform=ccrs.PlateCarree(),zorder=100)
    ax.plot(lon[0],lat[0],marker="p",c="k",transform=ccrs.PlateCarree(),zorder=100)

    for i in range(0,2500,25):
        cs = ax2.plot(time,alt[:,i],color="blue")

    plt.savefig(f"{timestr}:00_{firestr}.png")
