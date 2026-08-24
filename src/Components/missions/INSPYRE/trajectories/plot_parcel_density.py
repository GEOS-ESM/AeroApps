#!/usr/bin/env python3

from optparse import OptionParser   # Command-line args
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import os
from traj_make_plot import make_plot, kdeplot, add_satellite_track
from datetime import datetime, timedelta

satellite_colors = {"EarthCARE": "red", "NOAA-20": "magenta", "PACE": 'green'}

if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] filename",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        parcel = args[0]
    else:
        parser.error("must have 1 argument: filename")
        
    # Assuming parcel looks like: "parcel_traj.20260804_00.Bobcat_Lakes_2026-08-07T00:00:00"
    x = parcel.split(".")
    fcst_cycle = x[1]
    label = x[2]
    print(f"Forecast Cycle: {fcst_cycle}")
    print(f"Label: {label}")

    timestr = label[-19:-3]
    firestr = label[0:-20]

    print(parcel)
    ds = xr.open_mfdataset(parcel)
    lat = ds["lat"].values
    lon = ds["lon"].values
    alt = ds["PAlt"].values
    
    # Converts np.datetime64 -> datetime.datetime; note [us] is microseconds
    time = ds["time"].values.astype('datetime64[us]').astype('O')
    
    # Convert the extracted time string into a proper datetime object
    init_time = datetime.strptime(timestr, "%Y-%m-%dT%H:%M")
    base_date = datetime(init_time.year, init_time.month, init_time.day)
    print(f"Base Date: {base_date}")
    
    # Check the hour to determine when the first 21Z plot should occur
    if init_time.hour > 20:
        # If initialized late (21Z-23Z), the next 21Z is tomorrow
        offset = 1
    else:
        # If initialized early (0Z-20Z), the next 21Z is today
        offset = 0

    # Calculate forecast days dynamically based on the offset
    dp1  = base_date + timedelta(days=offset)
    dp2  = base_date + timedelta(days=offset + 1)
    dp3  = base_date + timedelta(days=offset + 2)
    
    ind = []
    valid_times = []  # Added to track the vtime string for the titles/filenames
    valid_dps = []    # Added to track the datetime for the satellite calculations
    
    for dp in [dp1,dp2,dp3]:
        vtime = dp.strftime("%Y-%m-%dT18:00:00")
        dstart = datetime(dp.year, dp.month, dp.day, 16, 59)
        dend   = datetime(dp.year, dp.month, dp.day, 18, 59)
        ind_   = [i for i,d in enumerate(time.tolist()) if (d>=dstart) & (d<=dend)]
        
        if ind_ != []:
            ind.append(ind_)
            valid_times.append(vtime)
            valid_dps.append(dp)
            
    print(ind)
    ds.close()

    lon[np.where(lon<0.)] = lon[np.where(lon<0.)]+360.

    # ---------------------------------------------------------
    # MAIN PLOTTING LOOP
    # ---------------------------------------------------------
    for i, it in enumerate(ind):
        
        current_vtime = valid_times[i]
        current_dp = valid_dps[i]
        
        ax, ax2 = make_plot(label, current_vtime)
        
        kdeplot(ax, it, lon, lat, cmap="Blues")
        
        for k in range(0, 2500, 25):
            ax.plot(lon[:,k], lat[:,k], c="blue", transform=ccrs.PlateCarree(), alpha=0.5, zorder=50)
            ax2.plot(time, alt[:,k], color="blue")
            
        # Draw the starting point
        ax.plot(lon[0], lat[0], marker="p", c="k", transform=ccrs.PlateCarree(), zorder=100)

        # 4. Add the satellites (using the current_dp for accurate tracks!)
        legend_handles = []
        legend_labels = []
        for satellite_name in ["EarthCARE", "PACE"]:
            track, line_handle = add_satellite_track(
                ax=ax,
                satellite=satellite_name,
                start=current_dp,
                stop=current_dp + timedelta(days=1),
                color=satellite_colors[satellite_name],
                interval_minutes=1.0,
                marker_interval_minutes=5,
                min_solar_elevation=0.0,
                linewidth=3,
                marker_size=100,
                label_fontsize=25,
                label_times=True
            )
            legend_handles.append(line_handle)
            legend_labels.append(satellite_name)
        ax.legend(handles=legend_handles, labels=legend_labels, loc='upper right', fontsize=20, framealpha=1.0, facecolor='white').set_zorder(999)
        filename = f"fcst{fcst_cycle}_{timestr}:00+{current_vtime}_{firestr}.density.png"
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        
