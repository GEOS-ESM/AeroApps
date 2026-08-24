#!/usr/bin/env python3

from optparse import OptionParser   # Command-line args
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import os
import glob
import sys
from traj_make_plot import make_plot, kdeplot, add_satellite_track
from datetime import datetime, timedelta

# Added PACE to the satellite colors dictionary
satellite_colors = {"EarthCARE": "crimson", "NOAA-20": "magenta", "PACE": 'green'}

if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] filename",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()
    if len(args) == 1:
        parcel_ = args[0]
    else:
        parser.error("must have 1 argument: filename")

    # Assuming parcel_ looks like: "parcel_traj.20260804_00.Bobcat_Lakes_2026-08-07T00:00:00"
    x = os.path.basename(parcel_).split(".")
    
    fcst_cycle = x[1] 
    label = x[2] 
    
    print(f"Forecast Cycle: {fcst_cycle}")
    print(f"Label: {label}")
    
    timestr = label[-19:-3]
    firestr = label[0:-20]
    
    # DYNAMICALLY FIND ALTITUDES 
    file_pattern = f"{parcel_}_*.nc"
    matched_files = glob.glob(file_pattern)
    
    if not matched_files:
        print(f" [!] CRITICAL ERROR: No altitude files found matching '{file_pattern}'")
        sys.exit(1)
        
    alts = []
    for f in matched_files:
        base_f = os.path.basename(f)
        base_p = os.path.basename(parcel_)
        alt_str = base_f.replace(base_p + "_", "").replace(".nc", "")
        alts.append(alt_str)
        
    # Sort them numerically so smaller altitudes plot first
    try:
        alts.sort(key=lambda val: float(val.replace("km", "")))
    except ValueError:
        alts.sort() # Fallback to alphabetical if formatting is weird
        
    print(f"Plotting Altitudes: {alts}")
    
    cmap_list = ["YlOrRd", "Blues", "YlGn"]
    mcolor_list = ["maroon", "cornflowerblue", "lime"]
  
    
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

    for dp in [dp1, dp2,]:#dp3]:
        vtime = dp.strftime("%Y-%m-%dT21:00:00")
        dstart = datetime(dp.year, dp.month, dp.day, 19, 59)
        dend   = datetime(dp.year, dp.month, dp.day, 21, 59)
        
        ax, ax2 = make_plot(label, vtime)
        for j, h in enumerate(alts):
            cmap = cmap_list[j]
            current_mcolor = mcolor_list[j]
            zorder = 97 + j
                
            parcel = f"{parcel_}_{h}.nc"
            print(f"Processing: {parcel}")
            
            try:
                ds = xr.open_mfdataset(parcel)
            except OSError:
                print(f" -> File not found: {parcel}. Skipping this altitude.")
                continue
                
            lat = ds["lat"].values
            lon = ds["lon"].values
            alt = ds["PAlt"].values
#           This converts np.datetime64 -> datetime.datetime; note [us] is microseconds
            time = ds["time"].values.astype('datetime64[us]').astype('O')
            ds.close()
            
            lon[np.where(lon < 0.)] = lon[np.where(lon < 0.)] + 360.

            ind_ = [i for i, d in enumerate(time.tolist()) if (d >= dstart) & (d <= dend)]

            if ind_:
                kdeplot(ax, ind_, lon, lat, cmap, zorder=zorder)
                
            # Plot the starting point ONLY ONCE (during the first altitude loop)
            if j == 0:
                ax.plot(lon[0], lat[0], marker="p", c="k", transform=ccrs.PlateCarree(), zorder=100)
                
            for i in range(0, 2500, 25):
                ax.plot(lon[:, i], lat[:, i], c=current_mcolor, transform=ccrs.PlateCarree(), zorder=50, alpha=.5)
                ax2.plot(time, alt[:, i], color=current_mcolor)

        # SATELLITE TRACKS & LEGEND
        legend_handles = []
        legend_labels = []
        for satellite_name in ["EarthCARE", "PACE"]:
            track, line_handle = add_satellite_track(
                ax=ax,
                satellite=satellite_name,
                start=dp,
                stop=dp + timedelta(days=1),
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
            
        ax.legend(handles=legend_handles, labels=legend_labels, loc='upper right', 
                  fontsize=20, framealpha=1.0, facecolor='white').set_zorder(999)

        filename = f"fcst{fcst_cycle}_{timestr}:00+{vtime}_{firestr}.density.png"
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Successfully saved {filename}\n")
