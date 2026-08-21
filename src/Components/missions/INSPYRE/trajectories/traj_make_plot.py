#!/usr/bin/env python3.11

import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon
from cartopy.geodesic import Geodesic  # from geographiclib
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.stats import gaussian_kde
import os
from satellite_groundtrack import get_daytime_ground_track


os.environ['CARTOPY_USER_BACKGROUNDS'] = "/home/pcolarco/silo/python/"

satellite_colors = {"EarthCARE": "deepskyblue","NOAA-20": "magenta", "PACE": 'green'}

def add_satellite_track(ax, satellite, start, stop, color, 
                        interval_minutes=1.0, marker_interval_minutes=5,
                        min_solar_elevation=0.0, linewidth=2.5,
                        marker_size=100, label_fontsize=20,label_times=True,):
    """
    Add a satellite daytime ground track to a Cartopy axis.

    Parameters
    ----------
    ax : cartopy.mpl.geoaxes.GeoAxes
        Cartopy map axis.

    satellite : str
        Satellite name accepted by get_daytime_ground_track().

    start, stop : datetime-like
        Beginning and ending UTC times.

    color : str
        Line and marker color.

    interval_minutes : float
        Orbit propagation interval.

    marker_interval_minutes : int
        Time interval between displayed markers.

    min_solar_elevation : float
        Minimum solar elevation defining daytime.

    linewidth : float
        Ground-track line width.

    marker_size : float
        Marker size.

    label_fontsize : float
        Marker-label font size.

    label_times : bool
        Whether to annotate each marker with UTC time.

    Returns
    -------
    satellite_track : dict
        Ground-track result returned by get_daytime_ground_track().
    """

    # ----------------------------------------------------------
    # Generate the satellite ground track
    # ----------------------------------------------------------
    satellite_track = get_daytime_ground_track(
        satellite=satellite,
        start=start,
        stop=stop,
        interval_minutes=interval_minutes,
        min_solar_elevation=min_solar_elevation,
    )

    times = np.asarray(satellite_track["times"],dtype=object,)

    latitude = np.asarray(satellite_track["lat"],dtype=float)

    longitude = np.asarray(satellite_track["lon"],dtype=float,)

    daytime = np.asarray(satellite_track["daytime"],dtype=bool)

    # ----------------------------------------------------------
    # Plot daytime ground track
    # ----------------------------------------------------------
    daytime_lat = np.where(daytime,latitude,np.nan)

    daytime_lon = np.where(daytime,longitude,np.nan)

    # Break the line when crossing the international date line.
    longitude_jump = (np.abs(np.diff(daytime_lon)) > 180.0)

    daytime_lon[1:][longitude_jump] = np.nan
    daytime_lat[1:][longitude_jump] = np.nan

    satellite_line, = ax.plot(daytime_lon,daytime_lat,linestyle="-",linewidth=linewidth,color=color,transform=ccrs.PlateCarree(),
                              zorder=80,clip_on=True,label="{} daytime ground track".format(satellite))

    # ----------------------------------------------------------
    # Select daytime markers
    # ----------------------------------------------------------
    marker_time_mask = np.array([(current_time.minute% marker_interval_minutes== 0) and current_time.second == 0 for current_time in times],dtype=bool)

    marker_mask  = daytime & marker_time_mask
    marker_lon   = longitude[marker_mask]
    marker_lat   = latitude[marker_mask]
    marker_times = times[marker_mask]

    # ----------------------------------------------------------
    # Retain only markers inside the map extent
    # ----------------------------------------------------------
    lon_min, lon_max, lat_min, lat_max = ax.get_extent(crs=ccrs.PlateCarree())

    inside_map = (
        np.isfinite(marker_lon)
        & np.isfinite(marker_lat)
        & (marker_lon >= lon_min)
        & (marker_lon <= lon_max)
        & (marker_lat >= lat_min)
        & (marker_lat <= lat_max)
    )

    marker_lon = marker_lon[inside_map]
    marker_lat = marker_lat[inside_map]
    marker_times = marker_times[inside_map]

    # ----------------------------------------------------------
    # Plot markers
    # ----------------------------------------------------------
    ax.scatter(
        marker_lon,
        marker_lat,
        s=marker_size,
        marker="o",
        facecolor=color,
        edgecolor="black",
        linewidth=1.0,
        transform=ccrs.PlateCarree(),
        zorder=300,
        clip_on=True,
    )

    # ----------------------------------------------------------
    # Add UTC labels
    # ----------------------------------------------------------
    if label_times:

        for i, (lon, lat, current_time) in enumerate(
            zip(
                marker_lon,
                marker_lat,
                marker_times,
            )
        ):
            time_label = current_time.strftime("%H:%M")

            # Alternate label position to reduce overlap.
            if i % 2 == 0:
                offset = (7, 7)
                vertical_alignment = "bottom"
            else:
                offset = (7, -7)
                vertical_alignment = "top"

            annotation = ax.annotate(
                time_label,
                xy=(lon, lat),
                xycoords=(
                    ccrs.PlateCarree()
                    ._as_mpl_transform(ax)
                ),
                xytext=offset,
                textcoords="offset points",
                fontsize=label_fontsize,
                color=color,
                ha="left",
                va=vertical_alignment,
                zorder=310,
                annotation_clip=True,
                clip_on=True,
            )

            # Prevent text from extending outside the map axis.
            annotation.set_clip_path(ax.patch)

    return satellite_track, satellite_line
    

def kdeplot(ax,it,lon,lat,cmap,zorder=99):
#   KDE
#   This is inspired by https://www.iditect.com/faq/python/how-to-plot-a-density-map-in-python.html
#   ToDo: how do do this not explicitly specifying the times but using datetime to get
#   +/- around 21z each of the days in the output?
    xx = np.squeeze(lon[it,:])
    yy = np.squeeze(lat[it,:])
    x  = xx.flatten()
    y  = yy.flatten()
    k = gaussian_kde(np.vstack([x,y]))
#   Create a contour plat
    xi, yi = np.mgrid[x.min():x.max():100j,y.min():y.max():100j]
    zi_ = k(np.vstack([xi.flatten(),yi.flatten()]))
    zi = np.ma.masked_array(zi_, zi_<0.001)
    ax.contourf(xi,yi,zi.reshape(xi.shape),transform=ccrs.PlateCarree(), cmap=cmap,zorder=zorder)
    #print(np.min(zi),np.max(zi))
    return


def make_plot(parcel,vtime):
    timestr = parcel[-19:-3]
    firestr = parcel[0:-20]
    projLcc = ccrs.LambertConformal(central_longitude=-100, central_latitude=40)
    fig = plt.figure(figsize=(16,20))
    gs = GridSpec(2, 1, height_ratios=[3.5, 1])
    plt.subplots_adjust(left=0.05,bottom=0.05,right=0.95,top=0.95,hspace=0.02)
    
    ax  = fig.add_subplot(gs[0],projection=projLcc)
    ax.set_extent([-120,-70,22.5,60],crs=ccrs.PlateCarree())
#    ax.set_extent([240,290,22.5,70],crs=ccrs.PlateCarree())
    ax.coastlines(resolution="50m",zorder=100)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=2, color='brown')
    ax.add_feature(cfeature.BORDERS, edgecolor='black',linewidth=2,zorder=100)
    ax.add_feature(cfeature.STATES, linestyle='--', edgecolor='black', linewidth=1,zorder=100)
    provinc_bodr = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lines', scale="50m", facecolor='none', edgecolor='k',zorder=100)
    ax.add_feature(provinc_bodr, linestyle='--', linewidth=1, edgecolor="k", zorder=10)
    ax.set_title(f"Fire: {firestr}  Init Time: {timestr}:00  Valid Time: {vtime}",fontsize=24)
    ax.background_img(name='NE', resolution='high')

#   Some markers
    markernames = ['Great Falls','Boulder']
    x = [-111.3008,-105.118]
    y = [47.5053,39.909]
    labels = ['1','2']
    ax.plot(x,y,markersize=28, marker="o", color="k", transform=ccrs.PlateCarree(),linestyle='')
    ax.plot(x,y,markersize=24, marker="o", color="red", transform=ccrs.PlateCarree(),linestyle='')
    for i in np.arange(0,2):
        ax.text(x[i],y[i],labels[i],color="k", transform=ccrs.PlateCarree(),
                ha="center", va="center", size=16, fontweight='black')

#   Some range rings
#   ER-2
    for radius in [550,925,1295,1650,2000]:  #km
        n_samples = 80
        circles = Polygon(Geodesic().circle(x[0], y[0], radius*1000., n_samples=n_samples))
        feature = cfeature.ShapelyFeature([circles], ccrs.PlateCarree(), fc='None', ec="black", lw=2, linestyle="-",zorder=101)
        circle = ax.add_feature(feature)

#   GV
    for radius in [926,1389,1852]:  #km
        n_samples = 80
        circles = Polygon(Geodesic().circle(x[1], y[1], radius*1000., n_samples=n_samples))
        feature = cfeature.ShapelyFeature([circles], ccrs.PlateCarree(), fc='None', ec="white", lw=3, linestyle="-",zorder=101)
        circle = ax.add_feature(feature)
        
#   Add height axis
    ax2 = fig.add_subplot(gs[1])
    ax2.set_ylim(6,16)
    ax2.set_ylabel("Altitude [km]",fontsize=24)
    ax2.set_xlabel("Days Since Init",fontsize=24)
    ax2.set_title("Parcel Altitude Over Time",fontsize=26)
    ax2.tick_params(axis='x', labelsize=20)
    ax2.tick_params(axis='y', labelsize=20)
    
    return ax, ax2

if __name__ == "__main__":
    print("PETE")
