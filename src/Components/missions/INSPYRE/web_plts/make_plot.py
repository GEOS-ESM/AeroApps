#!/usr/bin/env python

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
import os
import matplotlib as mpl
mpl.use('Agg')

os.environ['CARTOPY_USER_BACKGROUNDS'] = "/home/pcolarco/silo/python/"

def make_plot():
    projLcc = ccrs.LambertConformal(central_longitude=-100, central_latitude=40)
    fig = plt.figure(figsize=(14,12))
    gs = GridSpec(1, 1)

    ax  = fig.add_subplot(gs[0],projection=projLcc)
    ax.set_extent([-130,-70,22.5,70],crs=ccrs.PlateCarree())
    ax.coastlines(resolution="50m")
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax.add_feature(cfeature.BORDERS, edgecolor='black',linewidth=2)
    ax.add_feature(cfeature.STATES, linestyle='--', edgecolor='black', linewidth=1)
    provinc_bodr = cfeature.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lines', scale="50m", facecolor='none', edgecolor='k')
    ax.add_feature(provinc_bodr, linestyle='--', linewidth=1, edgecolor="k", zorder=10)

#    ax.stock_img()
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
        feature = cfeature.ShapelyFeature([circles], ccrs.PlateCarree(), fc='None', ec="limegreen", lw=3, linestyle="-",zorder=101)
        circle = ax.add_feature(feature)
        
    return ax

def make_plot_mercator():
    proj = ccrs.epsg(3857)
    fig = plt.figure(figsize=(19,24))
    gs = GridSpec(1, 1)

    ax  = fig.add_subplot(gs[0],projection=proj)
    ax.set_extent([-130,-70,22.5,70],crs=ccrs.PlateCarree())
    ax.coastlines(resolution="50m")
    ax.gridlines(dms=True, x_inline=False, y_inline=False)
    ax.add_feature(cfeature.BORDERS, edgecolor='black',linewidth=1)
    ax.add_feature(cfeature.STATES, linestyle='--', edgecolor='black', linewidth=1)
    provinc_bodr = cfeature.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lines', scale="50m",
                   facecolor='none', edgecolor='k')
    ax.add_feature(provinc_bodr, linestyle='--', linewidth=1, edgecolor="k", zorder=10)
    ax.axis("off")
        
    return ax

if __name__ == "__main__":
    print("PETE")
