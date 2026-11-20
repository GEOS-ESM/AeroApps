#!/usr/bin/env python3

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
from pyobs.icartt import ICARTT
from optparse   import OptionParser   # Command-line args
from datetime import datetime, timedelta
import os


os.environ['CARTOPY_USER_BACKGROUNDS'] = "/home/pcolarco/silo/python/"

def make_plot():
    projLcc = ccrs.PlateCarree()
    fig = plt.figure(figsize=(14,12))
    plt.subplots_adjust(left=0.05,bottom=0.05,right=0.95,top=0.95,hspace=0.05)
    
    ax  = fig.add_subplot(projection=projLcc)
    ax.set_extent([-45,-10,0.,30],crs=ccrs.PlateCarree())
    ax.coastlines(resolution="50m",zorder=99)
    cg = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=1, color='grey')
    ax.add_feature(cfeature.BORDERS, edgecolor='black',linewidth=2,zorder=99)
    ax.add_feature(cfeature.STATES, linestyle='--', edgecolor='black', linewidth=1,zorder=99)
    provinc_bodr = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lines', scale="50m",
                                                facecolor='none', edgecolor='k',zorder=99)
    ax.add_feature(provinc_bodr, linestyle='--', linewidth=1, edgecolor="k", zorder=10)
    ax.background_img(name='NE', resolution='higher')
    cg.xlabel_style = {"size": 22}
    cg.ylabel_style = {"size": 22}

    return ax



if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] ictFile",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        ictDir = args[0]
    else:
        parser.error("must have 1 argument: directory")

    ax = make_plot()

    ictFiles = []
    for x in sorted(os.listdir(ictDir)):
        if x.endswith(".ict"):
            # Prints only text file present in My Folder
            print(x)
            ictFiles.append(ictDir+x)
    for ictFile in ictFiles:
        m = ICARTT(ictFile)
        alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
        cd = ax.plot(lon,lat,lw=3,transform=ccrs.PlateCarree(),zorder=99,label=ictFile[-15:-7])
    ax.legend(loc="lower left",fontsize=24)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(24)
    plt.savefig("CPEX-CV_tracks.png")
