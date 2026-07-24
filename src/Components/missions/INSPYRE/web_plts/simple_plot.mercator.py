#!/usr/bin/env python
from optparse   import OptionParser   # Command-line args
import os
import model
import numpy as np
import datetime
import numpy.ma as ma
from make_plot import make_plot_mercator

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon
from cartopy.geodesic import Geodesic  # from geographiclib


def plot(lev='0',varn='TOTEXTTAU',cbarmax=2.,title='Total AOT', 
         filename='/dev/null'):

    levz = -1
    print(varn,filename)

    ext,lon,lat,lev,rc = model.read1(filename,varn=varn)
    ext = np.squeeze(ext[np.where(lat>20),:])
    lat=lat[np.where(lat>20)]
    lonx, latsx = np.meshgrid(lon,lat)

    ax = make_plot_mercator()

#   Plot on log scale
    clevs = np.arange(-1.5,0.5,0.05)
    cf = plt.contourf(lonx, latsx, np.log10(ext), clevs, cmap='YlOrBr',
                      extend='both',transform=ccrs.PlateCarree())

    tokens = os.path.basename(filename).split('.')
    ax.set_position([0,0,1,1])
    plt.savefig('fp.%s.%s.mercator.png'%(varn,tokens[4]),pad_inches=0)


if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] filename",
                          version='omi_level2a-1.0.0' )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        filename = args[0]
    else:
        parser.error("must have 1 argument: filename")
    print('main:',filename)
    plot(filename=filename)
    
