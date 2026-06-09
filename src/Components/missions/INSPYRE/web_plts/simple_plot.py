#!/usr/bin/env python
from optparse   import OptionParser   # Command-line args
import os
import model
import numpy as np
import datetime
import numpy.ma as ma
from make_plot import make_plot

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

    ax = make_plot()

#    clevs = np.arange(0,cbarmax,cbarmax/40)
#   Plot on log scale
    clevs = np.arange(-1.5,0.5,0.05)
    cf = plt.contourf(lonx, latsx, np.log10(ext), clevs, cmap='YlOrBr',
                      extend='both',transform=ccrs.PlateCarree())

    cbar1=plt.colorbar(cf, ax=ax, orientation='vertical', 
                       pad=0.1, extend="max", ticks=np.log10([0.05,0.1,0.2,0.5,1,2]),
                       format=ticker.FixedFormatter(['0.05', '0.1', '0.2','0.5','1','2']))
    cbar1.ax.tick_params(labelsize=12)
    cbar1.set_label(label='%s'%(title),size=12)

    tokens = os.path.basename(filename).split('.')

    plt.title('GEOS FP %s %s'%(varn,tokens[4]))
    plt.savefig('fp.%s.%s.png'%(varn,tokens[4]))


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
    plot(varn='BREXTTAU',cbarmax=.5,title='Brown Carbon AOT',filename=filename)
    plot(varn='OCEXTTAU',cbarmax=.5,title='Organic Carbon AOT',filename=filename)
    plot(varn='SSEXTTAU',cbarmax=.2,title='Sea Salt AOT',filename=filename)
    plot(varn='DUEXTTAU',cbarmax=.2,title='Dust AOT',filename=filename)
    plot(varn='SUEXTTAU',cbarmax=.2,title='Sulfate AOT',filename=filename)
    
