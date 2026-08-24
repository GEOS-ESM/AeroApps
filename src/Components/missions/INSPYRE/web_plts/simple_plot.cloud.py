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
    tokens = os.path.basename(filename).split('.')
    cldfilename = os.path.dirname(filename) + '/GEOS.fp.fcst.tavg1_2d_rad_Nx.'+tokens[4][0:23]+'30.V01.nc4'
    if os.path.exists(cldfilename):
        print(cldfilename)
    else:
        print(f"File does not exist: {cldfilename}")
        return

    cldh,lon,lat,lev,rc = model.read1(cldfilename,varn='CLDHGH')
    cldh = np.squeeze(cldh[np.where(lat>20),:])

    cldl,lon,lat,lev,rc = model.read1(cldfilename,varn='CLDLOW')
    cldl = np.squeeze(cldl[np.where(lat>20),:])

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

    cf3 = plt.contourf(lonx, latsx, cldl, [.99,10],colors='0',
                       transform=ccrs.PlateCarree(), alpha=.2)
    cf2 = plt.contourf(lonx, latsx, cldh, [.99,10],colors='1',
                       transform=ccrs.PlateCarree(), alpha=.6)
    cbar1=plt.colorbar(cf, ax=ax, orientation='vertical', 
                       pad=0.1, extend="max", ticks=np.log10([0.05,0.1,0.2,0.5,1,2]),
                       format=ticker.FixedFormatter(['0.05', '0.1', '0.2','0.5','1','2']))
    cbar1.ax.tick_params(labelsize=18)
    cbar1.set_label(label='%s'%(title),size=24)

    tokens = os.path.basename(filename).split('.')

    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.99,top=0.9)
    plt.title('GEOS FP %s %s'%(varn,tokens[4]),size=24)
    plt.savefig('fp.%s_wcloud.%s.png'%(varn,tokens[4]))


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
    
