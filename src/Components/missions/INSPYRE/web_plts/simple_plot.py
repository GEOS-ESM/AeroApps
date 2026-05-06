#!/usr/bin/env python
from optparse   import OptionParser   # Command-line args
import os
import model
import numpy as np
import datetime
import numpy.ma as ma

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

    fig = plt.figure(figsize=(14,9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='black', linewidth=0.5)
    
    ax.set_xlim([-140,-80])
    ax.set_ylim([20,70])

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

#   Some markers
    markernames = ['Great Falls','Boulder']
    x = [-111.3008,-105.118]
    y = [47.5053,39.909]
    labels = ['1','2']
    ax.plot(x,y,markersize=14, marker="o", color="k", transform=ccrs.PlateCarree(),linestyle='')
    ax.plot(x,y,markersize=12, marker="o", color="red", transform=ccrs.PlateCarree(),linestyle='')
    for i in np.arange(0,2):
        ax.text(x[i],y[i],labels[i],color="k", transform=ccrs.PlateCarree(),
                ha="center", va="center", size=8, fontweight='black')


#   Some range rings
#   ER-2
    for radius in [550,925,1295,1650,2000]:  #km
        n_samples = 80

        circles = Polygon(Geodesic().circle(x[0], y[0], radius*1000., n_samples=n_samples))
        feature = cfeature.ShapelyFeature([circles], ccrs.PlateCarree(), fc='None', ec="black", lw=1, linestyle="-")

        circle = ax.add_feature(feature)

#   GV
    for radius in [926,1389,1852]:  #km
        n_samples = 80

        circles = Polygon(Geodesic().circle(x[1], y[1], radius*1000., n_samples=n_samples))
        feature = cfeature.ShapelyFeature([circles], ccrs.PlateCarree(), fc='None', ec="black", lw=1, linestyle="dotted")

        circle = ax.add_feature(feature)

#   Get the buoy locations from file
#    labels = ['J','K','L','M','N','O','P','Q','R']
    '''
    labels = ['J','L','M','N','O','P','Q','R']
    i = 0
    with open('buoys.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in spamreader:
            xb = float(row[2])
            yb = float(row[3])
            ax.plot(xb,yb,markersize=14, marker="o", color="w", 
                    transform=ccrs.PlateCarree(),linestyle='')
            ax.plot(xb,yb,markersize=12, marker="o", color="k", 
                    transform=ccrs.PlateCarree(),linestyle='')
            ax.text(xb,yb,labels[i],color="w", transform=ccrs.PlateCarree(),
                    ha="center", va="center", size=8, fontweight='extra bold')
            i = i+1
    '''
    
    tokens = os.path.basename(filename).split('.')

    plt.title('GEOS FP %s %s'%(varn,tokens[4]))
    fig.savefig('fp.%s.%s.png'%(varn,tokens[4]))
    plt.close(fig)


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
    
