#!/usr/bin/env python
from optparse   import OptionParser   # Command-line args
import os
import model
import numpy as np
import datetime
import numpy.ma as ma

import csv

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import matplotlib.ticker as ticker
import cartopy.crs as ccrs


def plot(lev='0',varn='TOTEXTTAU',cbarmax=1.,title='Total AOT', 
         filename='/dev/null'):

    levz = -1
    print(varn,filename)

    ext,lon,lat,lev,rc = model.read1(filename,varn=varn)
    ext = np.squeeze(ext[np.where(lat>45),:])
    lat=lat[np.where(lat>45)]
    lonx, latsx = np.meshgrid(lon,lat)

    fig = plt.figure(figsize=(12,9))
    globe = ccrs.Globe(ellipse=None,
                       semimajor_axis=6370000,
                       semiminor_axis=6370000)
    projection=ccrs.RotatedPole(pole_longitude=-30 - 180,
                      pole_latitude=36,
                      globe=globe)
    ax = plt.axes(projection=projection)
    ax.coastlines(resolution='50m')
    ax.gridlines()
#    ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    xs, ys, zs = projection.transform_points(ccrs.PlateCarree(),
                                         np.array([-70, 80.]),
                                         np.array([50, 80])).T
    ax.set_xlim(xs)
    ax.set_ylim(ys)

    clevs = np.arange(0,cbarmax,cbarmax/40)
    cf = plt.contourf(lonx, latsx, ext, clevs, cmap='Spectral_r',
                      extend='both',transform=ccrs.PlateCarree())

    cbar1=plt.colorbar(cf, ax=ax, orientation='vertical', 
                       pad=0.1, extend="max")
    cbar1.ax.tick_params(labelsize=12)
    cbar1.set_label(label='%s'%(title),size=12)

#   Some markers
    markernames = ['Alert','Baffin Bay','Belle','Eureka','Fram Strait',
                   'Kaffeklubben Island','Lisee','Northern Arctic Ocean',
                   'Pituffik','Resolute Bay','Station Nord','Summit Greenland',
                   'Svalbard Zeppelin']
    x = [-62.51,-67.40,87,-86.42,0,-31.2,28,-60,-68.75,-94.97,-16.6,-38.46,11.88]
    y = [82.45,74.60,87,80.05,78.92,83.62,87,88,76.53,74.71,81.60,72.58,78.90]
    labels = ['1','2','3','4','5','6','7','8','9','10','11','12','13']
    ax.plot(x,y,markersize=14, marker="o", color="k", transform=ccrs.PlateCarree(),linestyle='')
    ax.plot(x,y,markersize=12, marker="o", color="red", transform=ccrs.PlateCarree(),linestyle='')
    for i in np.arange(0,13):
        ax.text(x[i],y[i],labels[i],color="k", transform=ccrs.PlateCarree(),
                ha="center", va="center", size=8, fontweight='black')


#   Get the buoy locations from file
#    labels = ['J','K','L','M','N','O','P','Q','R']
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
    plot(varn='OCEXTTAU',cbarmax=.5,title='Organic Carbon AOT',filename=filename)
    plot(varn='SSEXTTAU',cbarmax=.2,title='Sea Salt AOT',filename=filename)
    plot(varn='DUEXTTAU',cbarmax=.2,title='Dust AOT',filename=filename)
    plot(varn='SUEXTTAU',cbarmax=.2,title='Sulfate AOT',filename=filename)
    
