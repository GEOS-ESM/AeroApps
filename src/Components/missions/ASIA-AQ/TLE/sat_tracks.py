import matplotlib.pyplot as plt

import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pyobs.tle import TLE
from datetime  import datetime, timedelta

TLEs = dict( Aqua = 'aqua.tle',
             NOAA20 = 'noaa20.tle',
             Sentinel5P = 's5p.tle',
             SNPP = 'snpp.tle',
             Terra = 'terra.tle'
            )


def trim(t, x, y, bbox,t_offset):

    I = (x>=bbox[0])&(x<=bbox[1])&(y>=bbox[2])&(y<=bbox[3])
    t[~I] = np.nan
    x[~I] = np.nan
    y[~I] = np.nan
    ts, xs, ys = [], [], []
    i = 0
    t_, x_, y_ = t[I], x[I], y[I]
    for i in range(len(t_)):
        ts.append((t_[i]+t_offset).isoformat().split('T')[1][0:5]+'L')
        xs.append(x_[i])
        ys.append(y_[i])
        
        
    return t, x, y, ts, xs, ys
    
def plot_traj(year,month,day,t_offset):
    
    fig = plt.figure(figsize=[15, 15])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    bbox = [90, 130, -6, 25] # Southeast Asia
    ax.set_extent(bbox, crs=ccrs.PlateCarree())

    # Put a background image on for nice sea rendering.
    ax.stock_img()

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)

    dp = 0.2 # delta pixel for text
    for sat in TLEs:
        tle = TLE(TLEs[sat])
        t1, t2, dt = datetime(year,month,day,0), datetime(year,month,day,12), timedelta(seconds=60)
        t, x, y = tle.getSubpoint(t1,t2,dt)
        t, x, y, ts, xs, ys = trim(t,x,y,bbox,t_offset)
        ax.plot(x,y,'-o',label=sat)
        for i in range(len(ts)):
           ax.text(xs[i]+dp,ys[i]-dp,ts[i])


    ax.legend(loc='lower left')
    ax.title.set_text('Satellite Tracks for %0d-%02d-%02d'%(year,month,day))            
    #plt.show()

    fname = 'sattrack_%d-%02d-%02d'%(year,month,day)+'.png'
    print('Saving',fname)
    plt.savefig(fname,bbox_inches='tight')
    
    return

if __name__ == '__main__':
    t_offset = timedelta(hours = 8) # Manila offset
    for day in range(5,16):
        plot_traj(2024,2,day,t_offset)
