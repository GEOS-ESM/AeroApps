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
from pyobs.icartt import ICARTT


os.environ['CARTOPY_USER_BACKGROUNDS'] = "/home/pcolarco/silo/python/"



# Get the model configuration
def get_model_configuration(ictFile,model='m21c',collection='aer_inst_3hr_glo_Nv'):
    modname   = model
    cs        = False
    do_optics = True
    campaign  = "CPEX-AW"

    # Default
    fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.prog.eta.ddf',
              './'+model+'.inst2d_hwl_x.ddf']
    if(model == "m21c"):
        fpdata = [f"{model}.{collection}.ddf",
                  f"{model}.asm_inst_3hr_glo_Nv.ddf",
                  f"{model}.aer_inst_1hr_glo_L1152x721_slv.ddf"]
        collname = collection
        modname = "MERRA-21C"
        config = "./g2g_pm25.yaml"
    if(model == 'fp'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf',
                  './'+model+'.inst1_2d_hwl_Nx.ddf']
        if(collection == "inst3_3d_aer_Nv"):
            collname = "inst3-3d-AER-Nv"
        if(collection == "tavg1_2d_lfo_Nx"):
            do_optics = False
            collname = "tavg1-2d-LFO-Nx"
        modname = "GEOS-FP"
        config = './fp_pm25.yaml'
    if(model == 'MERRA2'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf']
        config = './m2_pm25.yaml'
        collname = collection
        if(collection == 'tavg1_2d_aer_Nx'):
            fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.tavg1_2d_adg_Nx.ddf']
            do_optics = False
            collname = collection
    
#   Get the ICARTT file describing the trajectory
    if ictFile.find('DC8')  > 0:
        aircraft = 'DC8'
        i0 = ictFile.find('DC8')+4
    m = ICARTT(ictFile)
    yyyymmdd  = ictFile[i0:-4]
    dateout   = ictFile[i0:i0+4]+"-"+ictFile[i0+4:i0+6]+"-"+ictFile[i0+6:i0+8]
    return yyyymmdd, dateout, modname, aircraft, campaign, cs, do_optics, fpdata, config, collname






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
    projLcc = ccrs.LambertConformal(central_longitude=-20, central_latitude=20)
    fig = plt.figure(figsize=(16,26))
    gs = GridSpec(2, 1, height_ratios=[3.5, 1])
    plt.subplots_adjust(left=0.05,bottom=0.05,right=0.95,top=0.95,hspace=0.05)
    
    ax  = fig.add_subplot(gs[0],projection=projLcc)
    ax.set_extent([-40,10,0.,40],crs=ccrs.PlateCarree())
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
