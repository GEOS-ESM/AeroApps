#!/usr/bin/env python3

"""
This provides tools to read model and obs AOD fields and plot
simple three panel difference plot (A, B, B-A) or four panel
"closeness" plot.
"""


import numpy as np
import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

month = ["January","February","March","April","May","June","July",
         "August","September","October","November","December"]

def plotcloseness(yy, mm, lon, lat, baseline, perturb, obs, baselinen, perturbn, sat,
                  extent=[-180,180,-90,90], varn="TOTEXTTAU550", saveplot=True,box=None):
    """
    Plots a 4-panel figure comparing:
    Panel A: Baseline model aod
    Panel B: Observed AOD
    Panel C: Perturbed model aod
    Panel D: Closeness metric = |perturb - obs AOD| - |baseline - obs AOD|

    Parameters:
        lon, lat (2D or 1D arrays): Longitude and Latitude data for the global grid
        baseline (array): Baseline model data (AOD)
        obs (array): Observational data (AOD)
        perturb (array): Perturbed model data (AOD)
        baselinen (str): Baseline experiment name
        perturbn (str):  Perturbed experiment name
        sat (str):       Satellite ID (MOD04, MYD04)
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    
    # Compute Panel D: Closeness metric
    closeness = np.abs(perturb - obs) - np.abs(baseline - obs)
    
    # Create the figure with a 2x2 grid of panels
    fig, axes = plt.subplots(2, 2, figsize=(10, 11), subplot_kw={'projection': ccrs.PlateCarree()})

    axes[0,0].set_extent(extent)
    axes[0,1].set_extent(extent)
    axes[1,0].set_extent(extent)
    axes[1,1].set_extent(extent)

    # Common color settings for AOD plots
    vmin_aod = 0
    vmax_aod = 1
    vmin_closeness = -0.1
    vmax_closeness = 0.1
    cmap_aod = 'Spectral_r'
    cmap_closeness = 'RdBu_r'
    
    # --- Panel A: Baseline Model AOD ---
    im_a = axes[0, 0].pcolormesh(lon, lat, baseline, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[0, 0].set_title(baselinen, fontsize=12)
    axes[0, 0].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)

    # --- Panel B: Observed AOD ---
    im_b = axes[0, 1].pcolormesh(lon, lat, obs, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[0, 1].set_title(f"MODIS NNR ({sat})", fontsize=12)
    axes[0, 1].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)

    # --- Panel C: Perturbed Model AOD ---
    im_c = axes[1, 0].pcolormesh(lon, lat, perturb, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[1, 0].set_title(perturbn, fontsize=12)
    axes[1, 0].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)

    # --- Panel D: Closeness Analysis ---
    im_d = axes[1, 1].pcolormesh(lon, lat, closeness, cmap=cmap_closeness, vmin=vmin_closeness, vmax=vmax_closeness, transform=ccrs.PlateCarree())
    axes[1, 1].set_title(f"Closeness\n|{perturbn}-obs|-|{baselinen}-obs|", fontsize=12)
    axes[1, 1].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)
    
    # Add colorbars for each row
    cax_c = fig.add_axes([0.15, 0.05, 0.3, 0.02])  # Colorbar for AOD rows
    cbar_c = fig.colorbar(im_c, cax=cax_c, orientation='horizontal', extend='both')
    cbar_c.set_label("AOD at 550 nm")

    cax_d = fig.add_axes([0.575, 0.05, 0.3, 0.02])  # Colorbar for closeness row
    cbar_d = fig.colorbar(im_d, cax=cax_d, orientation='horizontal', extend='both')
    cbar_d.set_label(f"Blue: {perturbn} closer\nRed: {baselinen} closer")

    if box != None:
        from matplotlib.patches import Rectangle
        ll = (box[0],box[2])
        wd = box[1]-box[0]
        ht = box[3]-box[2]
        rect0 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        rect1 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        rect2 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        rect3 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        axes[0,0].add_patch(rect0)
        axes[0,1].add_patch(rect1)
        axes[1,0].add_patch(rect2)
        axes[1,1].add_patch(rect3)
    
    # Main title for the figure
    monlab=None
    if(len(mm) == 12):
        monlab = "(Annual)"
        montit = "annual"
    elif(len(mm) == 1):
        monlab = f"({month[mm[0]-1]})"
        montit = f"{mm[0]:02d}"
    plt.suptitle(f"{varn}: {yy:04d} {monlab}", fontsize=14, y=0.95)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    if saveplot:
        plt.savefig(f"{varn}_closeness.{baselinen}_{perturbn}.{sat}.{yy:04d}{montit}.png")
#        plt.show()
    else:
        plt.show()


def plotdifference(yy, mm, lon, lat, var1, var2, var1n, var2n,
                   extent=[-180,180,-90,90],box=None, 
                   cbartitle="AOD at 550 nm", varn="TOTEXTTAU550",
                   prange=[0,.8],drange=[-0.2,0.2],saveplot=True,proj="PlateCarree"):
    """
    Simple difference plot of two input maps (could be model versus satellite 
    or model versus model)
    """
    if proj == "nps":
        fig, axes = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-30)})
        extent = [-90,30,50,90]
        axes[0].set_extent(extent,ccrs.PlateCarree())
        axes[1].set_extent(extent,ccrs.PlateCarree())
        axes[2].set_extent(extent,ccrs.PlateCarree())
    else:
        fig, axes = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})
        axes[0].set_extent(extent,ccrs.PlateCarree())
        axes[1].set_extent(extent,ccrs.PlateCarree())
        axes[2].set_extent(extent,ccrs.PlateCarree())
    
    # Common settings
    vmin_aod  = prange[0]    # Min value for AOD
    vmax_aod  = prange[1]    # Max value for AOD
    vmin_diff = drange[0]    # Min value for difference
    vmax_diff = drange[1]    # Max value for difference
    cmap_aod  = 'Spectral_r' # Colormap for AOD
    cmap_diff = 'RdBu_r'     # Colormap for difference (bias)
    
    # Plot Baseline
    im1 = axes[0].pcolormesh(lon,lat,var1, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[0].set_title(var1n)
    axes[0].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2)
    axes[0].add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    axes[0].gridlines(draw_labels=True)
    
    # Plot Perturb
    im2 = axes[1].pcolormesh(lon,lat,var2, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[1].set_title(var2n)
    axes[1].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2)
    axes[1].add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    axes[1].gridlines(draw_labels=True)
    
    # Compute Bias 
    bias = var2-var1
    # Plot Bias 
    im3 = axes[2].pcolormesh(lon,lat,bias, cmap=cmap_diff, vmin=vmin_diff, vmax=vmax_diff, transform=ccrs.PlateCarree())
    axes[2].set_title(f"Difference ({var2n} - {var1n})")
    axes[2].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2)
    axes[2].add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    axes[2].gridlines(draw_labels=True)

    if box != None:
        from matplotlib.patches import Rectangle
        ll = (box[0],box[2])
        wd = box[1]-box[0]
        ht = box[3]-box[2]
        rect0 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        rect1 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        rect2 = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        axes[0].add_patch(rect0)
        axes[1].add_patch(rect1)
        axes[2].add_patch(rect2)
    
    
    # Add Colorbars for Each Plot
    cbar1 = fig.colorbar(im1, ax=axes[0], orientation='horizontal',
                         pad=0.05, aspect=30,extend='both',shrink=0.7,)
    cbar1.set_label(cbartitle)
    
    cbar2 = fig.colorbar(im2, ax=axes[1], orientation='horizontal',
                         pad=0.05, aspect=30,extend='both',shrink=0.7,)
    cbar2.set_label(cbartitle)
    
    cbar3 = fig.colorbar(im3, ax=axes[2], orientation='horizontal',
                         pad=0.05, aspect=30,extend='both',shrink=0.7,)
    cbar3.set_label("Difference")

    # Main title for the figure
    monlab=None
    if(len(mm) == 12):
        monlab = "(Annual)"
        montit = "annual"
    elif(len(mm) == 1):
        monlab = f"({month[mm[0]-1]})"
        montit = f"{mm[0]:02d}"
    plt.suptitle(f"{varn}: {yy:04d} {monlab}", fontsize=14, y=0.95)
    
    # Display the plot
    plt.tight_layout()
    if saveplot:
        plt.savefig(f"{varn}_difference.{var1n}_{var2n}.{yy:04d}{montit}.png")
#        plt.show()
    else:
        plt.show()


def plotone(lon, lat, var, varn,
            extent=[-180,180,-90,90],box=None,
            prange=[0,1],
            cbartitle=None, title=None,fname=None,
            var2=None,proj="PlateCarree"):
    """
    Simple plot of single input map (could be model or satellite)
    """
    if proj == "nps":
        fig, axes = plt.subplots(1, 1, figsize=(8, 6), subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-30)})
        extent = [-90,30,50,90]
        axes.set_extent(extent,ccrs.PlateCarree())
    else:
        fig, axes = plt.subplots(1, 1, figsize=(32, 24), subplot_kw={'projection': ccrs.PlateCarree()})
        axes.set_extent(extent,ccrs.PlateCarree())

    
    # Common settings
    vmin_aod = prange[0]  # Min value for AOD
    vmax_aod = prange[1]  # Max value for AOD
    cmap_aod = 'Spectral_r'  # Colormap for AOD
    cmap_aod = "RdYlGn_r"
    cmap_cnt = 'grey'
    
    # Plot Baseline
    im1 = axes.pcolormesh(lon,lat,var, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes.set_title(varn)
    axes.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2)
    axes.set_extent(extent)
    axes.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    axes.gridlines(draw_labels=True)
#    if var2[0] != None:
#    im2 = axes.contour(lon, lat, var2, [40], colors='red',linewidths=1,transform=ccrs.PlateCarree())
#    im3 = axes.pcolormesh(lon, lat, var2, cmap=cmap_cnt, vmin=39.9,vmax=40.,transform=ccrs.PlateCarree(),alpha=.2)
#    im2 = axes.contour(lon, lat, var2, [0.1,0.5], colors='red',linewidths=1,transform=ccrs.PlateCarree())
#    im3 = axes.pcolormesh(lon, lat, var2, cmap=cmap_cnt, vmin=0.01,vmax=.1,transform=ccrs.PlateCarree(),alpha=.2)
    if box != None:
        from matplotlib.patches import Rectangle
        ll = (box[0],box[2])
        wd = box[1]-box[0]
        ht = box[3]-box[2]
        rect = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
        axes.add_patch(rect)
    
    
    # Add Colorbars for Each Plot
    cbar1 = fig.colorbar(im1, ax=axes, orientation='horizontal',
                         pad=0.05, aspect=30,extend='both',shrink=0.7,)
    if cbartitle:
        cbar1.set_label(cbartitle)
    if title:
        plt.suptitle(title, fontsize=14, y=0.95)
    
    # Display the plot
    plt.tight_layout()
    if fname:
        plt.savefig(fname)
#        plt.show()
    else:
        plt.show()


def read_model(EXPID,EXPDIR,yy,mm,sat,varn="TOTEXTTAU550"):
    varout = []
    lat, lon = None, None
    filen = []
    
    for m in mm:
        filen_ = f"{EXPDIR}/holding/inst2d_hwl_x/{yy:04d}{m:02d}/{EXPID}.inst2d_hwl_x.monthly.{sat}.{yy:04d}{m:02d}.nc4"
        if varn == "DUCMASS":
            filen_ = f"{EXPDIR}/holding/tavg2d_aer_x/{yy:04d}{m:02d}/{EXPID}.tavg2d_aer_x.monthly.{yy:04d}{m:02d}.nc4"
        if not os.path.exists(filen_):
            print(f"File not found: {filen_}. Skipping...")
            continue
        filen.append(filen_)
#        print(filen)
    
    ds = xr.open_mfdataset(filen)
    varout = ds[varn].mean(dim='time')
    lat = ds['lat'].values
    lon = ds['lon'].values
        
    return varout, lon, lat

def read_sat(sat,yy,mm,varn="tau_",weighted=False):
    varout = []
    lat, lon = None, None
    filen = []
    
    for m in mm:
        filen_ = f"/home/pcolarco/geos_aerosols/pcolarco/{sat}/Level3/Y{yy:04d}/M{m:02d}/nnr_003.{sat}_L3a.blend.monthly.{yy:04d}{m:02d}.nc4"
        if weighted:
            filen_ = f"/home/pcolarco/geos_aerosols/pcolarco/{sat}/Level3/Y{yy:04d}/M{m:02d}/nnr_003.{sat}_L3a.blend.monthly.weighted.{yy:04d}{m:02d}.nc4"
        if not os.path.exists(filen_):
            print(f"File not found: {filen_}. Skipping...")
            continue
        filen.append(filen_)
#        print(filen)
    
    ds = xr.open_mfdataset(filen)
    varout = ds[varn].mean(dim='time')
    lat = ds['lat'].values
    lon = ds['lon'].values
        
    return varout, lon, lat
    

if __name__ == "__main__":

#   Example: annual mean differences

    baseline = "c180R_v11.8.0_develop"
    perturb  = "c180R_v11.8.0_newdust"
    sat      = "MYD04"
    yy       = 2019
    mm       = [1,2,3,4,5,6,7,8,9,10,11,12]
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
    pertdir  = f"/home/pcolarco/geos_aerosols/pcolarco/{perturb}/"

    baseaod, lon, lat  = read_model(baseline,basedir,yy,mm,sat)
    pertaod, lon, lat  = read_model(perturb ,pertdir,yy,mm,sat)
    obsaod, lono, lato = read_sat(sat,yy,mm)
#    plotdifference(yy, mm, lon, lat, baseaod.values,baseline,pertaod.values,perturb)
    obsaodi = obsaod.interp_like(baseaod,method="nearest",assume_sorted=True)
#    plotdifference(yy, mm, lon, lat, baseaod.values,baseline,np.squeeze(obsaodi.values),sat)
    plotcloseness(yy, mm, lon, lat, baseaod.values, pertaod.values, np.squeeze(obsaodi.values),
                  baseline,perturb,sat)
