#! /usr/bin/env python3
import os, sys
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import xarray as xr
import cartopy
import matplotlib.pyplot as plt
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib as mpl

if __name__ == '__main__':
    row = 85
    col = 220
    rootDir = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/SBG/SSD650'
    vFiles = sorted(glob(rootDir + '/LevelC/Y2006/M01/D16/ssd650-sbg-g5nr.lc.vlidort.20060116_173[0-9]z.nc4'))
    brdfFiles = sorted(glob(rootDir + '/LevelB/Y2006/M01/D16/ssd650-g5nr.lb2.ames_brdf.20060116_173[0-9]z.nc4'))
    aerFiles = sorted(glob(rootDir + '/LevelB/Y2006/M01/D16/ssd650-g5nr.lb2.aer_Nv.20060116_173[0-9]z.nc4'))
    locFiles = sorted(glob(rootDir + '/LevelB/Y2006/M01/D16/ssd650-g5nr.lb2.sbg.20060116_173[0-9]z.nc4'))

    vlid = xr.open_mfdataset(vFiles,chunks='auto')
    brdf = xr.open_mfdataset(brdfFiles,chunks='auto')
    aerm = xr.open_mfdataset(aerFiles,chunks='auto')
    swath = xr.open_mfdataset(locFiles)

    channels = brdf.nwav.values*1e3


    # plot the surface reflectance
    transform=cartopy.crs.PlateCarree()
    cmap = mpl.cm.jet

    slons = swath.longitude.squeeze().values
    slats = swath.latitude.squeeze().values    

    # surface reflctance
    vname = 'surf_reflectance'
    sdata = vlid[vname]
    ich = 75
    sdata = sdata.isel(ch=ich)    
    condition = sdata != -500
    # sets all -500 to nan
    sdata = sdata.where(condition)

    bounds = np.linspace(sdata.min(),sdata.max(),100)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

    # get one spectrum
    itime = 2200
    iacross = 100
    spec = vlid[vname]
    spec = spec.isel(time=itime,across=iacross)


    projection=cartopy.crs.Orthographic(central_longitude=slons.mean(),central_latitude=slats.mean())
    axes_class = (GeoAxes,
                  dict(projection=projection))
    fig = plt.figure(figsize=[8,10])
    axgr = AxesGrid(fig,  (0.1, 0.5, 0.8, 0.4) , axes_class=axes_class,
            nrows_ncols=(1, 1),
            axes_pad=0.6,
            cbar_location='right',
            cbar_mode='single',
            cbar_pad=0.2,
            cbar_size='3%',
            label_mode='keep')  

    ax = axgr[0]
    ax.pcolormesh(slons,slats,sdata,transform=transform,norm=norm,cmap=cmap)
    ax.set_title(vname+' '+str(channels[ich])+' nm')
    ax.coastlines()
    ax.set_global()

    # colorbar
    cb = axgr.cbar_axes[0].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label=vname,ticks=bounds[::10])
#                            cb.ax.set_yticklabels(ticks.astype(str))


    # zoom in 
    # ------------------------------------
    projection=cartopy.crs.Orthographic(central_longitude=slons.mean(),central_latitude=slats.mean())
    axes_class = (GeoAxes,
                  dict(projection=projection))
    axgr = AxesGrid(fig,  (0.01, 0.01, 0.4, 0.45) , axes_class=axes_class,
            nrows_ncols=(1, 1),
            axes_pad=0.6,
            cbar_location='right',
            cbar_mode='single',
            cbar_pad=0.2,
            cbar_size='4%',
            label_mode='keep')  # note the empty label_mode

    ax = axgr[0]
    ax.pcolormesh(slons,slats,sdata,transform=transform,norm=norm,cmap=cmap)
    ax.set_title(vname+' '+str(channels[ich])+' nm')
    ax.coastlines()
    ax.add_feature(cartopy.feature.STATES)

    # label point
    ax.plot(slons[itime,iacross],slats[itime,iacross],marker='*',color='tab:red',markersize=5,transform=transform)

    # colorbar
    cb = axgr.cbar_axes[0].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label=vname,ticks=bounds[::10])
#                            cb.ax.set_yticklabels(ticks.astype(str))


    # plot spectrum
    # -----------------
    ax = fig.add_axes([0.5, 0.05, 0.45, 0.40])
    ax.plot(channels,spec,label=vname)
    ax.plot(channels,brdf['Ki'].isel(time=itime,ncross=iacross,nalong=0),label='Riso')
    ax.plot(channels,brdf['Kg'].isel(time=itime,ncross=iacross,nalong=0),label='Rgeo')
    ax.plot(channels,brdf['Kv'].isel(time=itime,ncross=iacross,nalong=0),label='Rvol')
    ax.plot(channels,vlid['toa_reflectance'].isel(time=itime,across=iacross),label='TOA Ref')
    ax.set_xlabel('wavelength')
    ax.legend()

    plt.savefig("swath_{}_{}.png".format(vname,channels[ich]))
    sys.exit()

