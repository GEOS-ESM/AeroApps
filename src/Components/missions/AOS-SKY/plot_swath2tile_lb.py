#!/usr/bin/env python3

"""
    Plot the input and output of swath2tile
    Make sure it looks sane

    P. Castellanos 2024

"""

import os
import sys

import argparse
from   dateutil.parser import parse         as isoparser
from   datetime        import datetime, timedelta
import numpy           as np
from   MAPL.config     import Config
import xarray          as xr
from   glob            import glob
import cartopy
import matplotlib.pyplot as plt
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib as mpl

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    algo   = 'linear'
    dt_days = 1
    col = 'polar07'

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("tile",
                        help="SBG tile")

    parser.add_argument("orbit",
                        help="orbit name")

    parser.add_argument("swath_pcf",
                        help="swath filenames pcf input file")

    parser.add_argument("tile_pcf",
                        help="tile filenames pcf input file")

    parser.add_argument("var",
                        help="variable name: time, vza, vaa, sza, saa, or scatAngle")

    parser.add_argument("start_isotime",
                        help="start isotime")

    parser.add_argument("end_isotime",
                        help="end isotime")

    parser.add_argument("--col",default=col,
              help="collection name - only works with polar07 right now (default=%s)"\
                        %col)

    parser.add_argument("-d", "--dt_days", default=dt_days,type=int,
              help="Timesetp in hours for swath files (default=%s)"%dt_days )


    args = parser.parse_args()

    cf_swath = Config(args.swath_pcf,delim=' = ')
    cf_tile  = Config(args.tile_pcf,delim=' = ')

    # get lat/lon for tile 
    nlon,nlat = int(cf_tile('nlon')), int(cf_tile('nlat'))
    H,V = int(args.tile[1:3]),int(args.tile[4:])

    tlats = 60.0 - (V*6 + np.arange(nlat)*0.01) - 0.005
    tlons = -180.0 + (H*6 + np.arange(nlon)*0.01) + 0.005
    # Make sure longitudes in [-180,180]
    # ----------------------------------
    if tlons.max()>180.:
        tlons[tlons>180] = tlons[tlons>180] - 360.



    tlonmin,tlonmax = tlons.min(),tlons.max()
    tlatmin,tlatmax = tlats.min(),tlats.max()
    grid_x,grid_y = np.meshgrid(tlons,tlats)  #nlat,nlon

    # loop through start to end time
    sdate = isoparser(args.start_isotime)
    edate = isoparser(args.end_isotime)
    DT    = timedelta(days=args.dt_days)

    transform=cartopy.crs.PlateCarree()
    cmap = mpl.cm.jet

    while sdate < edate:
      
        year = sdate.strftime('%Y')
        month = sdate.strftime('%m')
        day   = sdate.strftime('%d')
        nymd  = year+month+day

        for hh in range(24):
            # Get Filenames 
            # -------------
            pdate = sdate + timedelta(hours=hh)

            # tile files
            tileDir = cf_tile('inDir').replace('%ORBITNAME',args.orbit.upper()).replace('%year',year).replace('%month',month).replace('%day',day)
            tileFile = tileDir + '/' + cf_tile('inFile').replace('%tile',args.tile.lower()).replace('%orbitname',args.orbit.lower()).replace('%col',args.col).replace('%nymd',nymd).replace('%hour',"{:02d}".format(hh))

            tilelist = sorted(glob(tileFile))

            if len(tilelist) > 0:
                swathDir = cf_swath('inDir').replace('%ORBITNAME',args.orbit.upper()).replace('%year',year).replace('%month',month).replace('%day',day)    
                swathFile = swathDir + '/' + cf_swath('inFile').replace('%orbitname',args.orbit.lower()).replace('%col',args.col).replace('%nymd',nymd).replace('%hour',"{:02d}".format(hh))

                swath = xr.open_dataset(swathFile)
                slons = swath.longitude.values
                slats = swath.latitude.values
                nalong = swath.dims['nalong']
                ncross = swath.dims['ncross']

                for tFile in tilelist:
                    tile = xr.open_dataset(tFile)

                    for ialong in range(nalong):              
                        tdata = tile.sel(along=ialong).isel(date=0).variables[args.var].values 
                        if args.var != 'time':
                            sdata = swath.sel(nalong=ialong).variables[args.var].values
                            label = args.var
                            bounds = np.linspace(sdata.min(),sdata.max(),100)
                        else:
                            sdata = swath.sel(nalong=ialong).variables['time_ss'].values
                            # copy time across swath
                            sdata.shape = sdata.shape + (1,)
                            sdata = np.repeat(sdata,ncross,axis=1)

                            dt = sdata - np.datetime64(pdate)
                            sdata = dt.astype(np.dtype('timedelta64[s]'))
                            sdata = sdata.astype(int)

                            dt = tdata - np.datetime64(pdate)
                            tdata = dt.astype(np.dtype('timedelta64[s]'))
                            tdata = np.ma.array(tdata.astype(int),mask=np.isnan(tdata))

                            label = 'seconds since {}'.format(pdate.isoformat())
                            bounds = np.linspace(-100,60*60,100)
                            sys.exit()

                        x_points = swath.sel(nalong=ialong).longitude.values
                        y_points = swath.sel(nalong=ialong).latitude.values
                        norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

                        # plot swath on sphere
                        # ------------------
                        projection=cartopy.crs.Orthographic(central_longitude=tlons.mean(),central_latitude=tlats.mean())
                        axes_class = (GeoAxes,
                                      dict(map_projection=projection))
                        fig = plt.figure(figsize=[8,10])
                        axgr = AxesGrid(fig,  (0.1, 0.5, 0.8, 0.4) , axes_class=axes_class,
                                nrows_ncols=(1, 1),
                                axes_pad=0.6,
                                cbar_location='right',
                                cbar_mode='single',
                                cbar_pad=0.2,
                                cbar_size='3%',
                                label_mode='')  # note the empty label_mode

                        ax = axgr[0]
                        ax.pcolormesh(x_points,y_points,sdata,transform=transform,norm=norm,cmap=cmap)
                        ax.set_title('swath path')
                        ax.coastlines()
                        # draw tile box
                        ax.plot([tlonmin,tlonmin],[tlatmin,tlatmax],'k-',transform=transform) # left
                        ax.plot([tlonmin,tlonmax],[tlatmin,tlatmin],'k-',transform=transform) # bottom
                        ax.plot([tlonmin,tlonmax],[tlatmax,tlatmax],'k-',transform=transform) # top
                        ax.plot([tlonmax,tlonmax],[tlatmin,tlatmax],'k-',transform=transform) # right
                        ax.set_global()

                        # colorbar
                        cb = axgr.cbar_axes[0].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label=label,ticks=bounds[::10])
#                            cb.ax.set_yticklabels(ticks.astype(str))

                        # zoom in on part that overlaps tile
                        # ------------------------------------
                        projection=cartopy.crs.Orthographic(central_longitude=tlons.mean(),central_latitude=tlats.mean())
                        axes_class = (GeoAxes,
                                      dict(map_projection=projection))
                        axgr = AxesGrid(fig,  (0.01, 0.01, 0.4, 0.4) , axes_class=axes_class,
                                nrows_ncols=(1, 1),
                                axes_pad=0.6,
                                cbar_location='right',
                                cbar_mode='single',
                                cbar_pad=0.2,
                                cbar_size='3%',
                                label_mode='')  # note the empty label_mode

                        ax = axgr[0]
                        ax.pcolormesh(x_points,y_points,sdata,transform=transform,norm=norm,cmap=cmap)
                        ax.set_title('swath')
                        ax.coastlines()
                        # draw tile box
                        ax.plot([tlonmin,tlonmin],[tlatmin,tlatmax],'k-',transform=transform) # left
                        ax.plot([tlonmin,tlonmax],[tlatmin,tlatmin],'k-',transform=transform) # bottom
                        ax.plot([tlonmin,tlonmax],[tlatmax,tlatmax],'k-',transform=transform) # top
                        ax.plot([tlonmax,tlonmax],[tlatmin,tlatmax],'k-',transform=transform) # right

                        # colorbar
                        cb = axgr.cbar_axes[0].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label=label,ticks=bounds[::10])
#                            cb.ax.set_yticklabels(ticks.astype(str))

                        # plot interpolation to tile
                        # ------------------------------------
                        projection=cartopy.crs.Orthographic(central_longitude=tlons.mean(),central_latitude=tlats.mean())
                        axes_class = (GeoAxes,
                                      dict(map_projection=projection))
                        axgr = AxesGrid(fig,  (0.5, 0.01, 0.4, 0.4) , axes_class=axes_class,
                                nrows_ncols=(1, 1),
                                axes_pad=0.6,
                                cbar_location='right',
                                cbar_mode='single',
                                cbar_pad=0.2,
                                cbar_size='3%',
                                label_mode='')  # note the empty label_mode

                        ax = axgr[0]
                        ax.pcolormesh(tlons,tlats,tdata,transform=transform,norm=norm,cmap=cmap)
                        ax.set_title('interpolation to tile')
                        ax.coastlines()
                        # draw tile box
                        ax.plot([tlonmin,tlonmin],[tlatmin,tlatmax],'k-',transform=transform) # left
                        ax.plot([tlonmin,tlonmax],[tlatmin,tlatmin],'k-',transform=transform) # bottom
                        ax.plot([tlonmin,tlonmax],[tlatmax,tlatmax],'k-',transform=transform) # top
                        ax.plot([tlonmax,tlonmax],[tlatmin,tlatmax],'k-',transform=transform) # right

                        # colorbar
                        cb = axgr.cbar_axes[0].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label=label,ticks=bounds[::10])
#                            cb.ax.set_yticklabels(ticks.astype(str))

                        plt.suptitle("{} {} ialong={}".format(pdate.isoformat(),args.var,ialong))
                        plt.savefig("{}_{}_{}.png".format(pdate.strftime("%Y%m%d_%H"),args.var,ialong))
                    

                        


#                # loop through variables in swath, and interpolate the 3-D (swath level) vars 
#                # to the tile
#                for var in list(swath.variables):
#                    if (swath.variables[var].ndim == 3) & (var not in ['latitude','longitude']):
                        
#                        for ialong in swath.nalong:
                            
#                            indata = swath.sel(nalong=ialong).variables[var].values
#                            x_points = swath.sel(nalong=ialong).longitude.values
#                            y_points = swath.sel(nalong=ialong).latitude.values
#                            sub    = check[:,ialong,:]
#
#                            interp = griddata((x_points[sub], y_points[sub]), indata[sub],(grid_x,grid_y))
#                            interp = np.ma.array(interp,mask=np.isnan(interp))
#                            if var in tile.variables:
#                                tile.__dict__[var].data.append(interp)
#                            else:
#                                tile.variables.append(var)
#                                tile.__dict__[var] = HOLDER(var)
#                                tile.__dict__[var].long_name = swath.variables[var].attrs['long_name']
#                                tile.__dict__[var].units = swath.variables[var].attrs['units']
#                                tile.__dict__[var].km = 0
#                                tile.__dict__[var].data = [interp]
                            
#                    # special handling for time variable
#                    elif var == 'time_ss':
#                        for ialong in swath.nalong:

#                            indata = swath.sel(nalong=ialong).variables[var].values
#                            # copy time along the cross track
#                            indata.shape = indata.shape + (1,)
#                            indata = np.repeat(indata,ncross,axis=1)
#                            dt = indata - np.datetime64(sdate) 
#                            x_points = swath.sel(nalong=ialong).longitude.values
#                            y_points = swath.sel(nalong=ialong).latitude.values
#                            sub    = check[:,ialong,:]

#                            interp = griddata((x_points[sub], y_points[sub]), dt[sub],(grid_x,grid_y))
#                            # griddata returns float, convert back to timedelta
#                            interp = np.ma.array(interp,mask=np.isnan(interp))
#                            dt     = interp.astype(np.dtype('timedelta64[ns]'))
#                            tyme  = np.datetime64(sdate) + dt

#                            if 'tyme' in tile.__dict__:
#                                tile.tyme.append(tyme)
#                                tile.dt.append((interp*1e-9).astype(int)) #seconds
#                            else:
#                                tile.dt = [(interp*1e-9).astype(int)]  #seconds
#                                tile.tyme = [tyme]
                    

        


    

        sdate += DT
