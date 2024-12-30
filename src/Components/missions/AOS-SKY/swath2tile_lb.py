#!/usr/bin/env python3

"""
    Utility to sample G5NR model output in a SBG tile

    Adapted from ACCP/lidar_sampler.py
    Trying to modernize to use xarray

    P. Castellanos Feb 2024

"""

import os
import sys

import argparse
from   dateutil.parser import parse         as isoparser
from   datetime        import datetime, timedelta
import numpy           as np
from   MAPL            import eta
from   MAPL.config     import Config
from   MAPL.constants  import *
import xarray          as xr
from   glob            import glob
from scipy.interpolate import griddata
from netCDF4 import Dataset
import pandas as pd

class HOLDER(object):
    def __init__(self,name):
        self.name = name

class TILE(object):
    """
    class to hold dataset on tile
    """
    def __init__(self,args,date,lon,lat,nalong):
        self.tile = args.tile
        self.isoTime = args.isoTime
        self.date = date
        self.lon = lon
        self.lat = lat
        self.nalong = nalong
        self.variables = []



    def writeNC (self,doAkBk=False,zlib=True):
        """
        Write a NetCDF file with re-sampled variables on a satellite tile
        described by (lon,lat).
        """

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w',format='NETCDF4')

        # Set global attributes
        # ---------------------
        nc.title = 'GEOS sampled on tile'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from GEOS-5 OSSE collections by swath2tile_lb.py'
        nc.references = 'n/a'
        nc.comment = 'This file contains GEOS-5 related parameters along a satellite or aircraft track.'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        
     
        # Create dimensions
        # -----------------
        d = nc.createDimension('date',None)
        n = nc.createDimension('along',self.nalong)
        y = nc.createDimension('lat',len(self.lat))
        x = nc.createDimension('lon',len(self.lon))
        ls = nc.createDimension('ls',19)

        # keep this in for now.  not used.
        km = 0
        if km>0:
            nz = nc.createDimension('lev',km)
            if doAkBk:
                ne = nc.createDimension('ne',km+1)


        # Coordinate variables
        # --------------------
        date = nc.createVariable('date','i4',('date',),zlib=zlib)
        date.long_name = 'Date'
        t0 = self.date
        date.units = 'seconds since %s'%t0.isoformat(' ')
        date[:] = [0]

        along = nc.createVariable('along','f4',('along',),zlib=zlib)
        along.long_name = 'Along Track View'
        along.units = 'None'
        along[:] = np.arange(self.nalong)

        lat = nc.createVariable('lat','f4',('lat',),zlib=zlib)
        lat.long_name = 'Tile Latitude'
        lat.units = 'degrees_north'
        lat[:] = self.lat

        lon = nc.createVariable('lon','f4',('lon',),zlib=zlib)
        lon.long_name = 'Tile Longitude'
        lon.units = 'degrees_east'
        lon[:] = self.lon


        if km > 0: # pressure level not supported yet
            lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
            lev.long_name = 'Vertical Level'
            lev.units = levUnits.strip()
            lev.positive = 'down'
            lev.axis = 'z'
            lev[:] = levs[:]

            if doAkBk:
                ae, be = eta.getEdge(km) # Coefficients for Hybrid coordinates
                ak = nc.createVariable('ak','f4',('ne',),zlib=zlib)
                ak.long_name = 'Eta coordinate coefficient ak (p = ak + bk * ps)'
                ak.units = 'Pa'
                ak = ae[:]
                bk = nc.createVariable('bk','f4',('ne',),zlib=zlib)
                bk.long_name = 'Eta coordinate coefficient bk (p = ak + bk * ps)'
                bk.units = '1'
                bk = be[:]
        
        
        # Trajectory time coordinates
        # ----------------------
        time = nc.createVariable('time','i4',('date','along','lat','lon',),zlib=zlib)
        time.long_name = 'Time'
        time.missing_value = 999999
        t0 = self.date
        time.units = 'seconds since %s'%t0.isoformat(' ')
        time[0,:] = np.ma.array(self.dt)

        
        # Time in ISO format if so desired
        # ---------------------------------
        if self.isoTime:
            isotime = nc.createVariable('isotime','S1',('date','along','lat','lon','ls'),zlib=zlib)
            isotime.long_name = 'Time (ISO Format)'
            for n in range(self.nalong):
                for i in range(len(self.lat)):        
                    isotmp = np.ma.masked_all((len(self.lon),19),dtype='S1')
                    for j in range(len(self.lon)):
                        if not self.tyme[n].mask[i,j]:
                            isotmp[j,:] = list(pd.to_datetime(self.tyme[n][i,j]).isoformat()[0:19])
                            isotmp.mask[i,:] = False
                    isotime[0,n,i,:,:] = isotmp

        # Loop over varaibles and write 
        # --------------------------------------------------
        for var in self.variables:
            Var = self.__dict__[var]
            if Var.km == 0:
                dim = ('date','along','lat','lon')
            else:
                dim = ('date','along','lat','lon','lev')
            if Var.name != 'time':
                this = nc.createVariable(Var.name,'f4',dim,zlib=zlib)
                this.long_name = Var.long_name
                this.missing_value = np.float32(MAPL_UNDEF)
                this.units = Var.units

                
                this[:] = Var.data
                
        # Close the file
        # --------------
        nc.close()


        if args.verbose:
            print(" <> wrote %s file %s"%(args.format,args.outFile))
        
    #---




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

    parser.add_argument("start_isotime",
                        help="start isotime")

    parser.add_argument("end_isotime",
                        help="end isotime")

    parser.add_argument("--col",default=col,
              help="collection name - only works with polar07 right now (default=%s)"\
                        %col)

    parser.add_argument("--isoTime",action="store_true",
              help="write isoTime to out file (default=False)")

    parser.add_argument("-a", "--algorithm", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_argument("-d", "--dt_days", default=dt_days,type=int,
              help="Timesetp in hours for swath files (default=%s)"%dt_days )

    parser.add_argument("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

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


    while sdate < edate:
      
        year = sdate.strftime('%Y')
        month = sdate.strftime('%m')
        day   = sdate.strftime('%d')
        nymd  = year+month+day

        # Get Filenames 
        # -------------
        inDir = cf_swath('inDir').replace('%ORBITNAME',args.orbit.upper()).replace('%year',year).replace('%month',month).replace('%day',day)
        for hh in range(24):
            
            swathFile = inDir + '/' + cf_swath('inFile').replace('%orbitname',args.orbit.lower()).replace('%col',args.col).replace('%nymd',nymd).replace('%hour',"{:02d}".format(hh))

            swath = xr.open_dataset(swathFile)
            slons = swath.longitude.values
            slats = swath.latitude.values
            nalong = swath.dims['nalong']
            ncross = swath.dims['ncross']

            # check to see if swath intersetcts with the tile
            # --------------------------------------------------
            check = (slons <= tlonmax) & (slons >= tlonmin) & (slats <= tlatmax) & (slats >= tlatmin)
            if check.any():
            
                # initialize tile class
                tile = TILE(args,sdate,tlons,tlats,nalong)
                
                # output filename
                outDir = cf_tile('inDir').replace('%ORBITNAME',args.orbit.upper()).replace('%year',year).replace('%month',month).replace('%day',day)
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                tileFile = outDir + '/' + cf_tile('inFile').replace('%tile',args.tile.lower()).replace('%orbitname',args.orbit.lower()).replace('%col',args.col).replace('%nymd',nymd).replace('%hour',"{:02d}".format(hh))

                tile.outFile = tileFile

                # loop through variables in swath, and interpolate the 3-D (swath level) vars 
                # to the tile
                for var in list(swath.variables):
                    if (swath.variables[var].ndim == 3) & (var not in ['latitude','longitude']):
                        
                        for ialong in swath.nalong:
                            
                            indata = swath.sel(nalong=ialong).variables[var].values
                            x_points = swath.sel(nalong=ialong).longitude.values
                            y_points = swath.sel(nalong=ialong).latitude.values
                            sub    = check[:,ialong,:]

                            interp = griddata((x_points[sub], y_points[sub]), indata[sub],(grid_x,grid_y))
                            interp = np.ma.array(interp,mask=np.isnan(interp))
                            if var in tile.variables:
                                tile.__dict__[var].data.append(interp)
                            else:
                                tile.variables.append(var)
                                tile.__dict__[var] = HOLDER(var)
                                tile.__dict__[var].long_name = swath.variables[var].attrs['long_name']
                                tile.__dict__[var].units = swath.variables[var].attrs['units']
                                tile.__dict__[var].km = 0
                                tile.__dict__[var].data = [interp]
                            
                    # special handling for time variable
                    elif var == 'time_ss':
                        for ialong in swath.nalong:

                            indata = swath.sel(nalong=ialong).variables[var].values
                            # copy time along the cross track
                            indata.shape = indata.shape + (1,)
                            indata = np.repeat(indata,ncross,axis=1)
                            dt = indata - np.datetime64(sdate) 
                            x_points = swath.sel(nalong=ialong).longitude.values
                            y_points = swath.sel(nalong=ialong).latitude.values
                            sub    = check[:,ialong,:]

                            interp = griddata((x_points[sub], y_points[sub]), dt[sub],(grid_x,grid_y))
                            # griddata returns float, convert back to timedelta
                            interp = np.ma.array(interp,mask=np.isnan(interp))
                            dt     = interp.astype(np.dtype('timedelta64[ns]'))
                            tyme  = np.datetime64(sdate) + dt

                            if 'tyme' in tile.__dict__:
                                tile.tyme.append(tyme)
                                tile.dt.append((interp*1e-9).astype(int)) #seconds
                            else:
                                tile.dt = [(interp*1e-9).astype(int)]  #seconds
                                tile.tyme = [tyme]
                    

        


                # Write output file
                # -----------------
                tile.writeNC()
    

        sdate += DT
