"""
Interpolates OMI data to leo trajectory
"""

import os, sys
from   datetime               import date, datetime, timedelta
from   dateutil.parser        import parse as isoparser
from   dateutil.relativedelta import relativedelta
import numpy                  as     np
from   glob                   import glob
import  mpl_toolkits.basemap
from   netCDF4                import Dataset


SDS = { 'SRFLER388'           : ('surface lambertian equivalent reflectance at 388 nm','None'),
        'SRFLER354'           : ('surface lambertian equivalent reflectance at 354 nm','None'),
      }

MISSING = -999.
#----
def _copyVar(ncIn,ncOut,name,dtype='f4',zlib=False,verbose=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    if verbose:
        print('copy variable ',name,x.dimensions)

    try:
        fill_value = x.fill_value
    except:
        fill_value = None      
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib,fill_value=fill_value)
    if hasattr(x,'long_name'): y.long_name = x.long_name
    if hasattr(x,'units'): y.units = x.units
    try:
        y.missing_value = x.missing_value
    except:
        pass
    rank = len(x.shape)

    if rank == 1:
        y[:] = x[:]
    elif rank == 2:
        y[:,:] = x[:,:]
    elif rank == 3:
        y[:,:,:] = x[:,:,:]
    else:
        raise ValueError("invalid rank of <%s>: %d"%(name,rank))

class LER(object):
    def __init__(self):
        pass

class OMI_LER(object):
    def __init__(self,inDir,verbose=False):
        self.verbose = verbose

        self.readFile(inDir)
        self.readGrid(inDir)


    def readFile(self,inDir):
        flist = sorted(glob('{}/OMI_V158_SurfLER_Data_*.nc4'.format(inDir)))
        for inFile in flist:
            nc = Dataset(inFile)

            for sds in SDS:

                if not hasattr(self,'fill_value'):
                    var = nc.variables[sds]
                    self.fill_value = var.missing_value

                # time, lev, lat, lon
                data = nc.variables[sds][:,0,:,:]
                

                if hasattr(self,sds):
                    self.__dict__[sds] = np.ma.append(self.__dict__[sds],data,axis=0)
                else:
                    self.__dict__[sds] = data

    def readGrid(self,inDir):
        flist = sorted(glob('{}/OMI_V158_SurfLER_Data_*.nc4'.format(inDir)))
        gridFile = flist[0]
        nc = Dataset(gridFile)
        self.nEW = len(nc.dimensions['lon'])
        self.nNS = len(nc.dimensions['lat'])
        self.lon = nc.variables['lon'][:]
        self.lat = nc.variables['lat'][:]
        lonmax = self.lon.max()
        lonmin = self.lon.min()
        latmin = self.lat.min()
        latmax = self.lat.max()

        self.dlon = (lonmax-lonmin)/self.nEW
        self.dlat = (latmax-latmin)/self.nNS


            
    #----
    def writenc(self,nctrj,outFile,verbose=False,zlib=True):
        """
        Write a NetCDF file with sampled MCD43C1 kernel weights on lidar trajectory
        """

        # Dimensions
        # ------------------
        ntime = 1
        nlev  = 1

        # Open NC file
        # ------------
        nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = "OMI LER PACE Sampler"
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from OMI LER dataset by omi_ler_sampler.py'
        nc.references = 'n/a'
        nc.comment = 'This file contains OMI LER sampled on a PACE granule'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'    

        # Create dimensions
        # -----------------
        t = nc.createDimension('time',ntime)
        l = nc.createDimension('lev',nlev)
        x = nc.createDimension('ccd_pixels',self.npixel)
        y = nc.createDimension('number_of_scans',self.nscan)


        # lev
        dim = ('lev',)
        this = nc.createVariable('lev','f4',dim,zlib=zlib)  
        this.long_name = 'vertical level'
        this.units = 'layer'
        this.positive = 'down'
        this[:] = 1

        # Save lon/lat
        # --------------------------
        _copyVar(nctrj,nc,'ccd_pixels',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'number_of_scans',dtype='f4',zlib=zlib,verbose=verbose)            
        _copyVar(nctrj,nc,'longitude',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'latitude',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'time', dtype='f4',zlib=zlib,verbose=verbose)


        # Loop over variables writing each dataset
        #---------------------------------------
        dim = ('time','lev','number_of_scans','ccd_pixels')
        for sds in SDS:
            data = self.ler.__dict__[sds]
            this = nc.createVariable(sds,'f4',dim,zlib=zlib,fill_value=MISSING)  

            this.long_name = SDS[sds][0]
            this.missing_value = MISSING
            this.unit = SDS[sds][1]
            this[0,0,:,:] = data


        nc.close()          

    def sample(self,inFile,outFile,Verbose=False):

        # Open lidar sampled file
        nctrj = Dataset(inFile)

        # Read in locations
        trjLat = nctrj.variables['latitude'][:]
        trjLon = nctrj.variables['longitude'][:]

        nscan,npixel = trjLon.shape 

        # Read in tymes
        scanStart = isoparser(nctrj.time_coverage_start)
        midTime   = nctrj.variables['time'][:]
        tyme       = np.array([scanStart + timedelta(seconds=int(t)) for t in midTime])   
        # Round dates to month
        dtyme = np.array([d.month for d in tyme])
        utyme = np.sort(np.unique(dtyme))

        tyme.shape = (nscan,1)
        tyme       = np.repeat(tyme,npixel,axis=1)
        dtyme.shape = (nscan,1)
        dtyme       = np.repeat(dtyme,npixel,axis=1)

        self.ler = LER()
        self.tyme = tyme
        self.trjLon = trjLon
        self.trjLat = trjLat
        self.nscan  = nscan
        self.npixel = npixel


        for sds in SDS:
            data = self.__dict__[sds]
            self.ler.__dict__[sds] = np.zeros((nscan,npixel))

            print('sampling sds',sds)
            templin      = np.ma.zeros((nscan,npixel))
            templin.mask = np.ones((nscan,npixel)).astype(bool)
            tempnn       = np.ma.zeros((nscan,npixel))
            tempnn.mask  = np.ones((nscan,npixel)).astype(bool)
            for ut in utyme:  
                
                if Verbose:
                    print('Working on '+ str(ut.date()))

                Ityme = dtyme == ut

                lat = trjLat[Ityme]
                lon = trjLon[Ityme]

                datach = data[ut-1,:,:]
                
                if trjLon[0,0] > trjLon[-1,-1]:
                    #center at dateline, lon goes from 0 to 360
                    i = lon < 0
                    lon[i] = lon[i] + 360.0

                    nlat,nlon = datach.shape
                    lon1 = self.lon[0:int(nlon*0.5)]
                    lon2 = self.lon[int(nlon*0.5):]
                    datalon = np.append(lon2,lon1+360.)

                    data1  = datach[:,0:int(nlon*0.5)]
                    data2  = datach[:,int(nlon*0.5):]
                    datach = np.ma.hstack((data2, data1))                    

                else:
                    datalon = self.lon

                """
                If datain is a masked array and order=1 (bilinear interpolation) is
                used, elements of dataout will be masked if any of the four surrounding
                points in datain are masked.  To avoid this, do the interpolation in two
                passes, first with order=1 (producing dataout1), then with order=0
                (producing dataout2).  Then replace all the masked values in dataout1
                with the corresponding elements in dataout2 (using numpy.where).
                This effectively uses nearest neighbor interpolation if any of the
                four surrounding points in datain are masked, and bilinear interpolation
                otherwise.
                """                    
            
                tempnn[Ityme] = mpl_toolkits.basemap.interp(datach, datalon, self.lat,lon,lat,order=0)
                templin[Ityme] = mpl_toolkits.basemap.interp(datach, datalon, self.lat,lon,lat,order=1)

                if np.sum(templin.mask) > 0:
                    templin[templin.mask] = tempnn[templin.mask]

            self.ler.__dict__[sds][:,:] = templin


        self.writenc(nctrj,outFile,verbose=Verbose) 
        nctrj.close()     

