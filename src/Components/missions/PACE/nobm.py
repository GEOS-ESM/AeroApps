"""
Interpolates NOBM data to leo trajectory
"""

import os, sys
from   datetime               import date, datetime, timedelta
from   dateutil.parser        import parse as isoparser
from   dateutil.relativedelta import relativedelta
import numpy                  as     np
from   glob                   import glob
import  mpl_toolkits.basemap
from   netCDF4                import Dataset


SDS = { 'lwn'                 : 'sun normalized water leaving radiance'}

#----
def _copyVar(ncIn,ncOut,name,dtype='f4',zlib=False,verbose=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    if verbose:
        print 'copy variable ',name,x.dimensions
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
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
        raise ValueError, "invalid rank of <%s>: %d"%(name,rank)

class SLEAVE(object):
    def __init__(self,nch,nscan,npixel):
        for sds in SDS:
            self.__dict__[sds] = np.zeros((nch,nscan,npixel))


class NOBM(object):
    def __init__(self,inFile,verbose=False):
        self.verbose = verbose

        self.readFile(inFile)


    def readFile(self,inFile):

        nc = Dataset(inFile)

        self.nEW = len(nc.dimensions['lon'])
        self.nNS = len(nc.dimensions['lat'])
        lonmax = nc.Easternmost_Longitude
        lonmin = nc.Westernmost_Longitude
        latmin = nc.Southernmost_Latitude
        latmax = nc.Northernmost_Latitude

        self.dlon = (lonmax-lonmin)/self.nEW
        self.dlat = (latmax-latmin)/self.nNS

        self.lon = np.linspace(lonmin + 0.5*self.dlon,lonmax - 0.5*self.dlon,self.nEW)
        self.lat = np.linspace(latmin + 0.5*self.dlat,latmax - 0.5*self.dlat,self.nNS)

        self.nch = len(nc.dimensions['wavelength'])
        self.channels = np.linspace(250,775,self.nch)

        for sds in SDS:
            self.__dict__[sds] = nc.variables[sds][:]


            
    #----
    def writenc(self,nctrj,outFile,verbose=False):
        """
        Write a NetCDF file with sampled MCD43C1 kernel weights on lidar trajectory
        """

        # Dimensions
        # ------------------
        ntime = 1

        # Open NC file
        # ------------
        nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = "MERRA-NOBM Trajectory Sampler"
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from MERRA-NOBM Water Leaving Radiance dataset by nobm_sampler.py'
        nc.references = 'n/a'
        nc.comment = 'This file contains hyperspectral water leaving radiance sampled on a satellite trajectory'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'    

        # Create dimensions
        # -----------------
        t = nc.createDimension('time',ntime)
        x = nc.createDimension('ccd_pixels',self.npixel)
        y = nc.createDimension('number_of_scans',self.nscan)
        w = nc.createDimension('wavelength',self.nch)

        # Save lon/lat
        # --------------------------
        _copyVar(nctrj,nc,u'ccd_pixels',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'number_of_scans',dtype='f4',zlib=False,verbose=verbose)            
        _copyVar(nctrj,nc,u'longitude',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'latitude',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'ev_mid_time', dtype='f4',zlib=False,verbose=verbose)

        # Wavelenghts
        dim = ('wavelength',)
        this = nc.createVariable('wavelength','f4',dim)  
        this.long_name = 'center wavelegnth'
        this.units = 'nm'
        this[:] = self.channels

        # Loop over Bands writing each dataset
        #---------------------------------------
        dim = ('time','wavelength','number_of_scans','ccd_pixels')
        for sds in SDS:
          this = nc.createVariable(sds,'f4',dim)  

          this.long_name = SDS[sds]
          this.missing_value = -999.0
          this.unit = 'none'  
          this[:] = self.sleave.__dict__[sds]


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
        midTime   = nctrj.variables['ev_mid_time'][:]
        tyme       = np.array([scanStart + timedelta(seconds=int(t)) for t in midTime])   
        # Round dates to month
        dtyme = np.array([d.month for d in tyme])
        utyme = np.sort(np.unique(dtyme))

        tyme.shape = (nscan,1)
        tyme       = np.repeat(tyme,npixel,axis=1)
        dtyme.shape = (nscan,1)
        dtyme       = np.repeat(dtyme,npixel,axis=1)

        self.sleave = SLEAVE(self.nch,nscan,npixel)
        self.tyme = tyme
        self.trjLon = trjLon
        self.trjLat = trjLat
        self.nscan  = nscan
        self.npixel = npixel


        
        for ich in range(self.nch):
            temp = np.zeros((nscan,npixel))
            for ut in utyme:  
                
                if Verbose:
                    print 'Working on '+ str(ut.date())

                data = self.lwn[ut,ich,:,:]
                nlat,nlon = data.shape
                data1 = data[:,0:nlon*0.5]
                data2 = data[:,nlon*0.5:]
                data  = np.hstack((data2, data1))
                data = np.flipud(data)

                Ityme = dtyme == ut

                lat = trjLat[Ityme]
                lon = trjLon[Ityme]

            
                temp[Ityme] = mpl_toolkits.basemap.interp(data, self.lon, self.lat,lon,lat,order=0)

            
            self.sleave.lwn[ich,:,:] = temp


        self.writenc(nctrj,outFile,verbose=Verbose) 
        nctrj.close()     

