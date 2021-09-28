"""
Interpolates MCD12C1 data to lidar trajectory
"""

import os, sys, subprocess
from   datetime               import date, datetime, timedelta
from   dateutil.parser        import parse as isoparser
from   dateutil.relativedelta import relativedelta
from   pyhdf.SD               import SD, HDF4Error
import numpy                  as     np
from   glob                   import glob
from   scipy.interpolate      import RegularGridInterpolator
from   netCDF4                import Dataset

nkernels = 3
SDS = { 'Majority_Land_Cover_Type_1'                 : 'IGBPlc'}   

LANDCOVER = {"water"                              : 0,
             "evergreen needleleaf forest"        : 1,
             "evergreen broadleaf forest"         : 2,
             "deciduous needleleaf forest"        : 3,
             "deciduous broadleaf forest"         : 4,
             "mixed forests"                      : 5,
             "closed shrubland"                   : 6,
             "open shrublands"                    : 7,
             "woody savannas"                     : 8,
             "savannas"                           : 9,
             "grasslands"                         : 10,
             "permanent wetlands"                 : 11,
             "croplands"                          : 12,
             "urban and built-up"                 : 13,
             "cropland/natural vegetation mosaic" : 14,
             "snow and ice"                       : 15,
             "barren or sparsely vegetated"       : 16}

#from Maignan et al. RSE, 2009 Table 1
COEF = {0 : -999,
        1 : 4.98,
        2 : 7.36,
        3 : 5.98,
        4 : 6.87,
        5 : 4.95,
        6 : 4.71,
        7 : 5.99,
        8 : 6.16,
        9 : 6.66,
        10: 6.11,
        11: -999,
        12: 7.58,
        13: -999,
        14: 6.86,
        15: 7.07,
        16: 7.29}

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

class LC(object):
    def __init__(self,nscan,npixel):
        for sds in SDS:
            self.__dict__[SDS[sds]] = np.zeros((nscan,npixel))



class MCD12C(object):
    def __init__(self,datadir,verbose=False):
        self.verbose = verbose
        self.dlon = 0.05
        self.dlat = 0.05

        self.nEW = 360./self.dlon
        self.nNS = 180./self.dlat

        self.lon = np.linspace(-180 + 0.5*self.dlon,180 - 0.5*self.dlon,self.nEW)
        self.lat = np.linspace(-90 + 0.5*self.dlat,90 - 0.5*self.dlat,self.nNS)

        self.fileName = datadir+'/MCD12C1.A%year001.051.hdf'
        self.SDS      = SDS

    def readFile(self,tyme):

        hfile = SD(self.fileName.replace('%year',str(tyme.year)))

        for sds in SDS:
            v = hfile.select(sds).get()
            a = hfile.select(sds).attributes()

            v = v.astype('float')
            fill_value = a['_FillValue']
            v[np.abs(v - fill_value)/fill_value < 0.01] = -999.
            self.__dict__[SDS[sds]] = np.flipud(v)

            
    #----
    def writenc(self,nctrj,outFile,verbose=False):
        """
        Write a NetCDF file with sampled MCD12C1 Land Cover Type
        and BPDF coefficient on lidar trajectory
        """

        # Dimensions
        # ------------------
        ntime = 1

        # Open NC file
        # ------------
        nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = "MCD12C1 Trajectory Sampler"
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from MCD12C1 v051 collections by mcd12c_sampler.py'
        nc.references = 'n/a'
        nc.comment = 'This file contains IGBP Land Cover and BPDF coefficient from Table 1 in Maignan et al. 2009 sampled on a satellite trajectory'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'  

        for lc in LANDCOVER:  
            nc.__dict__['TYPE'+str(LANDCOVER[lc])] = lc

        # Create dimensions
        # -----------------
        t = nc.createDimension('time',ntime) # one 5-minute granule per file for now
        x = nc.createDimension('ccd_pixels',self.npixel)  
        y = nc.createDimension('number_of_scans',self.nscan)

        # Save lon/lat
        # --------------------------
        _copyVar(nctrj,nc,u'ccd_pixels',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'number_of_scans',dtype='f4',zlib=False,verbose=verbose)            
        _copyVar(nctrj,nc,u'longitude',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'latitude',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'time', dtype='f4',zlib=False,verbose=verbose)

 
        # Write land cover to file
        #---------------------------------------
        dim = ('time','number_of_scans','ccd_pixels')
        for sds in SDS:
          this = nc.createVariable(SDS[sds],'i2',dim)  

          this.long_name = sds
          this.missing_value = -999
          this.units = 'none'  
          this[:] = self.lc.__dict__[SDS[sds]]


        this = nc.createVariable('BPDFcoef','f4',dim)
        this.long_name = 'BPDF coefficient from Maignan et al.'
        this.units = 'none'
        this[:] = self.BPDFcoef
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
        # Round dates to day
        dtyme = np.array([isoparser(str(d.date())) for d in tyme])
        utyme = np.sort(np.unique(dtyme))

        tyme.shape = (nscan,1)
        tyme       = np.repeat(tyme,npixel,axis=1)
        dtyme.shape = (nscan,1)
        dtyme       = np.repeat(dtyme,npixel,axis=1)

        self.lc = LC(nscan,npixel)
        
        self.tyme = tyme
        self.trjLon = trjLon
        self.trjLat = trjLat
        self.npixel = npixel
        self.nscan  = nscan
        for ut in utyme:
            if Verbose:
                print 'Working on '+ str(ut.date())
            self.readFile(ut)

            Ityme = dtyme == ut

            lat = trjLat[Ityme]
            lon = trjLon[Ityme]
            pts = []
            for LAT,LON in zip(lat,lon): pts.append([LAT,LON])

            for sds in SDS:
                interpFunc = RegularGridInterpolator((self.lat, self.lon), self.__dict__[SDS[sds]],
                                method='nearest',bounds_error=False,fill_value=None)
                self.lc.__dict__[SDS[sds]][Ityme] = interpFunc(pts)      


        self.BPDFcoef = np.array([COEF[lc] for lc in self.lc.IGBPlc.ravel()]).reshape(nscan,npixel)
        self.writenc(nctrj,outFile,verbose=Verbose) 
        nctrj.close()     

