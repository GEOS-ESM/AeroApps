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
from   copyvar                import _copyVar

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
        ntime = len(self.tyme)

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
        t = nc.createDimension('time',ntime)
        s = nc.createDimension('ls',19)
        x = nc.createDimension('x',1)  
        y = nc.createDimension('y',1)

        # Save lon/lat
        # --------------------------
        _copyVar(nctrj,nc,u'trjLon',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'trjLat',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=verbose)            

        # Write land cover to file
        #---------------------------------------
        dim = ('time',)
        for sds in SDS:
          this = nc.createVariable(SDS[sds],'i2',dim)  

          this.long_name = sds
          this.missing_value = -999
          this.units = 'none'  
          this[:] = self.__dict__['trj'+SDS[sds]]


        this = nc.createVariable('BPDFcoef','f4',dim)
        this.long_name = 'BPDF coefficient from Maignan et al.'
        this.units = 'none'
        this[:] = self.BPDFcoef
        nc.close()          

    def sample(self,inFile,outFile,Verbose=False):
        # Open lidar sampled file
        nctrj = Dataset(inFile)

        # Read in tymes
        iso = nctrj.variables['isotime'][:]
        tyme = []
        for isotime in iso:
            tyme.append(isoparser(''.join(isotime)))

        tyme = np.array(tyme)        

        # Read in locations
        trjLat = nctrj.variables['trjLat'][:]
        trjLon = nctrj.variables['trjLon'][:]

        # Round dates to year
        dtyme = np.array([datetime(d.year,1,1) for d in tyme])
        utyme = np.sort(np.unique(dtyme))

        for sds in SDS:
            self.__dict__['trj'+SDS[sds]] = np.zeros(len(tyme))
        
        self.tyme = tyme
        self.trjLon = trjLon
        self.trjLat = trjLat
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
                interpFunc = RegularGridInterpolator((self.lat, self.lon), self.__dict__[SDS[sds]],method='nearest')
                self.__dict__['trj'+SDS[sds]][Ityme] = interpFunc(pts)      


        self.BPDFcoef = np.array([COEF[lc] for lc in self.trjIGBPlc])
        self.writenc(nctrj,outFile,verbose=Verbose) 
        nctrj.close()     

