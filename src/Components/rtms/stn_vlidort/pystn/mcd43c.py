"""
Interpolates MCD43C data to stations
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
from   stn_vlidort_aux        import _copyVar

nkernels = 3
SDS = { 'BRDF_Albedo_Parameter1_vis'                 : 'Risovis',
        'BRDF_Albedo_Parameter1_nir'                 : 'Risonir',
        'BRDF_Albedo_Parameter2_vis'                 : 'Rvolvis',
        'BRDF_Albedo_Parameter2_nir'                 : 'Rvolnir',
        'BRDF_Albedo_Parameter3_vis'                 : 'Rgeovis',
        'BRDF_Albedo_Parameter3_nir'                 : 'Rgeonir',
        'BRDF_Albedo_Parameter1_Band1'               : 'Riso650',
        'BRDF_Albedo_Parameter1_Band3'               : 'Riso470',        
        'BRDF_Albedo_Parameter2_Band1'               : 'Rvol650',
        'BRDF_Albedo_Parameter2_Band3'               : 'Rvol470',
        'BRDF_Albedo_Parameter3_Band1'               : 'Rgeo650',
        'BRDF_Albedo_Parameter3_Band3'               : 'Rgeo470',                                
        'BRDF_Albedo_Parameter1_Band2'               : 'Riso850',
        'BRDF_Albedo_Parameter1_Band4'               : 'Riso550',
        'BRDF_Albedo_Parameter1_Band5'               : 'Riso1200',
        'BRDF_Albedo_Parameter1_Band6'               : 'Riso1600',
        'BRDF_Albedo_Parameter1_Band7'               : 'Riso2100',
        'BRDF_Albedo_Parameter2_Band2'               : 'Rvol850',
        'BRDF_Albedo_Parameter2_Band4'               : 'Rvol550',
        'BRDF_Albedo_Parameter2_Band5'               : 'Rvol1200',
        'BRDF_Albedo_Parameter2_Band6'               : 'Rvol1600',
        'BRDF_Albedo_Parameter2_Band7'               : 'Rvol2100',
        'BRDF_Albedo_Parameter3_Band2'               : 'Rgeo850',
        'BRDF_Albedo_Parameter3_Band4'               : 'Rgeo550',
        'BRDF_Albedo_Parameter3_Band5'               : 'Rgeo1200',
        'BRDF_Albedo_Parameter3_Band6'               : 'Rgeo1600',
        'BRDF_Albedo_Parameter3_Band7'               : 'Rgeo2100'}   


class BRDF(object):
    def __init__(self,nstations,nobs,omp=False):
        if omp:      
            import pymp      
            for sds in SDS:
                self.__dict__[SDS[sds]] = pymp.shared.array((nstations,nobs,))                
        else:
            for sds in SDS:
                self.__dict__[SDS[sds]] = np.zeros((nstations,nobs,))


class MCD43C(object):
    def __init__(self,datadir,verbose=False):
        self.verbose = verbose
        self.dlon = 0.05
        self.dlat = 0.05

        self.nEW = 360./self.dlon
        self.nNS = 180./self.dlat

        self.lon = np.linspace(-180 + 0.5*self.dlon,180 - 0.5*self.dlon,self.nEW)
        self.lat = np.linspace(-90 + 0.5*self.dlat,90 - 0.5*self.dlat,self.nNS)

        self.inDir  = datadir
        self.HTTP = 'http://e4ftl01.cr.usgs.gov//MODV6_Cmp_C/MOTA/MCD43C1.006/'
        self.command = 'wget --user patticastellanos --password Hernan11617! -r -nH -nd -np -R "index.html*" -R "*.xml" -P '
        self.SDS     = SDS

    def readFile(self,inFile):

        hfile = SD(inFile)

        for sds in SDS:
            v = hfile.select(sds).get()
            a = hfile.select(sds).attributes()

            v = v.astype('float')
            if a['scale_factor']!=1.0 or a['add_offset']!=0.0:
                v = v * float(a['scale_factor'])  + float(a['add_offset'])
                fill_value = float(a['_FillValue'])*float(a['scale_factor'])   + float(a['add_offset'])
                v[np.abs(v - fill_value)/fill_value < 0.01] = -999.
                self.__dict__[SDS[sds]] = np.flipud(v)

            

    def downloadFile(self,tyme):
        MM = str(tyme.month).zfill(2)
        DD = str(tyme.day).zfill(2)
        doy  = tyme - datetime(tyme.year,1,1) + timedelta(days=1)
        doy  = str(doy.days).zfill(3)

        inFileList = glob("{}/Y{}/M{}/*A{}{}*.hdf".format(self.inDir,tyme.year,MM,tyme.year,doy))

        if len(inFileList) != 1:
            Outdir = "{}/Y{}/M{}/".format(self.inDir,tyme.year,MM)
            dd = '{}.{}.{}'.format(tyme.year,str(tyme.month).zfill(2),str(tyme.day).zfill(2))
            print('Downloading '+dd)
            subprocess.call(self.command+Outdir+' '+self.HTTP+dd+'/',shell=True)        
            inFileList = glob("{}/Y{}/M{}/*A{}{}*.hdf".format(self.inDir,tyme.year,MM,tyme.year,doy))
            if len(inFileList) != 1:
                raise Exception('problem downloading '+ dd)
            
        
        return inFileList[0]

    #----
    def writenc(self,ncstn,outFile,verbose=False):
        """
        Write a NetCDF file with sampled MCD43C1 kernel weights on lidar trajectory
        """

        # Dimensions
        # ------------------
        ntime = len(self.tyme)
        nstations = len(self.stnLat)

        # Open NC file
        # ------------
        nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = "MCD43C1 Trajectory Sampler"
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from MCD43C1 v006 collections by mcd43c_sampler.py'
        nc.references = 'n/a'
        nc.comment = 'This file contains BRDF Kernels weights for the RTLS model for 8 MODIS bands sampled on a satellite trajectory'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'    
        nc.BAND1 = "620-670nm"
        nc.BAND2 = "841-875nm"
        nc.BAND3 = "459-479nm"
        nc.BAND4 = "545-565nm"
        nc.BAND5 = "1230-1250nm"
        nc.BAND6 = "1628-1652nm"
        nc.BAND7 = "2105-2155nm"

        # Create dimensions
        # -----------------
        t = nc.createDimension('time',ntime)
        s = nc.createDimension('ls',19)
        x = nc.createDimension('x',1)  
        y = nc.createDimension('y',1)
        station = nc.createDimension('station',nstations)

        # Save lon/lat
        # --------------------------
        _copyVar(ncstn,nc,'time', dtype='i4',zlib=False,verbose=verbose)
        _copyVar(ncstn,nc,'x',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(ncstn,nc,'y',dtype='f4',zlib=False,verbose=verbose)      
        _copyVar(ncstn,nc,'station',dtype='f4',zlib=False,verbose=verbose)            
        _copyVar(ncstn,nc,'stnLon',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(ncstn,nc,'stnLat',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(ncstn,nc,'isotime', dtype='S1',zlib=False,verbose=verbose)

        # Loop over Bands writing each dataset
        #---------------------------------------
        dim = ('station','time',)
        for sds in SDS:
          this = nc.createVariable(SDS[sds],'f4',dim)  

          this.long_name = SDS[sds][4:] + ' BRDF Kernel weight'
          this.missing_value = -999.0
          this.unit = 'none'  
          this[:] = self.brdf.__dict__[SDS[sds]]


        nc.close()          

    def sample(self,inFile,outFile,Verbose=False,omp=False):
        if omp:
            import pymp

        # Open lidar sampled file
        ncstn = Dataset(inFile)

        # Read in tymes
        iso = ncstn.variables['isotime'][:]
        tyme = []
        for isotime in iso:
            tyme.append(isoparser(''.join(isotime)))

        tyme = np.array(tyme)        

        # Read in locations
        stnLat = ncstn.variables['stnLat'][:]
        stnLon = ncstn.variables['stnLon'][:]
        nstations = len(stnLat)

        # Round dates to day
        dtyme = np.array([isoparser(str(d.date())) for d in tyme])
        utyme = np.sort(np.unique(dtyme))


        self.brdf = BRDF(nstations,len(tyme),omp=omp)
        self.tyme = tyme
        self.stnLon = stnLon
        self.stnLat = stnLat

        pts = []
        for LAT,LON in zip(stnLat,stnLon): pts.append([LAT,LON])
        if omp:
            # use openmp type parallel processing
            with pymp.Parallel(10) as p:
                for ut in p.iterate(utyme):
                    if Verbose:
                        print('Working on '+ str(ut.date()))
                    inFile = self.downloadFile(ut)
                    self.readFile(inFile)

                    Ityme = dtyme == ut
                    Ityme = np.arange(len(dtyme))[Ityme]

                    for sds in SDS:
                        interpFunc = RegularGridInterpolator((self.lat, self.lon), self.__dict__[SDS[sds]],
                                        method='nearest',bounds_error=False,fill_value=None)
                        for t in Ityme:
                            self.brdf.__dict__[SDS[sds]][:,t] = interpFunc(pts)
        else:
            for ut in utyme:
                if Verbose:
                    print('Working on '+ str(ut.date()))
                inFile = self.downloadFile(ut)
                self.readFile(inFile)

                Ityme = dtyme == ut
                Ityme = np.arange(len(dtyme))[Ityme]

                for sds in SDS:
                    interpFunc = RegularGridInterpolator((self.lat, self.lon), self.__dict__[SDS[sds]],
                                    method='nearest',bounds_error=False,fill_value=None)
                    for t in Ityme:
                        self.brdf.__dict__[SDS[sds]][:,t] = interpFunc(pts)      


        self.writenc(ncstn,outFile,verbose=Verbose) 
        ncstn.close()     

