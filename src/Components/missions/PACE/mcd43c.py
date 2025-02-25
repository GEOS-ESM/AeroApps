"""
Interpolates MCD43C data to lidar trajectory
"""

import os, sys, subprocess
from   datetime               import date, datetime, timedelta
from   dateutil.parser        import parse as isoparser
from   dateutil.relativedelta import relativedelta
from   pyhdf.SD               import SD, HDF4Error
import numpy                  as     np
from   glob                   import glob
from   scipy.interpolate      import RegularGridInterpolator
import  mpl_toolkits.basemap
from   netCDF4                import Dataset


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

#----
def _copyVar(ncIn,ncOut,name,dtype='f4',zlib=False,verbose=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    if verbose:
        print('copy variable ',name,x.dimensions)
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
        raise ValueError("invalid rank of <%s>: %d"%(name,rank))

class BRDF(object):
    def __init__(self,nscan,npixel):
        for sds in SDS:
            self.__dict__[SDS[sds]]      = np.ma.zeros((nscan,npixel))
            self.__dict__[SDS[sds]].mask = np.ones((nscan,npixel)).astype(bool)


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
                I = np.abs(v - fill_value)/fill_value < 0.01
                v = np.ma.array(v,mask=I)
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
    def writenc(self,nctrj,outFile,verbose=False,zlib=True):
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
        x = nc.createDimension('ccd_pixels',self.npixel)
        y = nc.createDimension('number_of_scans',self.nscan)

        # Save lon/lat
        # --------------------------
        _copyVar(nctrj,nc,'ccd_pixels',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'number_of_scans',dtype='f4',zlib=zlib,verbose=verbose)            
        _copyVar(nctrj,nc,'longitude',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'latitude',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'time', dtype='f4',zlib=zlib,verbose=verbose)

        # Loop over Bands writing each dataset
        #---------------------------------------
        dim = ('time','number_of_scans','ccd_pixels')
        for sds in SDS:
          this = nc.createVariable(SDS[sds],'f4',dim,zlib=zlib,fill_value=-999.0)  

          this.long_name = SDS[sds][4:] + ' BRDF Kernel weight'
          this.missing_value = -999.0
          this.unit = 'none'  
          this[:] = self.brdf.__dict__[SDS[sds]]


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

        self.brdf = BRDF(nscan,npixel)
        self.tyme = tyme
        self.trjLon = trjLon
        self.trjLat = trjLat
        self.nscan  = nscan
        self.npixel = npixel


        for ut in utyme:
            if Verbose:
                print('Working on '+ str(ut.date()))
            inFile = self.downloadFile(ut)
            self.readFile(inFile)

            Ityme = dtyme == ut

            for sds in SDS:
                lat = trjLat[Ityme]
                lon = trjLon[Ityme]

                data = self.__dict__[SDS[sds]]

                if trjLon[0,0] > trjLon[-1,-1]:
                    #center at dateline, lon goes from 0 to 360
                    i = lon < 0
                    lon[i] = lon[i] + 360.0

                    nlat,nlon = data.shape
                    lon1 = self.lon[0:int(nlon*0.5)]
                    lon2 = self.lon[int(nlon*0.5):]
                    datalon = np.append(lon2,lon1+360)

                    data1 = data[:,0:int(nlon*0.5)]
                    data2 = data[:,int(nlon*0.5):]
                    data  = np.ma.hstack((data2, data1))
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
            
                tempnn  = mpl_toolkits.basemap.interp(data, datalon, self.lat,lon,lat,order=0)
                templin = mpl_toolkits.basemap.interp(data, datalon, self.lat,lon,lat,order=1)

                if np.sum(templin.mask) > 0:
                    templin[templin.mask] = tempnn[templin.mask]

                self.brdf.__dict__[SDS[sds]][Ityme] = templin


        self.writenc(nctrj,outFile,verbose=Verbose) 
        nctrj.close()     

