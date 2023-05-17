"""
Interpolates MCD43C data to aeronet locations
"""

import os, sys, subprocess
from   datetime               import date, datetime, timedelta
from   dateutil.parser        import parse
from   dateutil.relativedelta import relativedelta
from   pyhdf.SD               import SD, HDF4Error
import numpy                  as     np
from   glob                   import glob
from   scipy.interpolate      import RegularGridInterpolator


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
    def __init__(self,nobs,omp=False):
        if omp:      
            import pymp      
            for sds in SDS:
                self.__dict__[SDS[sds]] = pymp.shared.array((nobs,))                
        else:
            for sds in SDS:
                self.__dict__[SDS[sds]] = np.zeros(nobs)


class MCD43C(object):
    def __init__(self,inDir=None,coll='061'):
        self.dlon = 0.05
        self.dlat = 0.05

        self.nEW = 360./self.dlon
        self.nNS = 180./self.dlat

        self.lon = np.linspace(-180 + 0.5*self.dlon,180 - 0.5*self.dlon,self.nEW)
        self.lat = np.linspace(-90 + 0.5*self.dlat,90 - 0.5*self.dlat,self.nNS)
        
        if inDir is None:
            self.inDir  = '/nobackup/3/pcastell/MODIS/MCD43C1/061'
        else:
            self.inDir = inDir

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

            

    def getFileName(self,tyme):
        MM = str(tyme.month).zfill(2)
        DD = str(tyme.day).zfill(2)
        doy  = tyme - datetime(tyme.year,1,1) + timedelta(days=1)
        doy  = str(doy.days).zfill(3)

        inFileList = glob("{}/Y{}/M{}/*A{}{}*.hdf".format(self.inDir,tyme.year,MM,tyme.year,doy))

        if len(inFileList) != 1:
            print('Missing Files on {}'.format(tyme.strftime('%Y-%m-%d')))
            inFileList = [None] 
        
        return inFileList[0]

    def sample(self,giant,Verbose=False,omp=False):
        if omp:
            import pymp
        # Round dates to day
        tyme = np.array([parse(str(d.date())) for d in giant.tyme])
        utyme = np.sort(np.unique(tyme))


        giant.brdf = BRDF(giant.nobs,omp=omp)

        if omp:
            # use openmp type parallel processing
            with pymp.Parallel(10) as p:
                for ut in p.iterate(utyme):
                    if Verbose:
                        print('Working on '+ str(ut.date()))
                    inFile = self.getFileName(ut)
                    if inFile is not None:
                        self.readFile(inFile)

                        Ityme = tyme == ut

                        lat = giant.lat[Ityme]
                        lon = giant.lon[Ityme]
                        pts = []
                        for LAT,LON in zip(lat,lon): pts.append([LAT,LON])

                        for sds in SDS:
                            interpFunc = RegularGridInterpolator((self.lat, self.lon), self.__dict__[SDS[sds]],method='nearest')
                            giant.brdf.__dict__[SDS[sds]][Ityme] = interpFunc(pts)
        else:
            for ut in utyme:
                if Verbose:
                    print('Working on '+ str(ut.date()))
                inFile = self.getFileName(ut)
                if inFile is not None:
                    self.readFile(inFile)

                    Ityme = tyme == ut

                    lat = giant.lat[Ityme]
                    lon = giant.lon[Ityme]
                    pts = []
                    for LAT,LON in zip(lat,lon): pts.append([LAT,LON])

                    for sds in SDS:
                        interpFunc = RegularGridInterpolator((self.lat, self.lon), self.__dict__[SDS[sds]],method='nearest')
                        giant.brdf.__dict__[SDS[sds]][Ityme] = interpFunc(pts)            

