#!/usr/bin/env python3

"""
    Interpolates GEOS profile variables to standard PACE trace gas atmosphere
    Also calculate a few extra variables (e.g. water vapor VMR and O3)
"""

import os
import sys
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d
from glob import glob
from optparse        import OptionParser
from dateutil.parser import parse         as isoparser
from datetime import datetime, timedelta

SDS_MET = ['PS','U10M','V10M','TQV']
SDS_AER = ['RH']
SDS_CHM = ['O3']
class REGRID(object):
    def __init__(self,t1,t2,options):
        dt = timedelta(minutes=25)

        # read alpha Table
        # get new vertical grid
        self.read_alpha(options.alphaTable)

        # get filenames
        while t1 <= t2:
            self.aerFile = None
            self.metFile = None
            self.chmFile = None
            dir = t1.strftime('/Y%Y/M%m/D%d/')
            fname = t1.strftime('*%Y%m%d_%H%M00.nc4')
            # get filenames
            col = 'aer_Nv'
            self.aerFile = glob(options.lbRoot+dir+'*'+col+fname)[0]
            col = 'met_Nv'
            self.metFile = glob(options.lbRoot+dir+'*'+col+fname)[0]
            col = 'chm_Nv'
            self.chmFile = glob(options.lbRoot+dir+'*'+col+fname)[0]
            # get model ze
            self.get_ze()
          
            # create new files
            # copy over all variables regridded to new vertical grid
            if 'aer_Nv' in options.cols:
                col = 'aer_Nv'
                fname = self.aerFile
                self.newfile(col,fname)
            if 'met_Nv' in options.cols:
                col = 'met_Nv'
                fname = self.metFile
                self.newfile(col,fname)
                if 'aer_Nv' not in options.cols:
                    # need to get RH
                    self.newfile('aer_Nv',self.aerFile,SDS=['RH'],readonly=True)                    
            if 'chm_Nv' in options.cols:
                col = 'chm_Nv'
                fname = self.chmFile
                self.newfile(col,fname)
                self.get_to3()

            #convert RH to water vapor VMR and calculate TQV
            if ('aer_Nv' in options.cols) or ('met_Nv' in options.cols):
                self.get_wv_vmr()

            # write wv_vmr to aer file
            if 'aer_Nv' in options.cols:
                outFile = self.aerFile.replace('pace-g5nr','pace-g5nr-std')
                nc = Dataset(outFile,'r+')
                var = nc.createVariable('WV_VMR','f4',('time','lev','number_of_scans', 'ccd_pixels',),zlib=True,fill_value=1e15)
                var.long_name = 'water vapor mixing ratio'
                var.standard_name = 'WV_VMR'
                var.missing_value = np.float32(1e15)
                var.units = 'ppm'
                var[:] = self.wv_vmr

            # replace TQV with new values
            if 'met_Nv' in options.cols:
                outFile = self.metFile.replace('pace-g5nr','pace-g5nr-std')
                nc = Dataset(outFile,'r+')
                var = nc.variables['TQV']
                var[:] = self.tqv
                nc.close()

            t1 += dt
#---
    def get_to3(self):
        """
        Calculate total column ozone
        """
        DU = 2.69e20  #1DU = # molcules/m2
        AIR_MW = 28.964*1e-3   # kg/mole
        O3_MW  = 48.0*1e-3     # kg/mole
        Na = 6.022e23 # Avogadro's number molecules/mole

        # convert mass mixing ratio to molecules/m3
        airdens = self.AIRDENS    # kg/m3
        airdens.shape = (1,self.km,1,1)
        o3_conc = self.O3*airdens  # kg/m3
        o3_conc = o3_conc*Na/O3_MW # molecules/m3

        # sum up column  molecules/m2
        ze = self.ZE
        dz = ze[:-1]-ze[1:]
        dz.shape = (1,self.km,1,1)

        o3_col = np.sum(o3_conc*dz,axis=1)

        self.TO3 = o3_col/DU

        # get O3 VMR
        rhom = airdens*Na/AIR_MW   # molecules/m3
        o3_vmr = o3_conc/rhom

        # write to file
        outFile = self.chmFile.replace('pace-g5nr','pace-g5nr-std')
        nc = Dataset(outFile,'r+')
        var = nc.createVariable('O3_VMR','f4',('time','lev','number_of_scans', 'ccd_pixels',),zlib=True,fill_value=1e15)
        var.long_name = 'ozone mixing ratio'
        var.standard_name = 'O3_VMR'
        var.missing_value = np.float32(1e15)
        var.units = 'mol mol-1'
        var[:] = o3_vmr

        var = nc.createVariable('TO3','f4',('time','number_of_scans', 'ccd_pixels',),zlib=True,fill_value=1e15)
        var.long_name = 'ozone total column'
        var.standard_name = 'TO3'
        var.missing_value = np.float32(1e15)
        var.units = 'DU'
        var[:] = self.TO3

        nc.close()
        
#---
    def get_wv_vmr(self):
        """
        Convert RH to water vapor vmr [ppm]
        """
        R = 8.3145  # gas constant [J mol-1 K-1]
        Na = 6.022e23 # Avogadro's number molecules/mole
        tm = self.T
        tm.shape = (1,self.km,1,1)
        AIR_MW = 28.964*1e-3   # kg/mole
        ze = self.ZE
        airdens = self.AIRDENS  # kg/m3
        dz = ze[:-1]-ze[1:]
        dz.shape = (1,self.km,1,1)
        rhom = airdens*Na/AIR_MW   # molecules/m3
        
        #figure out saturation vapor pressure
        T0 = 273.15  #K
        e0 = 611.    #Pa
        L  = 2.5e6   #J/kg latend heat of vaporization
        Rv = 462.    # J/kg K specific gas constant for water vapor
        MWV = 18.0    # molecular weight of water g/mole

        tt = (1/T0) - (1/tm)
        esat = e0*np.exp(L*tt/Rv)

        esat.shape = (1,self.km,1,1)
    
        # get water vapor partial pressure from RH
        e = self.RH*esat

        # convert partial pressure to VMR [ppm]
        pm = self.P
        pm.shape = (1,self.km,1,1)
        self.wv_vmr = e*1e6/pm

        # calculate total precipitable water
        h2ocol = np.sum(self.wv_vmr*1e-6*rhom*dz,axis=1)   # molecules/m2
        tqv    = h2ocol*MWV/Na           # g/m2
        self.tqv    = tqv*1e-3                # kg/m2


#---
    def get_ze(self):
        nc = Dataset(self.aerFile)
        AIRDENS = nc.variables['AIRDENS'][:]
        PS      = nc.variables['PS'][:]
        DELP    = nc.variables['DELP'][:]

        ptop    = 1.0 # Pa
        grav = 9.80616  # m/s2

        ntime,km,nalong,ncross = AIRDENS.shape
        ZE = np.zeros([ntime,km+1,nalong,ncross])
        for k in range(km-1,-1,-1):
            ZE[:,k,:,:] = ZE[:,k+1,:] + DELP[:,k,:,:]/(AIRDENS[:,k,:,:]*grav)

        dz = ZE[:,:-1,:,:] - ZE[:,1:,:,:]
        self.zm = ZE[:,:-1,:,:] - 0.5*dz
        self.dz = dz
        self.nalong = nalong
        self.ncross = ncross

# ---
    def newfile(self,col,fname,SDS=None,readonly=False):
        outFile = fname.replace('pace-g5nr','pace-g5nr-std')

        nci = Dataset(fname)
        if not readonly:
            nc = Dataset(outFile,'w')
            # copy attributes
            for att in nci.ncattrs():
                nc.setncattr(att,nci.getncattr(att))

            # dimensions
            for dim in nci.dimensions:
                if dim == 'lev':
                    nc.createDimension('lev',self.km)
                else:
                    nc.createDimension(dim,len(nci.dimensions[dim]))

        if col == 'met_Nv':
            SDS = SDS_MET
#        elif col == 'aer_Nv':
#            SDS = SDS_AER
#        elif col == 'chm_Nv':
#            SDS = SDS_CHM
        elif SDS is None:
            SDS = list(nci.variables.keys())

        for sds in SDS:
            var = nci.variables[sds]

            if not readonly:
                varo = nc.createVariable(sds,var.dtype,var.dimensions,zlib=True)
            for att in var.ncattrs():
                if not readonly:
                    if att == 'missing_value':
                        varo.setncattr(att,np.float32(var.getncattr(att)))
                    else:
                        varo.setncattr(att,var.getncattr(att))

            if 'lev' not in var.dimensions:
                # copy
                if not readonly:
                    if sds == 'PS':
                        varo[:] = self.PE[-1]
                    else:
                        varo[:] = var[:]
            else:
                # regrid
                # lev dimensions is always second
                #data = var[:]
                if sds in ['DELP','AIRDENS']:
                    if not readonly:
                        outdata = self.__dict__[sds]
                        outdata.shape = (1,self.km,1,1)
                        varo[:] = np.tile(outdata,[1,1,self.nalong,self.ncross])
                        outdata.shape = (self.km,)
                elif sds == 'lev':
                    print('sds',sds)
                    if not readonly:
                        varo[:] = np.arange(self.km)+1
                else:
                    print('sds',sds)
                    data = var[:]
                    datanew = np.zeros([1,self.km,self.nalong,self.ncross])
                    for i in range(self.nalong):
                        for j in range(self.ncross):
                            f = interp1d(self.zm[0,:,i,j],data[0,:,i,j],fill_value="extrapolate")
                            datanew[0,:,i,j] = f(self.Z)

                    if col in ['aer_Nv','chm_Nv']:
                        I = datanew < 0
                        datanew[I] = 0.0

                    if not readonly:
                        varo[:] = datanew

                    if sds == 'RH':
                        self.RH = datanew
                    elif sds =='O3':
                        self.O3 = datanew        
            
        nci.close()
        if not readonly:
            nc.close()

    def read_alpha(self,fname):
        nc = Dataset(fname)
        SDS = 'PE','TE','ZE','P','T','Z','AIRDENS','DELP'
        for sds in SDS:
            self.__dict__[sds] = np.array(nc.variables[sds][:])
        nc.close()

        self.km = len(self.P)

    def get_sds(self,fname):
        nc = Dataset(fname)
        sdslist = list(nc.variables.keys())
        levlist = []
        for sds in sdslist:
            var = nc.variables[sds]
            if 'lev' in var.dimensions:
                levlist.append(sds)
        nc.close()

        self.levlist = levlist

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    cols = 'aer_Nv','met_Nv','chm_Nv'

    title   = 'GEOS-5 PACE Sampler'
    format  = 'NETCDF4_CLASSIC'
    alphaTable  = 'alphaTable_v0/alpha_CK_Thuillier_o3.nc4'

    # PACE default
    # -------------------
    calculon = '/nobackup/PACE'
    nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    if os.path.exists(nccs):
        lbRoot = nccs + '/LevelB'
        outdir = lbRoot
    elif os.path.exists(calculon):
        lbRoot = calculon + '/LevelB'
        outdir = lbRoot

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] iso_t1 [iso_t2]",
                          version='1.0.0' )

    parser.add_option("-p", "--path", dest="lbRoot", default=lbRoot,
              help='PACE LevelB path (default=%s)'%lbRoot)

    parser.add_option("-o", "--outdir", dest="outdir", default=outdir,
              help="Output NetCDF file in (default=%s)"\
                          %outdir )

    parser.add_option("-a", "--alphaTable", dest="alphaTable", default=alphaTable,
              help="alpha table (default=%s)"\
                          %alphaTable )

    parser.add_option("-t", "--title", dest="title", default=title,
              help="Output file title, typically the collection name (default=%s)"\
                          %title )

    parser.add_option("-c", "--col", dest="cols", default=cols,
              help="collections (default={})".\
                          format(cols) )

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")   

    (options, args) = parser.parse_args()
    options.cols = cols

    if len(args) == 2:
        iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 1:
        iso_t1 = args[0]
        iso_t2 = None
        t1     = isoparser(iso_t1)
        t2     = isoparser(iso_t1)
    else:
        parser.error("must have 1 or 2 arguments: iso_t1 [iso_t2]")


    if type(options.cols) is str:
        options.cols = [options.cols]

    regrid = REGRID(t1,t2,options)
