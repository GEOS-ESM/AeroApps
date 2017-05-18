#!/usr/bin/env python

"""
    Calculates polarized TOA radiance for a multiangle polarimeter on a lidar track.
    Model fields have already been sampled using trj_sampler

    Adapted from ext_sampler.py
    Patricia Castellanos, May, 2017

"""

import os
import MieObs_
from   netCDF4 import Dataset
from   mieobs  import  getAOPvector, getEdgeVars
import numpy   as np

from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from MAPL.constants import *
import LidarAngles_    
import VLIDORT_LIDAR_ 
from copyvar  import _copyVar

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

META    = ['DELP','RH','AIRDENS','LONGITUDE','LATITUDE','isotime']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = ['CLDTOT','PS']
SDS_INV = ['FRLAND']

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}


nMom     = 300

VZAdic = {'POLDER': np.array([3.66, 11., 18.33, 25.66, 33, 40.33, 47.66, 55.0])}

HGTdic = {'LEO': 705,
          'ISS': 400}


SurfaceFuncs = {'MODIS_BRDF': 'readSampledMODISBRDF'}
LandAlbedos  = 'MODIS_BRDF',


MISSING = -1.e+20

class POLAR_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,outFile,rcFile,albedoFile,albedoType,channel,VZA,hgtss,verbose=False):
        self.SDS_AER = SDS_AER
        self.SDS_MET = SDS_MET
        self.SDS_INV = SDS_INV
        self.AERNAMES = AERNAMES
        self.inFile  = inFile
        self.outFile = outFile
        self.albedoFile = albedoFile
        self.albedoType = albedoType
        self.rcFile  = rcFile
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom


        # initialize empty lists
        for sds in self.SDS_AER+self.SDS_MET+self.SDS_INV:
            self.__dict__[sds] = []

        # Read in data model data
        self.readSampledGEOS()

        # Make lists into arrays
        for sds in self.SDS_AER+self.SDS_MET:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        # convert isotime to datetime
        self.tyme = []
        for isotime in self.isotime:
            self.tyme.append(isoparser(''.join(isotime)))

        self.tyme = np.array(self.tyme)

        # Start out with all good obs
        self.nobs  = len(self.tyme)
        self.iGood = np.ones([self.nobs]).astype(bool)


        # Read in surface data
        albedoReader = getattr(self,SurfaceFuncs[albedoType])
        albedoReader()


        # Calculate aerosol optical properties
        self.computeMie()

        # Calculate atmospheric profile properties needed for Rayleigh calc
        self.computeAtmos()

        # Calculate Scene Geometry
        self.VZA = VZA
        self.hgtss = hgtss
        self.calcAngles()

        if self.nobs > 0:
            # Land-Sea Mask
            self.LandSeaMask()


    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)   

        # Make lists into arrays
        for sds in self.SDS_INV:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])            

        if self.albedoType in LandAlbedos:
            iGood = self.FRLAND >= 0.99
        else:
            iGood = self.FRLAND < 0.99

        self.iGood = self.iGood & iGood
        self.nobs  = np.sum(self.iGood)

    # --
    def computeAtmos(self):

        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze 
        self.te = te

    # --
    def calcAngles(self):
        SZA   = []
        SAA   = []
        VAAf  = []
        VAAb  = []
        for i,tyme in enumerate(self.tyme):
            if i == len(self.tyme)-1:
                CLAT = self.LATITUDE[i-1]
                CLON = self.LONGITUDE[i-1]
            else:
                CLAT = self.LATITUDE[i+1]
                CLON = self.LONGITUDE[i+1]
            
            SLAT = self.LATITUDE[i]
            SLON = self.LONGITUDE[i]
            sat_angles = LidarAngles_.satangles(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                SLAT,SLON,
                                                0.0,
                                                self.hgtss)

            SZA.append(sat_angles[3][0])
            SAA.append(sat_angles[2][0])

            if i == len(self.tyme)-1:
                VAAb.append(sat_angles[0][0])
                vaaf = sat_angles[0][0] - 180.0
                if vaaf < 0:
                    vaaf = vaaf + 360.
                VAAf.append(vaaf)

            else:
                VAAf.append(sat_angles[0][0])
                vaab = sat_angles[0][0]+180
                if vaab > 360.:
                    vaab = vaab - 360.
                VAAb.append(vaab)


        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.VAAf = np.array(VAAf)
        self.VAAb = np.array(VAAb)

        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I] = saa[I] - 360.

        self.RAAb = self.VAAb - saa
        self.RAAf = self.VAAf - saa

        # Limit SZAs
        self.iGood = self.SZA < 80
        self.nobs = np.sum(self.iGood)     

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_AER:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)

        col = 'met_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)        
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_MET:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)

    # --- 
    def readSampledMODISBRDF(self):
        """
        Read in MODIS BRDF kernel weights
        that have already been sampled on lidar track
        """
        chs = str(int(self.channel))
        SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs

        if self.verbose:
            print 'opening abledo file ',self.albedoFile
        nc = Dataset(self.albedoFile)

        for sds in SDS:
            self.__dict__[sds] = nc.variables[sds][:]

        nc.close()

        
        nobs = len(self.__dict__[sds])
        for sds in SDS:
            self.__dict__[sds].shape = (1,1,nobs)

        # [nkernel,nch,nobs]
        self.kernel_wt = np.append(self.__dict__['Riso'+chs],self.__dict__['Rgeo'+chs],axis=0)
        self.kernel_wt = np.append(self.kernel_wt,self.__dict__['Rvol'+chs],axis=0)

        param1 = np.array([2]*nobs)
        param2 = np.array([1]*nobs)

        param1.shape = (1,1,nobs)
        param2.shape = (1,1,nobs)

        # [nparam,nch,nobs]
        self.param = np.append(param1,param2,axis=0)

    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        tau,ssa,g,pmom = getAOPvector(self,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=True,
                                 rcfile=self.rcFile,
                                 nMom=self.nMom)
        self.tau = tau  #(km,nch,nobs)
        self.ssa = ssa  #(km,nch,nobs)
        self.g   = g    #(km,nch,nobs)
        self.pmom = pmom  #(km,nch,nobs,nMom,nPol)

    # --
    def runExt(self):
        """
        run ext_sampler.py 
        """
        col = 'aer_Nv'
        if self.verbose: 
            print 'running ext_sampler on file',self.inFile.replace('%col',col)

        outDir = os.path.dirname(self.outFile)
        instname = os.path.basename(self.inFile).split('.')[0]
        date_ch   = os.path.basename(self.inFile).split('.')[-2]
        outFile = '{}/{}.lc2.ext.{}.nc'.format(outDir,instname,date_ch)
        Options =     " --input=" + self.inFile.replace('%col',col)      + \
                      " --output=" + outFile       + \
                      " --rc=" + self.rcFile      + \
                      " --format=NETCDF4_CLASSIC"      + \
                      " --channel=%d" %self.channel     
                      

        cmd = 'ext_sampler.py {} '.format(Options)  

        self.cmd = cmd      


    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """
        sza  = self.SZA[self.iGood]
        raaf = self.RAAf[self.iGood]
        raab = self.RAAb[self.iGood]

        tau = self.tau[:,:,self.iGood]
        ssa = self.ssa[:,:,self.iGood]
        pmom = self.pmom[:,:,self.iGood,:,:]
        pe   = self.pe[:,self.iGood]
        ze   = self.ze[:,self.iGood]
        te   = self.te[:,self.iGood]

        kernel_wt = self.kernel_wt[:,:,self.iGood]
        param     = self.param[:,:,self.iGood]

        nangles = len(self.VZA)
        #for i in range(nangles):
        i = 0

        vza  = [self.VZA[i]]*self.nobs


        # with BRDF surface
        if self.albedoType == 'MODIS_BRDF':
            
            result = VLIDORT_LIDAR_.vector_brdf_modis(self.channel,tau, ssa, pmom, 
                                                      pe, ze, te, 
                                                      kernel_wt, param, 
                                                      SZA, RAAf, vza, 
                                                      MISSING,
                                                      self.verbose)        

        
        self.writenc(result)
        
    #---
    def writeNC (self,result,zlib=True):
        """
        Write a NetCDF file vlidort output
        """
        radiance_VL_SURF,reflectance_VL_SURF, ROT, BR, Q, U, rc = result
        km = 72

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'VLIDORT Simulation of GEOS-5 multiangle polarized reflectance'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        ls = nc.createDimension('ls',19)
        nz = nc.createDimension('lev',km)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        na = nc.createDimension('view_angles',2*len(self.VZA))

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nctrj       = Dataset(self.inFile.replace('%col',col))        
        _copyVar(nctrj,nc,u'trjLon',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'trjLat',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=verbose)   
        nctrj.close()


        vza = nc.createVariable('seonsor_zenith','f4',('time','view_angles',),zlib=zlib)
        vza.long_name = "sensor viewing zenith angle (VZA)"
        vza.missing_value = MISSING
        vza.units = "degrees (positive forward view)"

        vaa = nc.createVariable('sensor_azimuth','f4',('time','view_angles',),zlib=zlib)
        vaa.long_name = "sensor viewing azimuth angle (VAA)"
        vaa.missing_value = MISSING
        vaa.units = "degrees clockwise from North"

        sza = nc.createVariable('seonsor_zenith','f4',('time',),zlib=zlib)
        sza.long_name = "solar zenith angle (SZA)"
        sza.missing_value = MISSING
        sza.units = "degrees"

        saa = nc.createVariable('solar_azimuth','f4',('time',),zlib=zlib)
        saa.long_name = "solar azimuth angle (SAA)"
        saa.missing_value = MISSING
        saa.units = "degrees clockwise from North"


        # Write VLIDORT Outputs
        # ---------------------
        ref = nc.createVariable('I','f4',('time','view_angles',),zlib=zlib)
        ref.standard_name = '%.2f nm TOA Reflectance' %self.channel
        ref.long_name     = '%.2f nm Top of the Atmosphere Reflectance' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        ref._FillValue    = MISSING

        q = nc.createVariable('I','f4',('time','view_angles',),zlib=zlib)
        q.standard_name = '%.2f nm TOA Reflectance' %self.channel
        q.long_name     = '%.2f nm Top of the Atmosphere Reflectance' %self.channel
        q.missing_value = MISSING
        q.units         = "None"
        q._FillValue    = MISSING        

        sref = nc.createVariable('surf_reaf','f4',('time','view_angles',),zlib=zlib)
        ref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        ref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        ref._FillValue    = MISSING

        rot = nc.createVariable('ROT','f4',('time',),zlib=zlib)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        rot._FillValue    = MISSING



        for n, name in enumerate(MieVarsNames):

            var = squeeze(MieVars[name])
            size = len(var.shape)
            if options.station:
                if size == 2:
                    var = asarray([var])
                    size = len(var.shape)
            if size == 3:
                dim = ('station','time','lev')
            if size == 2:
                dim = ('time','lev')
            if size == 1:
                dim = ('time')
            this = nc.createVariable(name,'f4',dim,zlib=zlib)
            this.standard_name = name
            this.units = MieVarsUnits[n]
            this.missing_value = MAPL_UNDEF
            if options.station:
                this[:] = transpose(var,(0,2,1))
            else:
                this[:] = transpose(var)

        # Close the file
        # --------------
        nc.close()

        if options.verbose:
            print " <> wrote %s file %s"%(options.format,options.outFile)
    

def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,01,01,00)
    nymd     = str(date.date()).replace('-','')
    hour     = str(date.hour).zfill(2)
    format   = 'NETCDF4_CLASSIC'

    inDir    = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/LevelB/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    inFile   = '{}/calipso-g5nr.lb2.%col.{}_{}z.nc4'.format(inDir,nymd,hour)
    albedoDir = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BRDF/MCD43C1/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    albedoFile   = '{}/calipso-g5nr.lb2.brdf.{}_{}z.nc4'.format(albedoDir,nymd,hour)
    albedoType = 'MODIS_BRDF'

    channel  = 470
    chd      = get_chd(channel)
    outDir    = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/LevelC2/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    outFile   = '{}/calipso-g5nr.vlidort.vector.MCD43C.{}_{}z_{}nm.nc4'.format(outDir,nymd,hour,chd)
    
    rcFile   = 'Aod_EOS.rc'
    polarname = 'POLDER'
    orbit     = 'LEO'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = POLAR_VLIDORT(inFile,outFile,rcFile,
                            albedoFile,albedoType,
                            channel,
                            VZAdic[polarname],
                            HGTdic[orbit],
                            verbose=verbose)

   
    # Run ext_sampler
    #vlidort.runExt()

    # Run VLIDORT
    # if vlidort.nobs > 0:
    #     radiance_VL_SURF,reflectance_VL_SURF, ROT, BR, Q, U, rc = vlidort.runVLIDORT()