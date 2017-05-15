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

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}


nMom     = 300

VZAdic = {'POLDER': np.array([3.66, 11., 18.33, 25.66, 33, 40.33, 47.66, 55.0])}

HGTdic = {'LEO': 705,
          'ISS': 400}


class POLAR_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,outFile,rcFile,channel,VZA,hgtss,verbose=False):
        self.SDS_AER = SDS_AER
        self.SDS_MET = SDS_MET
        self.AERNAMES = AERNAMES
        self.inFile  = inFile
        self.outFile = outFile
        self.rcFile  = rcFile
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom


        # initialize empty lists
        for sds in self.SDS_AER+self.SDS_MET:
            self.__dict__[sds] = []

        # Read in data from path
        self.readSampledGEOS()


        # Make lists into arrays
        for sds in self.SDS_AER+self.SDS_MET:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        # convert isotime to datetime
        self.tyme = []
        for isotime in self.isotime:
            self.tyme.append(isoparser(''.join(isotime)))

        self.tyme = np.array(self.tyme)


        # Calculate aerosol optical properties
        self.computeMie()

        # Calculate atmospheric profile properties needed for Rayleigh calc
        self.computeAtmos()

        # Calculate Scene Geometry
        self.VZA = VZA
        self.hgtss = hgtss
        self.calcAngles()

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
    def readSampledBRDF            

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


    def brdfVLIDORT(self):
        """
        Calls VLIDORT with BRDF surface
        """
        radiance_VL_SURF,reflectance_VL_SURF, ROT, BR, Q, U, rc = VLIDORT_LIDAR_.vector_brdf_modis(self.tau, self.ssa, self.pmom, 
                                                                                                    self.pe, self.ze, self.te, 
                                                                                                    self.kernel_wt, self.param, 
                                                                                                    self.SZA, self.RAAf, self.VZA, 
                                                                                                    MISSING,
                                                                                                    self.verbose)        




#---
def writeNC ( stations, lons, lats, tyme, isotimeIn, MieVars, MieVarsNames,
              MieVarsUnits, inFile, outFile, zlib=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables along the satellite track
    described by (lon,lat,tyme).
    """
    km = 72

    # Open NC file
    # ------------
    nc = Dataset(outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.title = 'GEOS-5 Sampled Aerosol Optical Properties File'
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from sampled GEOS-5 collections'
    nc.references = 'n/a'
    nc.comment = 'This file contains sampled GEOS-5 aerosol optical properties.'
    nc.contact = 'Ed Nowottnick <edward.p.nowottnick@nasa.gov>'
    nc.Conventions = 'CF'
    nc.inFile = inFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',len(tyme))
    if options.station:
        ns = nc.createDimension('station',len(stations))
    ls = nc.createDimension('ls',19)
    if km>0:
        nz = nc.createDimension('lev',km)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)

    if options.station:
        # Station names
        # -------------
        stnName_ = nc.createVariable('stnName','S1',('station','ls'),zlib=zlib)
        stnName_.long_name = 'Station Names'
        stnName_.axis = 'e'
        stnName_[:] = stations[:]   

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    isot0 = isotimeIn[0]
    date0 = ''.join(isot0[:10])
    time0 = ''.join(isot0[-8:])
    time.units = 'seconds since '+date0+' '+time0
    time[:] = tyme
    if km > 0: # pressure level not supported yet
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = 'km'
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = range(1,km+1)

    # Add fake dimensions for GrADS compatibility
    # -------------------------------------------
    x = nc.createVariable('x','f4',('x',),zlib=zlib)
    x.long_name = 'Fake Longitude for GrADS Compatibility'
    x.units = 'degrees_east'
    x[:] = zeros(1)
    y = nc.createVariable('y','f4',('y',),zlib=zlib)
    y.long_name = 'Fake Latitude for GrADS Compatibility'
    y.units = 'degrees_north'
    y[:] = zeros(1)
    if options.station:
        e = nc.createVariable('station','i4',('station',),zlib=zlib)
        e.long_name = 'Station Ensemble Dimension'
        e.axis = 'e'
        e.grads_dim = 'e'
        e[:] = range(len(stations))
    
    # Lat/Lon Coordinates
    # ----------------------
    if options.station:
        lon = nc.createVariable('longitude','f4',('station',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = lons[:]
        lat = nc.createVariable('latitude','f4',('station',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = lats[:]
    else:
        lon = nc.createVariable('longitude','f4',('time',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = lons[:]
        lat = nc.createVariable('latitude','f4',('time',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = lats[:]        
    
    # Time in ISO format if so desired
    # ---------------------------------
    isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
    isotime.long_name = 'Time (ISO Format)'
    isotime[:] = isotimeIn[:]

    # Write each variable
    # --------------------------------------------------
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
    
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,01,01,00)
    nymd     = str(date.date()).replace('-','')
    hour     = str(date.hour).zfill(2)
    format   = 'NETCDF4_CLASSIC'
    inFile   = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/LevelB/Y2006/M01/calipso-g5nr.lb2.%col.rc.{}_{}z.nc4'.format(nymd,hour)
    outFile  = 'polar_vlidort.nc'
    channel  = 470
    rcFile   = 'Aod_EOS.rc'
    polarname = 'POLDER'
    orbit     = 'LEO'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = POLAR_VLIDORT(inFile,outFile,rcFile,channel,VZAdic[polarname],HGTdic[orbit],verbose=verbose)

   
    # Read in Surface Dataset
    vlidort.readSampledBRDF()

    # Run VLIDORT
    vlidort.brdfVLIDORT()