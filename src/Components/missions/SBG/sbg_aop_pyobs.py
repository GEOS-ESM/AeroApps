#!/usr/bin/env python3

"""
    Calculates optical property inputs for VLIDORT

    Adapted from sbg_vlidort.py
    Patricia Castellanos, Feb 2025

    Adapted from accp_polar_vlidort.py
    Patricia Castellanos, Jul 2024

    Adapted from polar_vlidort.py and lidar_vlidort.py
    Patricia Castellanos, Jan 2020

"""
import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL.config     import Config
from   netCDF4 import Dataset
import numpy   as np
import xarray  as xr
from pyobs import mietable as mt
from pyobs.aop import G2GAOP
import yaml
from pyobs.constants import MAPL_GRAV as GRAV
from pyobs.constants import MAPL_RDRY as RGAS
from pyobs.constants import MAPL_KAPPA as KAPPA


from py_leo_vlidort.vlidort import VLIDORT, WrapperFuncs
from py_leo_vlidort.copyvar  import _copyVar
from multiprocessing import Pool

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = [] #[CLDTOT]
ncALIAS = {'LONGITUDE': 'longitude',
           'LATITUDE' : 'latitude'}

MISSING = np.float32(-1.e+20)


class OPTICS_VLIDORT(VLIDORT,G2GAOP):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on satellite track
    """
    def __init__(self,ich,inFile,outFile,rcFile,brdfFile,
                instname,dryrun,
                verbose=False):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.AERNAMES    = AERNAMES
        self.inFile      = inFile
        self.outFile     = outFile
        self.rcFile      = rcFile
        self.brdfFile    = brdfFile
        self.verbose     = verbose
        self.instname   = instname


        # load optics tables
        self.getMie()

        # get granule dimensions
        self.getDims()

        # do one scanline at a time
        for ityme in range(self.ntyme)[1:2]:
            print('ityme ',ityme,self.ntyme)

            # Read in model data
            self.readSampledGEOS(ityme)

            # Calculate atmospheric profile properties needed for Rayleigh calc
            self.getEdgeVars()    

            if not dryrun:
                # Run VLIDORT
                self.calcOptics(ityme,ich)
#                self.writeNC(ityme)
    #---
    def getEdgeVars(self):
        """
        Get altitude, pressure, and temperature and edge of layers
        """

        # Get layer thicnkness from DELP & AIRDENS
        # DELP: Pa = kg m-1 s-2
        # AIRDENS: kg m-3
        # GRAV: m s-2
        # -----------------------------------------
        rhodz = self.aer['DELP'] / GRAV
        dz = rhodz / self.aer['AIRDENS']       # column thickness in m

        # add up the thicknesses to get edge level altitudes and pressures
        npts, nlev = dz.shape
        ze = np.array([dz[:,i:].sum(axis=1) for i in range(nlev)])
        # append surface level, altitude = 0
        # [nlev+1,nacross]
        self.ze = np.append(ze,np.zeros([1,npts]),axis=0)

        ptop = 1. # Pa
        pe = np.zeros([self.nlev+1,self.nacross])
        pe[0,:] = ptop
        for ilev in range(self.nlev):
            pe[ilev+1,:] = pe[ilev,:] + self.aer['DELP'][:,ilev]

        self.pe = pe

        # get the mid-level pressures and temperatures
        self.pm = (self.pe[:-1,:] + self.pe[1:,:])*0.5
        self.tm = self.pm/(self.aer['AIRDENS'].T*RGAS)

        # get the edge level temperature
        self.te = np.zeros([self.nlev+1,self.nacross])
        self.te[0,:] = self.tm[0,:]  #isothermal a highest level
        for ilev in range(1,self.nlev):
            alpha = np.log(self.pe[ilev,:]/self.pm[ilev-1,:])/np.log(self.pm[ilev,:]/self.pm[ilev-1,:])
            self.te[ilev,:] = self.tm[ilev-1] + alpha*(self.tm[ilev,:] - self.tm[ilev-1,:])

        # dry adiabatic
        self.te[self.nlev,:] = self.tm[self.nlev-1,:]*(self.pe[self.nlev,:]/self.pm[self.nlev-1,:])**KAPPA

    #---
    def getMie(self):
        """
        Use pyobs utilities to load mietables
        """
        self.mieTable = yaml.safe_load(open(self.rcFile))
        self.vector = True

        # load mietables
        # ---------------------------
        self.p, self.m = 0,0
        for s in self.mieTable:
            m = self.mieTable[s]
            m['mie'] = mt.MIETABLE(m['monoFile'])

        
            dims = dict(self.mieTable[s]['mie'].ds.sizes)
            self.p = max(self.p,dims['p'])
            self.m = max(self.m,dims['m'])

    #--
    def getDims(self):
        """
        Get granule dimensions
        """
        col = 'aer_Nv'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col)) 
        self.ntyme,self.nlev,self.nacross = ds.sizes['time'],ds.sizes['lev'],ds.sizes['ncross']
        self.nobs = self.nacross

        if self.verbose:
            print('opening BRDF file ',self.brdfFile)       
        ds = xr.open_dataset(self.brdfFile) 

        self.channels = ds.nwav.values  # microns
        self.nch = len(self.channels)
        self.channels = self.channels*1e3  # nm

        # split this into bands
        self.sband = np.arange(0,self.nch,4)
        self.nband = list(self.sband[1:]-self.sband[:-1])
        self.nband.append(self.nch - self.sband[-1] + 1)

    #---
    def readSampledGEOS(self,ityme):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))

        inList = [self.inFile.replace('%col',col)]

        if len(self.SDS_MET) > 0:
            col = 'met_Nv'
            inList.append(self.inFile.replace('%col',col))
            if self.verbose:
                print('opening file',self.inFile.replace('%col',col))

        self.aer = xr.open_mfdataset(inList)
        # make arrays [nobs,nlev]
        self.aer = self.aer.squeeze()
        self.aer = self.aer.isel(time=ityme)
        self.aer = self.aer.transpose()
        self.aer = self.aer.chunk(chunks={'ncross':10,'lev':5})

    # ---
    def calcOptics(self,ityme,ich):
        """
        Calls VLIDORT 
        """

        # Initiate output arrays
        nobs   = self.nobs
        nlev   = self.nlev
        nch    = self.nch
        nacross = self.nacross
        self.channel = [self.channels[ich]]

        # Calculate ROT
        args = [self.channel, self.pe.astype('float64'), self.ze.astype('float64'), self.te.astype('float64'), MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args)  
        self.depol_ratio = depol_ratio         
        self.ROT = ROT 

        self.aop = self.getAOPrt(wavelength=self.channels[ich],vector=True)

        # need to reshape these to [nlev,nobs]
        self.TAU = self.aop.AOT.astype('float64').transpose().to_numpy()
        self.SSA = self.aop.SSA.astype('float64').transpose().to_numpy()
        self.PMOM = self.aop.PMOM.astype('float64').transpose('lev','ncross','m','p').to_numpy()

    #---
    def writeNC (self,ityme,zlib=True):
        """
        Write a NetCDF file vlidort output
        """

        if ityme == 0:
            if not os.path.exists(os.path.dirname(self.outFile)):
                os.makedirs(os.path.dirname(self.outFile))

            # Open NC file
            # ------------
            nc = Dataset(self.outFile,'w',format='NETCDF4_CLASSIC')

            # Set global attributes
            # ---------------------
            nc.title = 'VLIDORT-GEOS-SBG Simulator'
            nc.institution = 'NASA/Goddard Space Flight Center'
            nc.source = 'Global Model and Assimilation Office'
            nc.history = 'VLIDORT simulation run on sampled GEOS'
            nc.references = 'n/a'
            nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
            nc.Conventions = 'CF'
            nc.inFile = self.inFile
         
            # Create dimensions
            # -----------------
            nt = nc.createDimension('time',self.ntyme)
            nz = nc.createDimension('lev',self.nlev)
            nx = nc.createDimension('across',self.nacross)
            nh = nc.createDimension('ch',1)
            pp = nc.createDimension('npol',self.p)
            nm = nc.createDimension('nmom',self.m)

            # Coordinate variables
            # --------------------
            col = 'aer_Nv'
            if self.verbose: 
                print('opening file',self.inFile.replace('%col',col))
            nctrj       = Dataset(self.inFile.replace('%col',col))        
            _copyVar(nctrj,nc,'trjLon',dtype='f4',zlib=False,verbose=self.verbose)
            _copyVar(nctrj,nc,'trjLat',dtype='f4',zlib=False,verbose=self.verbose)
            _copyVar(nctrj,nc,'time', dtype='i4',zlib=False,verbose=self.verbose)
            _copyVar(nctrj,nc,'lev', dtype='f4',zlib=False,verbose=self.verbose)

            dim = ('ch',)
            c = nc.createVariable('ch','f4',dim,zlib=zlib)
            c.standard_name = 'channel'
            c.long_name = 'wavelength'
            c.units = "nm"
            c[:] = self.channel

            dim = ('npol',)
            p = nc.createVariable('npol','i4',dim,zlib=zlib)
            p.standard_name = 'npol'
            p.long_name = 'number of elements in the scattering matrix'
            p.units = "None"
            p[:] = np.arange(self.p)

            dim = ('nmom',)
            m = nc.createVariable('nmom','i4',dim,zlib=zlib)
            m.standard_name = 'nmom'
            m.long_name = 'number of expansion moments of the scattering matrix'
            m.units = "None"
            m[:] = np.arange(self.m)


            nctrj.close()

        else:
            # Open NC file
            # ------------
            nc = Dataset(self.outFile,'a')            

        # Write VLIDORT Outputs
        # ---------------------
        if ityme == 0:
            dim = ('time','ch','lev','across')
            tau = nc.createVariable('TAU','f4',dim,zlib=zlib,fill_value=MISSING)
            tau.standard_name = 'tau' 
            tau.long_name     = 'aerosol optical thickness' 
            tau.missing_value = MISSING
            tau.units         = "None"
        else:
            tau = nc.variables['TAU']
        tau[ityme,0,:,:]            = self.TAU

        if ityme == 0:
            dim = ('time','ch','lev','across')
            ssa = nc.createVariable('SSA','f4',dim,zlib=zlib,fill_value=MISSING)
            ssa.standard_name = 'ssa' 
            ssa.long_name     = 'single scattering albedo' 
            ssa.missing_value = MISSING
            ssa.units         = "None"
        else:
            ssa = nc.variables['SSA']
        ssa[ityme,0,:,:]            = self.SSA

        if ityme == 0:
            dim = ('time','ch','lev','across','nmom','npol')
            p = nc.createVariable('PMOM','f4',dim,zlib=zlib,fill_value=MISSING)
            p.standard_name = 'pmom' 
            p.long_name     = 'spherical expansion coefficient of the phase function' 
            p.missing_value = MISSING
            p.units         = "None"
        else:
            p = nc.variables['PMOM']
        p[ityme,0,:,:,:,:]            = self.PMOM

        if ityme == 0:
            dim = ('time','lev','across','ch')
            r = nc.createVariable('ROT','f4',dim,zlib=zlib,fill_value=MISSING)
            r.standard_name = 'ROT'
            r.long_name     = 'rayleigh optical thickness' 
            r.missing_value = MISSING
            r.units         = "None"
        else:
            r = nc.variables['ROT']
        r[ityme,:,:,:]            = self.ROT

        if ityme == 0:
            dim = ('ch',)
            d = nc.createVariable('depol_ratio','f4',dim,zlib=zlib,fill_value=MISSING)
            d.standard_name = 'depol_ratio' 
            d.long_name     = 'depolarization retio' 
            d.missing_value = MISSING
            d.units         = "None"
            d[:]            = self.depol_ratio

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print(" <> wrote %s %d %d"%(self.outFile,ityme,self.ntyme))

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours   = 1
    rcFile     = 'm2_aop.yaml'
    albedoType = None

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with track input file names")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("ich",type=int,
                        help="channel index")

    parser.add_argument("--rcfile",default=rcFile,
                        help="rcFile (default=%s)"%rcFile)

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)


    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    args = parser.parse_args()
    rcFile         = args.rcfile

    # Parse prep config
    # -----------------
    cf             = Config(args.inst_pcf,delim=' = ')
    instname       = cf('instname')

    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()

    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')
    outTemplate    = cf('outDir')    + '/' + cf('outFile')

    try:
        brdfTemplate = cf('brdfDir') + '/' + cf('brdfFile')
    except:
        brdfTemplate = None


    # Loop through dates, running VLIDORT
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)
        minute = str(date.minute).zfill(2)

        inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%instname',instname).replace('%band',str(args.ich).zfill(3))

        if brdfTemplate is None:
            brdfFile = None
        else:
            brdfFile = brdfTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)



        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>channel:    ',args.ich)
        print('>>>inFile:    ',inFile)
        print('>>>outFile:   ',outFile)
        print('>>>rcFile:    ',rcFile)
        print('>>>brdfFile:  ',brdfFile)
        print('>>>verbose:   ',args.verbose)
        print('++++End of arguments+++')
        
        vlidort = OPTICS_VLIDORT(args.ich,inFile,outFile,rcFile,brdfFile, 
                            instname,
                            args.dryrun,
                            verbose=args.verbose)


        date += Dt
