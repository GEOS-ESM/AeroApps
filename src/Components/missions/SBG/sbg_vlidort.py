#!/usr/bin/env python3

"""
    Calculates polarized TOA radiance for sbg imager
    Model fields have already been sampled using trj_sampler
    Uses POLAR_VLIDORT as parent class

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


from py_leo_vlidort.vlidort import VLIDORT, get_chd, WrapperFuncs, MODIS_BRDF_run
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
SDS_INV = ['FRLAND']
SDS_ANG = ['SZA','SAA','VZA','VAA']
ncALIAS = {'LONGITUDE': 'longitude',
           'LATITUDE' : 'latitude',
           'SZA'      : 'sza',
           'SAA'      : 'saa',
           'VZA'      : 'vza',
           'VAA'      : 'vaa'}

MISSING = np.float32(-1.e+20)


class ACCP_POLAR_VLIDORT(VLIDORT,G2GAOP):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on satellite track
    """
    def __init__(self,inFile,outFile,rcFile,albedoType,
                instname,dryrun,
                nstreams=12,
                plane_parallel=True,
                brdfFile=None,
                verbose=False):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.SDS_INV     = SDS_INV
        self.SDS_ANG     = SDS_ANG
        self.AERNAMES    = AERNAMES
        self.inFile      = inFile
        self.outFile     = outFile
        self.albedoType  = albedoType
        self.rcFile      = rcFile
        self.verbose     = verbose
        self.brdfFile    = brdfFile
        self.nstreams    = nstreams
        self.plane_parallel = plane_parallel
        self.instname   = instname


        # load optics tables
        self.getMie()

        # get granule dimensions
        self.getDims()

        # Start out with all good obs
        self.iGood = np.ones([self.nobs]).astype(bool)

        # do one scanline at a time
        for ityme in range(self.ntyme)[0:1]:
            print('ityme ',ityme,self.ntyme)
            # Read in precalculated Scene Geometry
            # limit iGood to sza < 80
            self.readAngles(ityme)

            if self.nobs > 0:

                # Read in model data
                self.readSampledGEOS(ityme)

                # Read in surface data
                self.readSampledAMESBRDF(ityme)

                # Calculate atmospheric profile properties needed for Rayleigh calc
                self.getEdgeVars()    

                # Read in precalculated Scene Geometry
                # limit iGood to sza < 80
                self.readAngles(ityme)

                # Land-Sea Mask
                # limit iGood to land pixels
                self.LandSeaMask(ityme)       


                if not dryrun:
                    # Run VLIDORT
                    self.runVLIDORT()
                    self.writeNC(ityme) 
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

    # ---
    def LandSeaMask(self,ityme):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = ds[sds_].isel(time=ityme).squeeze().values
            self.__dict__[sds] = var

        self.iLand = self.FRLAND >= 0.99
        self.iSea  = self.FRLAND < 0.99

        # self.iGood = self.iGood & iGood
        self.nobsLand  = np.sum(self.iGood & self.iLand)
        self.nobsSea   = np.sum(self.iGood & self.iSea)

        self.iGood = self.iGood & self.iLand
        self.nobs = np.sum(self.iGood)


    #---
    def readSampledAMESBRDF(self,ityme):
        """
        Read in AMES BRDF kernel weights
        that have already been sampled on swath
        """
        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
        ds = xr.open_dataset(self.brdfFile)

        self.channels = ds.nwav.values  # microns
        self.nch = len(self.channels)
        self.channels = self.channels*1e3  # nm
        self.Riso = ds.Ki.isel(time=ityme).squeeze().values
        self.Rgeo = ds.Kg.isel(time=ityme).squeeze().values
        self.Rvol = ds.Kv.isel(time=ityme).squeeze().values


        self.Riso.shape = (1,self.nch,self.nacross)
        self.Rgeo.shape = (1,self.nch,self.nacross)
        self.Rvol.shape = (1,self.nch,self.nacross)

        # [nkernel,nch,nacross]
        self.kernel_wt = np.append(self.Riso,self.Rgeo,axis=0)
        self.kernel_wt = np.append(self.kernel_wt,self.Rvol,axis=0)
        self.kernel_wt = self.kernel_wt.reshape(3,self.nch,self.nacross)

        param1 = np.array([2]*self.nch*self.nacross)
        param2 = np.array([1]*self.nch*self.nacross)

        param1.shape = (1,self.nch,self.nacross)
        param2.shape = (1,self.nch,self.nacross)

        # [nparam,nch,nobs]
        self.RTLSparam = np.append(param1,param2,axis=0)
        self.RTLSparam = self.RTLSparam.reshape(2,self.nch,self.nacross)

        # filter for nans
        iGood = ~np.isnan(self.kernel_wt[0,0,:])
        self.iGood = self.iGood & iGood
        self.nobs = np.sum(self.iGood)

    def readAngles(self,ityme):
        """
        Read in viewing and solar Geometry from angFile
        """

        col = self.instname
        if self.verbose: 
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_ANG:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = ds[sds_].isel(time=ityme).squeeze().values
            self.__dict__[sds] = var

        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I] = saa[I] - 360.

        RAA = self.VAA - saa
        RAA[RAA < 0] = RAA[RAA<0]+360.0
        self.RAA = RAA

        # Limit SZAs
        self.iGood = self.iGood & (self.SZA < 80)
        self.nobs = np.sum(self.iGood)         

    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """

        # Initiate output arrays
        iGood  = np.arange(self.nacross)[self.iGood]
        nobs   = self.nobs
        nlev   = self.nlev
        nch    = self.nch
        nacross = self.nacross
        self.I = np.ones([nacross,nch])*MISSING
        self.Q = np.ones([nacross,nch])*MISSING
        self.U = np.ones([nacross,nch])*MISSING
        self.reflectance = np.ones([nacross,nch])*MISSING
        self.surf_reflectance = np.ones([nacross,nch])*MISSING
        self.BR_Q = np.ones([nacross,nch])*MISSING
        self.BR_U = np.ones([nacross,nch])*MISSING

        # Calculate ROT
        args = [self.channels, self.pe.astype('float64'), self.ze.astype('float64'), self.te.astype('float64'), MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args)  
        self.depol_ratio = depol_ratio         
        self.ROT = ROT 

        p = Pool(120) 
        # loop through channels
        for ich in np.arange(self.nch):
            print('ich ',ich,self.nch)            
            self.channel = [self.channels[ich]]
            self.aop = self.getAOPrt(wavelength=self.channels[ich],vector=True)
            rot  = ROT[:,iGood,ich:ich+1]
            depol_ratio = [self.depol_ratio[ich]]

            # need to reshape these to [nlev,nch,nobs]
            tau = self.aop.AOT[iGood,:].astype('float64').expand_dims(dim={"ch": 1},axis=1).transpose().to_numpy()
            ssa = self.aop.SSA[iGood,:,].astype('float64').expand_dims(dim={"ch": 1},axis=1).transpose().to_numpy()
            pmom = self.aop.PMOM[iGood,:,:,:].astype('float64').transpose('lev','ncross','m','p').expand_dims(dim={"ch": 1},axis=1).to_numpy()

            pe   = self.pe[:,iGood].astype('float64')
            ze   = self.ze[:,iGood].astype('float64')
            te   = self.te[:,iGood].astype('float64')


            # Get surface data
            param     = self.RTLSparam[:,ich:ich+1,iGood].astype('float64')
            kernel_wt = self.kernel_wt[:,ich:ich+1,iGood].astype('float64')

            vza = self.VZA[iGood].astype('float64')
            sza = self.SZA[iGood].astype('float64')
            raa = self.RAA[iGood].astype('float64')
            
            I = []
            Q = []
            U = []
            reflectance = []
            surf_reflectance = []
            BR_Q = []
            BR_U = []
            args = [(self.channel, self.nstreams, self.plane_parallel, rot[:,i:i+1,:], depol_ratio, 
                    tau[:,:,i:i+1], ssa[:,:,i:i+1], pmom[:,:,i:i+1,:,:],
                    pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                    kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
                    sza[i:i+1], raa[i:i+1], vza[i:i+1],
                    MISSING,
                    self.verbose) for i in range(nobs)]
            result = p.map(MODIS_BRDF_run,args)
            
            for r in result:
                I_r,Q_r,U_r,reflectance_r,surf_reflectance_r,BR_Q_r,BR_U_r = r
                I.append(I_r)
                Q.append(Q_r)
                U.append(U_r)
                reflectance.append(reflectance_r)
                surf_reflectance.append(surf_reflectance_r)
                BR_Q.append(BR_Q_r)
                BR_U.append(BR_U_r)
            I = np.concatenate(I)
            Q = np.concatenate(Q)
            U = np.concatenate(U)
            reflectance = np.concatenate(reflectance)
            surf_reflectance = np.concatenate(surf_reflectance)
            BR_Q = np.concatenate(BR_Q)
            BR_U = np.concatenate(BR_U)

            self.I[iGood,ich] = np.squeeze(I)
            self.reflectance[iGood,ich] = np.squeeze(reflectance)
            self.surf_reflectance[iGood,ich] = np.squeeze(surf_reflectance)
            self.Q[iGood,ich] = np.squeeze(Q)
            self.U[iGood,ich] = np.squeeze(U) 
            self.BR_Q[iGood,ich] = np.squeeze(BR_Q)
            self.BR_U[iGood,ich] = np.squeeze(BR_U) 
        p.close()
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
            nh = nc.createDimension('ch',self.nch)
            nx = nc.createDimension('across',self.nacross)

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
            nctrj.close()

        else:
            # Open NC file
            # ------------
            nc = Dataset(self.outFile,'a')            

        # Write VLIDORT Outputs
        # ---------------------
        if ityme == 0:
            dim = ('time','across','ch')
            ref = nc.createVariable('toa_reflectance','f4',dim,zlib=zlib,fill_value=MISSING)
            ref.standard_name = 'TOA Reflectance' 
            ref.long_name     = 'reflectance at the top of the atmosphere' 
            ref.missing_value = MISSING
            ref.units         = "None"
        else:
            ref = nc.variables['toa_reflectance']
        ref[ityme,:,:]            = self.reflectance

        if ityme == 0:
            i = nc.createVariable('I','f4',dim,zlib=zlib,fill_value=MISSING)
            i.standard_name = 'TOA I' 
            i.long_name     = 'intensity at the top of the atmosphere' 
            i.missing_value = MISSING
            i.units         = "W m-2 sr-1 nm-1"
        else:
            i = nc.variables['I']
        i[ityme,:,:]            = self.I

        if ityme == 0:
            q = nc.createVariable('Q','f4',dim,zlib=zlib,fill_value=MISSING)
            q.standard_name = 'TOA Q' 
            q.long_name     = 'Q-component of the stokes vector at the top of the atmopshere' 
            q.missing_value = MISSING
            q.units         = "W m-2 sr-1 nm-1"
        else:
            q = nc.variables['Q']
        q[ityme,:,:]            = self.Q    

        if ityme == 0:
            u = nc.createVariable('U','f4',dim,zlib=zlib,fill_value=MISSING)
            u.standard_name = 'TOA U'
            u.long_name     = 'U-component of the stokes vector at the top of the atmopshere' 
            u.missing_value = MISSING
            u.units         = "W m-2 sr-1 nm-1"
        else:
            u = nc.variables['U']
        u[ityme,:,:]            = self.U

        if ityme == 0:
            sref = nc.createVariable('surf_reflectance','f4',dim,zlib=zlib,fill_value=MISSING)
            sref.standard_name = 'Surface Reflectance' 
            sref.long_name     = 'Bi-Directional Surface Reflectance' 
            sref.missing_value = MISSING
            sref.units         = "None"
        else:
            sref = nc.variables['surf_reflectance']
        sref[ityme,:,:]            = self.surf_reflectance

        if ityme == 0:
            sref = nc.createVariable('surf_reflectance_Q','f4',dim,zlib=zlib,fill_value=MISSING)
            sref.standard_name = 'Surface Reflectance Q' 
            sref.long_name     = 'Bi-Directional Surface Reflectance Q'
            sref.missing_value = MISSING
            sref.units         = "None"
        else:
            sref = nc.variables['surf_reflectance_Q']
        sref[ityme,:,:]            = self.BR_Q

        if ityme == 0:
            sref = nc.createVariable('surf_reflectance_U','f4',dim,zlib=zlib,fill_value=MISSING)
            sref.standard_name = 'Surface Reflectance U' 
            sref.long_name     = 'Bi-Directional Surface Reflectance U' 
            sref.missing_value = MISSING
            sref.units         = "None"
        else:
            sref = nc.variables['surf_reflectance_U']
        sref[ityme,:,:]            = self.BR_U

        if ityme == 0:
            rot = nc.createVariable('ROT','f4',('time','lev','across','ch',),zlib=zlib,fill_value=MISSING)
            rot.long_name = 'Rayleigh Optical Thickness'
            rot.missing_value = MISSING
            rot.units         = "None"
        else:
            rot = nc.variables['ROT']
        rot[ityme,:,:,:]            = self.ROT

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

    parser.add_argument("-a","--albedotype", default=albedoType,
                        help="albedo type keyword. default is to figure out according to channel")

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
    albedoType     = args.albedotype

    # figure out albedoType keyword
    if albedoType is None:
        albedoType = 'AMES_BRDF'

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
        outFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%instname',instname)

        if brdfTemplate is None:
            brdfFile = None
        else:
            brdfFile = brdfTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)



        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>inFile:    ',inFile)
        print('>>>outFile:   ',outFile)
        print('>>>rcFile:    ',rcFile)
        print('>>>albedoType:',albedoType)
        print('>>>brdfFile:  ',brdfFile)
        print('>>>verbose:   ',args.verbose)
        print('++++End of arguments+++')
        
        vlidort = ACCP_POLAR_VLIDORT(inFile,outFile,rcFile,
                            albedoType, 
                            instname,
                            args.dryrun,
                            brdfFile=brdfFile,
                            verbose=args.verbose)


        date += Dt
