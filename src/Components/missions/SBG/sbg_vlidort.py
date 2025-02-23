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
                instname,ich,dryrun,
                nstreams=12,
                plane_parallel=True,
                brdfFile=None,
                verbose=False,
                nproc=125):
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
        self.instname    = instname
        self.ich         = ich
        self.nproc       = nproc

        # get granule dimensions
        self.getDims()

        # Start out with all good obs
        self.iGood = np.ones([self.nobs]).astype(bool)

        # Read in precalculated Scene Geometry
        # limit iGood to sza < 80
        self.readAngles()

        # Land-Sea Mask
        # limit iGood to land pixels
        self.LandSeaMask()

        if self.nobs > 0:

            # Read in model data
            self.readSampledGEOS()

            # Read in surface data
            self.readSampledAMESBRDF()

            # Set up AOP rcFile
            self.getrcFile()

            if not dryrun:
                # Run VLIDORT
                print('++++ Running VLIDORT ++++')
                self.runVLIDORT()
                self.writeNC() 

    #---
    def getrcFile(self):
        """
        Update the rcFile if needed
        """

        # Copy over rc and edit
        source = open('Aod_EOS.rc','r')
        destination = open(self.rcFile,'w')
        a = self.channels[self.ich]*1e-3
        for line in source:
            if (line[0:11] == 'r_channels:'):
                destination.write('r_channels: '+'{:0.10f}e-6'.format(a)+'\n')
            else:
                destination.write(line)
        source.close()
        destination.close()

    #---
    def getEdgeVars(self,iobs):
        """
        Get altitude, pressure, and temperature and edge of layers
        """

        # Get layer thicnkness from DELP & AIRDENS
        # DELP: Pa = kg m-1 s-2
        # AIRDENS: kg m-3
        # GRAV: m s-2
        # -----------------------------------------
        DELP = self.aer['DELP'].isel(nobs=iobs).values
        AIRDENS = self.aer['AIRDENS'].isel(nobs=iobs).values
        rhodz = DELP / GRAV
        dz = rhodz / AIRDENS       # column thickness in m

        # add up the thicknesses to get edge level altitudes and pressures
        npts, nlev = dz.shape
        ze = np.array([dz[:,i:].sum(axis=1) for i in range(nlev)])
        # append surface level, altitude = 0
        # [nlev+1,nacross]
        self.ze = np.append(ze,np.zeros([1,npts]),axis=0)

        ptop = 1. # Pa
        pe = np.zeros([self.nlev+1,npts])
        pe[0,:] = ptop
        for ilev in range(self.nlev):
            pe[ilev+1,:] = pe[ilev,:] + DELP[:,ilev]

        self.pe = pe

        # get the mid-level pressures and temperatures
        self.pm = (self.pe[:-1,:] + self.pe[1:,:])*0.5
        self.tm = self.pm/(AIRDENS.T*RGAS)

        # get the edge level temperature
        self.te = np.zeros([self.nlev+1,npts])
        self.te[0,:] = self.tm[0,:]  #isothermal a highest level
        for ilev in range(1,self.nlev):
            alpha = np.log(self.pe[ilev,:]/self.pm[ilev-1,:])/np.log(self.pm[ilev,:]/self.pm[ilev-1,:])
            self.te[ilev,:] = self.tm[ilev-1] + alpha*(self.tm[ilev,:] - self.tm[ilev-1,:])

        # dry adiabatic
        self.te[self.nlev,:] = self.tm[self.nlev-1,:]*(self.pe[self.nlev,:]/self.pm[self.nlev-1,:])**KAPPA

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
        self.nobs = self.nacross*self.ntyme
        self.nbatch = self.nproc
        self.nMom = 300 

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        Use xarray to lazy load
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

        self.aer = xr.open_mfdataset(inList,chunks="auto")
        # make arrays [nobs,nlev]
        self.aer = self.aer.squeeze()
        self.aer = self.aer.stack(nobs=("time","ncross"))
        self.aer = self.aer.transpose("nobs","lev")


    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col),chunks="auto")

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = ds[sds_].squeeze().stack(nobs=("time","ncross"))
            self.__dict__[sds] = var

        iGood = self.FRLAND >= 0.99
        self.iLand = iGood.values
        iGood = self.FRLAND < 0.99
        self.iSea  = iGood.values

        # self.iGood = self.iGood & iGood
        self.nobsLand  = np.sum(self.iGood & self.iLand)
        self.nobsSea   = np.sum(self.iGood & self.iSea)

        self.iGood = self.iGood & self.iLand
        self.nobs = np.sum(self.iGood)


    #---
    def readSampledAMESBRDF(self):
        """
        Read in AMES BRDF kernel weights
        that have already been sampled on swath
        """
        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
        ds = xr.open_dataset(self.brdfFile,chunks="auto")

        self.channels = ds.nwav.values  # microns
        self.nch = len(self.channels)
        self.channels = self.channels*1e3  # nm

        # concatenate kernel weights into one data array
        da = xr.concat([ds.Ki,ds.Kg,ds.Kv],dim="nalong")
        da = da.rename("kernel_wt")

        # reshape to [nkernel,nch,nobs]
        self.kernel_wt = da.squeeze().stack(nobs=("time","ncross")).transpose("nalong","nwav","nobs")
        self.kernel_wt = self.kernel_wt.isel(nwav=slice(self.ich,self.ich+1))

        # make an array of the params, just big enough for a batch
        param1 = np.array([2]*self.nbatch)
        param2 = np.array([1]*self.nbatch)

        #[nparam,nch,nobs]
        param1.shape = (1,1,self.nbatch)
        param2.shape = (1,1,self.nbatch)

        # [nparam,nch,nobs]
        self.RTLSparam = np.append(param1,param2,axis=0)
        self.RTLSparam = self.RTLSparam.reshape(2,1,self.nbatch)

        # filter for nans
        iGood = ~np.isnan(self.kernel_wt[0,0,:])
        self.iGood = self.iGood & iGood.values
        self.nobs = np.sum(self.iGood)

    def readAngles(self):
        """
        Read in viewing and solar Geometry from angFile
        """

        col = self.instname
        if self.verbose: 
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col),chunks="auto")

        for sds in self.SDS_ANG:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = ds[sds_].squeeze().stack(nobs=("time","ncross"))
            self.__dict__[sds] = var

        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I.compute()] = saa[I.compute()] - 360.

        RAA = self.VAA - saa
        I = RAA < 0
        RAA[I.compute()] = RAA[I.compute()]+360.0
        self.RAA = RAA

        # Limit SZAs
        iGood = self.SZA < 80
        self.iGood = self.iGood & iGood.values
        self.nobs = np.sum(self.iGood)         

    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """

        # Initiate output arrays
        iGood = np.arange(len(self.iGood))[self.iGood]
        nobs   = self.nobs
        nlev   = self.nlev
        nch    = 1
        self.I = np.ones([nobs,nch])*MISSING
        self.Q = np.ones([nobs,nch])*MISSING
        self.U = np.ones([nobs,nch])*MISSING
        self.reflectance = np.ones([nobs,nch])*MISSING
        self.surf_reflectance = np.ones([nobs,nch])*MISSING
        self.BR_Q = np.ones([nobs,nch])*MISSING
        self.BR_U = np.ones([nobs,nch])*MISSING
        self.ROT  = np.ones([nlev,nobs,nch])*MISSING

        p = Pool(self.nbatch)
        # loop through nobs in batches
        for sob in range(0,self.nobs,self.nbatch):
            print('sob, nobs',sob, self.nobs)
            eob = min([self.nobs,sob + self.nbatch])
            iobs = iGood[sob:eob]
            self.channel = [self.channels[self.ich]]

            # Calculate ROT
            self.getEdgeVars(iobs)
            args = [self.channel, self.pe.astype('float64'), self.ze.astype('float64'), self.te.astype('float64'), MISSING, self.verbose]
            vlidortWrapper = WrapperFuncs['ROT_CALC']
            # [lev,nobs,ch], [ch]
            rot, depol_ratio, rc = vlidortWrapper(*args) 
            pe   = self.pe.astype('float64')
            ze   = self.ze.astype('float64')
            te   = self.te.astype('float64') 

            # calculate AOPs
            for sds in self.AERNAMES + ['PS','DELP','RH']:
                self.__dict__[sds] = self.aer[sds].isel(nobs=iobs)

            self.computeMie()
            tau = self.tau.astype('float64')
            ssa = self.ssa.astype('float64')
            pmom = self.pmom.astype('float64')


            # Get surface data
            npts = eob - sob
            param     = self.RTLSparam[:,:,0:npts]
            kernel_wt = self.kernel_wt.isel(nobs=iobs).astype('float64')

            # get angles
            vza = self.VZA.isel(nobs=iobs).astype('float64')
            sza = self.SZA.isel(nobs=iobs).astype('float64')
            raa = self.RAA.isel(nobs=iobs).astype('float64')
            
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
                    self.verbose) for i in range(npts)]
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

            self.I[sob:eob,:] = I
            self.reflectance[sob:eob,:] = reflectance
            self.surf_reflectance[sob:eob,:] = surf_reflectance
            self.Q[sob:eob,:] = Q
            self.U[sob:eob,:] = U 
            self.BR_Q[sob:eob,:] = BR_Q
            self.BR_U[sob:eob,:] = BR_U
            self.ROT[:,sob:eob,:] = rot
            self.depol_ratio = depol_ratio
        p.close()
    #---
    def writeNC (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """

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
        nch = 1
        nt = nc.createDimension('time',self.ntyme)
        nz = nc.createDimension('lev',self.nlev)
        nh = nc.createDimension('ch',nch)
        nx = nc.createDimension('across',self.nacross)

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print('opening file',self.inFile.replace('%col',col))
        nctrj       = Dataset(self.inFile.replace('%col',col))        
        _copyVar(nctrj,nc,'trjLon',dtype='f4',zlib=zlib,verbose=self.verbose)
        _copyVar(nctrj,nc,'trjLat',dtype='f4',zlib=zlib,verbose=self.verbose)
        _copyVar(nctrj,nc,'time', dtype='i4',zlib=zlib,verbose=self.verbose)
        _copyVar(nctrj,nc,'lev', dtype='f4',zlib=zlib,verbose=self.verbose)
        nctrj.close()


        # Write VLIDORT Outputs
        # ---------------------
        dim = ('time','across','ch')
        ref = nc.createVariable('toa_reflectance','f4',dim,zlib=zlib,fill_value=MISSING)
        ref.standard_name = 'TOA Reflectance' 
        ref.long_name     = 'reflectance at the top of the atmosphere' 
        ref.missing_value = MISSING
        ref.units         = "None"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.reflectance
        ref[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        i = nc.createVariable('I','f4',dim,zlib=zlib,fill_value=MISSING)
        i.standard_name = 'TOA I' 
        i.long_name     = 'intensity at the top of the atmosphere' 
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.I
        i[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        q = nc.createVariable('Q','f4',dim,zlib=zlib,fill_value=MISSING)
        q.standard_name = 'TOA Q' 
        q.long_name     = 'Q-component of the stokes vector at the top of the atmopshere' 
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.Q
        q[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        u = nc.createVariable('U','f4',dim,zlib=zlib,fill_value=MISSING)
        u.standard_name = 'TOA U'
        u.long_name     = 'U-component of the stokes vector at the top of the atmopshere' 
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.U
        u[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        sref = nc.createVariable('surf_reflectance','f4',dim,zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'Surface Reflectance' 
        sref.long_name     = 'Bi-Directional Surface Reflectance' 
        sref.missing_value = MISSING
        sref.units         = "None"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.surf_reflectance
        sref[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        sref = nc.createVariable('surf_reflectance_Q','f4',dim,zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'Surface Reflectance Q' 
        sref.long_name     = 'Bi-Directional Surface Reflectance Q'
        sref.missing_value = MISSING
        sref.units         = "None"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.BR_Q
        sref[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        sref = nc.createVariable('surf_reflectance_U','f4',dim,zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'Surface Reflectance U' 
        sref.long_name     = 'Bi-Directional Surface Reflectance U' 
        sref.missing_value = MISSING
        sref.units         = "None"
        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        temp[self.iGood,:] = self.BR_U
        sref[:]            = temp.reshape([self.ntyme,self.nacross,nch])

        rot = nc.createVariable('ROT','f4',('lev','time','across','ch',),zlib=zlib,fill_value=MISSING)
        rot.long_name = 'Rayleigh Optical Thickness'
        rot.missing_value = MISSING
        rot.units         = "None"
        temp = np.ones([self.nlev,self.ntyme*self.nacross,nch])*MISSING
        temp[:,self.iGood,:] = self.ROT
        rot[:]            = temp.reshape([self.nlev,self.ntyme,self.nacross,nch])

        d = nc.createVariable('depol_ratio','f4',('ch',),zlib=zlib,fill_value=MISSING)
        d.long_name = 'depolarization ratio'
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
    DT_mins   = 1
    nproc     = 125
    rcFile     = 'rc/Aod_EOS_%ich.rc'
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

    parser.add_argument("-a","--albedotype", default=albedoType,
                        help="albedo type keyword. default is to figure out according to channel")

    parser.add_argument("--rcFile",default=rcFile,
                        help="rcFile (default=%s)"%rcFile)

    parser.add_argument("-D","--DT_mins", default=DT_mins, type=int,
                        help="Timestep in minutes for each file (default=%i)"%DT_mins)


    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    parser.add_argument("-n", "--nproc",default=nproc,
                        help="Number of processors (default=%i)"%nproc)

    args = parser.parse_args()
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
    Dt        = timedelta(minutes=args.DT_mins)

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


        rcFile = args.rcFile.replace('%ich',str(args.ich).zfill(3))

        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>channel:   ',args.ich)
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
                            args.ich,
                            args.dryrun,
                            brdfFile=brdfFile,
                            verbose=args.verbose,
                            nproc=args.nproc)


        date += Dt
