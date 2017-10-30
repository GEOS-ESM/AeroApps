#!/usr/bin/env python

"""
    Calculates polarized BOA radiance at ground stations.
    Model fields have already been sampled using stn_sampler

    Adapted from polar_vlidort.py
    Patricia Castellanos, June, 2017

"""

import os
import MieObs_
from   netCDF4 import Dataset
from   mieobs  import  getAOPvector, getEdgeVars
import numpy   as np

from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from MAPL.constants import *
import stnAngles_    
import VLIDORT_STN_ 
from copyvar  import _copyVar
from scipy import interpolate
import multiprocessing

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime','stnName']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = []#['CLDTOT']
SDS_INV = ['FRLAND']

ncALIAS = {'LONGITUDE': 'stnLon',
           'LATITUDE' : 'stnLat'}

nMom     = 300

SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER'}

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_STN_.vector_brdf_modis,
                'MODIS_BRDF_BPDF': VLIDORT_STN_.vector_brdf_modis_bpdf,
                'LAMBERTIAN'     : VLIDORT_STN_.vector_lambert}   

LandAlbedos  = 'MODIS_BRDF','MODIS_BRDF_BPDF','LAMBERTIAN'


MISSING = -1.e+20


def extrap1d(interpolator):
    """ extrapolator wrapper for an interpolator"""
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike


def unwrap_self_doMie(arg, **kwarg):
    return STN_VLIDORT.doMie(*arg, **kwarg)
 

class MieVARS(object):
    """
    container for mie vars calculations
    """
    pass

class STN_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,invFile,outFile,rcFile,albedoType,
                channel,
                brdfFile=None,
                ndviFile=None,
                lcFile=None,
                lerFile=None,
                verbose=False,
                extOnly=False):
        self.SDS_AER = SDS_AER
        self.SDS_MET = SDS_MET
        self.SDS_INV = SDS_INV
        self.AERNAMES = AERNAMES
        self.inFile  = inFile
        self.invFile = invFile
        self.outFile = outFile
        self.albedoType = albedoType
        self.rcFile  = rcFile
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom
        self.brdfFile = brdfFile        
        self.lcFile = lcFile
        self.ndviFile = ndviFile
        self.lerFile  = lerFile

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
        self.nstations = len(self.LONGITUDE)
        self.ntyme    = len(self.tyme)
        self.nobs  = [self.nstations]*self.ntyme
        self.iGood = [np.ones([self.nstations]).astype(bool)]*self.ntyme

        if not extOnly:
            # Read in surface data
            # Intensity
            if (self.channel < 470) & ("MODIS_BRDF" in albedoType):
                albedoReader = getattr(self,'readHybridMODISBRDF')
            else:
                albedoReader = getattr(self,SurfaceFuncs[albedoType])
            albedoReader()

            # Polarization
            if 'BPDF' in albedoType:
                self.BPDFinputs()


            # Calculate aerosol optical properties
            self.computeMie()

            # Calculate atmospheric profile properties needed for Rayleigh calc
            self.computeAtmos()

        # Calculate Scene Geometry
        self.calcAngles()

        if any(self.nobs):
            #Land-Sea Mask
            self.LandSeaMask()

    # --
    def BPDFinputs(self):
        """
        Read in NDVI and Landuse Coefficient- should have been sampled already
        """
        if self.verbose:
            print 'opening file',self.ndviFile

        nc = Dataset(self.ndviFile)
        NDVI = np.array(nc.variables['NDVI'][:])
        missing_value = nc.variables['NDVI'].missing_value
        I = NDVI < -900
        NDVI[I] = MISSING
        I = NDVI == missing_value
        NDVI[I] = MISSING
        nc.close()

        if self.verbose:
            print 'opening file',self.lcFile
        nc = Dataset(self.lcFile)
        BPDFcoef = nc.variables['BPDFcoef'][:]
        I = BPDFcoef < -900
        BPDFcoef[I] = MISSING
        nc.close()

        for s in range(self.ntyme):
            self.iGood[s] = self.iGood[s] & (NDVI[:,s] != MISSING)
            self.nobs[s]  = np.sum(self.iGood[s])

        #BPDFparam(nparam,nch,nobs)
        self.BPDFparam = [np.zeros([3,1,self.nstations])]*self.ntyme
        for t in range(self.ntyme):
            self.BPDFparam[t][0,0,:] = 1.5
            self.BPDFparam[t][1,0,:] = NDVI[:,t]
            self.BPDFparam[t][2,0,:] = BPDFcoef[:,t]        

    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose: 
            print 'opening file',self.invFile.replace('%col',col)
        nc       = Dataset(self.invFile.replace('%col',col))

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)   

        # Make lists into arrays
        for sds in self.SDS_INV:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])            

        for t in range(self.ntyme):
            if self.albedoType in LandAlbedos:
                iGood = self.FRLAND[:,0] >= 0.99
            else:
                iGood = self.FRLAND[:,0] < 0.99

            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t]  = np.sum(self.iGood[t])

    # --
    def computeAtmos(self):
        self.pe = []
        self.ze = []
        self.te = []
        NAMES = ['AIRDENS','DELP','PS']
        for t in range(self.ntyme):
            inVars = MieVARS()
            for sds in NAMES:
                var = self.__dict__[sds]
                if len(var.shape) == 3:
                    inVars.__dict__[sds] = var[:,t,:]
                elif len(var.shape) == 2:
                    inVars.__dict__[sds] = var[:,t]

            pe, ze, te = getEdgeVars(inVars)

            self.pe.append(pe) # (km,nobs)
            self.ze.append(ze) 
            self.te.append(te)

    # --
    def calcAngles(self):
        self.SZA   = []
        self.SAA   = []
        for t in self.tyme:
            SZA = []
            SAA = []
            for s in range(self.nstations):
                CLAT = self.LATITUDE[s]
                CLON = self.LONGITUDE[s]
                
                year  = t.year 
                month = t.month 
                day   = t.day 
                hour  = t.hour 
                minute = t.minute 
                second = t.second             

                saa, sza = stnAngles_.sunangles(year,month,day,hour,minute,second,
                                                  CLAT,CLON,
                                                  0.0)
                SZA.append(sza[0])
                SAA.append(saa[0])

            self.SZA.append(np.array(SZA))
            self.SAA.append(np.array(SAA))


        # RAA is zero when pointing at the sun
        self.RAA = [np.zeros([self.nstations])]*self.ntyme

        # Limit SZAs
        for t in range(self.ntyme):
            iGood = self.SZA[t] < 80
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])     

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

        # col = 'met_Nv'
        # if self.verbose: 
        #     print 'opening file',self.inFile.replace('%col',col)        
        # nc       = Dataset(self.inFile.replace('%col',col))

        # for sds in self.SDS_MET:
        #     sds_ = sds
        #     if sds in ncALIAS:
        #         sds_ = ncALIAS[sds]
        #     var = nc.variables[sds_][:]
        #     self.__dict__[sds].append(var)

    # --- 
    def readSampledMODISBRDF(self):
        """
        Read in MODIS BRDF kernel weights
        that have already been sampled on lidar track
        """
        MODIS_channels = np.array(['470','550','650','850','1200','1600','2100'])
        chs = str(int(self.channel))
        if chs in MODIS_channels:
            SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        else:
            dch = self.channel - np.array(MODIS_channels).astype('int')
            chmin = np.argmax(dch[dch<0])
            chmin = MODIS_channels[dch<0][chmin]
            chmax = np.argmin(dch[dch>0])
            chmax = MODIS_channels[dch>0][chmax]

            SDS = 'Riso'+chmin,'Rgeo'+chmin,'Rvol'+chmin,'Riso'+chmax,'Rgeo'+chmax,'Rvol'+chmax

        if self.verbose:
            print 'opening BRDF file ',self.brdfFile
        nc = Dataset(self.brdfFile)

        for sds in SDS:
            self.__dict__[sds] = np.array(nc.variables[sds][:])

        missing_value = nc.variables[sds].missing_value
        nc.close()
        
        # Interpolate if necessary
        if chs not in MODIS_channels:
            X = np.array([chmin,chmax]).astype('int')
            for R in ['Riso','Rgeo','Rvol']:
                sds = R + chs
                self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
                for s in range(self.nstations):
                    for i in range(self.ntyme):
                        Y = np.array([self.__dict__[R+chmin][s,i],self.__dict__[R+chmax][s,i]])
                        if missing_value in Y:
                            self.__dict__[sds][s,i] = missing_value
                        else:
                            f = interpolate.interp1d(X, Y)
                            self.__dict__[sds][s,i] = f([int(chs)])

            SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        
        # Check for missing kernel weights
        Riso = self.__dict__['Riso'+chs]
        Rgeo = self.__dict__['Rgeo'+chs]
        Rvol = self.__dict__['Rvol'+chs]
        for t in range(self.ntyme):
            iGood = (Riso[:,t] != missing_value) & (Rgeo[:,t] != missing_value) & (Rvol[:,t] != missing_value)
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])

        for sds in SDS:
            self.__dict__[sds].shape = (1,1,self.nstations,self.ntyme)

        # [nkernel,nch,nobs]
        self.kernel_wt = []
        for t in range(self.ntyme):
            kernel_wt = np.append(self.__dict__['Riso'+chs][:,:,:,t],self.__dict__['Rgeo'+chs][:,:,:,t],axis=0)
            kernel_wt = np.append(kernel_wt,self.__dict__['Rvol'+chs][:,:,:,t],axis=0)
            self.kernel_wt.append(kernel_wt)

        param1 = np.array([2]*self.nstations)
        param2 = np.array([1]*self.nstations)

        param1.shape = (1,1,self.nstations)
        param2.shape = (1,1,self.nstations)

        # [nparam,nch,nobs]
        self.RTLSparam = [np.append(param1,param2,axis=0)]*self.ntyme

    # --- 
    def readHybridMODISBRDF(self):
        """
        Read in MODIS BRDF kernel weights
        that have already been sampled on lidar track
        for wavelngths between 388 and 470
        """
        MODIS_channels = np.array(['470','550','650','850','1200','1600','2100'])
        LER_channels   = np.array(['354','388'])
        MODISchs = '470' 
        mSDS = 'Riso'+MODISchs,'Rgeo'+MODISchs,'Rvol'+MODISchs
        LERchs = '388'
        lSDS = 'SRFLER' + LERchs
        chs = str(int(self.channel))

        if self.verbose:
            print 'opening BRDF file ',self.brdfFile
        nc = Dataset(self.brdfFile)

        for sds in mSDS:
            self.__dict__[sds] = nc.variables[sds][:]

        missing_value = nc.variables[sds].missing_value
        nc.close()

        if self.verbose:
            print 'opening LER albedo file ',self.lerFile
        nc = Dataset(self.lerFile)

        self.__dict__[lSDS] = np.array(np.squeeze(nc.variables[lSDS][:]))
        nc.close()

        # Interpolate 
        X = np.array([388,470])
        #Riso
        sds = 'Riso' + chs
        R   = 'Riso' + MODISchs
        self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
        for s in range(self.nstations):
            for i in range(self.ntyme):
                Y = np.array([self.__dict__[lSDS][s,i],self.__dict__[R][s,i]])
                if missing_value in Y:
                    self.__dict__[sds][s,i] = missing_value
                else:
                    f = interpolate.interp1d(X, Y)
                    self.__dict__[sds][s,i] = f([int(chs)])

        #Rgeo and Rvol
        for R in 'Rgeo','Rvol':
            sds = R + chs
            self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
            for s in range(self.nstations):
                for i in range(self.ntyme):
                    Y = np.array([0.0,self.__dict__[R+MODISchs][s,i]])
                    if missing_value in Y:
                        self.__dict__[sds][s,i] = missing_value
                    else:
                        f = interpolate.interp1d(X, Y)
                        self.__dict__[sds][s,i] = f([int(chs)])

        SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        
        # Check for missing kernel weights
        Riso = self.__dict__['Riso'+chs]
        Rgeo = self.__dict__['Rgeo'+chs]
        Rvol = self.__dict__['Rvol'+chs]
        for t in range(self.ntyme):
            iGood = (Riso[:,t] != missing_value) & (Rgeo[:,t] != missing_value) & (Rvol[:,t] != missing_value)
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])

        for sds in SDS:
            self.__dict__[sds].shape = (1,1,self.nstations,self.ntyme)

        # [nkernel,nch,nobs]
        self.kernel_wt = []
        for t in range(self.ntyme):
            kernel_wt = np.append(self.__dict__['Riso'+chs][:,:,:,t],self.__dict__['Rgeo'+chs][:,:,:,t],axis=0)
            kernel_wt = np.append(kernel_wt,self.__dict__['Rvol'+chs][:,:,:,t],axis=0)
            self.kernel_wt.append(kernel_wt)

        param1 = np.array([2]*self.nstations)
        param2 = np.array([1]*self.nstations)

        param1.shape = (1,1,self.nstations)
        param2.shape = (1,1,self.nstations)

        # [nparam,nch,nobs]
        self.RTLSparam = [np.append(param1,param2,axis=0)]*self.ntyme

    def readSampledLER(self):
        """
        Read in sampler LER files.
        Fill albedo attribute
        """
        LER_channels   = np.array(['354','388'])   
        chs = str(int(self.channel))

        SDS = []
        for lchs in LER_channels:
            SDS.append('SRFLER'+lchs)

        if self.verbose:
            print 'opening LER albedo file ',self.lerFile
        nc = Dataset(self.lerFile)

        for sds in SDS:
            self.__dict__[sds] = np.array(np.squeeze(nc.variables[sds][:]))

        missing_value = nc.variables[sds].missing_value
        nc.close()

        sds = 'SRFLER'+chs
        # interpolate if needed
        if chs not in LER_channels:
            dch   = self.channel - LER_channels.astype('int')
            try:
                chmin = np.argmax(dch[dch<0])
                chmin = LER_channels[dch<0][chmin]
            except:
                chmin = str(LER_channels.astype('int').min())
            try:
                chmax = np.argmin(dch[dch>0])
                chmax = LER_channels[dch>0][chmax]            
            except:                
                chmax = str(LER_channels.astype('int').max())

            X = np.array([chmin,chmax]).astype('int')

            self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
            for s in range(self.nstations):
                for i in range(self.ntyme):
                    Y = np.array([self.__dict__['SRFLER'+chmin][s,i],self.__dict__['SRFLER'+chmax][s,i]])
                    if missing_value in Y:                        
                        self.__dict__[sds][s,i] = missing_value
                    else:
                        try:
                            f = interpolate.interp1d(X, Y,fill_value='extrapolate')
                            self.__dict__[sds][s,i] = f([int(chs)]) 
                        except:
                            f = extrap1d(interpolate.interp1d(X, Y))
                            self.__dict__[sds][s,i] = f([int(chs)]) 


        # Check for missing values
        for t in range(self.ntyme):
            iGood = (self.__dict__[sds][:,t] != missing_value) 
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])

        # (nobs,nch)
        self.__dict__[sds].shape = (self.nstations,self.ntyme,1) 

        self.albedo = []
        for t in range(self.ntyme):
            self.albedo.append(self.__dict__[sds][:,t,:])

    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        self.tau = []
        self.ssa = []
        self.g   = []
        self.pmom = []

        pool = multiprocessing.Pool(int(multiprocessing.cpu_count()*0.5))     
        args = zip([self]*self.ntyme,range(self.ntyme))
        result = pool.map(unwrap_self_doMie,args)
        
        for r in result:
            tau, ssa, g, pmom = r
            self.tau.append(tau)  #(km,nch,nobs)
            self.ssa.append(ssa)  #(km,nch,nobs)
            self.g.append(g)    #(km,nch,nobs)
            self.pmom.append(pmom)  #(km,nch,nobs,nMom,nPol)
        

    def doMie(self,t):
        NAMES = self.AERNAMES + ['PS','DELP','RH','AIRDENS']
        inVars = MieVARS()
        for sds in NAMES:
            var = self.__dict__[sds]
            if len(var.shape)==3:
                inVars.__dict__[sds] = var[:,t,:]
            elif len(var.shape) ==2:
                inVars.__dict__[sds] = var[:,t]

        tau,ssa,g,pmom = getAOPvector(inVars,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=False,
                                 rcfile=self.rcFile,
                                 nMom=self.nMom)
        return tau,ssa,g,pmom

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
        outFile = '{}/{}.lc2.ext.{}_{}nm.nc'.format(outDir,instname,date_ch,get_chd(self.channel))
        Options =     " --input=" + self.inFile.replace('%col',col)      + \
                      " --output=" + outFile       + \
                      " --rc=" + self.rcFile      + \
                      " --format=NETCDF4_CLASSIC"      + \
                      " --channel=%d" %self.channel + \
                      " --intensive"   +\
                      " --stn"    
                      

        if not os.path.exists(os.path.dirname(outFile)):
            os.makedirs(os.path.dirname(outFile))

        cmd = './ext_sampler.py {} '.format(Options)  
        print cmd
        if os.system(cmd):
            raise ValueError, "ext_sampler.py failed for %s "%(self.inFile.replace('%col',col))       


    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """

        # index of stations with good obs
        self.nobs = np.array(self.nobs)
        tymes     = np.arange(self.ntyme)[self.nobs>0]

        self.I = []
        self.Q = []
        self.U = []
        self.reflectance = []
        self.surf_reflectance = []
        self.ROT = []
        for t in tymes:            
            sza  = self.SZA[t][self.iGood[t]]
            vza  = sza
            raa  = self.RAA[t][self.iGood[t]]

            tau = self.tau[t][:,:,self.iGood[t]]
            ssa = self.ssa[t][:,:,self.iGood[t]]
            pmom = self.pmom[t][:,:,self.iGood[t],:,:]
            pe   = self.pe[t][:,self.iGood[t]]
            ze   = self.ze[t][:,self.iGood[t]]
            te   = self.te[t][:,self.iGood[t]]

            # Initiate output arrays
            nlev    = tau.shape[0]
            self.I_ = np.ones([self.nstations])*MISSING
            self.Q_ = np.ones([self.nstations])*MISSING
            self.U_ = np.ones([self.nstations])*MISSING
            self.reflectance_ = np.ones([self.nstations])*MISSING
            self.surf_reflectance_ = np.ones([self.nstations])*MISSING
            self.ROT_ = np.ones([self.nstations,nlev])*MISSING
            
            # Get VLIDORT wrapper name from dictionary
            vlidortWrapper = WrapperFuncs[self.albedoType]

            # run VLIDORT                
            # get args list for each surface model
            if self.albedoType == 'MODIS_BRDF':
                kernel_wt = self.kernel_wt[t][:,:,self.iGood[t]]
                param     = self.RTLSparam[t][:,:,self.iGood[t]]                
                
                args = [self.channel,tau, ssa, pmom, 
                        pe, ze, te, 
                        kernel_wt, param, 
                        sza, raa, vza, 
                        MISSING,
                        self.verbose]

                # Call VLIDORT wrapper function
                I, reflectance, ROT, surf_reflectance, Q, U, rc = vlidortWrapper(*args)                        
                
            elif self.albedoType == 'MODIS_BRDF_BPDF':
                # For albedo
                kernel_wt = self.kernel_wt[t][:,:,self.iGood[t]]
                RTLSparam = self.RTLSparam[t][:,:,self.iGood[t]] 
                RTLSparam = np.append(RTLSparam,np.zeros([1,1,self.nobs[t]]),axis=0) 

                # For BPDF
                BPDFparam = self.BPDFparam[t][:,:,self.iGood[t]]

                # Loop through one by one
                # Some land covers do not have polarization (i.e. urban)
                I = np.zeros([self.nobs[t],1])
                Q = np.zeros([self.nobs[t],1])
                U = np.zeros([self.nobs[t],1])
                reflectance = np.zeros([self.nobs[t],1])
                surf_reflectance = np.zeros([self.nobs[t],1])
                ROT = np.zeros([nlev,self.nobs[t],1])

                for p in range(self.nobs[t]):
                    if BPDFparam[2,0,p] == MISSING: 
                        args = [self.channel,
                                tau[:,:,p:p+1], 
                                ssa[:,:,p:p+1], 
                                pmom[:,:,p:p+1,:,:], 
                                pe[:,p:p+1], 
                                ze[:,p:p+1], 
                                te[:,p:p+1], 
                                kernel_wt[:,:,p:p+1], 
                                RTLSparam[:,:,p:p+1], 
                                sza[p:p+1], 
                                raa[p:p+1], 
                                vza[p:p+1], 
                                MISSING,
                                self.verbose]

                        BRDFvlidortWrapper = WrapperFuncs['MODIS_BRDF']
                        # Call VLIDORT wrapper function
                        I_, reflectance_, ROT_, surf_reflectance_, Q_, U_, rc = BRDFvlidortWrapper(*args)                         

                    else:
                        args = [self.channel,
                                tau[:,:,p:p+1], 
                                ssa[:,:,p:p+1], 
                                pmom[:,:,p:p+1,:,:], 
                                pe[:,p:p+1], 
                                ze[:,p:p+1], 
                                te[:,p:p+1], 
                                kernel_wt[:,:,p:p+1], 
                                RTLSparam[:,:,p:p+1], 
                                BPDFparam[:,:,p:p+1],
                                sza[p:p+1], 
                                raa[p:p+1], 
                                vza[p:p+1], 
                                MISSING,
                                self.verbose]

                        # Call VLIDORT wrapper function
                        I_, reflectance_, ROT_, surf_reflectance_, Q_, U_, rc = vlidortWrapper(*args)
            
                    I[p:p+1,:] = I_
                    Q[p:p+1,:] = Q_
                    U[p:p+1,:] = U_
                    reflectance[p:p+1,:] = reflectance_
                    surf_reflectance[p:p+1,:] = surf_reflectance_
                    ROT[:,p:p+1,:] = ROT_
            
            elif self.albedoType == 'LAMBERTIAN':
                albedo = self.albedo[t][self.iGood[t],:]
                
                args = [self.channel,tau, ssa, pmom, 
                        pe, ze, te, 
                        albedo, 
                        sza, raa, vza, 
                        MISSING,
                        self.verbose]

                # Call VLIDORT wrapper function
                I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
                surf_reflectance = albedo

            # Store values in initialized arrays
            self.ROT_[self.iGood[t],:] = np.squeeze(ROT).T
            self.I_[self.iGood[t]] = np.squeeze(I)
            self.reflectance_[self.iGood[t]] = np.squeeze(reflectance)
            self.surf_reflectance_[self.iGood[t]] = np.squeeze(surf_reflectance)
            self.Q_[self.iGood[t]] = np.squeeze(Q)
            self.U_[self.iGood[t]] = np.squeeze(U) 

            self.I.append(self.I_)
            self.Q.append(self.Q_)
            self.U.append(self.U_)
            self.reflectance.append(self.reflectance_)
            self.surf_reflectance.append(self.surf_reflectance_)
            self.ROT.append(self.ROT_)

        self.writeNC()
        
    #---
    def writeNC (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

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
        nc.note = 'Sun normalized output.  Solar flux set equal to 1.'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        ns = nc.createDimension('station',self.nstations)
        nz = nc.createDimension('lev',km)
        ls = nc.createDimension('ls',19)

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nctrj       = Dataset(self.inFile.replace('%col',col))     
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)   
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)  
        _copyVar(nctrj,nc,u'station',dtype='f4',zlib=False,verbose=self.verbose)     
        _copyVar(nctrj,nc,u'lev',dtype='f4',zlib=False,verbose=self.verbose)       
        _copyVar(nctrj,nc,u'stnLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnName',dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)

        nctrj.close()

        sza = nc.createVariable('solar_zenith','f4',('station','time',),zlib=zlib)
        sza.long_name     = "solar zenith angle (SZA)"
        sza.missing_value = MISSING
        sza.units         = "degrees"
        sza[:]            = np.array(self.SZA).T

        saa = nc.createVariable('solar_azimuth','f4',('station','time',),zlib=zlib)
        saa.long_name     = "solar azimuth angle (SAA)"
        saa.missing_value = MISSING
        saa.units         = "degrees clockwise from North"
        saa[:]            = np.array(self.SAA).T


        # Write VLIDORT Outputs
        # ---------------------
        ref = nc.createVariable('I_normalized','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm BOA normalized I' %self.channel
        ref.long_name     = '%.2f nm normalized intensity at the bottom of the atmosphere (pi*I/cosSZA*SOLAR_FLUX)' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        ref[:]            = np.array(self.reflectance).T

        i = nc.createVariable('I','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm BOA I' %self.channel
        i.long_name     = '%.2f nm intensity at the bottom of the atmosphere' %self.channel
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        i[:]            = np.array(self.I).T

        q = nc.createVariable('Q','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm BOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the bottom of the atmopshere' %self.channel
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = np.array(self.Q).T

        u = nc.createVariable('U','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm BOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the bottom of the atmopshere' %self.channel
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = np.array(self.U).T

        sref = nc.createVariable('surf_reflectance','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = np.array(self.surf_reflectance).T

        rot = nc.createVariable('ROT','f4',('station','time','lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        b = np.empty([self.nstations,self.ntyme,km])
        for t in range(self.ntyme):
            b[:,t,:] = self.ROT[t]
        rot[:]            = b

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFile)
    

def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,01,01,00)
    nymd     = str(date.date()).replace('-','')
    hour     = str(date.hour).zfill(2)
    format   = 'NETCDF4_CLASSIC'

    inDir        = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    inFile       = '{}/aeronet-g5nr.lb2.%col.{}_{}z.nc4'.format(inDir,nymd,hour)
    invDir       = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/invariant'
    invFile      = '{}/aeronet-g5nr.lb2.%col.nc4'.format(invDir)
    brdfDir      = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/surface/BRDF/MCD43C1/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    brdfFile     = '{}/aeronet-g5nr.lb2.brdf.{}_{}z.nc4'.format(brdfDir,nymd,hour)
    ndviDir      = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/surface/BPDF/NDVI/MYD13C2/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    ndviFile     = '{}/aeronet-g5nr.lb2.ndvi.{}_{}z.nc4'.format(ndviDir,nymd,hour)
    lcDir        = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/surface/BPDF/LAND_COVER/MCD12C1/051/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    lcFile       = '{}/aeronet-g5nr.lb2.land_cover.{}_{}z.nc4'.format(lcDir,nymd,hour)    
    albedoType   = 'MODIS_BRDF_BPDF'

    channel  = 470
    chd      = get_chd(channel)
    outDir    = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelC/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    outFile   = '{}/aeronet-g5nr.vlidort.vector.MCD43C.{}_{}z_{}nm.nc4'.format(outDir,nymd,hour,chd)
    
    rcFile   = 'Aod_EOS.rc'
    verbose  = True
    extOnly  = False

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = STN_VLIDORT(inFile,invFile,outFile,rcFile,
                            albedoType,
                            channel,
                            brdfFile=brdfFile,
                            ndviFile=ndviFile,
                            lcFile=lcFile,
                            verbose=verbose,
                            extOnly=extOnly)

   
    # Run ext_sampler
    vlidort.runExt()

    # Run VLIDORT
    if any(vlidort.nobs):
        vlidort.runVLIDORT()