#!/usr/bin/env python

"""
    Utility to compute optical properties for a PACE granule

    Adapted from ext_sampler
    P. Castellanos April 2019

"""

import os
import MieObs_
from netCDF4           import Dataset
from mieobs            import getAOPext, getAOPscalar, getAOPint, getAOPvector

from datetime          import datetime, timedelta
from dateutil.parser   import parse         as isoparser

from MAPL              import Config, eta
from MAPL.constants    import *
import argparse
from   pace            import granules
import numpy           as np

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']

META    = ['DELP','PS','RH','AIRDENS']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU

MieVars = {'ext'    : ['total aerosol extinction','km-1'],
           'scatext': ['scattering extinction','km-1'],
           'depol'  : ['depolarization ratio','unitless']}

IntVars = {'refi' : ['imaginary refractive index','unitless'],
           'refr' : ['real refractive index','unitless'],
           'reff' : ['aerosol effective raidus','m']}


class EXTCALC(object):
    """
    Deals with everything to do with aerosol extinction and intensive variables sampling
    """
    def __init__(self,args,AERNAMES=AERNAMES,MieVars=MieVars,IntVars=IntVars):
        self.AERNAMES  = AERNAMES
        self.MieVars   = MieVars
        self.IntVars   = IntVars
        self.SDS       = AERNAMES + META

        self.rootdir = args.rootdir
        self.Date    = isoparser(args.iso_t1)
        self.rcDir   = args.rcDir
        self.verbose = args.verbose
        self.iscan   = int(args.iscan)
        self.iccd    = int(args.iccd)

        # get PACE channels
        self.get_channels()        

        # Read in what you need from aerFile
        self.getVars()

        # Call Mie Calculator
        self.computeMie()
    #---
    def get_channels(self):
        pdate = isoparser(self.Date.strftime('2020-%m-%dT%H:%M:00'))
        inFile = granules(self.rootdir+'/L1B',pdate,pdate)
        if len(inFile) == 0:
            raise Exception('No PACE Infiles found, nothing to do, exiting.')
        
        nc = Dataset(inFile[0])
        group = 'sensor_band_parameters'
        SDS = ['blue_wavelength','red_wavelength','SWIR_wavelength']
        for sds in SDS:
            self.__dict__[sds] = np.array(nc.groups[group].variables[sds][:])

    #---
    def getVars(self):
        """
        Parse input file, create variable dictionary
        """
        YMDdir    = self.Date.strftime('Y%Y/M%m/D%d')
        nymd = self.Date.strftime('%Y%m%d')
        hms  = self.Date.strftime('%H%M00')

        LbDir     = '{}/LevelB/{}'.format(self.rootdir,YMDdir)
        AER_file  = '{}/pace-g5nr.lb.aer_Nv.{}_{}.nc4'.format(LbDir,nymd,hms)

        nc    = Dataset(AER_file)
        nccd  = len(nc.dimensions['ccd_pixels'])
        nscan = len(nc.dimensions['number_of_scans'])
        nlev  = len(nc.dimensions['lev'])

        for sds in self.SDS:
            var = np.squeeze(nc.variables[sds][:])
            rank = len(var.shape)
            if rank == 2:
                # [nscan,nccd]
                #var = var.reshape(nscan*nccd)
                var = var[self.iscan:self.iscan+1,self.iccd]
            else:
                # [nlev,nscan,nccd]
                #var = var.reshape([nlev,nscan*nccd])
                var = var[:,self.iscan:self.iscan+1,self.iccd]

            self.__dict__[sds] = var
        
        nc.close()
        self.nccd  = nccd
        self.nscan = nscan
        self.nlev  = nlev
        self.nobs  = nccd*nscan

    #---
    def computeMie(self):
        """
        Computes optical quantities 
        """

        for band in ['blue','red','SWIR']:
            channel = self.__dict__[band + '_wavelength']
            tau_ = []
            ssa_ = []
            pmom_ = []

            for ch in channel:
                rcFile = self.rcDir + '/Aod_EOS_{}.rc'.format(int(ch))
                tau,ssa,g,pmom = getAOPvector(self,[ch],vnames=self.AERNAMES,Verbose=self.verbose,rcfile=rcFile,nMom=300)
                tau_.append(tau)
                ssa_.append(ssa)
                pmom_.append(pmom)

            self.__dict__['tau_'+band] = np.concatenate(tau_,axis=1)
            self.__dict__['ssa_'+band] = np.concatenate(ssa_,axis=1)
            self.__dict__['pmom_'+band] = np.concatenate(pmom_,axis=1)

#     ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        #     ext2back = ones(backscat.shape)*MAPL_UNDEF
        #     I = backscat > 0
        #     ext2back[I] = ext[I]/backscat[I]
        #     MieVars = {"ext":[ext],"scatext":[sca],"backscat":[backscat],"aback_sfc":[aback_sfc],"aback_toa":[aback_toa],"depol":[depol],"ext2back":[ext2back],"tau":[tau],"ssa":[ssa],"g":[g]}

        #     if options.intensive:
        #         vol, area, refr, refi, reff = getAOPint(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        #         MieVars['vol']  = [vol]
        #         MieVars['area'] = [area]
        #         MieVars['refr'] = [refr]
        #         MieVars['refi'] = [refi]
        #         MieVars['reff'] = [reff]
        #         else:
        #             tau,ssa,g = getAOPscalar(VarsIn,channel,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        #             ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        #             ext2back = ones(backscat.shape)*MAPL_UNDEF
        #             I = backscat > 0
        #             ext2back[I] = ext[I]/backscat[I]
        #             MieVars['ext'].append(ext)
        #             MieVars['scatext'].append(sca)
        #             MieVars['backscat'].append(backscat)
        #             MieVars['aback_sfc'].append(aback_sfc)
        #             MieVars['aback_toa'].append(aback_toa)
        #             MieVars['depol'].append(depol) 
        #             MieVars['ext2back'].append(ext2back)
        #             MieVars['tau'].append(tau)
        #             MieVars['ssa'].append(ssa)
        #             MieVars['g'].append(g)

        #             if options.intensive:
        #                 vol, area, refr, refi, reff = getAOPint(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)
        #                 MieVars['vol'].append(vol)
        #                 MieVars['area'].append(area)
        #                 MieVars['refr'].append(refr)
        #                 MieVars['refi'].append(refi)
        #                 MieVars['reff'].append(reff)                    


    #---
    def writeNC ( self, outFile, zlib=False):
        """
        Write a NetCDF file with sampled GEOS-5 variables along the satellite track
        described by (lon,lat,tyme).
        """
        km = 72

        # Open NC file
        # ------------
        nc = Dataset(outFile,'w')

        # Set global attributes
        # ---------------------
        nc.title = 'GEOS-5 Sampled Aerosol Optical Properties File'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from sampled GEOS-5 collections'
        nc.references = 'n/a'
        nc.comment = 'This file contains sampled GEOS-5 aerosol optical properties.'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',1)
        nz = nc.createDimension('lev',km)
        nb = nc.createDimension('blue_wavelength',len(self.blue_wavelength))
        nr = nc.createDimension('red_wavelength',len(self.red_wavelength))
        ns = nc.createDimension('SWIR_wavelength',len(self.SWIR_wavelength))
        nm = nc.createDimension('nMoments',300)
        nP = nc.createDimension('nPol',6)
        # Coordinate variables
        # --------------------
        if km > 0: # pressure level not supported yet
            lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
            lev.long_name = 'Vertical Level'
            lev.units = 'km'
            lev.positive = 'down'
            lev.axis = 'z'
            lev[:] = range(1,km+1)

        
        # Write each variable
        # --------------------------------------------------
        for band in ['blue','red','SWIR']:
            for n, name in enumerate(['tau','ssa','pmom']):
                print name
                var = self.__dict__[name + '_' + band]
                if name == 'pmom':
                    dim = ('lev',band+'_wavelength','time','nMoments','nPol')
                else:
                    dim = ('lev',band+'_wavelength','time') 

                this = nc.createVariable(name+'_'+band,'f4',dim,zlib=zlib)
                this.standard_name = name
                this.units = 'None'
                this.missing_value = np.float32(MAPL_UNDEF)
                this[:] = var

        # Close the file
        # --------------
        nc.close()

        print " <> wrote file %s"%(outFile)
    
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format   = 'NETCDF4_CLASSIC'
    rootdir  = '/nobackup/PACE/'
    outFile  = 'pace_ext.nc4'
    rcDir    = './rc/'
    iscan    = "600"
    iccd     = "900"

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')

    parser.add_argument("--outFile",default=outFile,
                           help="out file name (default=%s)."%outFile)

    parser.add_argument("--iscan",default=iscan,
                           help="scan index (default=%s)."%iscan)

    parser.add_argument("--iccd",default=iccd,
                                       help="ccd index (default=%s)."%iccd)

    parser.add_argument("--rootdir",default=rootdir,
               help="root directory for PACE data (default=%s)."%rootdir)      

    parser.add_argument("-r", "--rcDir", default=rcDir,
              help="Directory for Resource files pointing to optical tables")

    parser.add_argument("--no_intensive",action="store_true",
                      help="do not return intensive variables")

    parser.add_argument("-v", "--verbose",
                      action="store_true", 
                      help="Verbose mode")

    parser.add_argument("--du",
                      action="store_true", 
                      help="Dust Only")

    parser.add_argument("--ss",
                      action="store_true",
                      help="Seasalt Only")

    parser.add_argument("--su",
                      action="store_true",
                      help="Sulfate Only")

    parser.add_argument("--bc",
                      action="store_true",
                      help="Black Carbon Only")

    parser.add_argument("--oc",
                      action="store_true",
                      help="Organic Carbon Only")

    args = parser.parse_args()
        
    # Get Variables
    # --------------------------
    ext = EXTCALC(args)
    ext.writeNC(args.outFile)
