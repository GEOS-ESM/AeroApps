#!/usr/bin/env python3

"""
    script to compare performance of different AOP calculators
"""
import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL.config     import Config
import numpy   as np
import xarray  as xr
from pyobs import mietable as mt
from pyobs.aop import G2GAOP
from py_leo_vlidort.vlidort import VLIDORT
import yaml
from pyobs.constants import MAPL_GRAV as GRAV
from pyobs.constants import MAPL_RDRY as RGAS
from pyobs.constants import MAPL_KAPPA as KAPPA
import time
import matplotlib.pyplot as plt
from scipy.special.orthogonal import legendre
from math import sqrt, cos, factorial
import eval_gsfun

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
           'LATITUDE' : 'latitude'}

MISSING = np.float32(-1.e+20)


class OPTICS_VLIDORT(VLIDORT,G2GAOP):
    """
    Calculate AOPS
    GEOS-5 has already been sampled on satellite track
    """
    def __init__(self,inFile,channels,miercFile,pyobsrcFile,
                instname,dryrun,
                verbose=False):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.AERNAMES    = AERNAMES
        self.inFile      = inFile
        self.miercFile   = miercFile
        self.pyobsrcFile = pyobsrcFile
        self.verbose     = verbose
        self.instname    = instname
        self.nch = len(channels)
        self.channels = channels


        # Read in model data
        self.readSampledGEOS()

        # calc AOP - MieMod
        ts = time.time()
        self.nMom = 300
        self.mietau = []
        self.miessa = []
        self.miepmom = []
        for sds in self.AERNAMES + ['PS','DELP','RH']:
            self.__dict__[sds] = self.aer[sds]
        for ich in range(self.nch): 
            self.getmiercFile(ich)
            self.getmieAOP(ich)
        te = time.time()
    
        print('MieMod %f seconds'%(te-ts))

        # calc AOP - pyobs
        ts = time.time()
        self.pyobstau = []
        self.pyobsssa = []
        self.pyobspmom = []
        self.getpyobsMT()
        for ich in range(self.nch):
            self.getpyobsAOP(ich)
        te = time.time()
        print('pyobs %f seconds'%(te-ts))

    #---
    def getpyobsAOP(self,ich):
        """
        use pyobs utilities to get AOP
        """
        aop = self.getAOPrt(wavelength=self.channels[ich],vector=True,m=self.nMom)

        # need to reshape these to [nlev,nobs]
        self.pyobstau.append(aop.AOT.astype('float64').transpose().to_numpy())
        self.pyobsssa.append(aop.SSA.astype('float64').transpose().to_numpy())
        self.pyobspmom.append(aop.PMOM.astype('float64').transpose('lev','nobs','m','p').to_numpy())


    #---
    def getpyobsMT(self):
        """
        Use pyobs utilities to load mietables
        """
        self.mieTable = yaml.safe_load(open(self.pyobsrcFile))
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

    #---
    def getmiercFile(self,ich):
        """
        Update the rcFile if needed
        """
        rcDir = os.path.dirname(self.miercFile)
        if not os.path.exists(rcDir):
            os.makedirs(rcDir)

        # Copy over rc and edit
        rcFile = self.miercFile.replace('ich',str(ich).zfill(3))
        source = open('Aod_EOS.rc','r')
        destination = open(rcFile,'w')
        a = self.channels[ich]*1e-3
        for line in source:
            if (line[0:11] == 'r_channels:'):
                destination.write('r_channels: '+'{:0.10f}e-6'.format(a)+'\n')
            else:
                destination.write(line)
        source.close()
        destination.close()


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

        # just do one pixel 
        self.aer = self.aer.isel(nobs=[0])

    def getmieAOP(self,ich):
        """
        gets AOPs using MieMod 
        """
        self.rcFile = self.miercFile.replace('ich',str(ich).zfill(3))
        self.channel = [self.channels[ich]]

        # calculate AOPs
        self.computeMie()
        self.mietau.append(self.tau.astype('float64'))  #(km,nch,nobs)
        self.miessa.append(self.ssa.astype('float64'))  #(km,nch,nobs)
        self.miepmom.append(self.pmom.astype('float64')) #(km,nch,nobs,nMom,nPol)


    #---
    def comparePMOM(self):
        """
        get back the full phase function from the moments
        """

        # compare the coefficients directly
        do_direct = False
        if do_direct:
            for i in np.arange(self.p):
                x_data,y_data = self.pyobspmom[0].squeeze()[0,:,i],self.miepmom[0].squeeze()[0,:,i]
                fig, ax = plt.subplots()
                ax.scatter(x_data, y_data)

                # Add the 1:1 line
                ax.axline((0, 0), slope=1, color='r', linestyle='--', label='1:1 line')

                # Set axis limits to make sure the 1:1 line spans the plot
                min_val = min(min(x_data), min(y_data))
                max_val = max(max(x_data), max(y_data))
                ax.set_xlim(min_val, max_val)
                ax.set_ylim(min_val, max_val)
                ax.set_title(i)

                plt.show()

        # Define Angles
        angle  = np.linspace(0,180,int(180/0.1)+1)
        theta  = np.radians(angle)
        mu    = np.cos(theta)
        ntheta = len(theta)


        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))
        # Reconstruct P11
        # ordered over nPol as P11, P12, P33, P34, P22, P44
        do_p11 = True
        if do_p11:
            ipol = 0
            p11pyobs = np.zeros(ntheta)
            p11mie   = np.zeros(ntheta)
            
            miepmom = self.miepmom[0].squeeze()[0,:,ipol]
            pyobspmom = self.pyobspmom[0].squeeze()[0,:,ipol]
            for s in range(self.nMom):
                P = legendre(s)
                p11pyobs += pyobspmom[s]* P(mu)
                p11mie += miepmom[s]* P(mu)

            ax = axes[0,0]
            ax.plot(angle, p11pyobs,label='pyobs') 
            ax.plot(angle,p11mie,label='mieobs')
            ax.legend()
            ax.set_title('P11=P1')

        # Reconstruct P12
        # ordered over nPol as P11, P12, P33, P34, P22, P44
        do_P12 = True
        if do_P12:
            ipol = 1
            m = 0
            n = 2
            gsf_coef = self.miepmom[0].squeeze()[0,:,ipol]
            leg_coef = self.pyobspmom[0].squeeze()[0,:,ipol]

            p12 = np.zeros(ntheta)
            for i,t in enumerate(theta):
                d0, dneg1 = eval_gsfun.d_mn(m,n,0,t)
                pfunc = d0/(1j**(n-m)).real
                p12[i] = gsf_coef[0]*pfunc
                d1, d0 = eval_gsfun.d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
                pfunc = d1/(1j**(n-m)).real
                p12[i] = p12[i] + gsf_coef[1]*pfunc

                dm2 = d0
                dm1 = d1
                for s in range(2,self.nMom):
                    dfunc, dm2 = eval_gsfun.d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                    pfunc = dfunc/(1j**(n-m)).real

                    p12[i] = p12[i] + gsf_coef[s]*pfunc

                    dm1 = dfunc
            p12mie = p12

            p12pyobs = np.zeros(ntheta)
            for s in range(self.nMom):
                P = legendre(s)
                p12pyobs += leg_coef[s]* P(mu)

            ax = axes[0,1] 
            ax.plot(angle,p12pyobs/p11pyobs,label='pyobs')
            ax.plot(angle,p12mie/p11mie,label='mie')
            ax.legend()
            ax.set_title('P12 = P2')

        # Reconstruct P33 & P22
        # ordered over nPol as P11, P12, P33, P34, P22, P44
        do_P33_P22 = True
        if do_P33_P22:
            # first get a2 + a3
            m = 2
            n = 2
            ipol = 4
            gsf_coef22 = self.miepmom[0].squeeze()[0,:,ipol]
            leg_coef22 = self.pyobspmom[0].squeeze()[0,:,ipol]

            ipol = 2
            gsf_coef33 = self.miepmom[0].squeeze()[0,:,ipol]
            leg_coef33 = self.pyobspmom[0].squeeze()[0,:,ipol]

            a2p3 = np.zeros(ntheta)
            for i,t in enumerate(theta):
                d0, dneg1 = eval_gsfun.d_mn(m,n,0,t)
                pfunc = d0/(1j**(n-m)).real
                a2p3[i] = (gsf_coef22[0]+gsf_coef33[0])*pfunc
                d1, d0 = eval_gsfun.d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
                pfunc = d1/(1j**(n-m)).real
                a2p3[i] = a2p3[i] + (gsf_coef22[1]+gsf_coef33[1])*pfunc

                dm2 = d0
                dm1 = d1
                for s in range(2,self.nMom):
                    dfunc, dm2 = eval_gsfun.d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                    pfunc = dfunc/(1j**(n-m)).real

                    a2p3[i] = a2p3[i] + (gsf_coef22[s]+gsf_coef33[s])*pfunc

                    dm1 = dfunc

            # next get a2 - a3
            m = 2
            n = -2
            a2m3 = np.zeros(ntheta)
            for i,t in enumerate(theta):
                d0, dneg1 = eval_gsfun.d_mn(m,n,0,t)
                pfunc = d0/(1j**(n-m)).real
                a2m3[i] = (gsf_coef22[0]-gsf_coef33[0])*pfunc
                d1, d0 = eval_gsfun.d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
                pfunc = d1/(1j**(n-m)).real
                a2m3[i] = a2m3[i] + (gsf_coef22[1]-gsf_coef33[1])*pfunc

                dm2 = d0
                dm1 = d1
                for s in range(2,self.nMom):
                    dfunc, dm2 = eval_gsfun.d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                    pfunc = dfunc/(1j**(n-m)).real

                    a2m3[i] = a2m3[i] + (gsf_coef22[s]-gsf_coef33[s])*pfunc

                    dm1 = dfunc

            # now combine a2 + a3 = a2p3 with a2 - a3 = a2m3
            p22mie = 0.5*(a2p3 + a2m3)
            p33mie = a2p3 - p22mie

            p22pyobs = np.zeros(ntheta)
            p33pyobs = np.zeros(ntheta)
            for s in range(self.nMom):
                P = legendre(s)
                p22pyobs += leg_coef22[s]* P(mu)
                p33pyobs += leg_coef33[s]* P(mu)

            ax = axes[1,1]
            ax.plot(angle, p22pyobs/p11pyobs,label='pyobs')
            ax.plot(angle,p22mie/p11mie,label='mieobs')
            ax.legend()
            ax.set_title('P22 = P5')

            ax = axes[0,2]
            ax.plot(angle, p33pyobs/p11pyobs,label='pyobs')
            ax.plot(angle,p33mie/p11mie,label='mieobs')
            ax.legend()
            ax.set_title('33 = P3')

        # P34
        # ordered over nPol as P11, P12, P33, P34, P22, P44
        do_P34 = True
        if do_P34:
            ipol = 3
            m = 0
            n = 2
            gsf_coef = self.miepmom[0].squeeze()[0,:,ipol]
            leg_coef = self.pyobspmom[0].squeeze()[0,:,ipol]

            p34 = np.zeros(ntheta)
            for i,t in enumerate(theta):
                d0, dneg1 = eval_gsfun.d_mn(m,n,0,t)
                pfunc = d0/(1j**(n-m)).real
                p34[i] = gsf_coef[0]*pfunc
                d1, d0 = eval_gsfun.d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
                pfunc = d1/(1j**(n-m)).real
                p34[i] = p34[i] + gsf_coef[1]*pfunc

                dm2 = d0
                dm1 = d1
                for s in range(2,self.nMom):
                    dfunc, dm2 = eval_gsfun.d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                    pfunc = dfunc/(1j**(n-m)).real

                    p34[i] = p34[i] + gsf_coef[s]*pfunc

                    dm1 = dfunc

            p34mie = p34

            p34pyobs = np.zeros(ntheta)
            for s in range(self.nMom):
                P = legendre(s)
                p34pyobs += leg_coef[s]* P(mu)

            ax = axes[1,0]
            ax.plot(angle, p34pyobs/p11pyobs,label='pyobs')
            ax.plot(angle,p34mie/p11mie,label='mieobs')
            ax.legend()
            ax.set_title('34 = P4')


        # Resonctruct P44
        do_p44 = True
        if do_p44:
            ipol = 5
            p44pyobs = np.zeros(ntheta)
            p44mie   = np.zeros(ntheta)

            miepmom = self.miepmom[0].squeeze()[0,:,ipol]
            pyobspmom = self.pyobspmom[0].squeeze()[0,:,ipol]
            for s in range(self.nMom):
                P = legendre(s)
                p44pyobs += pyobspmom[s]* P(mu)
                p44mie += miepmom[s]* P(mu)

            ax = axes[1,2]
            ax.plot(angle, p44pyobs/p11pyobs,label='pyobs')
            ax.plot(angle,p44mie/p11mie,label='mieobs')
            ax.legend()
            ax.set_title('P44 = P6')

        # Adjust layout and display
        plt.tight_layout()
        plt.show()
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    miercFile     = 'rc/Aod_EOS_ich.rc'
    pyobsrcFile   = 'm2_aop.yaml'
    channels      = '360:380:5'
#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t",
                        help="iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with track input file names")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("--channels",default=channels,
                        help="channels (detault={})".format(channels))

    parser.add_argument("--miercFile",default=miercFile,
                        help="miercFile (default={})".format(miercFile))

    parser.add_argument("--pyobsrcFile",default=pyobsrcFile,
                        help="pyobsrcFile (default={})".format(pyobsrcFile))

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf             = Config(args.inst_pcf,delim=' = ')
    instname       = cf('instname')

    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()

    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')


    # get input file names, calculate AOPS
    # ------------------------------------
    date      = isoparser(args.iso_t)

    nymd  = str(date.date()).replace('-','')
    year  = str(date.year)
    month = str(date.month).zfill(2)
    day   = str(date.day).zfill(2)
    hour  = str(date.hour).zfill(2)
    minute = str(date.minute).zfill(2)

    inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

    chs = args.channels.split(':')
    channels = np.arange(float(chs[0]),float(chs[1]),int(chs[2]))
    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    print('++++Running VLIDORT with the following arguments+++')
    print('>>>channel:   ',args.channels)
    print('>>>inFile:    ',inFile)
    print('>>>rcFile:    ',args.miercFile)
    print('>>>rcFile:    ',args.pyobsrcFile)
    print('>>>verbose:   ',args.verbose)
    print('++++End of arguments+++')
    
    vlidort = OPTICS_VLIDORT(inFile,channels,args.miercFile,args.pyobsrcFile,
                        instname,
                        args.dryrun,
                        verbose=args.verbose)
    vlidort.comparePMOM()
