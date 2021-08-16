#!/usr/bin/env python

"""
    Wrapper to loop through channels and submit pace_lc.j jobs to sbatch for one pace granule
"""

import os
import subprocess
import shutil
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import argparse
import numpy           as np
import time
from   pace            import granules
from   netCDF4         import Dataset
from   glob            import glob
from   MAPL.ShaveMantissa_ import shave32
from   MAPL.constants import MAPL_UNDEF
from   pyhdf.SD import SD, SDC
from   scipy.interpolate import interp1d
from collections import OrderedDict

MISSING = np.float32(1e20)

nlev = 72

SDS_RT = OrderedDict()
SDS_RT['surf_ref_I_blue'] = ['surface reflectance for blue CCD','None',None]
SDS_RT['surf_ref_I_red']  = ['surface reflectance for red CCD','None',None]
SDS_RT['surf_ref_I_SWIR'] = ['surface reflectance for SWIR bands','None',None]
SDS_RT['ref_blue'] = ['TOA reflectance for blue CCD','None',None]
SDS_RT['ref_red']  = ['TOA reflectance for red CCD','None',None]
SDS_RT['ref_SWIR'] = ['TOA reflectance for SWIR bands','None',None]


SDS_ADD =  OrderedDict()
SDS_ADD['SLEAVE_blue'] =  ['water leaving radiance for blue bands','W m-2 sr-1 nm-1',None]
SDS_ADD['SLEAVE_red']  =  ['water leaving radiance for red bands','W m-2 sr-1 nm-1',None]
SDS_ADD['SLEAVE_SWIR'] =  ['water leaving radiance for SWIR bands','W m-2 sr-1 nm-1',None]
SDS_ADD['IRR_blue']    =  ['solar irradiance radiance for blue bands','W m-2 sr-1 nm-1',None]
SDS_ADD['IRR_red']     =  ['solar irradiance for red bands','W m-2 sr-1 nm-1',None]
SDS_ADD['IRR_SWIR']    =  ['solar irradiance for SWIR bands','W m-2 sr-1 nm-1',None]

#         'surf_ref_Q_blue': ['surface polarized reflectance Q for blue CCD','None',None],
#         'surf_ref_Q_red':  ['surface polarized reflectance Q for red CCD','None',None],
#         'surf_ref_Q_SWIR': ['surface polarized reflectance Q for SWIR bands','None',None],
#         'surf_ref_U_blue': ['surface polarized reflectance U for blue CCD','None',None],
#         'surf_ref_U_red':  ['surface polarized reflectance U for red CCD','None',None],
#         'surf_ref_U_SWIR': ['surface polarized reflectance U for SWIR bands','None',None],
#         'Q_blue': ['TOA polarized radiance Q for blue CCD','None',None],
#         'Q_red':  ['TOA polarized radiance Q for red CCD','None',None],
#         'Q_SWIR': ['TOA polarized radiance Q for SWIR bands','None',None],
#         'U_blue': ['TOA polarized radiance U for blue CCD','None',None],
#         'U_red':  ['TOA polarized radiance U for red CCD','None',None],
#         'U_SWIR': ['TOA polarized radiance U for SWIR bands','None',None]         }

SDS_GAS =  OrderedDict()
SDS_GAS['ROD_blue']         = ['rayleigh optical depth for blue CCD','None',None]
SDS_GAS['ROD_red']          = ['rayleigh optical depth for red CCD','None',None]
SDS_GAS['ROD_SWIR']         = ['rayleigh optical depth for SWIR bands','None',None]
SDS_GAS['TRANS_RAY_blue']   = ['rayleigh transmittance for blue CCD','None',None]
SDS_GAS['TRANS_RAY_red']    = ['rayleigh transmittance for red CCD','None',None]
SDS_GAS['TRANS_RAY_SWIR']   = ['rayleigh transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_O2_blue']    = ['o2 transmittance for blue CCD','None',None]
SDS_GAS['TRANS_O2_red']     = ['o2 transmittance for red CCD','None',None]
SDS_GAS['TRANS_O2_SWIR']    = ['o2 transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_N2O_blue']   = ['n2o transmittance for blue CCD','None',None]
SDS_GAS['TRANS_N2O_red']    = ['n2o transmittance for red CCD','None',None]
SDS_GAS['TRANS_N2O_SWIR']   = ['n2o transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_CH4_blue']   = ['ch4 transmittance for blue CCD','None',None]
SDS_GAS['TRANS_CH4_red']    = ['ch4 transmittance for red CCD','None',None]
SDS_GAS['TRANS_CH4_SWIR']   = ['ch4 transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_O3_blue']    = ['o3 transmittance for blue CCD','None',None]
SDS_GAS['TRANS_O3_red']     = ['o3 transmittance for red CCD','None',None]
SDS_GAS['TRANS_O3_SWIR']    = ['o3 transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_H2O_blue']   = ['h2o transmittance for blue CCD','None',None]
SDS_GAS['TRANS_H2O_red']    = ['h2o transmittance for red CCD','None',None]
SDS_GAS['TRANS_H2O_SWIR']   = ['h2o transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_CO_blue']    = ['co transmittance for blue CCD','None',None]
SDS_GAS['TRANS_CO_red']     = ['co transmittance for red CCD','None',None]
SDS_GAS['TRANS_CO_SWIR']    = ['co transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_CO2_blue']   = ['co2 transmittance for blue CCD','None',None]
SDS_GAS['TRANS_CO2_red']    = ['co2 transmittance for red CCD','None',None]
SDS_GAS['TRANS_CO2_SWIR']   = ['co2 transmittance for SWIR bands','None',None]
SDS_GAS['TRANS_TOTAL_blue'] = ['total gas transmittance for blue CCD','None',None]
SDS_GAS['TRANS_TOTAL_red']  = ['total transmittance for red CCD','None',None]
SDS_GAS['TRANS_TOTAL_SWIR'] = ['total transmittance for SWIR bands','None',None]

SDS_CLD =  OrderedDict()
SDS_CLD['lcod_blue']    = ['liquid cloud optical depth for blue CCD','None',None]
SDS_CLD['lc_ssa_blue']  = ['liquid cloud single scattering albedo for blue CCD','None',None]
# SDS_CLD['lc_g_blue'] =  ['liquid cloud assymetry parameter for blue CCD','None',nlev]
SDS_CLD['icod_blue']    = ['ice cloud optical depth for blue CCD','None',None]
SDS_CLD['ic_ssa_blue']  =  ['ice  cloud single scattering albedo for blue CCD','None',None]
# SDS_CLD['ic_g_blue'] =  ['ice  cloud assymetry parameter for blue CCD','None',nlev]
SDS_CLD['lcod_red']     = ['liquid cloud optical depth for red CCD','None',None]
SDS_CLD['lc_ssa_red']   = ['liquid cloud single scattering albedo for red CCD','None',None]
# SDS_CLD['lc_g_red'] =  ['liquid cloud assymetry parameter for red CCD','None',nlev]
SDS_CLD['icod_red']     = ['ice  cloud optical depth for red CCD','None',None]
SDS_CLD['ic_ssa_red']   = ['ice  cloud single scattering albedo for red CCD','None',None]
# SDS_CLD['ic_g_red'] =  ['ice  cloud assymetry parameter for red CCD','None',nlev]
SDS_CLD['lcod_SWIR']    = ['liquid cloud optical depth for SWIR bands','None',None]
SDS_CLD['lc_ssa_SWIR']  = ['liquid cloud single scattering albedo for SWIR bands','None',None]
# SDS_CLD['lc_g_SWIR'] =  ['liquid cloud assymetry parameter for SWIR bands','None',nlev]
SDS_CLD['icod_SWIR']    = ['ice  cloud optical depth for SWIR bands','None',None]
SDS_CLD['ic_ssa_SWIR']  =  ['ice  cloud single scattering albedo for SWIR bands','None',None]
# SDS_CLD['ic_g_SWIR'] =  ['ice  cloud assymetry parameter for SWIR bands','None',nlev]
          
SDS_AER =  OrderedDict()
SDS_AER['aod_blue'] = ['aerosol optical depth for blue CCD','None',None]
SDS_AER['ssa_blue'] = ['aerosol single scattering albedo for blue CCD','None',None]
SDS_AER['g_blue']   = ['aerosol assymetry parameter for blue CCD','None',None]
SDS_AER['aod_red']  = ['aerosol optical depth for red CCD','None',None]
SDS_AER['ssa_red']  = ['aerosol single scattering albedo for red CCD','None',None]
SDS_AER['g_red']    = ['aerosol assymetry parameter for red CCD','None',None]
SDS_AER['aod_SWIR'] = ['aerosol optical depth for SWIR bands','None',None]
SDS_AER['ssa_SWIR'] = ['aerosol single scattering albedo for SWIR bands','None',None]
SDS_AER['g_SWIR']   = ['aerosol assymetry parameter for SWIR bands','None',None]
         

#---
def shave(q,undef=MAPL_UNDEF,has_undef=1,nbits=12):
    """
    Shave variable. On input, nbits is the number of mantissa bits to keep
    out of maximum of 24.
    """

    # Determine shaving parameters
    # ----------------------------
    xbits = 24 - nbits
    shp = q.shape
    rank = len(shp)
    if rank == 2:  # yx
        chunksize = shp[0]*shp[1] 
    elif rank == 3: # zyx
        chunksize = shp[1]*shp[2]
    else:
        raise ValueError, "invalid rank=%d"%rank

    # Shave it
    # --------
    qs, rc = shave32(q.ravel(),xbits,has_undef,undef,chunksize)
    if rc:
        raise ValueError, "error on return from shave32, rc=%d"%rc

    return qs.reshape(shp)


class LEVELC(object):
    """
    LevelC Object for populating L1B file
    """
    def __init__(self,args):
        self.Date    = isoparser(args.iso_t1)


        self.cwd           = os.getcwd()
        self.rootdir       = args.rootdir
        self.run_name      = args.run_name
        self.write_add     = args.write_add
        self.write_aer     = args.write_aer
        self.write_cld     = args.write_cld
        self.rm_ch         = args.rm_ch
        self.do_sleave_adjust = args.do_sleave_adjust
        self.do_single_xtrack = args.do_single_xtrack
        self.ALPHA_file    = args.ALPHA_file
        self.RSR_file      = args.RSR_file
        self.IRR_file      = args.IRR_file
        self.version       = args.version

        self.outfilelist = None
        self.addfilelist = None
        self.aerfilelist = None
        self.cldfilelist = None

        YMDdir    = self.Date.strftime('Y%Y/M%m/D%d')
        pYMDdir   = self.Date.strftime('Y2020/M%m/D%d')
        self.LcDir     = '{}/LevelC/{}/v{}'.format(self.rootdir,YMDdir,self.version)
        self.L1bDir    = '{}/L1B/{}'.format(self.rootdir,pYMDdir)
        self.Lc2Dir    = '{}/LevelC2/{}/v{}'.format(rootdir,YMDdir,self.version)
        self.hms       = self.Date.strftime('%H%M00')
        self.nymd      = self.Date.strftime('%Y%m%d')

        pDate = isoparser(self.Date.strftime('2020-%m-%dT%H:%M:00'))
        self.nyj = pDate.strftime('%Y%j')
        self.L1B_file = '{}/OCI{}{}.L1B_PACE.nc'.format(self.L1bDir,self.nyj,self.hms)

        # read in LBL and OCI channels
        self.get_channels()

        # read in the RSR
        self.get_RSR()

        # get list of file & merge them
        self.outfilelist = self.get_outfilelist('vlidort')
        self.rsr_weighted_radiance(self.outfilelist)
        self.populate_L1B()
        self.addfilelist = self.get_outfilelist('add')
        self.condense_LC()


        if self.write_add:
            self.condense_GAS('gas',self.addfilelist,SDS_GAS)

        if self.write_aer:
            self.aerfilelist = self.get_outfilelist('aerosol')
            self.condense_AER_CLD('aerosol',self.aerfilelist,SDS_AER)

        if self.write_cld:
            self.cldfilelist = self.get_outfilelist('cloud')
            self.condense_AER_CLD('cloud',self.cldfilelist,SDS_CLD)
# ---
    def rod_rsr_weighted(self,filelist,varname):
        """
        weight rod by RSR
        this is only for ck wavs
        """

        # get granule dimensions
        nc = Dataset(filelist[0])
        nccd = len(nc.dimensions['ccd_pixels'])
        nscan = len(nc.dimensions['number_of_scans'])
        nc.close()

        # read in data
        var = []
        for filename in filelist:
            nc = Dataset(filename)
            chname = self.get_chname(filename)
            pvar = nc.variables[varname+'_'+chname]
            rank = len(pvar.shape)
            if rank == 4:
                var.append(pvar[0,0:1,:,:])
            else:
                var.append(pvar[0:1,:,:])

            nc.close()

        # split into UV and CK wavelengths
        var = np.concatenate(var,axis=0)
        varck = var[self.ick,:,:]
        self.varuv = var[self.iuv,:,:]
        var = None

        # regrid CK so it's a consistent resolution
        # integrate over RSR
        # but only for OCI channels covered by CK bins
        wav_ck_grid = np.arange(self.wav_ck.min(),self.wav_ck.max()+0.5,0.5)
        istart = 53
        ck_smooth = np.ones([self.noci-istart,nscan,nccd])*MISSING

        rsr_f = interp1d(self.wav_rsr,self.rsr,kind='linear',fill_value=0.0,bounds_error=False)
        rsr_int = rsr_f(wav_ck_grid)
        norms = np.zeros(self.noci)
        for ich in range(self.noci):
            norms[ich] = np.trapz(rsr_int[ich,:],wav_ck_grid)

        if self.do_single_xtrack:
            iscan = 600
            for iccd in range(nccd):
                ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                varck_grid = ck_f(wav_ck_grid)
                for ich in range(istart,self.noci):
                    norm = norms[ich]
                    ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:],wav_ck_grid)/norm
        else:
            for iscan in range(nscan):
                for iccd in range(nccd):
                    ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                    varck_grid = ck_f(wav_ck_grid)
                    for ich in range(istart,self.noci):
                        norm = norms[ich]
                        ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:],wav_ck_grid)/norm


        self.varck = ck_smooth

# ---
    def trans_rsr_weighted(self,filelist,varname):
        """
        weight transmittance by RSR
        this is only for ck wavs
        """

        # get granule dimensions
        nc = Dataset(filelist[0])
        nccd = len(nc.dimensions['ccd_pixels'])
        nscan = len(nc.dimensions['number_of_scans'])
        nc.close()

        # read in data
        var = []
        F0  = []
        for filename in filelist:
            nc = Dataset(filename)
            chname = self.get_chname(filename)
            pvar = nc.variables[varname+'_'+chname]
            rank = len(pvar.shape)
            if rank == 4:
                var.append(pvar[0,0:1,:,:])
            else:
                var.append(pvar[0:1,:,:])

            pvar = nc.variables['IRR_'+chname]
            F0.append(pvar[0:1,600,0])
            nc.close()

        # split into UV and CK wavelengths
        var = np.concatenate(var,axis=0)
        F0  = np.concatenate(F0,axis=0)
        varck = var[self.ick,:,:]
        F0ck  = F0[self.ick]
        self.varuv = var[self.iuv,:,:]
        var = None
        F0  = None

        # regrid CK so it's a consistent resolution
        # integrate over RSR
        # but only for OCI channels covered by CK bins
        wav_ck_grid = np.arange(self.wav_ck.min(),self.wav_ck.max()+0.5,0.5)
        istart = 53
        ck_smooth = np.ones([self.noci-istart,nscan,nccd])*MISSING

        F0_f = interp1d(self.wav_ck,F0ck,kind='linear')
        F0_int = F0_f(wav_ck_grid)

        rsr_f = interp1d(self.wav_rsr,self.rsr,kind='linear',fill_value=0.0,bounds_error=False)
        rsr_int = rsr_f(wav_ck_grid)
        norms = np.zeros(self.noci)
        for ich in range(self.noci):
            norms[ich] = np.trapz(rsr_int[ich,:]*F0_int,wav_ck_grid)

        if self.do_single_xtrack:
            iscan = 600
            for iccd in range(nccd):
                ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                varck_grid = ck_f(wav_ck_grid)
                for ich in range(istart,self.noci):
                    norm = norms[ich]
                    ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:]*F0_int,wav_ck_grid)/norm
        else:
            for iscan in range(nscan):
                for iccd in range(nccd):
                    ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                    varck_grid = ck_f(wav_ck_grid)
                    for ich in range(istart,self.noci):
                        norm = norms[ich]
                        ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:]*F0_int,wav_ck_grid)/norm


        self.varck = ck_smooth
# ---
    def rsr_weighted_radiance(self,filelist):
        """
        weight get radiance from RSR weighted reflectance
        this is only for ck wavs
        """

        # get granule dimensions
        nc = Dataset(filelist[0])
        nccd = len(nc.dimensions['ccd_pixels'])
        nscan = len(nc.dimensions['number_of_scans'])
        nc.close()

        # read in reflectance data
        varname = 'ref'
        var = []
        for filename in filelist:
            nc = Dataset(filename)
            chname = self.get_chname(filename)
            pvar = nc.variables[varname+'_'+chname]
            rank = len(pvar.shape)
            if rank == 4:
                var.append(pvar[0,0:1,:,:])
            else:
                var.append(pvar[0:1,:,:])
            nc.close()

        # Read in RSR weighted solar irradiance
        nc = Dataset(self.IRR_file)
        irr     = nc.variables['F0_rsr'][:]
        # Convert to W/m2/nm
        irr = irr*1e-2
        irr.shape = (len(irr),1,1)
        nc.close()

        # read in sza
        nc = Dataset(self.L1B_file)
        sza = nc.groups['geolocation_data'].variables['solar_zenith'][:]
        csza = np.cos(np.radians(sza))
        nc.close()        

        # split into UV and CK wavelengths
        var = np.concatenate(var,axis=0)
        varck = var[self.ick,:,:]
        self.varuv = var[self.iuv,:,:]
        var = None

        # regrid CK so it's a consistent resolution
        # integrate over RSR
        # but only for OCI channels covered by CK bins
        wav_ck_grid = np.arange(self.wav_ck.min(),self.wav_ck.max()+0.5,0.5)
        istart = 53
        ck_smooth = np.ones([self.noci-istart,nscan,nccd])*MISSING
        rsr_f = interp1d(self.wav_rsr,self.rsr,kind='linear',fill_value=0.0,bounds_error=False)
        rsr_int = rsr_f(wav_ck_grid)
        norms = np.zeros(self.noci)
        for ich in range(self.noci):
            norms[ich] = np.trapz(rsr_int[ich,:],wav_ck_grid)

        if self.do_single_xtrack:
            iscan = 600
            for iccd in range(nccd):
                ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                varck_grid = ck_f(wav_ck_grid)
                for ich in range(istart,self.noci):
                    norm = norms[ich]
                    ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:],wav_ck_grid)/norm
        else:
            for iscan in range(nscan):
                for iccd in range(nccd):
                    ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                    varck_grid = ck_f(wav_ck_grid)
                    for ich in range(istart,self.noci):
                        norm = norms[ich]
                        ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:],wav_ck_grid)/norm

        self.varck = ck_smooth


        # go from reflectance to radiance
        self.varuv = self.varuv*irr[:istart,:,:]*csza/np.pi        
        self.varck = self.varck*irr[istart:,:,:]*csza/np.pi

# ---
    def rsr_weighted(self,filelist,varname):
        """
        weight parameter by RSR
        this is only for ck wavs
        """

        # get granule dimensions
        nc = Dataset(filelist[0])
        nccd = len(nc.dimensions['ccd_pixels'])
        nscan = len(nc.dimensions['number_of_scans'])
        nc.close()

        # read in data
        var = []
        for filename in filelist:
            nc = Dataset(filename)
            chname = self.get_chname(filename)                   
            pvar = nc.variables[varname+'_'+chname] 
            rank = len(pvar.shape)
            if rank == 4:
                var.append(pvar[0,0:1,:,:])
            else:
                var.append(pvar[0:1,:,:])
            nc.close()

        # split into UV and CK wavelengths
        var = np.concatenate(var,axis=0)
        varck = var[self.ick,:,:]
        self.varuv = var[self.iuv,:,:]
        var = None

        # regrid CK so it's a consistent resolution
        # integrate over RSR
        # but only for OCI channels covered by CK bins
        wav_ck_grid = np.arange(self.wav_ck.min(),self.wav_ck.max()+0.5,0.5)
        istart = 53
        ck_smooth = np.ones([self.noci-istart,nscan,nccd])*MISSING
        rsr_f = interp1d(self.wav_rsr,self.rsr,kind='linear',fill_value=0.0,bounds_error=False)
        rsr_int = rsr_f(wav_ck_grid)        
        norms = np.zeros(self.noci)
        for ich in range(self.noci):
            norms[ich] = np.trapz(rsr_int[ich,:],wav_ck_grid)

        if self.do_single_xtrack:
            iscan = 600
            for iccd in range(nccd):
                ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                varck_grid = ck_f(wav_ck_grid)
                for ich in range(istart,self.noci):
                    norm = norms[ich]
                    ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:],wav_ck_grid)/norm
        else:
            for iscan in range(nscan):
                for iccd in range(nccd):
                    ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                    varck_grid = ck_f(wav_ck_grid)
                    for ich in range(istart,self.noci):
                        norm = norms[ich]
                        ck_smooth[ich-istart,iscan,iccd] = np.trapz(varck_grid*rsr_int[ich,:],wav_ck_grid)/norm
            
            
        self.varck = ck_smooth
# ---
    def center_interpolated(self,filelist,varname):
        """
        parameter interpolated to OCI band centers
        this is only for ck wavs
        """

        # get granule dimensions
        nc = Dataset(filelist[0])
        nccd = len(nc.dimensions['ccd_pixels'])
        nscan = len(nc.dimensions['number_of_scans'])
        nc.close()

        # read in data
        var = []
        for filename in filelist:
            nc = Dataset(filename)
            chname = self.get_chname(filename)
            pvar = nc.variables[varname+'_'+chname]
            rank = len(pvar.shape)
            if rank == 4:
                var.append(pvar[0,0:1,:,:])
            else:
                var.append(pvar[0:1,:,:])
            nc.close()

        # split into UV and CK wavelengths
        var = np.concatenate(var,axis=0)
        varck = var[self.ick,:,:]
        self.varuv = var[self.iuv,:,:]
        var = None

        # interpolate to OCI center wavelengths
        # but only for OCI channels covered by CK bins
        istart = 53
        ich = range(istart,self.noci)
        ck_center = np.ones([self.noci-istart,nscan,nccd])*MISSING
        if self.do_single_xtrack:
            iscan = 600
            for iccd in range(nccd):
                ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                ck_center[:,iscan,iccd] = ck_f(self.wav_cen[ich])
        else:
            for iscan in range(nscan):
                for iccd in range(nccd):
                        ck_f = interp1d(self.wav_ck,varck[:,iscan,iccd],kind='linear')
                        ck_center[:,iscan,iccd] = ck_f(self.wav_cen[ich])


        self.varck = ck_center

        
# ---
    def get_chname(self,filename):
        fname = os.path.basename(filename)
        lchname = fname.split('.')[-2]
        chname = lchname.split('-')[-1].replace('d','.')

        return chname
# ---
    def get_RSR(self):
        """
        Read in OCI RSR File
        """
        hdf = SD(self.RSR_file, SDC.READ)
        rsr = hdf.select('RSR')[:]
        wav_rsr = hdf.select('rsrwave')[:]
        wav_oci = hdf.select('wave')[:]
        hdf.end()

        #re-alias 1038 to 1040 in L1b file
        ich = wav_oci == 1038
        wav_oci[ich] = 1040

        # right now just copy this over.
        # will have to be read for newer RSR
        wav_cen = wav_oci

        self.rsr = rsr.T
        self.wav_rsr = wav_rsr
        self.wav_oci = wav_oci
        self.wav_cen = wav_cen
        self.noci = len(wav_oci) 
# ---
    def get_outfilelist(self,name):
         # Create file if you need to
        if self.run_name is None:
            outfile = '{}/channel/pace-g5nr.lc.{}.{}_{}.*.nc4'.format(self.LcDir,name,self.nymd,self.hms)
        else:
            outfile = '{}/channel/pace-g5nr.lc.{}.{}.{}_{}.*.nc4'.format(self.LcDir,name,self.run_name,self.nymd,self.hms)

        return sorted(glob(outfile))

# ---
    def get_channels(self):
        nc = Dataset(self.ALPHA_file)
        alpha_channels  = nc.variables['channels'][:]
        bins            = nc.variables['nbins'][:]
        iuv = bins == 1
        self.wav_uv = np.array(alpha_channels[iuv])
        self.iuv    = iuv
        ick = bins > 1
        self.wav_ck = np.array(alpha_channels[ick])
        self.ick = ick
        nc.close()

# ---        

    def populate_L1B(self):
        """
        insert TOA reflectance to L1B file
        """
        if self.run_name is None:
            outfile = '{}/OCI{}{}.L1B_PACE.nc'.format(self.Lc2Dir,self.nyj,self.hms)
        else:
            outfile = '{}/OCI{}{}.L1B_PACE.{}.nc'.format(self.Lc2Dir,self.nyj,self.hms,self.run_name)
        shutil.copyfile(self.L1B_file,outfile)


        ncmerge = Dataset(outfile,mode='r+')
        SDS = {'blue_wavelength' : 'Lt_blue',
               'red_wavelength'  : 'Lt_red',
               'SWIR_wavelength' : 'Lt_SWIR'}
        for sds in SDS:
            pchannels = np.array(ncmerge.groups['sensor_band_parameters'].variables[sds][:])
            pvar      = ncmerge.groups['observation_data'].variables[SDS[sds]]

            # first do UV channels
            for ich,ch in enumerate(self.wav_uv):
                for i,pch in enumerate(pchannels):
                    if float(ch) == pch:
                        print 'inserting ',ch
                        data = self.varuv[ich,:,:]
                        ii = data == -500.0
                        if np.sum(ii) > 0:
                            if (not hasattr(data,'mask')):
                                data = np.ma.array(data)
                                data.mask = ii
                            elif (type(data.mask) is np.bool_):
                                data.mask = ii
                            else:
                                data.mask[ii] = True
                        pvar[i,:,:] = data

            # now do other channels
            istart = 53
            for ich, ch in enumerate(self.wav_oci[istart:]):
                for i,pch in enumerate(pchannels):
                    if float(ch) == pch:
                        print 'inserting ',ch
                        data = self.varck[ich,:,:]
                        ii = data == -500.0
                        if np.sum(ii) > 0:
                            if (not hasattr(data,'mask')):
                                data = np.ma.array(data)
                                data.mask = ii
                            elif (type(data.mask) is np.bool_):
                                data.mask = ii
                            else:
                                data.mask[ii] = True
                        pvar[i,:,:] = data                 

        ncmerge.close()
# ---
    def condense_LC(self):

        if self.run_name is None:
            outfile = '{}/pace-g5nr.lc.vlidort.{}_{}.nc4'.format(self.LcDir,self.nymd,self.hms)
        else:
            outfile = '{}/pace-g5nr.lc.vlidort.{}.{}_{}.nc4'.format(self.LcDir,self.run_name,self.nymd,self.hms)

        # create new outfile
        if self.do_sleave_adjust:
            SDS_RT['ADJUSTED_SLEAVE_red']= ['transmittance adjusted water leaving radiance for red bands','W m-2 sr-1 nm-1',None]
            SDS_RT['ADJUSTED_SLEAVE_blue']= ['transmittance adjusted water leaving radiance for blue bands','W m-2 sr-1 nm-1',None]
            SDS_RT['ADJUSTED_SLEAVE_SWIR']= ['transmittance adjusted water leaving radiance for SWIR bands','W m-2 sr-1 nm-1',None]

        create_condenseFile(self.L1B_file,outfile,self.Date,dict(SDS_RT, **SDS_ADD))

        # Condense RT stuff
        SDS = []
        for sds in SDS_RT:
            sds_ = '_'.join(sds.split('_')[:-1])
            SDS.append(sds_)

        SDS = np.unique(np.array(SDS))
        for sds in SDS:
            if 'ref' in sds:
                self.center_interpolated(self.outfilelist,sds)
            else:
                self.rsr_weighted(self.outfilelist,sds)
            self.insert_condenseVar(outfile,sds)

        # Condense ADD stuff
        SDS = []
        for sds in SDS_ADD:
            sds_ = '_'.join(sds.split('_')[:-1])
            SDS.append(sds_)

        SDS = np.unique(np.array(SDS))
        for sds in SDS:
            if 'ref' in sds:
                self.center_interpolated(self.addfilelist,sds)
            else:
                self.rsr_weighted(self.addfilelist,sds)
            self.insert_condenseVar(outfile,sds)


        #Compress
        cmd = "$BASEDIR/Linux/bin/ncks --hst --no_abc -O -4 -L 2 --cnk_plc=g2d --cnk_dmn ccd_pixels,91 --cnk_dmn number_of_scans,144 --cnk_dmn lev,1 {} {}"
        newcmd = cmd.format(outfile,outfile)

        devnull = open(os.devnull, 'w')
        stat = subprocess.call(newcmd, shell=True, stdout=devnull)
        devnull.close()

        # Remove channel files if desired and compression was successful
        if self.rm_ch and (stat ==0):
            for ofile in self.outfilelist:
                os.remove(ofile)


    def insert_condenseVar(self,outfile,sds_):
        # insert data into correct place in outfile
        SDS = []
        for ch in ['_blue','_red','_SWIR']:
            SDS.append(sds_+ch)
 
        for sds in SDS:
            if 'blue' in sds:
                chname = 'blue_wavelength'
            elif 'red' in sds:
                chname = 'red_wavelength'
            elif 'SWIR' in sds:
                chname = 'SWIR_wavelength'

            ncmerge = Dataset(outfile,mode='r')
            pchannels = np.array(ncmerge.variables[chname][:])
            ncmerge.close()

            if '_' in sds:
                oname = '_'.join(sds.split('_')[:-1])
            else:
                oname = sds

            ncmerge = Dataset(outfile,mode='r+')
            pvar      = ncmerge.variables[sds]
            # first do UV channels
            for ich,ch in enumerate(self.wav_uv):
                for i,pch in enumerate(pchannels):
                    if float(ch) == pch:
                        print 'inserting ',ch, ' to ',sds
                        data = self.varuv[ich,:,:]
                        pvar[:,:,i] = data

            # now do other channels
            istart = 53
            for ich, ch in enumerate(self.wav_oci[istart:]):
                for i,pch in enumerate(pchannels):
                    if float(ch) == pch:
                        print 'inserting ',ch, ' to ',sds
                        data = self.varck[ich,:,:]
                        pvar[:,:,i] = data

            ncmerge.close()
# ---                        
    
    def condense_GAS(self,name,filelist,SDS):
        
        if self.run_name is None:
            outfile = '{}/pace-g5nr.lc.{}.{}_{}.nc4'.format(self.LcDir,name,self.nymd,self.hms)
        else:
            outfile = '{}/pace-g5nr.lc.{}.{}.{}_{}.nc4'.format(self.LcDir,name,self.run_name,self.nymd,self.hms)

        create_condenseFile(self.L1B_file,outfile,self.Date,SDS)

        # Condense stuff
        SDS_ = []
        for sds in SDS:
            sds_ = '_'.join(sds.split('_')[:-1])
            SDS_.append(sds_)

        SDS_ = np.unique(np.array(SDS_))
        for sds in SDS_:
            if 'ROD' in sds:
                self.rod_rsr_weighted(filelist,sds)
                self.insert_condenseVar(outfile,sds)
            else:
                self.trans_rsr_weighted(filelist,sds)
                self.insert_condenseVar(outfile,sds)


        #Compress
        cmd = "$BASEDIR/Linux/bin/ncks --hst --no_abc -O -4 -L 2 --cnk_plc=g2d --cnk_dmn ccd_pixels,91 --cnk_dmn number_of_scans,144 --cnk_dmn lev,1 {} {}"
        newcmd = cmd.format(outfile,outfile)

        devnull = open(os.devnull, 'w')
        stat = subprocess.call(newcmd, shell=True, stdout=devnull)
        devnull.close()

        # Remove channel files if desired and compression was successful
        if self.rm_ch and (stat == 0):
            for ofile in filelist:
                os.remove(ofile)
# ---
    def condense_AER_CLD(self,name,filelist,SDS):

        if self.run_name is None:
            outfile = '{}/pace-g5nr.lc.{}.{}_{}.nc4'.format(self.LcDir,name,self.nymd,self.hms)
        else:
            outfile = '{}/pace-g5nr.lc.{}.{}.{}_{}.nc4'.format(self.LcDir,name,self.run_name,self.nymd,self.hms)

        create_condenseFile(self.L1B_file,outfile,self.Date,SDS)

        # Condense stuff
        SDS_ = []
        for sds in SDS:
            sds_ = '_'.join(sds.split('_')[:-1])
            SDS_.append(sds_)

        SDS_ = np.unique(np.array(SDS_))
        for sds in SDS_:
            self.center_interpolated(filelist,sds)
            self.insert_condenseVar(outfile,sds)

        #Compress
        cmd = "$BASEDIR/Linux/bin/ncks --hst --no_abc -O -4 -L 2 --cnk_plc=g2d --cnk_dmn ccd_pixels,91 --cnk_dmn number_of_scans,144 --cnk_dmn lev,1 {} {}"
        newcmd = cmd.format(outfile,outfile)

        devnull = open(os.devnull, 'w')
        stat = subprocess.call(newcmd, shell=True, stdout=devnull)
        devnull.close()

        # Remove channel files if desired and compression was successful
        if self.rm_ch and (stat == 0):
            for ofile in filelist:
                os.remove(ofile)

# ---
def create_condenseFile(L1B_file,outfile,Date,SDS):

    # Open NC file
    # ------------
    print 'create outfile',outfile
    nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    print 'copy from',L1B_file
    nctrj = Dataset(L1B_file)

    # Get Dimensions
    # ----------------------
    npixel = len(nctrj.dimensions['ccd_pixels'])
    nscan  = len(nctrj.dimensions['number_of_scans'])
    nblue  = len(nctrj.dimensions['blue_bands'])
    nred   = len(nctrj.dimensions['red_bands'])
    nswir  = len(nctrj.dimensions['SWIR_bands'])


    # Set global attributes
    # ---------------------
    nc.title = "RT inputs for VLIDORT Simulation of GEOS-5 PACE Sampled Data"
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Condensed outputs from channel files'
    nc.references = 'n/a'
    nc.comment = 'n/a'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.Conventions = 'CF'    

    # Create dimensions
    # -----------------
    l = nc.createDimension('lev',nlev)
    le = nc.createDimension('leve',nlev+1)
    x = nc.createDimension('ccd_pixels',npixel)
    xs = nc.createDimension('SWIR_pixels',npixel)
    y = nc.createDimension('number_of_scans',nscan)
    b = nc.createDimension('blue_bands',nblue)
    r = nc.createDimension('red_bands',nred)
    s = nc.createDimension('SWIR_bands',nswir)

    # Save lon/lat
    # --------------------------
    _copyVar(nctrj,nc,u'longitude','geolocation_data',dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'latitude','geolocation_data',dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'ev_mid_time','scan_line_attributes', dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'blue_wavelength','sensor_band_parameters', dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'red_wavelength','sensor_band_parameters', dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'SWIR_wavelength','sensor_band_parameters', dtype='f4',zlib=True,verbose=True)


    # Loop over SDS creating each dataset
    #---------------------------------------
    for sds in SDS:
        lname,unit,levs = SDS[sds]
        if levs == nlev:
            if 'red' in sds:
                dim = ('lev','number_of_scans','ccd_pixels','red_bands',)
            elif 'blue' in sds:
                dim = ('lev','number_of_scans','ccd_pixels','blue_bands',)
            elif 'SWIR' in sds:
                dim = ('lev','number_of_scans','ccd_pixels','SWIR_bands',)
            else:
                dim = ('lev','number_of_scans','ccd_pixels')

        elif levs == nlev+1:
            if 'red' in sds:
                dim = ('leve','number_of_scans','ccd_pixels','red_bands',)
            elif 'blue' in sds:
                dim = ('leve','number_of_scans','ccd_pixels','blue_bands',)
            elif 'SWIR' in sds:
                dim = ('leve','number_of_scans','ccd_pixels','SWIR_bands',)
            else:
                dim = ('leve','number_of_scans','ccd_pixels')

        else:
            if 'red' in sds:
                dim = ('number_of_scans','ccd_pixels','red_bands',)
            elif 'blue' in sds:
                dim = ('number_of_scans','ccd_pixels','blue_bands',)
            elif 'SWIR' in sds:
                dim = ('number_of_scans','ccd_pixels','SWIR_bands',)
            else:
                dim = ('number_of_scans','ccd_pixels')


        this = nc.createVariable(sds,'f4',dim,zlib=True,fill_value=-999.0)  

        this.long_name = lname
        this.missing_value = -999.0
        this.unit = unit  


    nc.close()
    nctrj.close()          


#----
def _copyVar(ncIn,ncOut,name,group,dtype='f4',zlib=False,verbose=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.groups[group].variables[name]
    if verbose:
        print 'copy variable ',name,x.dimensions
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
    if hasattr(x,'long_name'): y.long_name = x.long_name
    if hasattr(x,'units'): y.units = x.units
    if hasattr(x,'missing_value'): y.missing_value = x.missing_value
    rank = len(x.shape)

    if rank == 1:
        y[:] = x[:]
    elif rank == 2:
        if hasattr(x,'missing_value'):
            y[:,:] = shave(x[:,:],undef=x.missing_value)
        else:
            y[:,:] = shave(x[:,:],has_undef=0)
    elif rank == 3:
        if hasattr(x,'missing_value'):
            y[:,:,:] = shave(x[:,:,:],undef=x.missing_value)
        else:
            y[:,:,:] = shave(x[:,:,:],has_undef=0)
    else:
        raise ValueError, "invalid rank of <%s>: %d"%(name,rank)


if __name__ == '__main__':
    
    #Defaults
    rootdir    = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    ALPHA_file = 'alphaTable_v0/alpha_CK_Thuillier_o3.nc4'
    RSR_file   = 'rsrFile/OCI_RSR_v0.hdf'
    IRR_file   = 'irrTable_v0/F0_rsr_weighted_V0.nc'
    version    = 2.0


    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')

    parser.add_argument("--rootdir",default=rootdir,
                        help="root directory for PACE data (default=%s)."%rootdir)       

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument('--version',default=version,
                        help="version numbering (default={})".format(version))

    parser.add_argument('--run_name',default=None,
                        help="Optional version naming for outfiles (default=None)")    

    parser.add_argument('--ALPHA_file',default=ALPHA_file,
                        help="gas cross section file (default=%s)"%ALPHA_file)

    parser.add_argument('--RSR_file',default=RSR_file,
                        help="OCI spectral response function (default=%s)"%RSR_file)

    parser.add_argument('--IRR_file',default=IRR_file,
                        help="RSR weighted solar irradiance file (default=%s)"%IRR_file)

    parser.add_argument("--do_single_xtrack",action="store_true",
                        help="for testing - mask all but 1 xtrack (default=False).")

    parser.add_argument("--do_sleave_adjust",action="store_true",
                        help="Do transmittance adjustment of surface leaving radiance (default=False).")

    parser.add_argument("--no_add",action="store_true",
                        help="Do NOT write additional output in VLIDORT call (default=False).")  

    parser.add_argument("--no_write_aer",action="store_true",
                        help="Do NOT aerosol output in VLIDORT call (default=False).")  

    parser.add_argument("--no_write_cld",action="store_true",
                        help="Do NOT cloud output in VLIDORT call (default=False).")                                                    
    parser.add_argument("--force",action="store_true",
                        help="Overwrite existing L1B file in LevelC directory (default=False).")    

    parser.add_argument("--rm_ch",action="store_true",
                        help="Remove channel specific files after condensing to one file (default=False).")


    args = parser.parse_args()
    args.write_add = True
    if args.no_add:
        args.write_add = False
        
    args.write_cld = True
    if args.no_write_cld:
        args.write_cld = False

    args.write_aer = True
    if args.no_write_aer:
        args.write_aer = False     

    lc = LEVELC(args)


