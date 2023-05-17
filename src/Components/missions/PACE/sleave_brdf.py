#!/usr/bin/env python3

"""
    Calculates correction to nobm water leaving radiance
    for sun-satellite geometry
    Interpolate morel table for correction factor
    adapted from get_brdf.py from Amir Ibrahim and Sean Bailey
"""

import os
from  optparse import OptionParser
import numpy as np
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator as rgi
from nobm import _copyVar, shave

class SLEAVE_BRDF(object):
    def __init__(self,t1,t2,options):
        dt = timedelta(minutes=5)

        # get filenames
        self.morelFile = options.morelFile
        while t1 <= t2:
            dirname = t1.strftime('/Y%Y/M%m/D%d/')
            fname = t1.strftime('%Y%m%d_%H%M00.nc4')
            self.nobmFile = options.lbRoot + '/surface/SLEAVE/NOBM' + dirname + 'pace-g5nr.lb.sleave.' + fname

            self.outFile = options.lbRoot + '/surface/SLEAVE/NOBM' + dirname + 'pace-g5nr.lb.fq.' + fname

            dirname = t1.strftime('/Y2020/M%m/D%d/')
            tj    = isoparser(t1.strftime('2020-%m-%dT%H:%M'))
            jdate = tj.strftime('%Y%j%H%M00') 
            self.l1bFile = options.l1bRoot + dirname + 'OCI' + jdate + '.L1B_PACE.nc'

            # read in geometry
            self.read_geom()

            # read chlorophyll and wavelengths
            self.read_chl()

            # set up morel file interpolator
            self.fq_table_interp()

            # interpolate morel table
            self.get_fq()

            # write to file
            self.writenc()

            t1 += dt
# ---
    def fq_table_interp(self):
        nc = Dataset(self.morelFile)
        data = nc['foq'][:]

        grdwvl = nc['wave'][:]
        grdchl = np.log10(nc['chl'][:])
        grdsolz = nc['solz'][:]
        grdviewz = nc['senz'][:]
        grdrelaz = nc['phi'][:]

        interp_func = rgi((grdwvl, grdsolz, grdchl, grdviewz, grdrelaz), data,\
                      fill_value=None,bounds_error=False)

        self.interp_func = interp_func

# ---
    def get_fq(self):
        """
        purpose: estimate f/Q for a given suite of observations
        references:
        1) Morel and Gentilli, Applied Optics, 35(24), 1996
        2) Morel et al., Applied Optics, 41(3), 2002
        input:
        1) wavelength 
        2) solar zenith angle [degrees]
        3) chl concentration float [log10(chl) in mg/m3]
        4) nadir zenith angle float [degrees]
        5) azimuth angle float [degrees]
        output: brdf correction factor
        """
        # loop through wavelengths
        npts = len(self.chl.ravel())
        brdf = []
        for wav in self.wav:
            # input array to interpolator is dimensions [npts,nvar]
            wavs = np.ones([npts])*wav
            inputs = [wavs,self.sza.ravel(),self.chl.ravel(),self.thetap.ravel(),self.raa.ravel()]
            inputs = np.array(inputs).T
            fqA = self.interp_func(inputs)

            zeros = np.zeros([npts])
            inputs = [wavs,zeros,self.chl.ravel(),zeros,zeros]
            inputs = np.array(inputs).T
            fq0 = self.interp_func(inputs)
            fq = fq0/fqA
            fq = fq.reshape(self.chl.shape)
            fq.shape = (1,) + fq.shape
            brdf.append(fq)

        self.fq = np.concatenate(brdf)


# ---
    def read_chl(self):
        nc = Dataset(self.nobmFile)
        self.chl = np.squeeze(nc.variables['montot'][:])
        self.chl = np.log10(self.chl)
        self.wav = nc.variables['wavelength'][:]
        nc.close()    
# ---
    def read_geom(self):
        nc = Dataset(self.l1bFile)
        self.sza = nc.groups['geolocation_data'].variables['solar_zenith'][:]
        self.saa = nc.groups['geolocation_data'].variables['solar_azimuth'][:] 
        self.vza = nc.groups['geolocation_data'].variables['sensor_zenith'][:]
        self.vaa = nc.groups['geolocation_data'].variables['sensor_azimuth'][:]
        nc.close()

        # calc expected raa
        self.raa = self.vaa - 180.0 - self.saa
        Il = self.raa < -180.
        Ia = self.raa  > 180.
        self.raa[Il] = self.raa[Il] + 360.0
        self.raa[Ia] = self.raa[Ia] - 360.0
        # fq is symmetric in azimuth
        self.raa = np.abs(self.raa)

        # get thetap
        h2o = 1.34
        self.thetap = np.degrees(np.arcsin(np.sin(np.radians(self.vza)) / h2o))

# ---        
    def writenc(self,zlib=True,verbose=True):
        """
        Write netCDF file of fq values
        """
        # Dimensions
        # ------------------
        ntime = 1

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = "water leaving radiance BRDF correction factors"
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'created from sleave_brdf.py'
        nc.references = 'Morel and Gentilli, Applied Optics, 35(24), 1996, Morel et al., Applied Optics, 41(3), 2002'
        nc.comment = ''
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'

        # Create dimensions
        # -----------------
        nch,nscan,npixel = self.fq.shape
        t = nc.createDimension('time',ntime)
        x = nc.createDimension('ccd_pixels',npixel)
        y = nc.createDimension('number_of_scans',nscan)
        w = nc.createDimension('wavelength',nch)

        nctrj = Dataset(self.nobmFile)

        # Save lon/lat
        # --------------------------
        _copyVar(nctrj,nc,'ccd_pixels',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'number_of_scans',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'longitude',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'latitude',dtype='f4',zlib=zlib,verbose=verbose)
        _copyVar(nctrj,nc,'time', dtype='f4',zlib=zlib,verbose=verbose)

        nctrj.close()

        # Wavelenghts
        dim = ('wavelength',)
        this = nc.createVariable('wavelength','f4',dim,zlib=zlib)
        this.long_name = 'center wavelegnth'
        this.units = 'nm'
        this[:] = self.wav

        # fq
        dim = ('time','wavelength','number_of_scans','ccd_pixels')
        this = nc.createVariable('fq','f4',dim,zlib=zlib)
        this.long_name = 'brdf correction factor'
        this.units = 'None'
        this[:] = self.fq

        nc.close()
if __name__ == "__main__":

    # PACE default
    # -------------------
    morelFile = 'oci_tables/morel_fq.nc' 

    calculon = '/nobackup/PACE'
    nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    if os.path.exists(nccs):
        lbRoot = nccs + '/LevelB'
        l1bRoot = nccs + '/L1B'
        outdir = lbRoot
    elif os.path.exists(calculon):
        lbRoot = calculon + '/LevelB'
        l1bRoot = calculon + '/L1B'
        outdir = lbRoot

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] iso_t1 [iso_t2]",
                          version='1.0.0' )

    parser.add_option("--levelB", dest="lbRoot", default=lbRoot,
              help='PACE LevelB path (default=%s)'%lbRoot)

    parser.add_option("--level1B", dest="l1bRoot", default=l1bRoot,
              help='PACE Level1B path (default=%s)'%l1bRoot)

    parser.add_option("--morel", dest="morelFile", default=morelFile,
                          help='file with morel fq values (default=%s)'%morelFile)

    parser.add_option("-o", "--outdir", dest="outdir", default=outdir,
              help="Output NetCDF file in (default=%s)"\
                          %outdir )

    (options, args) = parser.parse_args()

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


    brdf = SLEAVE_BRDF(t1,t2,options)
