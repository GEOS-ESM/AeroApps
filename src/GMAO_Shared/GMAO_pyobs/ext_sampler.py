#!/usr/bin/env python

"""
    Utility to compute optical properties along a sampled file.

    Ed Nowottnick, January, 2015.

"""

import os
import MieObs_
from types import *
from netCDF4 import Dataset
from mieobs import VNAMES, getAOPext
from numpy import zeros, arange, array, ones, zeros, interp, isnan, ma, NaN, squeeze
from math import pi, sin, cos, asin, acos

from types           import *
from optparse        import OptionParser
from datetime        import *
from dateutil.parser import parse         as isoparser
from csv             import DictReader

from MAPL           import Config, eta
from MAPL.constants import *
from pyobs          import ICARTT

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS004']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless']

class MieCalc(object):
    pass
    """                                                                                  
    Generic container for Variables
    """

#---
def getVars(inFile):
    """
#    Parse input file, create variable dictionary
#    """

    Vars       = MieCalc()
    levs       = MieCalc()
    levUnits   = MieCalc()
    file       = Dataset(inFile)
    names      = file.variables.keys()
    lon        = file.variables['trjLon'][:]
    lat        = file.variables['trjLat'][:]
    tyme       = file.variables['time'][:]
    MIENAMES = names
    for n, name in enumerate(MIENAMES):
        var = file.variables[name]
        size = len(var.shape)
        if size == 2:
            setattr(Vars,name,var[:,:])
        if size == 1:
            setattr(Vars,name,var[:])

    return Vars        

#---
def computeMie(Vars, channel, varnames, rcFile):
    """
#    Computes optical quantities and combines into a dictionary
#   """
    ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(Vars,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=True,rcfile=rcFile)

    MieVars = dict()
    MieVars['ext'] = ext
    MieVars['scatext'] = sca
    MieVars['backscat'] = backscat
    MieVars['aback_sfc'] = aback_sfc
    MieVars['aback_toa'] = aback_toa
    MieVars['depol'] = depol

    return MieVars

#---
def writeNC ( lons, lats, tyme, isotimeIn, MieVars, MieVarsNames,
              MieVarsUnits, inFile, outFile, zlib=False):
#, levs, levUnits, trjFile, options,
#              title='GEOS-5 Trajectory Sampler',
#              doAkBk=False, zlib=False):
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
    ls = nc.createDimension('ls',19)
    if km>0:
        nz = nc.createDimension('lev',km)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since tstart'
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
    
    # Trajectory coordinates
    # ----------------------
    lon = nc.createVariable('trjLon','f4',('time',),zlib=zlib)
    lon.long_name = 'Trajectory Longitude'
    lon.units = 'degrees_east'
    lon[:] = lons[:]
    lat = nc.createVariable('trjLat','f4',('time',),zlib=zlib)
    lat.long_name = 'Trajectory Latitude'
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
        if size == 2:
            dim = ('time','lev')
        if size == 1:
            dim = ('time')
        this = nc.createVariable(name,'f4',dim,zlib=zlib)
        this.standard_name = name
        this.units = MieVarsUnits[n]
        this.missing_value = MAPL_UNDEF
        this[:] = var

    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    inFile  = 'trj_sampler.nc'
    outFile = 'ext_sampler.nc'
    channel = (532,)
    rcFile = 'Aod3d_532nm.rc'

#   Parse command line options
#   --------------------------
    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inFile", default=inFile,
              help="Sampled input file")

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output file containing optical properties")

    parser.add_option("-r", "--rc", dest="rcFile", default=rcFile,
              help="Resource file pointing to optical tables")

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT or EXCEL (default=%s)"%format )

    parser.add_option("-c", "--channel", dest="channel", default=channel,
              help="Channel for Mie calculation")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode")

    parser.add_option("--du",
                      action="store_true", dest="dust",
                      help="Dust Only")

    parser.add_option("--ss",
                      action="store_true", dest="seasalt",
                      help="Seasalt Only")

    parser.add_option("--su",
                      action="store_true", dest="sulfate",
                      help="Sulfate Only")

    parser.add_option("--bc",
                      action="store_true", dest="bcarbon",
                      help="Black Carbon Only")

    parser.add_option("--oc",
                      action="store_true", dest="ocarbon",
                      help="Organic Carbon Only")

    (options, args) = parser.parse_args()

         
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if ext.upper() == '.XLS':
        options.format = 'EXCEL'
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    elif 'EXCEL' in options.format:
        options.outFile = name + '.xls'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
 

    # Get Variables
    # --------------------------
    Vars = getVars(options.inFile)

    # Run Mie Calculator and Write Output Files
    # --------------------------
    channelIn = float(options.channel)
    outFile = 'sampled.aerosol.optical.properties.'+options.channel+'nm.total.nc'
    MieVars = computeMie(Vars,channelIn,VNAMES,rcFile)
    writeNC(Vars.trjLon,Vars.trjLat,Vars.time,Vars.isotime,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.dust:
        outFile = 'sampled.aerosol.optical.properties.'+options.channel+'nm.dust.nc'
        MieVars = computeMie(Vars,channelIn,VNAMES_DU,rcFile)
        writeNC(Vars.trjLon,Vars.trjLat,Vars.time,Vars.isotime,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.seasalt:
        outFile = 'sampled.aerosol.optical.properties.'+options.channel+'nm.ss.nc'
        MieVars = computeMie(Vars,channelIn,VNAMES_SS,rcFile)
        writeNC(Vars.trjLon,Vars.trjLat,Vars.time,Vars.isotime,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.sulfate:
        outFile = 'sampled.aerosol.optical.properties.'+options.channel+'nm.su.nc'
        MieVars = computeMie(Vars,channelIn,VNAMES_SU,rcFile)
        writeNC(Vars.trjLon,Vars.trjLat,Vars.time,Vars.isotime,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.bcarbon:
        outFile = 'sampled.aerosol.optical.properties.'+options.channel+'nm.bc.nc'
        MieVars = computeMie(Vars,channelIn,VNAMES_BC,rcFile)
        writeNC(Vars.trjLon,Vars.trjLat,Vars.time,Vars.isotime,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)

    if options.ocarbon:
        outFile = 'sampled.aerosol.optical.properties.'+options.channel+'nm.oc.nc'
        MieVars = computeMie(Vars,channelIn,VNAMES_OC,rcFile)
        writeNC(Vars.trjLon,Vars.trjLat,Vars.time,Vars.isotime,MieVars,MieVarsNames,MieVarsUnits,options.inFile,outFile)
   
