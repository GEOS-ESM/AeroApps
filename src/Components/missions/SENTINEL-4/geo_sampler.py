#!/usr/bin/env python3

"""
    Utility to GEOS-5 Collections on TEMPO's domain.
    Based on trj_sampler.py.

    Arlindo da Silva, November 2014.

    Adapted for GOES-R 
    P. Castellanos, Sept 2015

    Adapted for GEMS
    P. Castellanos, Nov 2015

    Adapted for SENTINEL-4
    P. Castellanos, Dec 2015    
"""

import os

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from numpy import zeros, ones, arange, array, tile, empty, ma, where

from netCDF4 import Dataset

from MAPL           import Config, eta
from MAPL.constants import *

from pyobs.nc4ctl   import NC4ctl  

from MAPL.ShaveMantissa_ import shave32

MISSING = 1E15

class TempoVar(object):
    """                                                                                  
    Generic container for Variables
    """
    def __init__(self,name):
        self.name = name

class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation

#---
def Open(filename,doNC4ctl=True):
    """
    Uses GFIO or GFIOctl to open either a NetCDF-4 or a control file.
    Very heuristic.

    TO DO: Remove GFIO dependencies, using NC4ctl instead.

    """
    from gfio import GFIOurl, GFIOctl
    name, ext = os.path.splitext(filename)

    # Open the GFIO dataset
    # ---------------------
    f = None
    if 'HTTP://' == name[:7].upper():
        f = GFIOurl(filename)
        f.lower = True # force variable names to be lower case when sampling.
    elif ext.upper() in ('.NC4','.NC','.HDF','.H5'):
        f = GFIOurl(filename)
        f.lower = False
    else:
        f = GFIOctl(filename)
        f.lower = False

    # Create variable dictionary
    # --------------------------
    Vars = dict()
    for i in range(len(f.vname)):
        if f.lower:
            v = f.vname[i].upper()
        else:
            v = f.vname[i]
        var = TempoVar(v)
        #var.title = f.vtitle[i]
        var.km = f.kmvar[i]
        if var.km>0:
            var.levunits = f.levunits[:]
            var.levs = f.levs[:]    
        try:
            var.units = f.vunits[i]
        except:
            var.units = 'unknown'  # opendap currently lacks units
        Vars[v] = var

    f.Vars = Vars

    if doNC4ctl:
        f.nc4 = NC4ctl_(filename)

    return f

#---
def TEMPOscanTimes(nNS,nEW,tBeg):
    """
    Calculate TEMPO scantimes. 
    """
    tScan = array([tBeg + i * timedelta(seconds=2.85) for i in arange(nEW) ])
    tScan = tScan[-1::-1] # flip times
    return tile(tScan,(nNS,1))

#---
def SENTINEL4scanTimes(nNS,nEW,tBeg):
    """
    Calculate SENTINEL-4 scantimes. 
    """
    s     = 30*60./nEW
    tScan = array([tBeg + i * timedelta(seconds=s) for i in arange(nEW) ])
    tScan = tScan[-1::-1] # flip times
   
    return tile(tScan,(nNS,1))

#---
def GEMSscanTimes(nNS,nEW,tBeg):
    """
    Calculate GEMS scantimes. 
    """
    s     = 30*60./nEW
    tScan = array([tBeg + i * timedelta(seconds=s) for i in arange(nEW) ])
    tScan[int(0.5*nEW):] = tScan[int(0.5*nEW):] + timedelta(seconds=30.*60)
    tScan = tScan[-1::-1] # flip times
   
    return tile(tScan,(nNS,1))

#---
def GOESRscanTimes(scanTimeFile,nNS,nEW,tBeg):
    """
    Read in GOES-R scantimes - have already been synchronized with TEMPO. 
    """
    nc = Dataset(scanTimeFile)   

    seconds = nc.variables['scanTime'][:]
    MISSING = nc.variables['scanTime'].missing_value

    seconds[where(abs(seconds-MISSING)/MISSING < 0.1)] = ma.masked

    tScan = empty([nNS, nEW],dtype='object')
    for j in arange(nNS):
        tScan[j,:] = array([tBeg + timedelta(seconds=s) for s in seconds[j,:]])
    return tScan

#---
def getCoords(geoFile,verbose):
    """
    Get coordinates from geo-location file.
    Assumes the file has albed been tighten in the E-W domain
    """
    if verbose:
       print(" <> Getting GEO coordinates from ", geoFile)
    nc = Dataset(geoFile)
    lon = nc.variables['clon'][:,:]
    lat = nc.variables['clat'][:,:]
    missing = nc.variables['clon'].missing_value
    return (nc,lon,lat,missing)

def getVars(rcFile):
    """
    Parse reource file, create variable dictionary with relevant
    metadata.
    """
    cf = Config(rcFile)
    Vars = dict()
    levUnits = 'none'
    levs = []
    for V in list(cf.keys()):
        path = cf(V)
        f = Open(path)
        varList = []
        for v in V.split(','):
            var = f.Vars[v.strip()]
            if var.km>0:
                levUnits = var.levunits
                levs = var.levs
            varList += [var,]
        Vars[path] = varList
        
    return (Vars, levs, levUnits)

#----
def _copyVar(ncIn,ncOut,name,dtype='f4',layout=None,zlib=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    print(x.dimensions)
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
    y.long_name = x.long_name
    y.units = x.units 
    try:
        y.missing_value = x.missing_value
    except:
        pass
    rank = len(x.shape)

    if layout is not None:
        nX = int(layout[0])
        nY = int(layout[1])
        i  = int(layout[2])

        i_ = i%nX
        j_ = int(i/nX) 
        print('_copyVar',nX,nY,i,i_,j_)


    if rank == 1:
        if layout is None:
            y[:] = x[:]
        else:
            nEW, = x.shape
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            y[:] = x[Xstart:Xend]
    elif rank == 2:
        if layout is None:
            y[:,:] = x[:,:]
        else:
            nNS, nEW = x.shape 
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            Ystart = j_*nNS/nY
            Yend   = Ystart + nNS/nY

            y[:,:] = x[Ystart:Yend,Xstart:Xend]
    elif rank == 3:
        if layout is None:
            y[:,:,:] = x[:,:,:]
        else:
            nNS, nEW = x.shape 
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            Ystart = j_*nNS/nY
            Yend   = Ystart + nNS/nY

            y[:,:,:] = x[:,Ystart:Yend,Xstart:Xend]
    else:
        raise ValueError("invalid rank of <%s>: %d"%(name,rank))

#---
def shave(q,options,undef=MISSING,has_undef=1,nbits=12):
    """
    Shave variable. On input, nbits is the number of mantissa bits to keep
    out of maximum of 24.
    """
    # no compression, no shave
    # ------------------------
    if not options.zlib:
        return q

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
        raise ValueError("invalid rank=%d"%rank)

    # Shave it
    # --------
    qs, rc = shave32(q.ravel(),xbits,has_undef,undef,chunksize)
    if rc:
        raise ValueError("error on return from shave32, rc=%d"%rc)

    return qs.reshape(shp)

#---
def writeNC ( ncGeo, clon, clat, ctyme, tBeg, Vars, levs, levUnits, options, 
              layout=None, xchunk=150, ychunk=200, zchunk=1,
              doAkBk=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables pn geostationary grid
    described by (clon,clat,tyme).
    """

    # Gridded Dimensions
    # ------------------
    km = len(levs)
    nNS, nEW = clon.shape

    # Open NC file
    # ------------
    if layout is None:
        nc = Dataset(options.outFile,'w',format=options.format)
    else:
        nc = Dataset(options.outFile[0:-4]+'_'+layout+'.nc4','w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.title = options.title
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from GEOS-5 standard collections by tempo_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains GEOS-5 parameters sampled on a Geostationary grid'
    nc.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    nc.Conventions = 'CF'
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',1) # one time per file for now
    if km>0:
        nz = nc.createDimension('lev',km)
        if doAkBk:
            ne = nc.createDimension('ne',km+1)
    x = nc.createDimension('ew',nEW)
    y = nc.createDimension('ns',nNS)

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=False)
    time.long_name = 'Initial Time of Scan'
    time.units = 'seconds since %s'%tBeg.isoformat(' ')
    time[0] = 0
    if km > 0: # pressure level not supported yet
        lev = nc.createVariable('lev','f4',('lev',),zlib=False)
        lev.long_name = 'Vertical Level'
        lev.units = levUnits.strip()
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = levs[:]

        if doAkBk:
            ae, be = eta.getEdge(km) # Coefficients for Hybrid coordinates
            ak = nc.createVariable('ak','f4',('ne',),zlib=False)
            ak.long_name = 'Eta coordinate coefficient ak (p = ak + bk * ps)'
            ak.units = 'Pa'
            ak = ae[:]
            bk = nc.createVariable('bk','f4',('ne',),zlib=False)
            bk.long_name = 'Eta coordinate coefficient bk (p = ak + bk * ps)'
            bk.units = '1'
            bk = be[:]
    
    # Add pseudo dimensions for GrADS compatibility
    # -------------------------------------------
    if layout is not None:
        nX = int(layout[0])
        nY = int(layout[1])
        i  = int(layout[2])

        i_ = i%nX
        j_ = int(i/nX) 
        _copyVar(ncGeo,nc,'ew',dtype='f4',zlib=False, layout = layout[0]+'1'+str(i_))
        _copyVar(ncGeo,nc,'ns',dtype='f4',zlib=False, layout = layout[1]+'1'+str(j_))
    else:
         _copyVar(ncGeo,nc,'ew',dtype='f4',zlib=False)
         _copyVar(ncGeo,nc,'ns',dtype='f4',zlib=False)
    dt = nc.createVariable('scanTime','f4',('ew',),zlib=False)
    dt.long_name = 'Time of Scan'
    dt.units = 'seconds since %s'%tBeg.isoformat(' ')
    DT = ctyme[0,:] - tBeg
    dt[:] = array([ds.total_seconds() for ds in DT])

    
    # Save lon/lat if so desired
    # --------------------------
    if options.coords:
        _copyVar(ncGeo,nc,'clon',dtype='f4',zlib=False, layout=layout)
        _copyVar(ncGeo,nc,'clat',dtype='f4',zlib=False, layout=layout)            
        
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    for path in Vars:
        if options.verbose:
            print(" <> opening "+path)
        g = Open(path)
        for var in Vars[path]:
            if var.km == 0:
                dim = ('time','ns','ew')
                chunks = (1, ychunk, xchunk)
                W = MISSING * ones((nNS,nEW))
            else:
                dim = ('time','lev','ns','ew')
                chunks = (1,zchunk,ychunk, xchunk)
                W = MISSING * ones((var.km,nNS,nEW))
            rank = len(dim)
            this = nc.createVariable(var.name,'f4',dim,
                                     zlib=options.zlib,
                                     chunksizes=chunks)

            #this.standard_name = var.title
            this.standard_name = var.name
            #this.long_name = var.title.replace('_',' ')
            this.long_name = ''
            this.missing_value = MAPL_UNDEF
            this.units = var.units
            if g.lower:
                name = var.name.lower() # GDS always uses lower case
            else:
                name = var.name
            if options.verbose:
                print(" [] Interpolating <%s>"%name.upper())

            # Use NC4ctl for linear interpolation
            # -----------------------------------
            I = (clon<0.1*MISSING)&(clat<0.1*MISSING)
            Z = g.nc4.sample(name,clon[I],clat[I],ctyme[I],
                             Transpose=False,squeeze=True,Verbose=options.verbose)
            if options.verbose: print(" <> Writing <%s> "%name)
            if rank == 3:
               W[I] = Z
               this[0,:,:] = shave(W[:,:],options)
            elif rank == 4:
               W[:,I] = Z
               this[0,:,:,:] = shave(W[:,:,:],options)

    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print(" <> wrote %s file %s"%(options.format,options.outFile))
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    inst   = 'TEMPO'
    title  = 'GEOS-5 TEMPO Sampler'
    format = 'NETCDF4_CLASSIC'
    rcFile  = 'tempo_sampler.rc'
    outFile = 'tempo_sampler.nc4'

    inst   = 'GOES-R'
    title  = 'GEOS-5 GOES-R Sampler'
    format = 'NETCDF4_CLASSIC'
    rcFile  = 'geo_sampler_omi_ler.rc'
    outFile = 'goes-r_sampler.nc4'
    scanTimeFile = 'goes-r.lg1.scantime.nc4'
    layout       = '41'

    inst   = 'GEMS'
    title  = 'GEOS-5 GEMS Sampler'
    format = 'NETCDF4_CLASSIC'
    rcFile  = 'geo_sampler_omi_ler.rc'
    outFile = 'gems_sampler.nc4'
    layout       = None   

    inst   = 'SENTINEL-4'
    title  = 'GEOS-5 SENTINEL-4 Sampler'
    format = 'NETCDF4_CLASSIC'
    rcFile  = 'geo_sampler_omi_ler.rc'
    outFile = 'sentinel-4_sampler.nc4'
    layout  = None       

    # Geolocation default
    # -------------------
    calculon = '/nobackup/'+inst+'/LevelG/invariant/'
    nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'+inst+'/DATA/LevelG/invariant/'
    if os.path.exists(calculon): 
        geoFile      = calculon + inst.lower()+'.lg1.invariant.nc4'
        scanTimeFile = calculon + inst.lower()+'.lg1.scantime.nc4'
    elif os.path.exists(nccs): 
        geoFile      = nccs  +inst.lower()+'.lg1.invariant.nc4'
        scanTimeFile = nccs + inst.lower()+'.lg1.scantime.nc4'
    else:
        geoFile      = inst.lower()+'.lg1.invariant.nc4'
        scanTimeFile = inst.lower()+'.lg1.scantime.nc4'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] iso_t1 [iso_t2]",
                          version='1.0.0' )

    parser.add_option("-C", "--coords",
                      action="store_true", dest="coords",
                      help="Save 2D lon/lat coordinates on output file.")

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-g", "--geolocation", dest="geoFile", default=geoFile,
              help='Level G1 geo-location file (default=%s)'%geoFile)

    parser.add_option("-N", "--nozip",
                      action="store_true", dest="nozip",
                      help="Do not shave/compress output files.")

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-r", "--rcFile", dest="rcFile", default=rcFile,
              help="Resource file defining parameters to sample (default=%s)"\
                          %rcFile )

    parser.add_option("-t", "--title", dest="title", default=title,
              help="Output file title, typically the collection name (default=%s)"\
                          %title )

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    parser.add_option("-i", "--instrument",
                  dest="inst", default=inst, help="Instrument name (default=%s)"\
                          %inst)

    parser.add_option("-T", "--timefile",
                  dest="scanTimeFile", default=scanTimeFile, 
                  help="Scan times file name (default=None)") 
    parser.add_option("-l", "--layout",
              dest="layout", default=None, 
              help="layout to break up domain into tiles (default=None)")    

    (options, args) = parser.parse_args()
    
    options.layout = layout
    if len(args) == 2:
        iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 1:
        iso_t1, iso_t2 = (args[0], None)
        t1, t2 = (isoparser(iso_t1), None)
    else:
        parser.error("must have 1 or 2 arguments: iso_t1 [iso_t2]")

    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError('invalid extension <%s>'%ext)
    options.zlib = not options.nozip

    # Get Variables and Metadata
    # --------------------------
    Vars, levs, levUnits = getVars(options.rcFile)

    # Create (x,y,t) coordinates
    # --------------------------
    ncGeo, clon, clat, MISSING = getCoords(options.geoFile,options.verbose)
    nNS, nEW = clon.shape
    
    if (inst == 'TEMPO'):
        ctyme = TEMPOscanTimes(nNS,nEW,t1)
    if (inst == 'GEMS'):
        ctyme = GEMSscanTimes(nNS,nEW,t1)  
    if (inst == 'SENTINEL-4'):
        ctyme = SENTINEL4scanTimes(nNS,nEW,t1)                
    if (inst == 'GOES-R'):
        ctyme = GOESRscanTimes(options.scanTimeFile,nNS,nEW,t1)


    if options.layout is not None:
        nX = int(options.layout[0])
        nY = int(options.layout[1])
        options.coords = True
        for i in arange(nX*nY):
            # Write output file
            # -----------------
            i_ = i%nX
            j_ = int(i/nX) 

            print(i_, j_)
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            Ystart = j_*nNS/nY
            Yend   = Ystart + nNS/nY
            print(Xstart, Xend, Ystart, Yend)
            writeNC ( ncGeo, clon[Ystart:Yend,Xstart:Xend], clat[Ystart:Yend,Xstart:Xend], 
                      ctyme[Ystart:Yend,Xstart:Xend], t1, Vars, levs, levUnits, options, 
                      layout = options.layout+str(i), doAkBk=False )
    else:
        writeNC ( ncGeo, clon, clat, ctyme, t1, Vars, levs, levUnits, options,
                  doAkBk=False )

    # All done
    # --------
    ncGeo.close()
