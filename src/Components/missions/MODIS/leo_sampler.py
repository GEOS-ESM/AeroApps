#!/usr/bin/env python

"""
    Utility to create GEOS-5 Collections on MODIS Level 2 granule.
    Based on trj_sampler.py and geo_sampler.py

    Arlindo da Silva, November 2014.
    P. Castellanos, Aug 2016 - Modified from GEO to LEO
"""

import os

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from   numpy import zeros, ones, arange, array, tile
import numpy as np

from netCDF4 import Dataset

from MAPL           import Config, eta
from MAPL.constants import *

from pyobs.nc4ctl   import NC4ctl  

from MAPL.ShaveMantissa_ import shave32
from glob                import glob

from pyhdf.SD            import SD, HDF4Error
MISSING = 1E15

SDS = ('Longitude', 'Latitude', 'Scan_Start_Time' )
DATE_START = datetime(1993,1,1,0,0,0)

class SampleVar(object):
    """                                                                                  
    Generic container for Variables
    """
    def __init__(self,name):
        self.name = name

class MODIS(object):
    """
    Generic container for MODIS SDS
    """

    def __init__(self, Path,verb=False):

        # Read each granule, appending them to the list
        # ---------------------------------------------
        if type(Path) is list:
            if len(Path) == 0:
               self.nobs = 0
               print "WARNING: Empty MxD04_L2 object created"
               return
            else:
                self.nobs = len(Path)
        else:
           Path = [Path, ]
           self.nobs = 1

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        self.verb = verb
        self.SDS  = SDS

        for name in self.SDS:
           self.__dict__[name] = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        self._readList(Path)

        # Create corresponding python time
        # --------------------------------
        self.tyme = []        
        for start_time in self.Scan_Start_Time:
            tyme = np.ma.masked_array([DATE_START+timedelta(seconds=s) for s in array(start_time).ravel()]).reshape(start_time.shape)
            tyme.mask = start_time.mask
            self.tyme.append(tyme)
                    

    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readGranule(item)
            else:
                print "%s is not a valid file or directory, ignoring it"%item

#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print "%s is not a valid file or directory, ignoring it"%item

#---
    def _readGranule(self,filename):
        """Reads one MOD04/MYD04 granule with Level 2 aerosol data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 2:
                print "- %s: not recognized as an HDF file"%filename
            return 

        # Read select variables (do not reshape)
        # --------------------------------------
        for sds_ in self.SDS:
            sds = sds_
            try:
                v = hfile.select(sds).get()
            except:
                if sds in NEW_SDS: # cope with new names in Coll. 6
                    sds = NEW_SDS[sds_]
                    v = hfile.select(sds).get()
            a = hfile.select(sds).attributes()
            fill = a['_FillValue']
            if sds == 'Scan_Start_Time':
                fill = -9999  #error in MODIS files
            v = np.ma.masked_array(v,fill_value=-999)
            v.mask = np.abs(v-fill) < 0.001            
            if a['scale_factor']!=1.0 or a['add_offset']!=0.0:
                v = a['scale_factor'] * v + a['add_offset']

            self.__dict__[sds_].append(v) # Keep Collection 5 names!

            

class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation


# ---
def granules ( path, prod, t1, t2, coll='006'):
    """
    Returns a list of MxD04 granules for a given product at given time.
    On input,

    path      ---  mounting point for the MxD04 Level 2 files
    prod      ---  either MOD04 or MYD04
    t1        ---  starting time (timedate format)
    t2        ---  ending time (timedate format)

    coll      ---  collection: 005, 051 (optional)

    """

    # Find MODIS granules in the time range
    # ------------------------------------------
    dt = timedelta(minutes=5)
    t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
    Granules = []
    while t < t2:
        if t >= t1:
            doy = t.timetuple()[7]
            basen = "%s/%s/%04d/%03d/%s_L2.A%04d%03d.%02d%02d.%s.*.hdf"\
                     %(path,prod,t.year,doy,prod,t.year,doy,t.hour,t.minute,coll)
            
            try:
                filen = glob(basen)[0]
                Granules += [filen,]
#               print " [x] Found "+filen
            except:
                pass
        t += dt

    if len(Granules) == 0:
        print "WARNING: no %s collection %s granules found for %s through %s"%(prod,coll, str(t1), str(t2))

    return Granules


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
        var = SampleVar(v)
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

#----
def getVars(rcFile):
    """
    Parse reource file, create variable dictionary with relevant
    metadata.
    """
    cf = Config(rcFile)
    Vars = dict()
    levUnits = 'none'
    levs = []
    for V in cf.keys():
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
def _copyVar(ncIn,ncOut,name,dtype='f4',zlib=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
    y.long_name = x.long_name
    y.units = x.units 
    try:
        y.missing_value = x.missing_value
    except:
        pass
    rank = len(x.shape)
    if rank == 1:
        y[:] = x[:]
    elif rank == 2:
        y[:,:] = x[:,:]
    elif rank == 3:
        y[:,:,:] = x[:,:,:]
    else:
        raise ValueError, "invalid rank of <%s>: %d"%(name,rank)

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
        raise ValueError, "invalid rank=%d"%rank

    # Shave it
    # --------
    qs, rc = shave32(q.ravel(),xbits,has_undef,undef,chunksize)
    if rc:
        raise ValueError, "error on return from shave32, rc=%d"%rc

    return qs.reshape(shp)

#---
def writeNC ( ncGeo, clon, clat, ctyme, tBeg, Vars, levs, levUnits, options,
              xchunk=150, ychunk=200, zchunk=1,
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
    nc = Dataset(options.outFile,'w',format=options.format)

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
    _copyVar(ncGeo,nc,u'ew',dtype='f4',zlib=False)
    _copyVar(ncGeo,nc,u'ns',dtype='f4',zlib=False)
    dt = nc.createVariable('scanTime','f4',('ew',),zlib=False)
    dt.long_name = 'Time of Scan'
    dt.units = 'seconds since %s'%tBeg.isoformat(' ')
    DT = ctyme[0,:] - tBeg
    dt[:] = array([ds.total_seconds() for ds in DT])

    
    # Save lon/lat if so desired
    # --------------------------
    if options.coords:
        _copyVar(ncGeo,nc,u'clon',dtype='f4',zlib=False)
        _copyVar(ncGeo,nc,u'clat',dtype='f4',zlib=False)
        
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    for path in Vars:
        if options.verbose:
            print " <> opening "+path
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
                print " [] Interpolating <%s>"%name.upper()

            # Use NC4ctl for linear interpolation
            # -----------------------------------
            I = (clon<0.1*MISSING)&(clat<0.1*MISSING)
            Z = g.nc4.sample(name,clon[I],clat[I],ctyme[I],
                             Transpose=False,squeeze=True,Verbose=options.verbose)
            if options.verbose: print " <> Writing <%s> "%name
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
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    title   = 'GEOS-5 MODIS Sampler'
    format  = 'NETCDF4_CLASSIC'
    rcFile  = 'modis_sampler.rc'
    outFile = 'modis_sampler.nc4'
    coll    = '006'

    # MODIS Level 2 default
    # -------------------
    calculon = '/nobackup/MODIS/{}/Level2/'.format(coll)
    nccs = '/discover/nobackup/projects/gmao/iesa/aerosol/Data/MODIS/Level2/{}/'.format(coll)
    if os.path.exists(calculon): 
        L2Root = calculon
    elif os.path.exists(nccs): 
        L2Root = nccs
    else:
        L2Root = './Level2/'

    dt_hours = 24 # hard-wire one-day for now

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] prod iso_t1 [iso_t2]",
                          version='1.0.0' )

    parser.add_option("-C", "--coords",
                      action="store_true", dest="coords",
                      help="Save 2D lon/lat coordinates on output file.")

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-p", "--path", dest="L2Root", default=L2Root,
              help='MODIS Level2 path (default=%s)'%L2Root)

    parser.add_option("-c", "--coll", dest="coll", default=coll,
              help='MODIS collection code (default=%s)'%coll)    

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

    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        prod, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 2:
        prod, iso_t1 = args
        iso_t2 = None
        t1     = isoparser(iso_t1)
        t2     = t1 + timedelta(hours=dt_hours)
    else:
        parser.error("must have 2 or 3 arguments: prod iso_t1 [iso_t2]")

    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
    options.zlib = not options.nozip

    # Get Granules to work on
    # -----------------------
    options.granules = granules ( L2Root, prod, t1, t2, coll=options.coll)

    # Create (x,y,t) coordinates
    # --------------------------
    mxd = MODIS(options.granules)

    # # Get Variables and Metadata
    # # --------------------------
    # Vars, levs, levUnits = getVars(options.rcFile)


    # # Write output file
    # # -----------------
    # writeNC ( mxd, Vars, levs, levUnits, options,
    #           doAkBk=False)

    # # All done
    # # --------
    # ncGeo.close()
