#!/usr/bin/env python

"""
    Utility to create GEOS-5 Collections on PACE L1B granule.
    Based on trj_sampler.py and geo_sampler.py

    Arlindo da Silva, November 2014.
    P. Castellanos, Aug 2016 - Modified from GEO to LEO
    P. Castellanos, 2018 - Modified from MODIS to PACE
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


SDS = {'longitude'  : 'geolocation_data',
       'latitude'   : 'geolocation_data',
       'ev_mid_time': 'scan_line_attributes' }


ALIAS = {'longitude': 'lon',
         'latitude' : 'lat',
         'ev_mid_time': 'midTime'}


class SampleVar(object):
    """                                                                                  
    Generic container for Variables
    """
    def __init__(self,name):
        self.name = name

class PACE(object):
    """
    Generic container for PACE SDS
    """

    def __init__(self, Path,verb=False):

        # Read each granule, appending them to the list
        # ---------------------------------------------
        if type(Path) is list:
            if len(Path) == 0:
               self.nobs = 0
               print "WARNING: Empty PACE_L1B object created"
               return
            else:
                self.nobs = len(Path)
        else:
           Path = [Path, ]
           self.nobs = 1

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        self.granules = Path
        self.verb = verb
        self.SDS  = SDS.keys()
        self.SDSg = SDS

        for name in self.SDS:
           self.__dict__[name] = []

        self.scanStart = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        self._readList(Path)

        # Alias
        for sds in self.SDS:
            self.__dict__[ALIAS[sds]] = self.__dict__[sds]

        # Create corresponding python time
        # --------------------------------
        nscan,npixel = self.lon[0].shape
        self.tyme = []        
        for scanStart,midTime,lon in zip(self.scanStart,self.midTime,self.lon):
            scanStart  = isoparser(scanStart.strftime('2006-%m-%dT00:00:00'))
            tyme       = np.array([scanStart + timedelta(seconds=t) for t in midTime])    
            tyme.shape = (nscan,1)
            tyme       = np.repeat(tyme,npixel,axis=1)
            tyme       = np.ma.array(tyme)
            tyme.mask  = lon.mask
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
        """Reads one PACE granule."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            nc = Dataset(filename)
        except:
            if self.verb > 2:
                print "- %s: not recognized as an netCDF file"%filename
            return 

        # Read select variables (do not reshape)
        # --------------------------------------
        self.scanStart.append(isoparser(nc.time_coverage_start))
        for sds in self.SDS:            
            group = self.SDSg[sds]
            v = nc.groups[group].variables[sds][:]
            if not hasattr(v,'mask'):
                v = np.ma.array(v)
                v.mask = np.zeros(v.shape).astype('bool')
            self.__dict__[sds].append(v) 

            

class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation


# ---
def granules ( path, t1, t2):
    """
    Returns a list of PACE granules for a given product at given time.
    On input,

    path      ---  mounting point for the PACE files
    t1        ---  starting time (timedate format)
    t2        ---  ending time (timedate format)
    """

    # Find MODIS granules in the time range
    # ------------------------------------------
    dt = timedelta(minutes=5)
    t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
    Granules = []
    while t < t2:
        if t >= t1:
            doy = t.timetuple()[7]
            basen = "%s/Y%04d/M%02d/D%02d/OCI%04d%03d%02d%02d00.L1B_PACE.nc"\
                     %(path,t.year,t.month,t.day,t.year,doy,t.hour,t.minute)
            
            try:
                filen = glob(basen)[0]
                Granules += [filen,]
#               print " [x] Found "+filen
            except:
                pass
        t += dt

    if len(Granules) == 0:
        print "WARNING: no PACE granules found for %s through %s"%(str(t1), str(t2))

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
def _copyVar(ncIn,ncOut,name,group,dtype='f4',zlib=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.group[group[name]].variables[name]
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
def shave(q,options,undef=MAPL_UNDEF,has_undef=1,nbits=12):
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
def writeNC ( pace, Vars, levs, levUnits, options,
              xchunk=50, ychunk=100, zchunk=1,
              doAkBk=False):


    """
    Write a NetCDF file with sampled GEOS-5 variables on PACE grid
    described by (clon,clat,tyme).
    """
    for i, path in enumerate(pace.granules):
        # Gridded Dimensions
        # ------------------
        km = len(levs)
        nAtrack, nXtrack = pace.longitude[i].shape

        scanStart = pace.scanStart[i]
        scanStart  = isoparser(scanStart.strftime('2006-%m-%dT%H:%M:00'))
        year = scanStart.year
        doy  = scanStart.strftime('%j')
        month = scanStart.strftime('%m')
        day = scanStart.strftime('%d')
        hhmmss = scanStart.strftime('%H%M00')

        # Root name for outfile
        # -------------------------
        outdir = options.outdir + '/Y{}/M{}/D{}'.format(year,month,day)
        if not os.path.exists(outdir):
            os.makedirs(outdir)


        coll = options.rcFile[:-3]
        outFile = '{}/pace-g5nr.lb.{}{}_{}.{}'.format(outdir,coll,year,doy,hhmmss,options.ext)

        # Open NC file
        # ------------
        nc = Dataset(outFile,'w',format=options.format)
        ncIn = Dataset(path)

        # Set global attributes
        # ---------------------
        nc.title = options.title
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from GEOS-5 standard collections by leo_sampler.py'
        nc.references = 'n/a'
        nc.comment = 'This file contains GEOS-5 parameters sampled on a PACE granule'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.time_coverage_start = scanStart.isoformat()
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',1) # one 5-minute granule per file for now
        if km>0:
            nz = nc.createDimension('lev',km)
            if doAkBk:
                ne = nc.createDimension('ne',km+1)
        x = nc.createDimension('ccd_pixels',nXtrack)
        y = nc.createDimension('number_of_scans',nAtrack)

        # Save lon/lat if so desired
        # --------------------------
        if options.coords:
            for sds in pace.SDS:
                _copyVar(ncIn,nc,sds,pace.SDSg,dtype='f4',zlib=False)

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
        ew = nc.createVariable('ew','f4',('ccd_pixels',),
                                fill_value=MAPL_UNDEF,zlib=False)
        ew.long_name    = 'pseudo longitude'
        ew.units        = 'degrees_east'
        ew[:]           = pace.longitude[i][int(nAtrack*0.5),:]

        ns = nc.createVariable('ns','f4',('number_of_scans',),
                                fill_value=MAPL_UNDEF,zlib=False)
        ns.long_name    = 'pseudo latitude'
        ns.units        = 'degrees_north'
        ns[:]           = pace.latitude[i][:,int(nXtrack*0.5)]


        # Loop over datasets, sample and write each variable
        # --------------------------------------------------
        for path in Vars:
            if options.verbose:
                print " <> opening "+path
            g = Open(path)
            for var in Vars[path]:
                if var.km == 0:
                    dim = ('time','number_of_scans','ccd_pixels')
                    chunks = (1, ychunk, xchunk)
                    W = MAPL_UNDEF * ones((nAtrack,nXtrack))
                else:
                    dim = ('time','lev','number_of_scans','ccd_pixels')
                    chunks = (1,zchunk,ychunk, xchunk)
                    W = MAPL_UNDEF * ones((var.km,nAtrack,nXtrack))
                rank = len(dim)
                this = nc.createVariable(var.name,'f4',dim,
                                         zlib=options.zlib,
                                         chunksizes=chunks, fill_value=MAPL_UNDEF)

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
                I = (~pace.longitude[i].mask)&(~pace.latitude[i].mask)&(~pace.tyme[i].mask)
                Z = g.nc4.sample(name,array(pace.longitude[i][I]),array(pace.latitude[i][I]),array(pace.tyme[i][I]),
                                 Transpose=False,squeeze=True,Verbose=options.verbose)
                if options.verbose: print " <> Writing <%s> "%name
                if rank == 3:
                   W[I] = Z
                   W = np.ma.masked_array(shave(W[:,:],options))
                   W.mask = W>0.1*MAPL_UNDEF
                   this[0,:,:] = W
                elif rank == 4:
                   W[:,I] = Z
                   W = np.ma.masked_array(shave(W[:,:,:],options))
                   W.mask = W>0.1*MAPL_UNDEF
                   this[0,:,:,:] = W

        # Close the file
        # --------------
        nc.close()

        if options.verbose:
            print " <> wrote %s file %s"%(options.format,outFile)
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    title   = 'GEOS-5 PACE Sampler'
    format  = 'NETCDF4_CLASSIC'
    rcFile  = 'leo_sampler.rc'

    # MODIS Level 2 default
    # -------------------
    calculon = '/nobackup/3/pcastell/PACE/L1B'
    nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    if os.path.exists(nccs): 
        L1Root = nccs + '/L1B'
        outdir = nccs + '/LevelB'
    elif os.path.exists(calculon): 
        L1Root = calculon + '/L1B'
        outdir = calculon + '/LevelB'


    dt_minutes = 5 #  5 minutes for now

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] iso_t1 [iso_t2]",
                          version='1.0.0' )

    parser.add_option("-C", "--coords",
                      action="store_true", dest="coords",
                      help="Save 2D lon/lat coordinates on output file.")

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-p", "--path", dest="L1Root", default=L1Root,
              help='PACE L1B path (default=%s)'%L1Root)    

    parser.add_option("-N", "--nozip",
                      action="store_true", dest="nozip",
                      help="Do not shave/compress output files.")

    parser.add_option("-o", "--outdir", dest="outdir", default=outdir,
              help="Output NetCDF file in (default=%s)"\
                          %outdir )

    parser.add_option("-r", "--rcFile", dest="rcFile", default=rcFile,
              help="Resource file defining parameters to sample (default=%s)"\
                          %rcFile )

    parser.add_option("--dt_minutes", dest="dt_minutes", default=dt_minutes,
              help="if no iso_t2, time from iso_t1 (default=%d)"\
                          %dt_minutes )    

    parser.add_option("-t", "--title", dest="title", default=title,
              help="Output file title, typically the collection name (default=%s)"\
                          %title )

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 1:
        iso_t1 = args[0]
        iso_t2 = None
        t1     = isoparser(iso_t1)
        t2     = t1 + timedelta(minutes=options.dt_minutes)
    else:
        parser.error("must have 1 or 2 arguments: iso_t1 [iso_t2]")

    # Create consistent file name extension
    # -------------------------------------
    if 'NETCDF4' in options.format:
        options.ext = 'nc4'
    elif 'NETCDF3' in options.format:
        options.ext = 'nc'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
    options.zlib = not options.nozip

    if not os.path.exists(options.outdir): 
        print 'outdir does not exist:',options.outdir
        raise ValueError, 'check path. exiting'


    # Get Granules to work on
    # -----------------------
    options.granules = granules (options.L1Root, t1, t2)

    # Create (x,y,t) coordinates
    # --------------------------
    pace = PACE(options.granules)

    # Get Variables and Metadata
    # --------------------------
    Vars, levs, levUnits = getVars(options.rcFile)

    # Write output file
    # -----------------
    writeNC ( pace, Vars, levs, levUnits, options, doAkBk=False)
