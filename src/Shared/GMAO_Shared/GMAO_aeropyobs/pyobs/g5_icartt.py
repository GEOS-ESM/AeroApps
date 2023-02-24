#!/usr/bin/env python
"""
     Produces ICARTT files at specific locations.
     A. da Silva <arlindo.dasilva@nasa.gov>, July 2011.
"""

import os
import string

from types            import *
from numpy            import ones, zeros, ndarray, interp, array
from datetime         import datetime, timedelta
from matplotlib.dates import num2date, date2num

from MAPL  import Config
from grads import GaLab

import icartt

COLLECTIONS = [ 'inst1_2d_hwl_Nx',  # Default collections
                'inst3_2d_asm_Nx',
                'tavg1_2d_flx_Nx',
                'tavg1_2d_lnd_Nx',
                'tavg1_2d_ocn_Nx',
                'tavg1_2d_rad_Nx',
                'tavg1_2d_slv_Nx',
                'tavg3_2d_adg_Nx',
                'tavg3_2d_aer_Nx',
                'tavg3_2d_chm_Nx',
                'tavg3_2d_tag_Nx',
                ]

Months = ('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC')

MISSING = -99999
#...........................................................................................
class G5_ICARTT(object):
    """
    Implements GEOS-5 sampler in terms of ICARTT files.
    This functionality is needed for supporting NASA's
    Airborne Campaigns.
    """

    def __init__ (self, rc_Sites,
                  coll_names = None,
                  coll_dir = './Collections',
                  top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim',
                  template='discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict',
                  outTopDir='./Archive',
                  PI='DASILVA.ARLINDO',
                  Verbose=False,
                  ):
        """
        Instantiate a G5_ICARTT object. On input,

        rc_Sites    ---  resource file defining sites to interpolate to

        coll_names  ---  list with names of GEOS-5 collections to sample, e.g.,
                                inst1_2d_hwl_Nx
                                inst3_2d_asm_Nx
                                etc.
                         If not specified, an internal defined set of
                         collections (defined in global variable COLLECTIONS)
                         will be used.

        coll_dir    ---  directory where to find the collection resource
                         file; these resource files are assumed to have
                         the same name as the collection and extension ".rc",
                         e.g., "inst1_2d_hwl_Nx.rc"

        top_url     ---  OPeNDAP URL where to find collections

        template    ---  template for ICARTT file.
        PI          ---  PI's last time as register on mission data site
        """

        # Defaults
        # --------
        if coll_names is None:
            coll_names = COLLECTIONS

        # Parse sites and collection information
        # --------------------------------------
        self.Sites = _parseSites(rc_Sites)
        self.Collections = _parseColls(coll_names, dirn=coll_dir, top_url=top_url)

        # Save this for later
        # -------------------
        self.template = template
        self.outTopDir = outTopDir
        self.PI = PI
        self.Verbose = Verbose
        
        # Place holders
        # -------------
        self.Collection = None
        self.Table = None
        self.date = None

#---
    def Sample_N_Write (self, date):
        """
        Sample all collections on this date, and write ICARTT files for all sites.
        Use methods sample() and write() to perform these operations one collection
        at a time.

        Each ICARTT file will contain all the time steps available within the
        collection on this date.
        """
        for c in self.Collections:
            self.sample(c,date)
            self.write()

    __call__ = Sample_N_Write
    
#---
    def sample (self, coll_name, date):
        """
        Interpolate collection named *coll_name* on this date (a datetime object)
        to all sites.
        """

        self.Collection = self.Collections[coll_name]
        self.date = date
        self.Table = _getTable( self.Sites, self.Collection, date, Verbose=self.Verbose)

#---
    def write (self):
        """
        Write ICARTT files for each individual site. You must use method
        sample() before doing a write().
        """
        if self.Verbose:
            print("")
        s = 0
        for site in self.Sites:
            _writeICARTT ( self.Sites[site], self.Collection, self.Table[:,:,s], self.date,
                           self.template, self.outTopDir, self.PI, self.Verbose)
            s += 1
    
#...........................................................................................

class G5_AIRCRAFT(object):
    """
    Implements GEOS-5 sampler using space-time coordinates from a glight path given
    bt an ICARTT file.
    This functionality is needed for supporting NASA's Airborne Campaigns.
    """

    def __init__ (self, ict_flight,
                  coll_names = None,
                  coll_dir = './Collections',
                  top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim',
                  template='discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict',
                  outTopDir='./Archive',
                  PI='DASILVA.ARLINDO',
                  Verbose=False,
                  ):
        """
        Instantiate a G5_AIRCRAFT object given flight navigation information (in ICARTT
        file flight_ict) and a set of GEOS-5 collections. On input,
        
        ict_flight  ---  ICARTT file name with navigation information

        coll_names  ---  list with names of GEOS-5 collections to sample, e.g.,
                                inst1_2d_hwl_Nx
                                inst3_2d_asm_Nx
                                etc.
                         If not specified, an internal defined set of
                         collections (defined in global variable COLLECTIONS)
                         will be used.

        coll_dir    ---  directory where to find the collection resource
                         file; these resource files are assumed to have
                         the same name as the collection and extension ".rc",
                         e.g., "inst1_2d_hwl_Nx.rc"

        top_url     ---  OPeNDAP URL where to find collections

        template    ---  template for ICARTT file.
        PI          ---  PI's last time as register on mission data site

        """

        # Defaults
        # --------
        if coll_names is None:
            coll_names = COLLECTIONS

        # Navigation information
        # ----------------------
        self._getNav(ict_flight)
        self.when = [ self.tyme[0], self.tyme[-1] ] # time range
        self.date = self.when[0].replace(hour=0,minute=0,second=0,microsecond=0)

        # "Site" information
        # ------------------
        self.aircraft = ict_flight.split('_')[1] # aircraft name
        Site = dict ( name = self.aircraft, 
                      lon = self.lon,
                      lat = self.lat,
                      elevation = self.alt )
        self.Sites = { self.aircraft: Site }
        self.Site = Site

        # Parse collection information
        # ----------------------------
        self.Collections = _parseColls(coll_names, dirn=coll_dir, top_url=top_url)

        # Save this for later
        # -------------------
        self.template = template
        self.outTopDir = outTopDir
        self.PI = PI
        self.Verbose = Verbose
        
        # Place holders
        # -------------
        self.Collection = None
        self.Table = None

#---
    def sample (self, coll_name, AltVar='h', AltColl_name=None):
        """
        Interpolate collection named *coll_name* to navigation trajectory.
        """

        if AltColl_name is not None:
            AltColl = self.Collections[AltColl_name]
        else:
            AltColl = None
            
        #  Get Table at model times
        #  ------------------------
        self.Collection = self.Collections[coll_name]
        Table_ = _getTable( self.Sites, self.Collection, self.when,
                            Ghosted = True,
                            AltVar=AltVar, AltColl=AltColl,
                            Verbose=self.Verbose)

        nt, nv, N = Table_.shape

        # Interpolate to aircraft time
        # ----------------------------
        tnum = date2num(self.tyme)
        Table = icartt.MISSING * ones((nv,N)) # will transpose later
        for n in range(nt-1):
            t1 = Table_[n,0,0]
            t2 = Table_[n+1,0,0]
            I = (tnum>=t1) & (tnum<t2)
            a = (tnum[I] - t1) / ( t2 - t1 )
            for j in range(nv):
                Table[j,I] = (1-a) * Table_[n,j,I] + a * Table_[n+1,j,I] 

        # Express time as seconds from beginning of day
        # ---------------------------------------------
        for n in range(N):
            Table[0,n] = (self.tyme[n]-self.date).total_seconds()

        # Transpose to a form extected by _writeICARTT
        # --------------------------------------------
        self.Table = Table.T # (N times, variables)
        
#---
    def write (self):
        """
        Write ICARTT file with each aircraft time. You must use method
        sample() before doing a write().
        """

        _writeICARTT ( self.Site, self.Collection, self.Table, self.date,
                       self.template, self.outTopDir, self.PI, self.Verbose)
    
#---

    def _getNav(self,ict_flight):
        """
        Retrieve navigation information.
        """
        
        Nav = icartt.ICARTT(ict_flight).Nav
        if Nav['Time']     ==None or \
               Nav['Longitude']==None or \
               Nav['Latitude'] ==None or \
               Nav['Altitude'] ==None:
            raise ValueError('Invalid or imcomplet navigation in '+ict_flight) 
       
        self.tyme = Nav['Time']
        self.lon = Nav['Longitude']
        self.lat = Nav['Latitude']
        self.alt = Nav['Altitude']
        self.lev = Nav['Pressure']

def _parseSites ( rc_Sites ):
    """
    Get Coordinates and archival info from resource file.
    Returns a dictionary that can be address like this:

      Sites[site]['lat']
      Sites[site]['lon']
      
    """
    cf = Config(rc_Sites)
    Sites = {}
    for site in list(cf.keys()):
        Info = {}        
        for token in cf(site).replace(' ','').split(';'):
            att, value = token.split('=')
            Info[att] = value
        Info['name'] = site
        Sites[site] = Info

    return Sites

#---
def _parseColl ( rc_Collection, name, url ):
    """
    Parse collection resource file.
       rc_Collection --- resource file defining collection
       url           --- url/filename for opening the collection
    """
    cf = Config(rc_Collection)
    Shorts = list(cf.keys())
    del Shorts[Shorts.index('__TITLE__')]
    Variables = {}
    for short in Shorts:
        long, units = cf(short).split(';')
        units = units.replace(' ','').replace('1','none')
        Variables[short] = dict(long=long, units=units)

    Collection = dict (  name = name,
                        title = cf('__TITLE__'),
                         vars = Variables,
                           rc = rc_Collection,
                          url = url )
    return Collection

#---
def _parseColls ( Names, dirn='./Collections',
                  top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim' ):
    """
    Returns a dictionary of collections.
    """
    Collections = {}
    for name in Names:
        rc = dirn + '/' + name + '.rc'
        url = top_url + '/' + name
        Collections[name] = _parseColl(rc,name,url) 
    return Collections

#---
def _gatime2dt(gatime):
    """
    Convert gatime to datetime.
    """
    time, date = gatime.upper().split('Z')
    if time.count(':') > 0:
        h, m = time.split(":")
    else:
        h = time
        m = '0'
    mmm = date[-7:-4]
    dd, yy = date.split(mmm)
    mm = Months.index(mmm) + 1
    dt = datetime(int(yy),int(mm),int(dd),int(h),int(m))
    return dt

#---
def _dt2gatime(t):
    """Convert datetime to grads time."""
    gatime = "%d:%dZ%d%s%d"%(t.hour,t.minute,t.day,Months[t.month-1],t.year)
    return gatime
    

#---
def _getTable ( Sites, Collection, when,
                Ghosted=False,
                AltVar='h', AltColl=None, zrange=None,
                Verbose=False, Echo=False ):
    """
    Generate data table for ICARTT files with GEOS-5 output.
    Data is returned at models discrete times, with only
    lat/lon interpolation.

    Sites      ---  a *Sites* dictionary (see _parseSites);
                    Note: "aircraft mode" is assumed if there 
                    is only one site, and the coordinates for
                    this site are arrays rather than single values.
    Collection ---  GEOS-5 collection (see _parseCollections)
    when       ---  single date or range
    Ghosted    ---  whether to include ghosting for time interp
    AltVar     ---  altitude variable on GEOS-5 file
    AltColl    ---  Collection with altitude variable         
    zrange     ---  vertical level range (in GrADS index space);
                    default: all levels on file.
    """

    # Type of Sites determin whether we are in aircraft mode or not
    # -------------------------------------------------------------
    if len(Sites) == 1 and isinstance(list(Sites.values())[0]['lon'],ndarray):
        Aircraft_mode = True
        aircraft = list(Sites.keys())[0]
        Site = Sites[aircraft]
    else:
        Aircraft_mode = False

    if Verbose:
        print("")
        print("             GEOS-5 ICARTT Sampler")
        print("             ---------------------")
        print("")
        print("Collection: ", Collection['name'])
        print("            ", Collection['title'])
        if type(when) in (ListType,TupleType):
            print("      When: ", when[0], ' to ', when[1])
        else:
            print("      When: ", when.date())
        print("     Sites: ")
        for s in Sites:
            print("            ", s)
        print("")

    # Guess ftype
    # -----------
    if 'http' in Collection['url']:
        ftype = 'sdf'
    else:
        ftype = 'ctl'

    # For robustness, start a fresh GrADS process
    # -------------------------------------------
    ga = GaLab(Window=False,Echo=Echo)

    # Open the collection
    # -------------------
    fh = ga.open(Collection['url'],ftype=ftype)

    # If in Aircraft mode, handle altitude variable
    # ---------------------------------------------
    if Aircraft_mode:
        if AltColl is not None:
            fhh = ga.open(AltColl['url'],ftype=ftype)
            AltVar = AltVar + '.2'
        else:
            fhh = fh
            AltVar = AltVar + '.1'
        try:
            nlevs = fhh.var_levs[fhh.vars.index(AltVar.split('.')[0])]
            if nlevs != fhh.nz:
                raise ValueError("Strange, altitude not defined in all levels")
        except:
            raise ValueError("Altitude variable no present")

    # otherwise, cannot yet handle 3D variables
    # ------------------------------------------
    else:
        for var in Collection['vars']:
            if fh.var_levs[fh.vars.index(var.lower())]>0:
                raise ValueError("cannot yet handle 3D variables unless in aircraft mode")

    # Serialize lat/lon coordinates
    # -----------------------------
    if Aircraft_mode:
        lons, lats, alts = (Site['lon'], Site['lat'], Site['elevation']) 
    else:
        lats = ones(len(Sites))
        lons = ones(len(Sites))
        alts = ones(len(Sites))
        i = 0
        for s in Sites:
            lats[i] = float(Sites[s]['lat'])
            lons[i] = float(Sites[s]['lon'])
            if type(Sites[s]['elevation']) == type("abc"):
                alts[i] = float(Sites[s]['elevation'].replace('m',''))
            else:
                alts[i] = float(Sites[s]['elevation'])
            i += 1
        
    # Tight bounding box to speed up I/O
    # ----------------------------------
    ga('set lon %f %f'%(lons.min(),lons.max()))
    ga('set lat %f %f'%(lats.min(),lats.max()))
    qh = ga.query("dims")
    x1, x2 = (qh.xi[0]-1,qh.xi[1]+1)
    y1, y2 = (max(1,qh.yi[0]-1),min(fh.ny,qh.yi[1]+1)) # in [1,ny]
    ga('set x %d %d'%(x1,x2),Quiet=True) # make sure x range is int
    ga('set y %d %d'%(y1,y2),Quiet=True) # make sure y range is int

    # Set vertical coordinate
    # -----------------------
    if zrange == None:
        ga("set z %d %d"%(1,fh.nz))
    else:
        ga("set z %d %d"%zrange)

    # Find time indices
    # -----------------
    if type(when) is ListType: # generic time range
        ga('set time %s'%_dt2gatime(when[0]))
        qh = ga.query("dims")
        tbeg =int(qh.t[0])
        if when[0] < _gatime2dt(qh.time[0]):
            tbeg = tbeg - 1
        ga('set time %s'%_dt2gatime(when[1]))
        qh = ga.query("dims")
        tend = int(qh.t[0])
        if when[1] > _gatime2dt(qh.time[0]):
            tend = tend + 1
        t0 = when[0].replace(hour=0,minute=0,second=0,microsecond=0)
        Ghosted = True # implicit above
        
    else: # single date
        ga('set time %d%s%d'%(when.day,when.ctime()[4:7],when.year))
        qh = ga.query("dims")
        tbeg =int(qh.t[0])
        next = when + timedelta(seconds=86400)
        ga('set time %d%s%d'%(next.day,next.ctime()[4:7],next.year))
        qh = ga.query("dims")
        tend = int(qh.t[0])-1
        t0 = when.replace(hour=0,minute=0,second=0,microsecond=0)
        if Ghosted:
            tbeg, tend = tbeg-1,tend+1 # needed for time interpolation

    T = list(range(tbeg,tend+1)) 

    # Sample variables at site locations
    # ----------------------------------
    ns, nt, nv = (len(lons), len(T), len(Collection['vars']))
    Table = zeros((nt,nv+1,ns)).astype('float') # to hold numeric values
    i = 0 # time index

    # Loop over time on file...
    # -------------------------
    for t in T:

       # Set time
       # --------
       ga('set t %d'%t)
       ga('query time')
       time_ = _gatime2dt(ga.rword(1,3))
       dt = time_ - t0
       if Verbose:
          print("[] Sampling", time_, '-->', dt.total_seconds(), 'secs')

       # Fill in time variable (by convention the first one in Table)
       # ------------------------------------------------------------
       if Ghosted:
           Table[i,0,:] = date2num(time_) * ones(ns) # Use this for time interp
       else:
           Table[i,0,:] = dt.total_seconds() * ones(ns) # first var is time, for icartt

       # In aircraft mode, get altitude variable
       # ---------------------------------------
       if Aircraft_mode:
           h, levs = ga.interp(AltVar,lons,lats)
           if h[0,0]>h[0,-1]: # make sure data is in ascending height order
               h[:,:] = h[:,-1::-1]
               hflip = True
           else:
               hflip = False

       # Loop over variables in Collection
       # ---------------------------------
       j = 1 # variable index
       for var in list(Collection['vars'].keys()):
          if Verbose>1:
             print("   - Doing <%s>"%var)
             
          # Interpolate variable to la/lon
          # ------------------------------
          try:
              values, levs = ga.interp(var,lons,lats)
          except:
              ga.flush()
              values = MISSING * ones(size(lons))

          # If needed, interpolate to altitude
          # (This only happens in aircraft mode)
          # ------------------------------------
          if len(lons) == 1:
              values = values * ones(1)
          if len(values.shape) == 1:
              Table[i,j,:] = values
          else:
              if hflip:
                  values[:,:] = values[:,-1::-1]
              for k in range(ns):
                  Table[i,j,k] = interp(alts[k],h[k,:],values[k,:]) # vertical interp
                    
          j += 1 # bump variable index
       i += 1    # bump time     index

    return Table

#---
def _writeICARTT ( Site, Collection, Table, date,
                   template='discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict',
                   outTopDir='./Archive', PI='DASILVA.ARLINDO', Verbose=False):
    """
    Writes the ICARTT file from a template.
    """

    # Special handling for aircraft mode
    # ----------------------------------
    if isinstance(Site['lon'],ndarray):
        Aircraft_mode = True
    else:
        Aircraft_mode = False

    product = Collection['name'].replace('_','-')
    date_ = str(date.year * 10000 + date.month*100 + date.day)
    now = date.utcnow()

    # Output file name
    # ----------------
    outDir = outTopDir + '/' + Site['name'] + '/' + PI
    filename = template.replace('PRODUCT',product).\
                        replace('SITE',Site['name']).\
                        replace('DATE',date_).\
                        replace('REVISION','RA').\
                        replace('.1001','')

    # Header template
    # ---------------
    Header = open(template).readlines()
    

    # Information for filling in the header
    # -------------------------------------
    nv = len(Collection['vars'])
    if Aircraft_mode:
        nv = nv + 3 # also include lon, lat, alt
    nh = len(Header) + nv - 1
    Product = Collection['title']
    first_date = '%d, %d, %d'%(date.year,date.month,date.day)
    rev_date = '%d, %d, %d'%(now.year,now.month,now.day)
    dt = 0 # most people use 0, not sure why
    scale = str(list(ones(nv).astype('int')))\
                .replace('[','')\
                .replace(']','')
    missing = str(list((icartt.MISSING*ones(nv)).astype('int')))\
                  .replace('[','')\
                  .replace(']','')

    if Aircraft_mode:
        Coords = ['UTC_start','Longitude','Latitude','Altitude',]
    else:
        Coords = ['UTC_start',]

    vnames = str(Coords + list(Collection['vars'].keys()))\
                 .replace('[','')\
                 .replace(']','')\
                 .replace("'","")\
                 .replace(' ','')

    if Aircraft_mode:
        lon_, lat_, elev_ = ("N/A", "N/A", "N/A")
    else:
        lon_, lat_, elev_ = (Site['lon'], Site['lat'], Site['elevation'])

    # Update Header
    # -------------
    Info = dict ( NUM_HEADER = nh,
                  PRODUCT = Product,
                  FIRST_DATA_DATE = first_date,
                  REV_DATE = rev_date,
                  DT_seconds = dt,
                  NUM_VARS = nv,
                  SCALE_FACTORS = scale,
                  MISSING_VALUES = missing,
                  lon = lon_,
                  lat = lat_,
                  elevation = elev_,
                  VAR_NAMES = vnames,
                )
        
    if Verbose:
        print("<> Writing <%s>"%filename)

    # Create directory and write out the header
    # -----------------------------------------
    os.system("mkdir -p "+outDir)
    f = open(outDir+'/'+filename,'w')
    for line in Header:
        if line.count('__VARIABLE_METADATA__')>0:
            if Aircraft_mode:
                f.write('Longitude, degrees_east, Aircraft Longitude\n')
                f.write('Latitude,  degrees_north, Aircraft Latitude\n')
                f.write('Altitude,   meters, Aircraft Altitude\n')
            for var in Collection['vars']:
                f.write('%s, %s, %s\n'%(var,Collection['vars'][var]['units'],
                                          Collection['vars'][var]['long']))
        else:
            f.write(string.Template(line).safe_substitute(Info))                

    # Write the data
    # --------------
    if Aircraft_mode:
        lon, lat, alt = (Site['lon'], Site['lat'], Site['elevation'])
    nt = Table.shape[0]
    for n in range(nt):
        f.write("%f,"%Table[n,0]) # time
        if Aircraft_mode:
            f.write("%f,%f,%f,"%(lon[n],lat[n],alt[n]))
        line = str(list(Table[n,1:]))\
                 .replace('[','')\
                 .replace(']','')\
                 .replace(' ','')\
                 .replace('.0,',',')\
               + '\n'
        f.write(line)
    f.close()
     
#---
def _g5_icartt ( Sites, Collection, date,
                template='discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict',
                outTopDir='./Archive',
                PI='DASILVA.ARLINDO',
                Verbose=False, Echo=False ):
    """
    Produces ICARTT files with GEOS-5 output.

    filename ---  grADS filename or OPeNDAP URL
    
    """

    # Interpolate Data to site location
    # ---------------------------------
    Table = _getTable( Sites,Collection,date,Verbose=Verbose,Echo=Echo)

    # Write ICARTT file for each site
    # -------------------------------
    s = 0
    for site in Sites:
        _writeICARTT ( Sites[site], Collection, Table[:,:,s], date,
                       template, outTopDir, PI, Verbose)
        s += 1
        
#    return Table

#--
def _ut():
    
    Verbose = 1
    
    rc_Sites = './d-aq_sites.rc'
    dir_colls = './Collections'
    
    Sites = _parseSites(rc_Sites)
    Collections = _parseColls ([ 'inst1_2d_hwl_Nx',
                                 'inst3_2d_asm_Nx',
                                 'tavg1_2d_flx_Nx',
                                 'tavg1_2d_lnd_Nx',
                                 'tavg1_2d_ocn_Nx',
                                 'tavg1_2d_rad_Nx',
                                 'tavg1_2d_slv_Nx',
                                 'tavg3_2d_adg_Nx',
                                 'tavg3_2d_aer_Nx',
                                 'tavg3_2d_chm_Nx',
                                 'tavg3_2d_tag_Nx',
                                 ])
    Collections = _parseColls ([ 'inst1_2d_hwl_Nx', ])

    # Sample date
    # -----------
    dt = timedelta(seconds=86400)
    t     = datetime(2011,0o7,0o1)
    t_max = datetime(2011,0o7,17)
    while t <= t_max:
        for c in Collections:
            _g5_icartt ( Sites, Collections[c], 
                        t, Verbose=Verbose)
        t = t + dt

#---

if __name__ == "__main__":

    # Sample collection
    # -----------------
    t = datetime(2011,0o7,0o1)
    g = G5_ICARTT('./d-aq_sites.rc', Verbose=True)
    g.sample(COLLECTIONS[0],t)
    g.write()

    # All collections
    # ---------------
    t = datetime(2011,0o7,0o2)
    g.Sample_N_Write(t)
    
