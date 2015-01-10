"""
   Classes for reading CSV files from the AERONET Dragon subset during DISCOVER-AQ.
"""

MISSING = -999.

import os
import string

from numpy import loadtxt, ones, savez, pi, log, interp, sort

from datetime         import datetime
from matplotlib.dates import num2date, date2num

import icartt
from   g5_icartt      import _getTable, _parseColls

Months = ('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC')

Alias = dict(          Long = 'lon',
                        Lat = 'lat',
                Elevation_m = 'elev',
              Temperature_C = 'an_T2m',
                 WaterVapor = 'an_TQV',
                  TOTEXTTAU = 'g5_AOD_550',
                  DUEXTTAU  = 'g5_AODdu_550',
                  BCEXTTAU  = 'g5_AODbc_550',
                  OCEXTTAU  = 'g5_AODoc_550',
                  SSEXTTAU  = 'g5_AODss_550',
                  SUEXTTAU  = 'g5_AODsu_550',
                        SLP = 'g5_SLP',
                       PS   = 'g5_PS',
                       U10M = 'g5_U10m',
                       U2M  = 'g5_U2m',
                       V10M = 'g5_V10m',
                        V2M = 'g5_v2m',
                       T10M = 'g5_T10m',
                        T2M = 'g5_T2m',
                      QV10M = 'g5_QV10m',
                      QV2M = 'g5_QV2m',
                        TS = 'g5_TS',
                       DISPH = 'g5_DISPH',
                       CLDPRS = 'g5_CLDPRS',
                       CLDTMP = 'g5_CLDTMP',
                       TQV = 'g5_TQV',
                       TQI = 'g5_TQI',
                       TQL = 'g5_TQL',
                       TOX = 'g5_TOX',
                   )

iAlias = {}
for v in Alias:
    iAlias[Alias[v]] = v # inverse alias

#----
class DRAGON(object):
    """Base class for the DRAGON aeronet subset."""

    def __init__ (self,filename=None):
        """
        Creates a DRAGON objec from a file name. If no file name specified,
        an empty object is created.
        """

        if filename is None:
            return

        Vars = open(filename).readline().replace('\n','').\
                    replace('(','_').replace(')','').replace(' ','').\
                    split(',')

        self.filename = filename

        # Read relevant columns from MAPSS granule
        # ----------------------------------------
        formats = ()
        converters = {}
        i = 0
        for name in Vars:
            if name=='Date':
                formats = formats + ('S10',)
            elif name=='Time':
                formats = formats + ('S8',)
            elif name=='Site':
                formats = formats + ('S20',)
            else:
                converters[i] = lambda s: float(s or MISSING)
                formats += ('f4',)
            i += 1
                
        # Read the data
        # -------------
        data = loadtxt(filename, delimiter=',',
                       dtype={'names':Vars,'formats':formats},
                       converters = converters,
                       skiprows=1)
 
        N = len(data)
        
        self.N = N

#       Save data columns as attributes, with nicer names
#       -------------------------------------------------
        for i in range(len(Vars)):
            name = Vars[i]
            if formats[i]=='f4':
                v = ones(N)
                for j in range(N):
                    v[j] = data[j][i]
            else:
                v = []
                for j in range(N):
                    v.append(data[j][i])

            self.__dict__[name] = v

            if Alias.__contains__(name):
                self.__dict__[Alias[name]] = v


#       Interpolate AOD to 550 (raw Dragon files only)
#       ----------------------------------------------
        if 'an_AOD_550' not in self.__dict__:
            i = (self.AOD_500>0) & (self.AOD_675>0)
            alpha = MISSING * ones(N)
            self.AOD_550 =  MISSING * ones(N)
            if any(i):
                alpha = (log(0.550) - log(self.CWL_500[i]) ) / ( log(self.CWL_675[i]) - log(self.CWL_500[i]) )
                self.AOD_550[i] =  (1.-alpha) * self.AOD_500[i] + alpha * self.AOD_675[i]
            self.an_AOD_550 = self.AOD_550

#       Create datetime object
#       ----------------------
        self.tyme = []
        self.gatime = []
        for d, t in zip(self.Date,self.Time):
            dd, mm, yy = d.split(':')
            h, m, s = t.split(':')
            if yy=="1900":
                raise ValueError, 'Invalid date on file: '+d
            else:
                tyme = datetime(int(yy),int(mm),int(dd),int(h),int(m),int(s))
            self.tyme.append(tyme)
            self.gatime.append(_dt2gatime(tyme))

        self.tnum = date2num(self.tyme)

#--
    def getSites(self,PI='DASILVA.ARLINDO'):
        """
        Return dictionary with DRAGON data stratified by sites.
        """

        N = self.N
        
        # Register sites
        # --------------
        SITES = {}
        for i in range(N):
            name = self.Site[i]
            path = 'AERONET/'+name.replace('DRAGON_','')+'/'+PI
            SITES[name] = dict ( name = name,
                                 lon = self.lon[i],
                                 lat = self.lat[i],
                                 elevation = self.elev[i],
                                 path = path,
                               )
            
        # Subset at site
        # --------------
        for s in SITES:
            site = SITES[s]
            I = (ones(N)<0) # all false
            Date = []
            Time = []
            for i in range(N):
                if self.Site[i] == s:
                    I[i] = True
                    Date += (self.Date[i],)
                    Time += (self.Time[i],)
                    
            self_ = DRAGON()
            self_.I = I
            self_.Date = Date
            self_.Time = Time
            self_.tnum = self.tnum[I]
            self_.tyme = num2date(self.tnum[I])
            self_.an_AOD_550 = self.an_AOD_550[I] 
            self_.an_T2m = self.an_T2m[I]
            self_.an_TQV = self.an_TQV[I]
            for v in self.__dict__:
                if v[0:3] == 'g5_':
                    self_.__dict__[v] = self.__dict__[v][I]
            self_.N = len(self_.tyme)
            
            SITES[s]['data'] = self_

        return SITES

#---
    def sample_N_writeDragon (self, filename,
                      coll_names = ['inst1_2d_hwl_Nx','inst3_2d_asm_Nx',], 
                      coll_dir = './Collections/Dragon',
                      top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim',
                      outTopDir='./Archive',
                      PI='DASILVA.ARLINDO',
                      Verbose=False,
                      ):
        """
        Sample GEOS-5 collections at individual DRAGON sites writing to a file.
        """

        # Sample GEOS-5 collections
        # -------------------------
        Sites, Collections = self.sampleG5atSites ( coll_names = coll_names,
                                                      coll_dir = coll_dir,
                                                       top_url = top_url,
                                                            PI = PI,
                                                       Verbose = Verbose
                                                   )

        # Determine AERONET and GEOS-5 variables
        # --------------------------------------
        allVars = Sites.values()[0]['data'].__dict__.keys()
        anetVars = [ 'Date', 'Time', 'Site', 'Lat', 'Long', 'Elevation(m)',
                     'an_AOD_550','an_T2m', 'an_TQV' ]
        xVars = [ 'tnum', 'I', 'tyme', 'N' ]
        g5Vars = []
        for v in allVars:
            if v not in anetVars+xVars:
                g5Vars += [v,]
        header = str(anetVars+g5Vars).replace('[','').replace(']','').replace("'","").replace(' ','')

        # Write out Dragon-stype file
        # ---------------------------
        f = open(filename,'w')
        f.write(header+'\n')
        for s in Sites:
            d = Sites[s]['data']
            lat, lon, elev = Sites[s]['lat'], Sites[s]['lon'], Sites[s]['elevation'] 
            for n in range(d.N):
                f.write("%s,%s,%s,%f,%f,%f,%f,%f,%f"%(d.Date[n],d.Time[n],s,lat,lon,elev,
                             d.an_AOD_550[n],d.an_T2m[n],d.an_TQV[n]))
                for v in g5Vars:
                    f.write(",%f"%d.__dict__[v][n])
                f.write("\n")
        f.close()
        if Verbose:
            print "[] wrote %s"%filename

#---
    def sample_N_writeICARTT (self, 
                      coll_names = ['inst1_2d_hwl_Nx','inst3_2d_asm_Nx',], 
                      coll_dir = './Collections/Dragon',
                      top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim',
                      template='discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict',
                      outTopDir='./Archive',
                      PI='DASILVA.ARLINDO',
                      Verbose=False,
                      ):
        """
        Sample GEOS-5 collections at individual DRAGON sites writing to a file.
        """

        # Assumes single date
        # -------------------
        if self.tyme[0].date() != self.tyme[-1].date():
            raise ValueError, "time range span multiple days, cannot proceed"
        t0 = self.tyme[0]
        date = t0.replace(hour=0,minute=0,second=0,microsecond=0)

        # Sample GEOS-5 collections
        # -------------------------
        Sites, Collections = self.sampleG5atSites ( coll_names = coll_names,
                                                      coll_dir = coll_dir,
                                                       top_url = top_url,
                                                            PI = PI,
                                                       Verbose = Verbose
                                                   )

        # Determine AERONET and GEOS-5 variables
        # --------------------------------------
        allVars = Sites.values()[0]['data'].__dict__.keys()
        anetVars = [ 'an_AOD_550','an_T2m', 'an_TQV' ]
        xVars = [ 'Date', 'Time', 'Site', 'Lat', 'Long', 'Elevation(m)', 'tnum', 'I', 'tyme', 'N' ]
        g5Vars = []
        for v in allVars:
            if v not in anetVars+xVars:
                g5Vars += [v,]
        keepVars = anetVars + g5Vars
        
        # Variable description for ICARTT files
        # -------------------------------------
        Vars = {}

        # For AERONET, hardwire it
        # ------------------------
        Vars['an_AOD_550'] = dict ( long = 'Aeronet AOD at 550 nm (interpolated from 500/675 nm)' ,
                                    units = 'none' )
        Vars['an_T2m']     = dict ( long = 'Aeronet 2m Temperature',
                                    units = 'C' )
        Vars['an_TQV']     = dict ( long = 'Aeronet column water vapor' ,
                                    units = 'cm' )
        
        # For GEOS-5 variables, get this info from Collections
        # ----------------------------------------------------
        for v in g5Vars:
            v_ = iAlias[v] # original variable name in collection
            for c in Collections:
                if v_ in Collections[c]['vars']:
                    Vars[v] = Collections[c]['vars'][v_]

        # Write ICARTT file for each site
        # -------------------------------
        for s in Sites:

            # Write ICARTT file
            # -----------------
            _writeICARTT ( Sites[s], date, Vars, 
                            template = template,
                           outTopDir = outTopDir,
                                  PI = PI,
                             Verbose = Verbose)

#---

    def sampleG5atSites (self, 
                         coll_names = ['inst1_2d_hwl_Nx','inst3_2d_asm_Nx',], 
                         coll_dir = './Collections/Dragon',
                         top_url='http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/assim',
                         PI='DASILVA.ARLINDO',
                         Verbose=False,
                         ):
        """
        Sample GEOS-5 collections at individual DRAGON sites.

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

        PI          ---  PI's last time as register on mission data site

        Example:

          Sites = d.sampleG5atSites()

        """

        # Assumes single date
        # -------------------
        if self.tyme[0].date() != self.tyme[-1].date():
            raise ValueError, "time range span multiple days, cannot proceed"
        t0 = self.tyme[0]
        date = t0.replace(hour=0,minute=0,second=0,microsecond=0)

        # Parse sites and collection information
        # --------------------------------------
        Sites = self.getSites(PI=PI)
        Collections = _parseColls(coll_names, dirn=coll_dir, top_url=top_url)

        # Sample GEOS-5 for this date
        # ---------------------------
        for c in Collections:

            # Get data at GEOS-5 times bracketing the day
            # -------------------------------------------
            Table = _getTable( Sites, Collections[c], date,
                               Ghosted=True, Verbose=Verbose)

            # Time Interpolation for each site
            # --------------------------------
            i = 0
            for s in Sites:
                d = Sites[s]['data']
                j = 1 # first in table is time
                for v in Collections[c]['vars'].keys():
                    if v in Alias:
                        v = Alias[v]
                    d.__dict__[v] = interp(d.tnum,Table[:,0,i],Table[:,j,i])
                    j += 1
                i += 1
                                      
        return (Sites, Collections)

#---
    def addVar(self,ga,expr='tau',vname=None, clmYear=None):
        """
        Given a grads object having the correct file as default,
        creates an attribute name *vname* with values interpolated
        to (lon,lat) at the correct time window (time "interpolation"
        is nearest neighbours.)
        """

        raise ValueError, "method addVar not tested yet"

        N = self.N
        U = ones(N)
        U[:] = MISSING

        if vname == None:
            vname = expr

        # Tight bounding box to speed up I/O
        # ----------------------------------
        fh = ga.query("file")
        ga('set lon %f %f'%(self.lon.min(),self.lon.max()))
        ga('set lat %f %f'%(self.lat.min(),self.lat.max()))
        qh = ga.query("dims")
        x1, x2 = (qh.xi[0]-1,qh.xi[1]+1)
        y1, y2 = (max(1,qh.yi[0]-1),min(fh.ny,qh.yi[1]+1)) # in [1,ny]
        ga('set x %d %d'%(x1,x2),Quiet=True) # make sure x range is int
        ga('set y %d %d'%(y1,y2),Quiet=True) # make sure y range is int

        # Find grads time indices associated with each station time
        # ---------------------------------------------------------
        T = ones(N).astype('int')
        for i in range(N):
            ga("set time "+self.gatime[i])
            qh = ga.query("dims")
            T[i] = qh.t[0]

        ta, tb = T.min(), T.max()
        
        for t in range(ta,tb+1):

            I = (T==t) # find data in this time bracket

            if any(I):
                
                ga('set t %d'%t)
                
                lons = self.lon[I]
                lats = self.lat[I]

                U[I], levs = ga.interp(expr,lons,lats)
                
        self.__dict__[vname] = U

#..........................................................................................

def _writeICARTT ( Site, date, Vars, 
                   template='discoveraq-geos5das-PRODUCT_SITE_DATE_RA.1001.ict',
                   outTopDir='./Archive', PI='DASILVA.ARLINDO', Verbose=False):
    """
    Writes the ICARTT file from a template.
    """

    # Header template
    # ---------------
    Header = open(template).readlines()
    
    # Output file name
    # ----------------
    product = 'aeronet'
    date_ = str(date.year * 10000 + date.month*100 + date.day)
    now = date.utcnow()
 
    outDir = outTopDir + '/' + Site['path']
    filename = template.replace('PRODUCT',product).\
                        replace('SITE',Site['name']).\
                        replace('DATE',date_).\
                        replace('.1001','')

    # Information for filling in the header
    # -------------------------------------
    nv = len(Vars) # does not include leading time column
    nh = len(Header) + nv - 1
    Product = 'AERONET measuments and GEOS-5 fields interpolated to location, time'
    first_date = '%d, %d, %d'%(date.year,date.month,date.day)
    rev_date = '%d, %d, %d'%(now.year,now.month,now.day)
    dt = 0 # most people use 0, not sure why
    scale = str(list(ones(nv).astype('int')))\
                .replace('[','')\
                .replace(']','')
    missing = str(list((icartt.MISSING*ones(nv)).astype('int')))\
                  .replace('[','')\
                  .replace(']','')
    vnames = str(['UTC_start',] + list(sort(Vars.keys())))\
                 .replace('[','')\
                 .replace(']','')\
                 .replace("'","")\
                 .replace(' ','')
    

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
                  lon = Site['lon'],
                  lat = Site['lat'],
                  elevation = Site['elevation'], # perhaps we should use the model elevation there?
                  VAR_NAMES = vnames,
                )
        
    if Verbose:
        print "<> Writing <%s>"%filename

    # Create directory and write out the header
    # -----------------------------------------
    os.system("mkdir -p "+outDir)
    f = open(outDir+'/'+filename,'w')
    for line in Header:
        if line.count('__VARIABLE_METADATA__')>0:
            for v in sort(Vars.keys()):
                f.write('%s, %s, %s\n'%(v,Vars[v]['units'],
                                            Vars[v]['long']))
        else:
            f.write(string.Template(line).safe_substitute(Info))                

    # Write the data
    # --------------
    d = Site['data']
    t0 = d.tyme[0].replace(hour=0,minute=0,second=0,microsecond=0)
    for n in range(d.N):
        utc_secs = (d.tyme[n] - t0).total_seconds()
        f.write("%s"%utc_secs)
        for v in sort(Vars.keys()):
            f.write(",%f"%d.__dict__[v][n])
        f.write("\n")

    # All done
    # --------
    f.close()
     
#---
def _gatime2dt(gatime):
    """
    Convert grads time to datetime.
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
def _gatrange(ga,tyme):
    """
    Given a sequence of datetime objects, find the grads
    time indices (t1,t2) that brackets this sequence.
    """

    num = date2num(tyme)

    time = _dt2gatime(num2date(num.min())) 
    ga('set time '+time)
    qh = ga.query("dims")
    ta = qh.t[0] - 1

    time = _dt2gatime(num2date(num.max())) 
    ga('set time '+time)
    qh = ga.query("dims")
    tb = qh.t[0] + 1

    return (ta,tb)

#....................................................................

if __name__ == "__main__":

#    d = DRAGON('dragon_aod_01-JUL-2011.txt')
    d = DRAGON('/nobackup/AERONET/Level1.5/Dragon/dragon_aod_26-JUL-2011.txt')
    
