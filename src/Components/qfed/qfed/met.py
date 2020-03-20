"""

   Implements interface to Met Fields for QFED.

"""

from types import *
from pylab import find
from numpy import zeros, ones
from grads import GrADS
from datetime import date, timedelta
from eta   import getPm, getEdge

DAY = timedelta(seconds=60*60*24)

class MET(object):

    def __init__(self,filename,lon=None,lat=None,Verb=0,pm_top=None,Vars=None,ga=None):
        """
        Start GrADS, select variables and levels for subsequent reading.
        If pm_top (in Pascal) is specified, the cvertical levels will be reduced.
        """

#       Expects GrADS v2
#       ----------------
        if ga is None:
            ga = GrADS(Bin='grads',Window=False,Echo=False)

#       Open the file
#       -------------
        fh = ga.open(filename)
        ga('set lon -180 180') # good for interpolation

#       Either all varibles on file or user subset
#       ------------------------------------------
        vinfo = []
        if Vars == None:
            for v,k,d,l in fh.var_info:
                vinfo.append((v,k,l,fh.fid)) # notice additional fid element
        else:
            for v,k,d,l in fh.var_info:
                if v in Vars:
                    vinfo.append((v,k,l,fh.fid))

        if len(vinfo)==0:
            print >>sys.stderr, "IndexError: requested variables - ", Vars
            raise IndexError, "cannot find any matchig variable in file %f"\
                  %filename

#       Find the appropriate level range
#       --------------------------------
        if fh.nz>1: # special handle 2D files
            ak, bk = getEdge(fh.nz) # must be a standard number of levels known by set_eta
            pm = getPm(fh.nz)
        else:
            ak, bk = (ones(1), ones(1))
            pm = ones(1)
        if pm_top is None:
            ktop = 0
            lev1 = fh.nz
            lev2 = 1
        else:
            ktop = find(pm>=pm_top)[0]  
            lev1 = ktop + 1     # 1-offset index for GrADS
            lev2 = fh.nz 
            ak = ak[ktop:] 
            bk = bk[ktop:] 
            pm = pm[ktop:] 

#       Find number of timesteps for each day
#       -------------------------------------
        qh = ga.query("ctlinfo")
        ntd = _nint(24. * 60. / qh.dt)

#       Set the vertical levels
#       CAUTION:
#          Because GrADS has z=1 corresponding to the bottom,
#          we use "lev" here to set the vertical dimension
#       -----------------------------------------------------
###        ga('set z %d %d'%(z1,z2))
        ga('set lev %d %d'%(lev2,lev1))  # notice order switch

#       Save the metadata info
#       ----------------------
        self.ga = ga
        self.verb = Verb
        self.fields = {}
        self.fh = fh
        self.qh = qh
        self.ntd = ntd
        self.vinfo = vinfo
        self.lon = lon          # could be None
        self.lat = lat          # could be None
        self.lev = (lev1,lev2)
        self.ak = ak
        self.bk = bk
        self.lev = pm
        self.ptop = ak[0]
        self.ktop = lev1     # 1-offset, like GrADS/Fortran

#---
    def appendVar(self,filename,Vars=None):
        """
        Open filename and append variables to the list. File
        must be conformant to the one used during _init__().
        """
        fh = ga.open(filename)

        missing = True
        if Vars is None:
            for v,k,l in fh.var_info:
                missing = False
                self.vinfo.append((v,k,l,fh.fid))
        else:
            for v,k,l in fh.var_info:
                if v in Vars:
                    missing = False
                    self.vinfo.append((v,k,l,fh.fid))
        if missing:
            print >>sys.stderr, "IndexError: requested variables - ", Vars
            raise IndexError, "cannot find any matchig variable in file %f"\
                  %filename
        
        
#---
    def interp(self,yyyy,jjj,t,lon=None,lat=None,Cache=False):
        """
        Given a date (yyyy is year, jjj the Julian day within that year), and the
        timestep t within that day, interpolate the met fields to the fire location
        for that particular time.

        You must call method prepMet() to start GrADS, open the particular file
        and select the variables of interest.

        When *Cache* is true, the variables are saved locally to directory
        "__cache__". When the input file name is "__cache__", the variables are read
        from cache.

        """

        ntd = self.ntd
        if t < 1 or t > ntd:
            raise ValueError, "invalid time index %d, maximum value is %d"%(t,ntd)

        if Cache:
            raise NotImplementedError, "cache not implemented yet, sorry."

        ga = self.ga

#       Reset lon/lat here, if specified
#       --------------------------------
        if lon is not None:
            self.lon = lon
        if lat is not None:
            self.lat = lat
        if self.lon is None or self.lat is None:
            raise ValueError, 'missing lon/lat, cannot proceed'

#       Set the requested time
#       ----------------------
        ga('set time %s'%_gatime(yyyy,jjj))
        qh = ga.query("dims")
        t_ = qh.ti[0] + t - 1
        ga('set t %d'%t_)
        ga('q time')
        now = ga.rword(1,3)

        vinfo = self.vinfo

        if self.verb > 0:
            print ""
            print "                 Met Fields' Interpolation"
            print "                 -------------------------"
            print ""

#       Loop over each desired variable and interpolate it to
#        fire location at this time
#       -----------------------------------------------------
        self.fields = {}
        nobs = self.lon.size
        for v,nlevs,l,fid in vinfo:
            if self.verb > 0:
                print "- Interpolating %6s to %5d fire locations on %s"%('<'+v+'>',nobs,now)
            ga('set dfile %d'%fid)
            self.fields[v], levs = ga.interp(v,self.lon,self.lat) # interp 
            self.fields[v] = self.fields[v].data

#............................................................................

def _gatime(yyyy,jjj,hh=None,nn=None):
    """Forms a GrADS time string from year, julian day, hour and minute"""
    t = date((yyyy,1,1)) + (jjj - 1)*DAY
    if hh is None:
        return '%02d%s%4d'%(t.day,t.ctime()[4:7],yyyy)
    else:
        return '%02d:%02dZ%02d%s%4d'%(hh,nn,t.day,t.ctime()[4:7],yyyy)

def _nint(str):
    """Nearest integer, internal use."""
    x = float(str)
    if x >= 0: return int(x+0.5)
    else:      return int(x-0.5)

