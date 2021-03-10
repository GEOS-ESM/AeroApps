"""
 The *pyods* module implements a Python interface to Observation Data 
Streams (ODS). It supports both traditional HDF-4/NetCDF ODS files as
well as GSI binary diag output.

"""
from types   import *
from pyods_  import *   # Fortran extension using f2py
from odsmeta import *   # useful constants
from numpy   import zeros, ones, arange, shape

__VERSION__ = '1.0.3'

ODS_UNDEF  = 1.0e15

class ODS(object):

   def __init__(self,filename=None,nymd=None,nhms=None,
                nobs=None,kx=None,kt=None,only_good=False):
      """
      Initializes the ODS class. If the date/time parameters (nymd,nhms)
      are specified the data for that particular time is read in. Otherwise,
      the file is simply opened. You must specify (nymd,nhms) for reading 
      GSI diag output.

      When *filename*is omitted and *nobs* is specified, an allocated ods object
      is returned. In this case, optional constants kt, kx can be specified for these
      attributes.
      
      The usual ODS attributes are defined:
      
           nobs --- nuber of observations for this date/time
           kid  --- observtion ID (sequential)
           lat  --- latitude (degrees)
           lon  --- longitude (degrees)
           lev  --- vertical level (usually in hPa)
           kx   --- data source index
           kt   --- data type index
           ks   --- sounding index
           xm   --- metadata index
           time --- time relative to synoptic time
           obs  --- observed value (units depend on kt)
           omf  --- innovation (obs minus background residual)
           oma  --- O-A (obs minus background residual)
           xvec --- tradiationally, PSAS intermediate solution
           qcx  --- QC exclusion flag (0 is good)
           qch  --- QC history flag

      When only_good=True only those observations with qcx=0 are returned.

      Examples:
      
            import pyods

      1) Read ODS from a file:
      
            misr = pyods.ODS('misr.ods',20080629,120000)
            print misr.nobs
            
      2) Create an ODS with space for 20,000 observations

            ods = pyods.ODS(nobs=20000)

      """

      self.fid = None
      self.kid = None
      self.nymd = None # to be set on read
      self.nhms = None # to be set on read

      if filename is None:
         if nobs != None:
            self._alloc(nobs,kx=kx,kt=kt)
         return
         
      if (nymd is None) or (nhms is None):
         self.fid, rc = pyods_open(filename)
         if rc:
            self.fid = None
            raise IOError, "cannot open ODS file" + filename
      else:
         self._Read(filename,nymd,nhms,only_good) 
         self.fid=None

   def __del__(self):
      if self.fid is not None:
         rc = pyods_close(fid)
         if rc:
            raise IOError, "could not close file"
      self.fid = None

   def nget(self,nymd,nhms,rc):
      """ Returns number of observations for this date/time."""
      if self.fid is None:
         raise IOError, "ODS file is not open"
      self.nobs, rc = pyods_nget(self.fid,nymd,nhms)
      if rc:
         raise IOError, "could not determine NOBS"
 
   def getInt(self,name,nymd,nhms):
      if self.fid is None:
         raise IOError, "ODS file is not open"
      self.nget(nymd,nhms)
      v, rc = pyods_getint(self.fid,name,nymd,nhms,self.nobs)
      if rc:
         raise IOError, "cannot get "+name
      return v

   def getFloat(self,name,nymd,nhms):
      if self.fid is None:
         raise IOError, "ODS file is not open"
      self.nget(nymd,nhms)
      v, rc = pyods_getfloat(self.fid,name,nymd,nhms,self.nobs)
      if rc:
         raise IOError, "cannot get "+name
      return v

   def select(self,kt=None,kx=None,qcx=None,qch=None,ks=None,lev=None):
      """Returns a copy of the ODS object restricting to the
         specified values of the arguments. For example,

                subset = ods.select(kt=(45,46))

         returns those observations having kt=45 or kt=46.
      """

      if (self.kid is None) or (self.nobs==0):
         return None

      I = (range(self.nobs) >  -1) # all True
      I = (I & _getMatchingIndices(self.kt,kt,  'kt'))
      I = (I & _getMatchingIndices(self.kx,kx,  'kx'))
      I = (I & _getMatchingIndices(self.ks,ks,  'ks'))
      I = (I & _getMatchingIndices(self.qcx,qcx,'qcx'))
      I = (I & _getMatchingIndices(self.qch,qch,'qch'))
      I = (I & _getMatchingIndices(self.lev,lev,'lev'))

      return self.__copy__(I)

#---
   def write(self,filename,nymd,nhms,nsyn=4,ftype='post_anal',only_good=False):
      """
      Writes an ODS file at a given date (nymd) and time (nhms).
      """

      nobs = self.nobs

      if only_good:
         I = (self.qcx[0:nobs] == 0)
      else:
         I = range(self.nobs)

      kid  = self.kid[I]
      lat  = self.lat[I]
      lon  = self.lon[I]
      lev  = self.lev[I]
      kx   = self.kx[I]
      kt   = self.kt[I]
      ks   = self.ks[I]
      xm   = self.xm[I]
      time = self.time[I]
      obs  = self.obs[I]
      omf  = self.omf[I]
      oma  = self.oma[I]
      xvec = self.xvec[I]
      qcx  = self.qcx[I]
      qch  = self.qch[I]

      rc = pyods_putall(filename,ftype,nymd,nhms,nsyn,
                        kid,lat,lon,lev,kx,kt,ks,xm,time,
                        obs,omf,oma,xvec,qcx,qch)

      if rc:
         print "On return from pyods_putall, rc = ", rc
         raise IOError, "cannot write ODS file "

#---
   def addVar(self,ga,expr='mag(u10m,v10m)',vname=None,clmYear=None,tight=True):
        """
        Given a grads object *ga* having the correct file as default,
        interpolates *var* to obs location and saves it as an attribute
        named *vname*.

        If *tight* is True, domain will be restricted conserve memory. This feature
        has proven somewhat unstable for reasons yet TBD.
        """

        U = ODS_UNDEF * ones(self.nobs)
        if vname is None:
            vname = expr

        # nearest time
        # ------------
        t = _gatime(self.nymd,self.nhms)
        if clmYear != None:
            t = t[:-4] + str(clmYear) # replace year
        ga('set time '+t,Quiet=True)

        # To conserve memory, restrict domain with 1 gridpoint halo
        # ---------------------------------------------------------
        if tight:
            fh = ga.query("file")
            x1, x2  = self.lon.min(),self.lon.max()
            y1, y2  = self.lat.min(),self.lat.max()
            ga('set lon %f %f'%(x1,x2),Quiet=True)
            ga('set lat %f %f'%(y1,y2),Quiet=True)
            qh = ga.query("dims")
            x1, x2 = (qh.xi[0]-1,qh.xi[1]+1)
            y1, y2 = (max(1,qh.yi[0]-1),min(fh.ny,qh.yi[1]+1)) # in [1,ny]
            ga('set x %d %d'%(x1,x2),Quiet=True) # make sure x range is int
            ga('set y %d %d'%(y1,y2),Quiet=True) # make sure y range is int
            expr_ = ga.exp(expr)
        else:
            expr_ = ga.expr(expr)
        u, levs = ga.interp(expr_, self.lon, self.lat )
        U = u.data
        if len(shape(U)) == 0:
             U = U * ones(1) # so that we can slice it later

        self.__dict__[vname] = U

 #---
   def _Read(self,filename,nymd,nhms,only_good=False):
      """
      Reads an ODS file at a given date (nymd) and time (nhms).
      This works for both traditional ODS as well as for GSI
      diag output.
      """
      nobs, kid, lat, lon, lev, \
      kx, kt, ks, xm, time, \
      obs, omf, oma, xvec, qcx, qch, rc \
      = pyods_getall(filename,nymd,nhms)
      if rc:
         print "On return from pyods_getall, rc = ", rc
         raise IOError, "cannot read ODS file "
 
      self.nobs = nobs
      self.nymd = nymd
      self.nhms = nhms

      if only_good:
         I = (qcx[0:nobs] == 0)
      else:
         I = range(self.nobs)

      self.kid  = kid[I]
      self.lat  = lat[I]
      self.lon  = lon[I]
      self.lev  = lev[I]
      self.kx   = kx[I]
      self.kt   = kt[I]
      self.ks   = ks[I]
      self.xm   = xm[I]
      self.time = time[I]
      self.obs  = obs[I]
      self.omf  = omf[I]
      self.oma  = oma[I]
      self.xvec = xvec[I]
      self.qcx  = qcx[I]
      self.qch  = qch[I]

#---

   def _alloc(self,nobs,kx=None,kt=None):
      """
      Create ODS object of size *nobs*. If (kx,kt) are specified the
      corresponding attributes are set to these constants. Othwerwise,
      except for *kid*, all attributes are set to zero.
      """
      if kx is None: kx = 0
      if kt is None: kt = 0
      self.nobs = nobs
      self.kid  = (1+arange(nobs)).astype('int')
      self.kx   = (kx * ones(nobs)).astype('int')
      self.kt   = (kt * ones(nobs)).astype('int')
      self.xm   = zeros(nobs)
      self.time = zeros(nobs).astype('int')
      self.omf  = zeros(nobs)
      self.oma  = zeros(nobs)
      self.xvec = zeros(nobs)
      self.qcx  = zeros(nobs).astype('int')
      self.qch  = zeros(nobs).astype('int')
      
      self.ks   = zeros(nobs).astype('int')
      self.lat  = zeros(nobs)
      self.lon  = zeros(nobs)
      self.lev  = zeros(nobs)
      self.obs  = zeros(nobs)
      
#---
   def __copy__(self,Indices=None):
      """Implements copy operation; if Indices are specified, a subset
         based on these indices are returned."""
 
      ods = ODS()
      ods.fid = self.fid
      ods.kid = self.kid
      ods.nymd = self.nymd
      ods.nhms = self.nhms

      if self.kid is None:
         return ods

      nobs = self.nobs

#     Copy data attributes
#     --------------------
      if nobs > 0:

         if Indices is None:
            I = range(nobs)
         else:
            I = Indices

         ods.kid  = self.kid[I]
         ods.nobs = ods.kid.size
         ods.lat  = self.lat[I]
         ods.lon  = self.lon[I]
         ods.lev  = self.lev[I]
         ods.kx   = self.kx[I]
         ods.kt   = self.kt[I]
         ods.ks   = self.ks[I]
         ods.xm   = self.xm[I]
         ods.time = self.time[I]
         ods.obs  = self.obs[I]
         ods.omf  = self.omf[I]
         ods.oma  = self.oma[I]
         ods.xvec = self.xvec[I]
         ods.qcx  = self.qcx[I]
         ods.qch  = self.qch[I]

#     All done
#     --------
      return ods

#---

def _gatime(nymd,nhms):
        Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
        cymd = "%08d"%int(nymd)
        chms = "%06d"%int(nhms)
        t = chms[0:2]+":"+chms[2:4]+"Z"+\
            cymd[6:8]+Months[int(cymd[4:6])-1]+cymd[0:4]
        return t

def _getMatchingIndices(att,kt,name):
   """ Return a tuple."""
   if kt is None:
      kt_list = ()
      J = ( range(att.size) >  -1) # all True
      return J
   if type(kt) is TupleType or type(kt) is ListType:
      kt_list = kt
   else:
      kt_list = (kt,)

   J = (range(att.size) <  0) # all False
   for k in kt_list:
      J = (J | (att==k) ) # True if matching kt
   return J

      

      

