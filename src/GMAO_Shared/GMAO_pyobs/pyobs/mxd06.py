"""
Reads Level 2 MOD04/MYD04 granules for a single day and returns a
single object with the relevant data.

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov
"""

import os
import sys
from types     import *
from numpy     import zeros, ones, sqrt, std, mean, unique,\
                      concatenate, where, array, linspace,\
                      shape, arange, interp
from datetime  import date, datetime, timedelta
from glob      import glob

from pyobs.npz import NPZ
from binObs_   import binobs2d, binobs3d, decimateswath

from pyhdf.SD import SD, HDF4Error

from bits import BITS

#---  

DATE_START = datetime(1993,1,1,0,0,0)

SDS = dict (
      META = ('Longitude', 'Latitude', 'Scan_Start_Time',
#             'Sensor_Zenith', 'Sensor_Azimuth',
#             'Solar_Zenith',  'Solar_Azimuth',
            ),
      DATA = ( 'Cloud_Top_Pressure', 
               'Cloud_Top_Temperature',
               'Brightness_Temperature',
#              'Cloud_Phase_Infrared',
               'Cloud_Height_Method',
               'Cloud_Fraction',
               'Cloud_Mask_5km',
               'Surface_Type',
            ),
      COP = ( 'Cloud_Effective_Radius',
              'Cloud_Optical_Thickness',
              'Cloud_Water_Path',
              'Cloud_Phase_Optical_Properties',
              'Cloud_Multi_Layer_Flag',
              'Cloud_Effective_Radius_Uncertainty',
 #            'Cloud_Mask_1km',
 #            'Quality_Assurance_1km',
            ),
      )

ALIAS = dict (              Longitude = 'lon',
                             Latitude = 'lat',
                        Sensor_Zenith = 'SensorZenith',
                       Sensor_Azimuth = 'SensorAzimuth',
                         Solar_Zenith = 'SolarZenith',
                        Solar_Azimuth = 'SolarAzimuth',
                     Scattering_Angle = 'ScatteringAngle',
                   Cloud_Top_Pressure = 'CTP',  
                Cloud_Top_Temperature = 'CTT',
               Brightness_Temperature = 'BT',
                       Cloud_Fraction = 'F',
               Cloud_Effective_Radius = 'RE',
              Cloud_Optical_Thickness = 'TAU',
                     Cloud_Water_Path = 'CWP',
       Cloud_Phase_Optical_Properties = 'PHASE',
               Cloud_Multi_Layer_Flag = 'MLAYER',
   Cloud_Effective_Radius_Uncertainty = 'eRE',
             )

BAD, MARGINAL, GOOD, BEST = ( 0, 1, 2, 3 ) # QA marks

MISSING = -999.999

#...........................................................................
class MxD06Handle(object):
    """
    Generic container for MxD06.
    """
    def __init__(self,name):
        self.name = name
        self.nobs = 0
        
class MxD06Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
  
class MxD06_L2(object):
    """
    This class implements the MODIS Level 2 CLOUDS products, usually
    referred to as MOD06 (TERRA satellite) and MYD06 (AQUA satellite).
    """

    def __init__ (self,Path,decimate='thinning',Verb=0,only_good=True):
       """
       Reads individual granules or a full day of Level 2 MOD06/MYD06 files
       present on a given *Path* and returns a single object with
       all data concatenated for 5 km and 1 km. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories.  Directories are
                 transversed recursively. If a non MOD06/MYD06 Level 2
                 file is encountered, it is simply ignored.

       Optional parameters:
         decimate  --- decimate 1 km data, 'binning' or 'thinning' it to 5 km
         only_good --- keep only *good* observations
         Verb      --- Verbose level:
                       0 - really quiet (default)
                       1 - Warns if invalid file is found
                       2 - Prints out non-zero number of aerosols in each file.

       """

#      Initially are lists of numpy arrays for each granule
#      ------------------------------------------------
       self.verb = Verb
       self.nobs = 0
       self.sat  = None # Satellite name
       self.col  = None # collection, e.g., 005
       self.SDS = SDS['META'] + SDS['DATA'] + SDS['COP']
       self.decimate = decimate
       self.iFilter = None # used to filter I/O

       # Resolution dependent containers
       # -------------------------------
       self.COP = MxD06Handle('1km')

       # Create empty lists for SDS to be read from file
       # -----------------------------------------------
       for name in self.SDS:
           self._setName(name,[])

       # Read each granule, appending them to the list
       # ---------------------------------------------
       if type(Path) is ListType:
           if len(Path) == 0:
               print "WARNING: Empty MxD06_L2 object created"
               return
       else:
           Path = [Path, ]
       self._readList(Path)

       # Make each attribute a single numpy array
       # ----------------------------------------
       for sds in self.SDS:
           try:
               self._catName(sds)
               if self.verb>1:
                        print "() Concatenated <%s>"%sds
           except:
               print "WARNING: Failed concatenating "+sds

       # If decimating, keep attributes flat
       # -----------------------------------
       if decimate in ( 'thinning', 'binning' ):
             for name in SDS['COP']:
                   self.__dict__[name] = self.COP.__dict__[name] 
                   if self.verb>1:
                        print "/\ Moved <%s> up"%name
       
       # Make aliases for compatibility with older code 
       # ----------------------------------------------
       Alias = ALIAS.keys()
       for sds in self.SDS:
           if sds in Alias:
               if self.verb>1:
                        print ">< Aliasing %16s --> %s"%(sds,ALIAS[sds])
               self._aliasName(sds)

       # Create corresponding python time
       # --------------------------------
       im, jm = self.Scan_Start_Time.shape
       T = self.Scan_Start_Time[:,0]
       Time = array([DATE_START+timedelta(seconds=s) for s in T])
       self.Time = Time.repeat(jm).reshape((im,jm)) # just like Scan_Start_Time

#---
    def _setName(self,name,value):
        if name in SDS['COP']:
              self.COP.__dict__[name] = value
        else:
              self.__dict__[name] = value
#---
    def _appendName(self,name,value):
        if name in SDS['COP']:
              self.COP.__dict__[name].append(value)
        else:
              self.__dict__[name].append(value)
#---
    def _catName(self,name):

        if name in SDS['COP']:
              b = self.COP
        else:
              b = self
              
        V = b.__dict__[name]
        if len(V[0].shape)==2:
              b.__dict__[name] = concatenate(V)
        elif len(V[0].shape)==3:
              w = None # will hold 3D array: (channel,along,cross)
              nc = V[0].shape[0]
              for c in range(nc):
                    U = []  # 2D list
                    for v in V:
                          U.append(v[c])
                    u = concatenate(U) # concatenated 2D array
                    if w is None:
                          w = ones((nc,) + u.shape) # only now know shape
                    w[c] = u # 3D array
              b.__dict__[name] = w
        else:
              raise MxD06Error, 'Invalid SDS rank=%d for <%s> '%(len(v.shape),name)
        

#---
    def _aliasName(self,name):
        if self.decimate is None and name in SDS['COP']:
              self.COP.__dict__[ALIAS[name]] = self.COP.__dict__[name] 
        else:
              self.__dict__[ALIAS[name]] = self.__dict__[name] 

#---
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
        """Reads one MOD06/MYD06 granule with Level 2 cloud data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print " [] Working on "+filename
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 2:
                print "- %s: not recognized as an HDF file"%filename
            return 

        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
        for sds in self.SDS:
            # print ' <> Doing %s'%sds
            v = hfile.select(sds).get()
            a = hfile.select(sds).attributes()
            I = (v!=a['_FillValue'])
            if a['scale_factor']!=1 or a['add_offset']!=0:
                  # print ' <> Scaling %s'%sds
                v[I] = a['scale_factor'] * (v[I] - a['add_offset']) # strange
            v[I==False] = MISSING
            self._appendName(sds,v) 

        # Decimate 1km data (Needs better Q/C)
        # ------------------------------------
        if self.decimate == 'binning':
            b = self.COP
            im, jm = self.Longitude[-1].shape
            idel = 5
            for name in SDS['COP']:
                  v = self.COP.__dict__[name][-1]
                  a, rc = decimateswath(im,jm,v,idel,MISSING)
                  if rc:
                        raise MxD06Error, 'ERROR: cannot decimate %s'%name
                  self.COP.__dict__[name][-1] = a
        elif self.decimate == 'thinning':
            im, jm = self.Longitude[-1].shape
            im1 = None
            for name in SDS['COP']:
                  v = self.COP.__dict__[name][-1]
                  if im1 is None:
                        im1, jm1 = v.shape
                        I, J = (range(2,im1,5),range(2,jm1,5))
                        I, J = I[:im], J[:jm] # trimm rough edges
                  a = v[I]
                  self.COP.__dict__[name][-1] = a[:,J] # TO DO: check for rank

        # Core Metadata
        # -------------
        cm = hfile.attributes()['CoreMetadata.0']

#       Satellite name
#       --------------
        if self.sat is None:
            sat = cm.split('ASSOCIATEDPLATFORMSHORTNAME')[1].split('\n')[3].split('=')[1]
            self.sat = sat.lstrip().replace('"','')
            
#       Collection
#       ----------
        if self.col is None:
            col = int(cm.split('COLLECTION')[1].split('VERSIONID')[1].split('\n')[2].split('=')[1])
            self.col = "%03d"%col
        
#---
    def write(self,filename, syn_time, iFilter=None, 
              refine=8,res=None, Verb=1):
       """
        Writes gridded MODIS measurements to file.

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a   4  x  5    grid
                       refine=2  produces a   2  x2.50   grid
                       refine=4  produces a   1  x1,25   grid
                       refine=8  produces a  0.50x0.625  grid
                       refine=16 produces a  0.25x0.3125 grid
        Alternatively, one can specify the grid resolution with a
        single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a   4  x  5    grid
                       res='b'  produces a   2  x2.50   grid
                       res='c'  produces a   1  x1,25   grid
                       res='d'  produces a  0.50x0.625  grid
                       res='e'  produces a  0.25x0.3125 grid

                   NOTE: *res*, if specified, supersedes *refine*.

         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.


       """
       from gfio import GFIO

       # Optional filter
       # ---------------
       self.iFilter = iFilter
       
#      Output grid resolution
#      ----------------------
       if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16

#      Lat lon grid
#      ------------
       dx = 5. / refine
       dy = 4. / refine
       im = int(360. / dx)
       jm = int(180. / dy + 1)

       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)

       nch = 7
       levs = linspace(1,7,7) # for now

       t = syn_time
       Y, M, D, h, m, s = (t.year,t.month,t.day,t.hour,t.minute,t.second)

       nymd = Y*10000+M*100+D
       nhms = h*10000+m*100+s

       vtitle = [ 'Cloud_Top_Pressure', 
                  'Cloud_Top_Temperature',
                  'Brightness_Temperature',
                  'Cloud_Fraction',
                  'Cloud_Effective_Radius',
                  'Cloud_Optical_Thickness',
                  'Cloud_Water_Path',
                 ]

       vname  = ['CTP', 'CTT', 'BT', 'F', 'RE',     'TAU', 'CWP'   ]
       vunits = [ 'Pa',   'K',   'K',  '1', 'micron', '1',   'g m-2' ]
       kmvar  = [  0,      0,    nch,   0,   0,        0,     0      ]

       title = 'Gridded MODIS Cloud Retrievals'
       source = 'MODIS PG06 gridded at NASA/GSFC/GMAO'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.sfc.%d_%02dz.nc4'%(dir,expid,self.nymd,self.nhms/10000)

       # Create the file
       # ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=levs, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # Grid variables and write to file
       # -------------------------------
       f.write ('CTP', nymd, nhms, self._binobs2d(self.CTP,im,jm) )
       f.write ('CTT', nymd, nhms, self._binobs2d(self.CTT,im,jm) )
       f.write ('BT',  nymd, nhms, self._binobs3d(self.BT, im,jm) )
       f.write ('F',   nymd, nhms, self._binobs2d(self.F,  im,jm) )
       f.write ('RE',  nymd, nhms, self._binobs2d(self.RE, im,jm) )
       f.write ('TAU', nymd, nhms, self._binobs2d(self.TAU,im,jm) )
       f.write ('CWP', nymd, nhms, self._binobs2d(self.CWP,im,jm) )

       f.close()

       # Reset filter
       # ------------
       self.iFilter = None

       if Verb >=1:
           print "[w] Wrote file "+filename

    def _binobs2d(self,q,im,jm):
        """Remove UNDEFs and bin it. Inout array is 2D."""
        I = (q!=MISSING)
        if self.iFilter is not None:
            I = I & self.iFilter 
        return binobs2d(self.lon[I],self.lat[I],q[I],im,jm,MISSING) 

    def _binobs3d(self,q,im,jm):
        """Remove UNDEFs and bin it. Inout array is 3D."""
        nc = q.shape[0]
        J = (q[0]==MISSING)
        for n in range(1,nc):
            J = J | (q[n]==MISSING)
        I = (J==False) # indices of good obs
        if self.iFilter is not None:
            I = I & self.iFilter 
        N = len(q[0][I])
        u = zeros((nc,N))
        for n in range(nc):
            u[n,:] = q[n][I]
        return binobs3d(self.lon[I],self.lat[I],u.T,im,jm,MISSING) 

#---
    def addVar(self,ctlfile,vname,**kwds):
        """
        Interpolates *var* to obs location and saves it as an attribute
        named *vname*.
        """
        from gfio import GFIOctl
        f = GFIOctl(ctlfile)
        self.__dict__[vname] = f.interpXY(vname,self.lon,self.lat,**kwds)
        f.close()
        
#............................................................................

def granules ( path, prod, syn_time, coll='051', nsyn=8 ):
    """
    Returns a list of MxD06 granules for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the MxD06 Level 2 files
    prod      ---  either MOD06 or MYD06
    syn_time  ---  synoptic time (timedate format)

    coll      ---  collection: 005, 051 (optional)
    nsyn      ---  number of synoptic times per day (optional)

    """

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
    t1, t2 = (syn_time-dt,syn_time+dt)

    # Find MODIS granules in synoptic time range
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
        print "WARNING: no %s collection %s granules found for"%(prod,coll), syn_time

    return Granules

#--

def print_stats(name,x=None):
    "Prints simple stats"
    from pylab import prctile
    if type(name) is not StringType:
        x = name
        name = 'mean,stdv,rms,min,25%,median,75%,max: '
    if name == '__header__':
        print ''
        n = (80 - len(x))/2
        print n * ' ' + x
        print n * ' ' + len(x) * '-'
        print ''
        print '   Name       mean      stdv      rms      min     25%    median     75%      max'
        print ' ---------  -------  -------  -------  -------  -------  -------  -------  -------'
    elif name == '__sep__':
        print ' ---------  -------  -------  -------  -------  -------  -------  -------  -------'
    elif name == '__footer__':
        print ' ---------  -------  -------  -------  -------  -------  -------  -------  -------'
        print ''
    else:
        ave = x.mean()
        std = x.std()
        rms = sqrt(ave*ave+std*std)
        prc = prctile(x)
        print '%10s  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  '%\
            (name,ave,std,rms,prc[0],prc[1],prc[2],prc[3],prc[4])

#--

def _gatime(nymd,nhms):
        Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
        cymd = "%08d"%int(nymd)
        chms = "%06d"%int(nhms)
        t = chms[0:2]+":"+chms[2:4]+"Z"+\
            cymd[6:8]+Months[int(cymd[4:6])-1]+cymd[0:4]
        return t
    
#............................................................................

if __name__ == "__main__":

      syn_time = datetime(2012,8,15,12,0,0)
      Files = granules('./Level2','MOD06',syn_time,coll='051')
      m = MxD06_L2(Files,Verb=1)

      aer_Nv = '/nobackup/fp/opendap/aer_Nv.ddf'
      #m.addVar(aer_Nv,'DU001',kbeg=36,kount=36)

def hold():
      m = MxD06_L2('/nobackup/MODIS/051/Level2/MOD06/2012/228',Verb=1)

      Files = sorted(glob('/nobackup/MODIS/051/Level2/MOD06/2012/228/*.hdf'))
      Files = '/nobackup/MODIS/051/Level2/MOD06/2012/228/MOD06_L2.A2012228.1200.051.2012228213410.hdf'
      m = MxD06_L2(Files,decimate='thinning',Verb=1)
      m = MxD06_L2(Files,decimate='binning',Verb=1)

    
