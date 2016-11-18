"""
Reads MODIS Level 2 MOD06/MYD06 Cloud Product granules and returns
a single object with the relevant data. Can read a single granule,
a directory of granules (such as for a whole day of data), or a list
of either or both of the above. Directory reads proceed recursively
to subdirectories.

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

# Time Reference used by MODIS granules
DATE_START = datetime(1993,1,1,0,0,0)

# Scientific Data Sets 
# COP = Cloud Optical Properties
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

BAD, MARGINAL, GOOD, BEST = ( 0, 1, 2, 3 )  # QA marks

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
    This class implements the MODIS Level 2 CLOUD product, usually
    referred to as MOD06 (TERRA satellite) and MYD06 (AQUA satellite).
    """

    def __init__ (self,Path,decimate='thinning',Verb=0,only_good=True):
       """
       Reads all Level 2 MOD06/MYD06 granules present on a given *Path*
       and returns a single object with all data concatenated for 5 km
       and 1 km. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories. Directories are transversed
                 recursively. If a non MOD06/MYD06 Level 2 file is
                 encountered, it is simply ignored.

       Optional parameters:
         decimate  --- decimate 1 km data, 'binning' or 'thinning' it to 5 km
         only_good --- keep only *good* observations
         Verb      --- Verbose level:
                       0 - really quiet (default)
                       1 - basic commentary, warns if invalid file is found
                       2 - detailed commentary
       """

       self.verb = Verb
       self.nobs = 0
       self.sat  = None # Satellite name
       self.coll = None # collection, e.g., 005
       self.SDS = SDS['META'] + SDS['DATA'] + SDS['COP']
       self.decimate = decimate
       self.iFilter = None  # used to filter I/O

       # Resolution dependent containers
       # -------------------------------
       self.COP = MxD06Handle('1km')

       # Create an empty list for each SDS
       # ---------------------------------
       for name in self.SDS:
           self._setName(name,[])

       # Read each granule, appending to the SDS lists
       # ---------------------------------------------
       if type(Path) is ListType:
           if len(Path) == 0:
               print "WARNING: Empty MxD06_L2 object created"
               return
       else:
           Path = [Path, ]
       self._readList(Path)

       # Concatenate each SDS list along swath to a single numpy array
       # -------------------------------------------------------------
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
                        print ">< Aliasing %s --> %s"%(sds,ALIAS[sds])
               self._aliasName(sds)

       # Create corresponding python time
       # --------------------------------
       im, jm = self.Scan_Start_Time.shape
       T = self.Scan_Start_Time[:,0]
       Time = array([DATE_START+timedelta(seconds=s) for s in T])
       self.Time = Time.repeat(jm).reshape((im,jm)) # just like Scan_Start_Time

#---
    def _setName(self,name,value):
        """Assign the value to SDS 'name'."""

        if name in SDS['COP']:
              self.COP.__dict__[name] = value
        else:
              self.__dict__[name] = value

#---
    def _appendName(self,name,value):
        """Append the value to SDS 'name'."""

        if name in SDS['COP']:
              self.COP.__dict__[name].append(value)
        else:
              self.__dict__[name].append(value)

#---
    def _catName(self,name):
        """Concatenate the granules in 'name' along swath to a single numpy array."""

        # reference the SDS 'name' as V
        # V is a list of the SDS for many granules
        if name in SDS['COP']:
              b = self.COP
        else:
              b = self
        V = b.__dict__[name]
 
        # find the *along* swath axis of the SDS
        rank = V[0].ndim
        if rank==2:
              # (along,across)
              along = 0
        elif rank==3:
              # ("channel",along,across)
              along = 1
        else:
              raise MxD06Error, 'Invalid SDS rank=%d for <%s> '%(rank,name)

        # concatenate along swath to a single numpy array
        b.__dict__[name] = concatenate(V,axis=along)

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
            if   os.path.isdir(item):    self._readDir(item)
            elif os.path.isfile(item):   self._readGranule(item)
            else:
                print "%s is not a valid file or directory, ignoring it" % item

#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if   os.path.isdir(path):    self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print "%s is not a valid file or directory, ignoring it" % item

#---
    def _readGranule(self,filename):
        """Reads one MOD06/MYD06 granule with Level 2 cloud data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb > 0:
                print " [] Working on " + filename
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 0:
                print "- %s: not recognized as an HDF file" % filename
            return 

        # Read selected variables
        # -----------------------
        for sds in self.SDS:
            if self.verb > 1:
                print ' <> Doing %s' % sds
            v = hfile.select(sds).get()
            a = hfile.select(sds).attributes()
            iOK = (v != a['_FillValue'])
            if a['scale_factor'] != 1 or a['add_offset'] != 0:
                if self.verb > 1:
                    print ' <> Scaling %s' % sds
                v[iOK] = a['scale_factor'] * (v[iOK] - a['add_offset'])
                # PS: I know that looks strange.
                #     The add_offset is what was used to get the HDF stored values within range.
                #     So in backing out the physical value, it must actually be subtracted!
            v[iOK==False] = MISSING
            self._appendName(sds,v) 

        # Decimate 1km data to lower 5km resolution (Needs better Q/C)
        # ------------------------------------------------------------
        # Note: [-1] below gets the current granule (that we just appended above)
        # TODO: Only 2D SDSs are currently implemented for decimation. This should
        #   be all we need, but otherwise add a 3D case with a "channel" loop.

        if self.decimate == 'binning':
            # average 5x5 1km pixel blocks ...

            im, jm = self.Longitude[-1].shape  # lower 5km resoln
            for name in SDS['COP']:
                  v = self.COP.__dict__[name][-1]
                  if v.ndim != 2:
                    raise MxD06Error, 'Only 2D SDSs implemented for decimation'
                  a, rc = decimateswath(im,jm,v,5,MISSING)
                  if rc: raise MxD06Error, 'ERROR: cannot decimate %s' % name
                  self.COP.__dict__[name][-1] = a

        elif self.decimate == 'thinning':
            # thin 5x5 1km pixel blocks to central pixel ...

            im, jm = self.Longitude[-1].shape
            im1 = None
            for name in SDS['COP']:
                  v = self.COP.__dict__[name][-1]
                  if v.ndim != 2:
                    raise MxD06Error, 'Only 2D SDSs implemented for decimation'
                  if im1 is None:
                        im1, jm1 = v.shape
                        I, J = (range(2,im1,5),range(2,jm1,5))
                        I, J = I[:im], J[:jm] # trim rough edges
                  a = v[I]
                  self.COP.__dict__[name][-1] = a[:,J]

        # Core Metadata
        # -------------
        cm = hfile.attributes()['CoreMetadata.0']

        # Satellite name
        # --------------
        sat = cm.split('ASSOCIATEDPLATFORMSHORTNAME')[1].split('\n')[3].split('=')[1]
        sat = sat.lstrip().replace('"','')
        if self.sat is None:
          self.sat = sat
        else:
          if sat != self.sat:
            raise MxD06Error, 'ERROR: Mixed satellites encountered'
 
        # Collection
        # ----------
        coll = int(cm.split('COLLECTION')[1].split('VERSIONID')[1].split('\n')[2].split('=')[1])
        coll = "%03d"%coll
        if self.coll is None:
            self.coll = coll
        else:
          if coll != self.coll:
            raise MxD06Error, 'ERROR: Mixed collections encountered'
        
#---
    def write(self, filename, syn_time,
              dir='.', iFilter=None, refine=8, res=None, Verb=1):
       """
        Writes gridded MODIS measurements to file.

         filename -- full path of filename to write to,
                       or None to use the standardized filename:
                           <dir>/<prod>_c<coll>.cld_Nx.yyyymmdd_hhmmz.nc4
                       In that case, provide the directory <dir> (default '.').
                       <prod> = [MOD06|MYD06] is obtained from the object's recorded satellite.
                       The collection <coll> is the object's recorded collection.

         syn_time -- python datetime appropriate for the MXD06_L2 object

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
                 0 - really quiet
                 1 - basic commentary (default)


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
       levs = linspace(1,7,7) # for now   ## <<<< needs fixing and levunits below

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

       vname  = ['CTP', 'CTT',  'BT', 'F', 'RE',     'TAU', 'CWP'   ]
       vunits = ['Pa',  'K',    'K',  '1', 'micron', '1',   'g m-2' ]
       kmvar  = [ 0,     0,     nch,   0,   0,        0,     0      ]

       if self.sat.lower() == 'terra':
         prod = 'MOD06'
       elif self.sat.lower() == 'aqua':
         prod = 'MYD06'
       else:
         raise MxD06Error, 'Unknown satellite <%s>' % self.sat

       if self.coll is None: raise MxD06Error, 'Missing collection'

       title = 'Gridded MODIS Cloud Retrievals'
       source = '%s_L2 collection %s gridded at NASA/GSFC/GMAO' % (prod, self.coll)
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = os.sep.join((dir,
             '%s_c%s.cld_Nx.%08d_%04dz.nc4' % (prod,self.coll,nymd,nhms/100)))

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

       if Verb > 0:
           print "[w] Wrote file "+filename

    def _binobs2d(self,q,im,jm):
        """Remove UNDEFs and bin it. Inout array is 2D."""
        iOK = (q!=MISSING)
        if self.iFilter is not None:
            iOK = iOK & self.iFilter 
        return binobs2d(self.lon[iOK],self.lat[iOK],q[iOK],im,jm,MISSING) 

    def _binobs3d(self,q,im,jm):
        """
        Remove UNDEFs and bin it. Inout array is 3D.
        Every channel must be defined for the the pixel to be gridded.
        """
        nc = q.shape[0]
        iBAD = (q[0]==MISSING)
        for n in range(1,nc):
            iBAD = iBAD | (q[n]==MISSING)
        iOK = (iBAD==False) # indices of good obs
        if self.iFilter is not None:
            iOK = iOK & self.iFilter 
        numOK = len(q[0][iOK])
        u = zeros((nc,numOK))
        for n in range(nc):
            u[n,:] = q[n][iOK]
        return binobs3d(self.lon[iOK],self.lat[iOK],u.T,im,jm,MISSING) 

#---
    def addVar(self,ctlfile,vname,**kwds):
        """
        Interpolates variable <vname> in <ctlfile> to pixel locations
        and saves it as an attribute named <vname>.
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
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)  # half of time window
    t1, t2 = (syn_time-dt,syn_time+dt)  # time window boundaries

    # Find MODIS granules in synoptic time range
    # ------------------------------------------
    dtGranule = timedelta(minutes=5)
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
        t += dtGranule

    if len(Granules) == 0:
        print "WARNING: no %s collection %s granules found for time" % (prod, coll), syn_time

    return Granules

#............................................................................

if __name__ == "__main__":

      Collection = '051'
      syn_time = datetime(2011,7,1,12,0,0)
      L2mount = os.sep.join((os.environ['NOBACKUP'],'MODIS',Collection,'Level2'))
      Files = granules(L2mount,'MOD06',syn_time,nsyn=24,coll=Collection)
#     m = MxD06_L2(Files,Verb=1)
      m = MxD06_L2(Files,Verb=1,decimate='binning')
      m.write(None, syn_time)

def hold():

      m = MxD06_L2('/nobackup/MODIS/051/Level2/MOD06/2012/228',Verb=1)

      Files = sorted(glob('/nobackup/MODIS/051/Level2/MOD06/2012/228/*.hdf'))
      Files = '/nobackup/MODIS/051/Level2/MOD06/2012/228/MOD06_L2.A2012228.1200.051.2012228213410.hdf'
      m = MxD06_L2(Files,decimate='thinning',Verb=1)
      m = MxD06_L2(Files,decimate='binning',Verb=1)

      aer_Nv = '/nobackup/fp/opendap/aer_Nv.ddf'
      m.addVar(aer_Nv,'DU001',kbeg=36,kount=36)

    
