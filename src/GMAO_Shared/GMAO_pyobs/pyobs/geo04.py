"""
Reads Level 2 Aerosol ABI/AHI granules produced by the MODIS Aerosol Group (Rob Levy) 
and returns a single object with the relevant data.

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov
"""

import os
import sys
from types    import *
from numpy    import zeros, ones, sqrt, std, mean, unique,\
                     concatenate, where, array, linspace,\
                     shape, arange, interp
from datetime import date, datetime, timedelta
from glob     import glob

from pyobs.npz import NPZ
from binObs_   import binobs2d, binobs3d

from netCDF4 import Dataset

try:
    from pyods import ODS # must be inported before pyhdf
except:
    pass

from bits import BITS

#---  

DATE_START = datetime(1993,1,1,0,0,0)

GROUPS = dict ( META = u'geolocation_data',
                LAND = u'geophysical_data',
               OCEAN = u'geophysical_data')

SDS = dict (
                     META = (u'longitude',
                             u'latitude',
                             u'solar_zenith_angle',
                             u'solar_azimuth_angle',
                             u'sensor_zenith_angle',
                             u'sensor_azimuth_angle',
                             u'Scattering_Angle',
                             u'Glint_Angle',
                             ),
                    LAND = ( u'Land_Ocean_Quality_Flag',
                           #  u'Cloud_Pixel_Distance_Land_Ocean',
                             u'Surface_Reflectance_Land',
                             u'Corrected_Optical_Depth_Land',
                             u'Mean_Reflectance_Land',
                             u'Aerosol_Cloud_Fraction_Land',
                             ),
                    OCEAN = (u'Land_Ocean_Quality_Flag',
                           #  u'Cloud_Pixel_Distance_Land_Ocean',
                             u'Effective_Optical_Depth_Average_Ocean',
                             u'Optical_Depth_Small_Average_Ocean',
                             u'Aerosol_Cloud_Fraction_Ocean',
                             u'Effective_Radius_Ocean',
                             u'Asymmetry_Factor_Average_Ocean',
                             u'Angstrom_Exponent_1_Ocean',
                             u'Angstrom_Exponent_2_Ocean',
                             u'Mean_Reflectance_Ocean',
                             ),
                     ALL  = (u'Land_Sea_Flag',
                             u'Aerosol_Cldmask_Land_Ocean',
                             u'Cloud_Pixel_Distance_Land_Ocean',
                             u'Land_Ocean_Quality_Flag',
                             u'Optical_Depth_Land_And_Ocean',
                             u'Image_Optical_Depth_Land_And_Ocean',
                             u'Aerosol_Type_Land',
                             u'Fitting_Error_Land',
                             u'Surface_Reflectance_Land',
                             u'Corrected_Optical_Depth_Land',
                             u'Optical_Depth_Ratio_Small_Land',
                             u'Number_Pixels_Used_Land',
                             u'Mean_Reflectance_Land',
                             u'STD_Reflectance_Land',
                             u'Mass_Concentration_Land',
                             u'Aerosol_Cloud_Fraction_Land',
                             u'Effective_Optical_Depth_Average_Ocean',
                             u'Optical_Depth_Small_Average_Ocean',
                             u'Optical_Depth_Large_Average_Ocean',
                             u'Mass_Concentration_Ocean',
                             u'Aerosol_Cloud_Fraction_Ocean',
                             u'Effective_Radius_Ocean',
                             u'PSML003_Ocean',
                             u'Asymmetry_Factor_Average_Ocean',
                             u'Backscattering_Ratio_Average_Ocean',
                             u'Angstrom_Exponent_1_Ocean',
                             u'Angstrom_Exponent_2_Ocean',
                             u'Least_Squares_Error_Ocean',
                             u'Optical_Depth_Ratio_Small_Ocean_0p55micron',
                             u'Optical_Depth_By_Models_Ocean',
                             u'Number_Pixels_Used_Ocean',
                             u'Mean_Reflectance_Ocean',
                             u'STD_Reflectance_Ocean',
                             u'Wind_Speed_Ncep_Ocean',
                             u'Topographic_Altitude_Land',
                             u'Error_Flag_Land_And_Ocean',
                             ),
               )
         

CHANNELS = dict (
                   LAND = ( 470, 550, 660 ),
                  OCEAN = ( 470, 550, 660, 870, 1200, 1600, 2100 ),
                   SREF = ( 470, 660, 2100 ),
                )

ALIAS = dict (  longitude = 'lon',
                latitude = 'lat',
                sensor_zenith_angle = 'SensorZenith',
                sensor_azimuth_angle = 'SensorAzimuth',
                solar_zenith_angle = 'SolarZenith',
                solar_azimuth_angle = 'SolarAzimuth',
                Scattering_Angle = 'ScatteringAngle',
                Glint_Angle = 'GlintAngle',
                Mean_Reflectance_Land = 'reflectance',
                Surface_Reflectance_Land = 'sfc_reflectance',
                Corrected_Optical_Depth_Land = 'aod',
                Optical_Depth_Small_Land = 'aod_fine',
                Cloud_Fraction_Land = 'cloud',
                Effective_Optical_Depth_Average_Ocean = 'aod',
                Optical_Depth_Small_Best_Ocean = 'aod_fine',
                Cloud_Fraction_Ocean = 'cloud',
                Mean_Reflectance_Ocean = 'reflectance',
             )

BAD, MARGINAL, GOOD, BEST = ( 0, 1, 2, 3 ) # QA marks

KX = dict ( G16_OCEAN = 327,
            G16_LAND  = 328,
            G17_OCEAN = 329,
            G17_LAND  = 330,
          )

KT = dict ( AOD = 45, )

IDENT = dict ( G16_OCEAN = 'g16o',
               G16_LAND  = 'g16l',
               G17_OCEAN = 'g17o',
               G17_LAND  = 'g17l',
          )

MISSING = 999.999

#...........................................................................

class GEO04_L2(object):
    """
    This class implements the MODIS Level 2 AEROSOL products, usually
    referred to as MOD04 (TERRA satellite) and MYD04 (AQUA satellite).
    """

    def __init__ (self,Path,Algo,syn_time=None,nsyn=8,Verb=0,
                  only_good=True,SDS=SDS,GROUPS=GROUPS,alias=None):
       """
       Reads individual granules or a full day of Level 2 MOD04/MYD04 files
       present on a given *Path* and returns a single object with
       all data concatenated for a given algorithm. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories.  Directories are
                 transversed recursively. If a non GEO Level 2
                 file is encountered, it is simply ignored.
         Algo -- Algorithm: LAND or OCEAN

       Optional parameters:
         syn_type  --- synoptic time
         nsyn      --- number of synoptic times per day
         only_good --- keep only *good* observations
         Verb      -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of aerosols in each file.
         SDS      --- Variables to be read from MODIS hdf files.  Must 
                      be a dictionary with keys 'META' and Algo
         ALIAS    --- dictionary of alises for SDSs

       """

       if Algo not in ('LAND', 'OCEAN'):
           raise ValueError, "invalid algorithm "+Algo+" --- must be LAND or OCEAN"

#      Initially are lists of numpy arrays for each granule
#      ------------------------------------------------
       self.verb = Verb
       self.sat  = None # Satellite name
       self.col  = None # collection, e.g., 005
       self.algo = Algo
       self.SDS  = SDS['META'] + SDS[Algo] + ('Time',)
       self.META = SDS['META']
       self.VARS = SDS[Algo]

       # Add/Substitute some aliases if given
       # ------------------------------------
       self.ALIAS = ALIAS.copy()
       if alias is not None:
           for a in alias: self.ALIAS[a] = alias[a]  
       

       # Create empty lists for SDS to be read from file
       # -----------------------------------------------
       for name in self.META+self.VARS:
           self.__dict__[name] = []
       self.__dict__['Time'] = []

       # Read each granule, appending them to the list
       # ---------------------------------------------
       if type(Path) is ListType:
           if len(Path) == 0:
               self.nobs = 0
               print "WARNING: Empty GEO04_L2 object created (1)"
               return
       else:
           Path = [Path, ]
       self._readList(Path)

       # Protect against empty GEO04 files
       # --------------------------------
       if len(self.longitude) == 0:
           self.nobs = 0
           print "WARNING: Empty GEO04_L2 object created (2)"
           return           

       # Make each attribute a single numpy array
       # ----------------------------------------
       for sds in self.SDS:
           try:
               self.__dict__[sds] = concatenate(self.__dict__[sds])
           except:
               print "Failed concatenating "+sds

       # QA flag is now readily available
       # --------------------------------
       self.qa_flag = self.Land_Ocean_Quality_Flag

       # Determine index of "good" observations
       # --------------------------------------
       if Algo == 'LAND':
           self.iGood = self.qa_flag==BEST
       elif Algo == 'OCEAN':
           self.iGood = self.qa_flag>BAD
       else:
           raise ValueError, 'invalid algorithm (very strange)'

       # Keep only "good" observations
       # -----------------------------
       if only_good:
           m = self.iGood
           for sds in self.SDS:
               rank = len(self.__dict__[sds].shape)
               if rank == 1:
                   self.__dict__[sds] = self.__dict__[sds][m]
               elif rank == 2:
                   self.__dict__[sds] = self.__dict__[sds][m,:]
               else:
                   raise IndexError, 'invalid rank=%d'%rank
           self.qa_flag = self.qa_flag[m]
           self.iGood = self.iGood[m]

       # Make aliases for compatibility with older code 
       # ----------------------------------------------
       for sds in self.SDS:
           if sds in self.ALIAS:
               self.__dict__[self.ALIAS[sds]] = self.__dict__[sds] 


       # ODS friendly attributes
       # -----------------------
       self.nobs = self.longitude.shape[0]
       self.kx = KX[self.sat+'_'+self.algo]
       self.ident = IDENT[self.sat+'_'+self.algo]
       self.Channels = CHANNELS["OCEAN"]   # all channels
       self.sChannels = CHANNELS["SREF"]   # LAND surface reflectivity (not the same as algo)
       self.channels = CHANNELS[self.algo]
       if syn_time == None:
           self.syn_time = None
           self.time = None
           self.nymd = None
           self.nhms = None
           self.nsyn = None
       else:
           Dt = [ t-syn_time for t in self.Time ]
           self.time = array([ (86400.*dt.days+dt.seconds)/60. for dt in Dt ])   # in minutes
           self.syn_time = syn_time
           self.nymd = 10000 * syn_time.year + 100*syn_time.month  + syn_time.day
           self.nhms = 10000 * syn_time.hour + 100*syn_time.minute + syn_time.second
           self.nsyn = nsyn # number of synoptic times per day

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
        """Reads one MOD04/MYD04 granule with Level 2 aerosol data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            nc = Dataset(filename)
        except:

            if self.verb > 2:
                print "- %s: cannot open"%filename
            return 

        # Read GEO location and geophysical variables
        # -------------------------------------------
        for (g, VARS) in ( (nc.groups[u'geolocation_data'], self.META),
                           (nc.groups[u'geophysical_data'], self.VARS), ):
            for sds in VARS:
                V = g.variables[sds]
                if len(V.shape) == 3:
                    v = V[:,:,:].data
                    i, j, k = v.shape
                    v = v.reshape((i*j,k))
                elif len(V.shape) == 2:
                    v = V[:,:]
                    v = v.ravel()
                    n = len(v)
                else:
                    raise IndexError, "invalid shape for SDS <%s>"%sds
                if self.verb > 2:
                    print "<> ", sds, v.shape
                if V.scale !=1. or V.add_offset != 0.:
                    v = V.scale * v + add_offset
                self.__dict__[sds].append(v)

        # Seeting time for this granule based on file name
        # ------------------------------------------------
        tyme = _gettyme(filename)
        Time = array([tyme,]*n )
        self.__dict__['Time'].append(Time)

#       Satellite name
#       --------------
        self.sat = 'G16' # hardwire this for now
            
#       Collection
#       ----------
        self.col = 0

#---

    def reduce(self,I):
        """
        Reduce observations according to index I. 
        """
        Nicknames = self.ALIAS.values()
        for name in self.__dict__:
            if name in Nicknames:
                continue # alias do not get reduced
            q = self.__dict__[name]
            if type(q) is type(self.lon):
                if len(q) == self.nobs:
                    # print "{} Reducing "+name
                    self.__dict__[name] = q[I]

        Alias = self.ALIAS.keys()
        for sds in self.SDS:
            if sds in Alias:
                self.__dict__[self.ALIAS[sds]] = self.__dict__[sds] # redefine aliases

            self.nobs = len(self.lon)

#---
    def write(self,filename=None,dir='.',expid=None,Verb=1):
        """
        Writes the un-gridded OMI object to a numpy npz file. 
        """

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        if expid == None:
            expid = self.ident

        if filename is None:
            filename = '%s/%s.omi.%d_%02dz.npz'%(dir,expid,self.nymd,self.nhms/10000)

        version = 1 # File format version
        meta = [self.nymd,self.nhms,self.nobs,self.nch,self.kx,version]
        savez(filename,
                            meta = meta,
                             lon = self.lon,
                             lat = self.lat,
                              ks = self.ks,
                        channels = self.channels,
                         qa_flag = self.qa_flag,
                     SolarZenith = self.SolarZenith,
                    SolarAzimuth = self.SolarAzimuth,
                    SensorZenith = self.SensorZenith,
                  SensorAzimuth = self.SensorAzimuth,
                 ScatteringAngle = self.ScatteringAngle,
                           cloud = self.cloud,
                             aod = self.aod,
                        aod_fine = self.aod_fine,
                     reflectance = self.reflectance)

        if Verb >=1:
            print "[w] Wrote file "+filename
#---

    def writeODS(self,filename=None,dir='.',expid=None,channels=None,
                 revised=False,nsyn=8,Verb=1):
        """
        Writes the un-gridded OMI object to an ODS file. If *revised*
        is True, the revid *aod_* parameter is written to file.
        """
        
        if self.syn_time == None:
            raise ValuError, "synoptic time missing, cannot write ODS"
            
        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        if expid == None:
            expid = self.ident

        if filename is None:
            filename = '%s/%s.obs.%d_%02dz.ods'%(dir,expid,self.nymd,self.nhms/10000)

        if channels is None:
            channels = self.channels

        # Create and populated ODS object
        # -------------------------------
        ns = self.nobs
        nobs = len(channels) * ns
        ods = ODS(nobs=nobs, kx=self.kx, kt=KT['AOD'])
        i = 0
        ks = arange(ns) + 1
        for ch in channels:
            I = range(i,i+ns)
            j = channels.index(ch)
            ods.ks[I]  = ks
            ods.lat[I] = self.lat[:]
            ods.lon[I] = self.lon[:]
            ods.time[I] = self.time[:].astype('int')
            ods.lev[I] = ch
            ods.qch[I] = self.qa_flag[:].astype('int')
            if revised:
                ods.obs[I]  = self.aod_[:,j].astype('float32')
                ods.xvec[I] = self.aod[:,j].astype('float32')
            else:
                ods.obs[I] = self.aod[:,j].astype('float32')
            i += ns

        # Handle corrupted coordinates
        # ----------------------------
        iBad = (ods.lon<-180) | (ods.lon>180.) | \
               (ods.lat<-90)  | (ods.lat>90.)  | \
               (abs(ods.time)>1440./self.nsyn)
        ods.lon[iBad] = 0.
        ods.lat[iBad] = -90.
        ods.time[iBad] = 0.
        ods.qcx[iBad] = 2

        # Exclusion flag
        # --------------
        iGood = (ods.qch>0) & (ods.obs<10.) & (ods.qcx==0)
        ods.qcx[:] = 1     # All bad...
        ods.qcx[iGood] = 0 # ... unless good

        ods_ = ods.select(qcx=0)
        if Verb >=1:
            print "[w] Writing file <"+filename+"> with %d observations"%ods_.nobs

        ods_.write(filename,self.nymd,self.nhms,nsyn=8,ftype='pre_anal')
        
#---
    def writeg(self,filename=None,dir='.',expid=None,refine=8,res=None,
               channels=None,Verb=1):
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
       
       # Stop here is no good obs available
       # ----------------------------------
       if self.nobs == 0:
           return # no data to work with
       if any(self.iGood) == False:
           return # no good data to work with

       if expid == None:
           expid = self.ident

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

       if channels is None:
           channels = self.channels
 
       levs = array(channels)

       nch = len(channels)
       nymd = self.nymd
       nhms = self.nhms

       vtitle = [ 'Aerosol Optical Depth',
                  'Aerosol Optical Depth (Revised)',
                  'Aerosol Optical Depth (Fine Mode)',
                  'Cloud Fraction' ]

       vname  = ['tau', 'tau_', 'tau_fine', 'cloud' ]
       vunits = [ '1',    '1',     '1',       '1',  ]
       kmvar  = [ nch,    nch,     nch,        0    ]

       title = 'Gridded MODIS Aerosol Retrievals'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
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

       # Subset AOD at specified channels
       # --------------------------------
       I = []
       for ch in channels:
           i = list(self.channels).index(ch)
           I = I + [i,]
       aod = self.aod[:,I]

       # Fine mode
       # ---------
       try:
           aod_fine = self.aod_fine[:,I]
       except:
           aod_fine = MISSING * ones(aod.shape) # will compress like a charm
       
       # The Revised AOD may not exist
       # -------------------------------
       try:
           aod_ = self.aod_[:,I]
       except:
           aod_ = MISSING * ones(aod.shape) # will compress like a charm

       # Grid variable and write to file
       # -------------------------------
       f.write('tau', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod,im,jm,MISSING) )
       f.write('tau_', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod_,im,jm,MISSING) )
       f.write('tau_fine', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod_fine,im,jm,MISSING) )
       f.write('cloud', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.cloud,im,jm,MISSING) )
           
#       try:
#           f.close()
#       except:
#           pass

       if Verb >=1:
           print "[w] Wrote file "+filename

#---
    def addVar(self,ga,expr='mag(u10m,v10m)',vname='wind',clmYear=None,tight=True):
        """
        Given a grads object *ga* having the correct MERRA file as default,
        interpolates *var* to obs location and saves it as an attribute
        named *vname*.

        If *tight* is True, domain will be restricted conserve memory. This feature
        has proven somewhat unstable for reasons yet TBD.
        """

        U = MISSING * ones(self.nobs)
        if vname == None:
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
    def getCoxMunk(self,filename='/nobackup/NNR/Misc/coxmunk_lut.npz',channel=550):
        """
        Returns ocean albedo.
        """
        
        # Get precomputed albedo LUT
        # --------------------------
        lut = NPZ(filename)
        
        # Trimmed wind speed
        # ------------------
        w10m = self.wind.copy()
        w10m[w10m<0] = 0
        w10m[w10m>50.] = 50.

        j = list(lut.channels).index(channel)

        # Interpolate albedo
        # ------------------
        albedo = zeros(len(w10m))
        albedo[:] = interp(w10m,lut.speed,lut.albedo[:,j])

        self.albedo = albedo

#............................................................................

def _gettyme(pathname):
    filename = os.path.basename(pathname)
    yyyyjjj = filename.split('_')[2].split('.')[0]
    hhmm = filename.split('.')[1]
    year = int(yyyyjjj[0:4])
    jjj = int(yyyyjjj[4:7])
    hour = int(hhmm[0:2])
    minute = int(hhmm[2:4])
    dt = timedelta(seconds=(jjj-1)*24*60*60)  
    return datetime(year,1,1,hour,minute)+dt

def granules ( path, syn_time, nsyn=8 ):
    """
    Returns a list of GEO04 granules for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the GEO04 Level 2 files
    syn_time  ---  synoptic time (timedate format)

    nsyn      ---  number of synoptic times per day (optional)

    """

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
    t1, t2 = (syn_time-dt,syn_time+dt)

    # Find MODIS granules in synoptic time range
    # ------------------------------------------
    dt = timedelta(minutes=15)
    t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
    Granules = []
    while t < t2:
        if t >= t1:
            doy = t.timetuple()[7]
            basen = "%s/%s/%s/ABI_L2_%04d%03d.%02d%02d.nc"\
                     %(path,t.year,doy,t.year,doy,t.hour,t.minute)
            try:
                print basen
                filen = glob(basen)[0]
                Granules += [filen,]
#               print " [x] Found "+filen
            except:
                pass
        t += dt

    if len(Granules) == 0:
        print "WARNING: no granules found for", syn_time

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

#     path = '/nobackup/1/GEO/ABI_Testing/Level2/2018/227'
#    files = [path+'/ABI_L2_2018227.1945.nc',
#              path+'/ABI_L2_2018227.2000.nc']

    
     path = '/nobackup/1/GEO/ABI_Testing/Level2'
     syn_time = datetime(2018,8,15,18,0)
     files = granules ( path, syn_time, nsyn=8 )

     go = GEO04_L2 (files,'OCEAN',syn_time=None,nsyn=8,Verb=3)
     gl = GEO04_L2 (files,'LAND',syn_time=None,nsyn=8,Verb=3)



