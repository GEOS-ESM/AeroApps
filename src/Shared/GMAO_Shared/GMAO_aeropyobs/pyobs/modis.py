"""
  MODIS class based on ODS.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
from pyods    import ODS
from binObs_  import binobs2d, binobs3d
from gfio     import GFIO
from numpy    import save, zeros, ones, linspace, any, \
                     shape, arange, savez, load, array

kxTERRA_OCEAN = 301
kxTERRA_LAND = 302
kxAQUA_OCEAN = 311
kxAQUA_LAND = 312
kxAQUA_DEEP = 320

kxMODO = kxTERRA_OCEAN
kxMODL = kxTERRA_LAND
kxMYDO = kxAQUA_OCEAN
kxMYDL = kxAQUA_LAND
kxDEEP = kxAQUA_DEEP

ktAOD = 45
ktANGSTROM = 48
ktREFLECTANCE = 48
ktSOLARZENITH = 54
ktSOLARAZIMUTH = 55
ktSENSORZENITH = 56
ktSENSORAZIMUTH = 56
ktSCATTERINGANGLE = 58
ktATYPE = 60
ktAOD_FINE = 61
ktFINE_FRACTION = 70
ktQA_FLAG = 74
ktCLOUD = 90

MISSING = 999.999

class MODIS(object):

    def __init__ (self,filename,nymd=None,nhms=None,kx=None,only_good=True):
        """
        Reads an ODS or NPZ file and creates a MODIS object with the following
        attributes:

             nobs           ---  number of profiles
             nch            ---  number of channels
             lon            ---  longitudes in degrees        (nobs)
             lat            ---  latitudes in degrees         (nobs)
             channels       ---  channels in nm               (nch)
             qa_flag        ---  quality assurance flag       (nobs)

             SolarZenith    ---  solar zenith angle           (nobs)
             SolarAzimuth   ---  solar azymuth angle          (nobs)
             SensorZenith   ---  sensor zenith angle          (nobs)
             SensorAzimuth  ---  sensor azymuth angle         (nobs)
             ScatteringAngle --  scattering angle             (nobs)

             aod            ---  Aerosol Optical Depth        (nobs,nch)
             aod_fine       ---  AOD of fine mode at 550nm    (nobs)
             reflectance    ---  reflectance                  (nobs,nch)

             cloud ---  cloud fraction               (nobs)
           
        """

#       Read in the data
#       ----------------
        if not os.path.exists(filename):
            raise IOError, "cannot find file <%s>"%filename
        self.nymd = -1
        self.nhms = -1
        ext = os.path.splitext(filename)[1] 
        if ext == '.ods':
            self.readODS (filename,nymd,nhms,kx,only_good)
        elif ext == '.npz':
            self._readNPZ (filename)
        else:
            raise IOError, "Unknown extension %s, it must be either .ods or .npz"%ext 

#---
    def readODS (self, filename, nymd, nhms, kx=None, only_good=True):
        """
        Reads an ODS file and creates a MODIS object.
        """
        
#       Load ODS file (unless we already have it)
#       -----------------------------------------
        try:
            if nymd==self.nymd and nhms==self.nhms:
               ods = self.ods
               kx = self.kx
            else:
                raise ValueError  # Different date, needs to read it
        except:
            ods = ODS(filename,nymd,nhms,only_good)

        if ods.nobs < 1:
            self.nobs = 0
            return # no obs, nothing to do
        
        if kx != None:
            ods = ods.select(kx=kx)  # restrict to kx, usually not needed

#       Consistency check
#       -----------------
        if any(ods.kx!=ods.kx[0]):
            raise ValueError, "more than one kx in ods; in this case must specify kx"
            
#       Keep this as ODS for now
#       ------------------------
        self.SolarZenith = ods.select(kt=ktSOLARZENITH)
        self.reflectance  = ods.select(kt=ktREFLECTANCE)
        self.aod  = ods.select(kt=ktAOD)
        self.aod_fine  = ods.select(kt=ktAOD_FINE)
        
#       Consistency checks: use reflectance and solar zenith as reference; 
#       these must be defined for all soundings
#       -----------------------------------------------------------
        nobs = self.SolarZenith.nobs
        nch  = self.reflectance.nobs / self.SolarZenith.nobs 
        
        if self.reflectance.nobs%self.SolarZenith.nobs != 0:
            raise ValueError, \
                  "reflectance.size=% is not a multiple of SolarZenith.size=%d"%\
             (self.reflectance.nobs,self.SolarZenith.nobs)
 
#       Pull out relevant arrays out of ODS object
#       ------------------------------------------
        self.nymd         = nymd
        self.nhms         = nhms   
        self.nobs         = nobs
        self.nch          = nch
        self.kx           = ods.kx[0]
        self.ks           = self.SolarZenith.ks
        self.time         = self.SolarZenith.time
        self.lon          = self.SolarZenith.lon
        self.lat          = self.SolarZenith.lat
        self.channels     = self.reflectance.lev.reshape(nobs,nch)[0,:]

        self.SolarAzimuth  = ods.select(kt=ktSOLARAZIMUTH).obs
        self.SolarZenith  = self.SolarZenith.obs
        self.SensorZenith  = ods.select(kt=ktSENSORZENITH).obs
        self.SensorAzimuth = ods.select(kt=ktSENSORAZIMUTH).obs
        self.ScatteringAngle = ods.select(kt=ktSCATTERINGANGLE).obs
        self.cloud = ods.select(kt=ktCLOUD).obs
        self.qa_flag = ods.select(kt=ktQA_FLAG).obs

#       Construct hash table with index of each ks
#       ------------------------------------------
        KS = {}
        for k in range(nobs):
            id_ = str(self.ks[k]) + str(self.time[k]) # Because of a bug, ks is not unique
            KS[id_] = k 
        if len(KS) != self.ks.size:
            raise ValueError, 'Strange, len(KS)=%d but len(ks)=%d'%(len(KS),self.ks.size)

#       Get obs, adding missing if sounding is missing
#       ----------------------------------------------
        self.aod = _getObsNch('aod',self.aod,KS,self.channels)
        self.aod_fine = _getObsNch('aod_fine',self.aod_fine,KS,self.channels)
        self.reflectance = _getObsNch('reflectance',\
                                      self.reflectance,KS,self.channels)

#---
    def _readNPZ (self, filename ):
        """
        Reads an NPZ file (written by method writeu) and creates an OMI object.        
        """
        npz = load(filename)

        self.nymd, self.nhms, self.nobs, self.nch, self.kx, version = npz['meta']

        self.lon = npz['lon']
        self.lat = npz['lat']
        self.ks = npz['ks']
        self.channels = npz['channels']
        self.qa_flag = npz['qa_flag']
        self.SolarZenith = npz['SolarZenith']
        self.SolarAzimuth = npz['SolarAzimuth']
        self.SensorZenith = npz['SensorZenith']
        self.SensorAzimuth = npz['SensorAzimuth']
        self.ScatteringAngle = npz['ScatteringAngle']
        self.cloud = npz['cloud']
        self.aod = npz['aod']
        self.aod_fine = npz['aod_fine']
        self.reflectance = npz['reflectance']

#---

    def write(self,filename=None,dir='.',expid='modis',Verb=1):
        """
        Writes the un-gridded OMI object to a numpy npz file. 
        """

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

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

    def writeODS(self,filename=None,dir='.',expid='modis',channels=None,
                 revised=False,Verb=1):
        """
        Writes the un-gridded OMI object to an ODS file. If *revised*
        is True, the revid *aod_* parameter is written to file.
        """

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        if filename is None:
           filename = '%s/%s.obs.%d_%02dz.ods'%(dir,expid,self.nymd,self.nhms/10000)

        if channels is None:
            channels = list(self.channels)

        # Create and populated ODS object
        # -------------------------------
        ns = self.nobs
        nobs = nobs=len(channels) * ns
        ods = ODS(nobs=nobs, kx=self.kx, kt=ktAOD)
        i = 0
        for ch in channels:
            I = range(i,i+ns)
            j = channels.index(ch)
            ods.ks[I]  = i+1
            ods.lat[I] = self.lat[:]
            ods.lon[I] = self.lon[:]
            ods.time[I] = self.time[:]
            ods.lev[I] = ch
            ods.qch[I] = self.qa_flag[:].astype('int')
            if revised:
                ods.obs[I]  = self.aod_[:,j]
                ods.xvec[I] = self.aod[:,j]
            else:
                ods.obs[I] = self.aod[:,j]
            i += ns

        # Exclusion flag
        # --------------
        iGood = (ods.qch>0) & (ods.obs<10.)
        ods.qcx[:] = 1     # All bad...
        ods.qcx[iGood] = 0 # ... unless good

        ods.write(filename,self.nymd,self.nhms)
        
        if Verb >=1:
            print "[w] Wrote file "+filename

#---
    def writeg(self,filename=None,dir='.',expid='modis',refine=8,res=None,
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

       # Stop here is no good obs available
       # ----------------------------------
       if self.nobs == 0:
           return # no data to work with
       if any(self.iGood) == False:
           return # no good data to work with

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

       title = 'Gridded MODIS Measurements'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.obs.%d_%02dz.nc4'%(dir,expid,self.nymd,self.nhms/10000)

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
       aod_fine = self.aod_fine[:,I]
       
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
    def addVar(self,ga,expr='mag(u10m,v10m)',vname='wind',clmYear=None):
        """
        Given a grads object *ga* having the correct MERRA file as default,
        interpolates *var* to obs location and saves it as an attribute
        named *vname*.
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

        # To conserve memory, restrict domain with 1 deg halo
        # ---------------------------------------------------
        x1, x2  = self.lon.min()-1.,self.lon.max()+1.
        y1, y2  = self.lat.min()-1.,self.lat.max()+1.
        ga('set lon %f %f'%(x1,x2),Quiet=True)
        ga('set lat %f %f'%(y1,y2),Quiet=True)
        expr_ = ga.expr(expr)
        u, levs = ga.interp(expr_, self.lon, self.lat )
        U = u.data
        if len(shape(U)) == 0:
             U = [U,] # so that we can slice it later

        self.__dict__[vname] = U

#--

def _gatime(nymd,nhms):
        Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
        cymd = "%08d"%int(nymd)
        chms = "%06d"%int(nhms)
        t = chms[0:2]+":"+chms[2:4]+"Z"+\
            cymd[6:8]+Months[int(cymd[4:6])-1]+cymd[0:4]
        return t
    
def _getObs(name,ods,KS):
    """
    Given an ODS structure, returns a 1D array with the "obs"
    attribute, with MISSING values used for the missing soundings.
    On input, *KS* is a hash with the index of each sounding index.
    """
    obs = MISSING * ones(len(KS))
    i = [ KS[str(s)+str(t)] for s,t in zip(ods.ks,ods.time) ] # index of each available sounding
    obs[i] = ods.obs
    return obs

def _getObsNch(name,ods,KS,channels):
    """
    Given an ODS structure, returns a 1D array with the "obs"
    attribute, with MISSING values used for the missing soundings.
    On input, *ns* is the number of sounds (number of pixels in the
    orbit), and *channels* are the channels.
    """
    ns  = len(KS)
    nch = channels.shape[0] 
    obs = MISSING * ones((ns,nch))
    j = 0
    for ch in channels:
        ods_ = ods.select(lev=ch)
        obs[:,j] = _getObs(name,ods_,KS)
        j = j + 1
    return obs
        
#...............................................................................

if __name__ == "__main__":

    """Simple unit testing."""

    
    filename = 'MOD04_L2_051.aero_tc8.ocean.20080630.ods'
    nymd = 20080630
    nhms = 120000
    
    modis = MODIS(filename,nymd,nhms)
    modis.writeg(channels=[870,660,550,470])
    modis.writeODS(channels=[870,660,550,470])

    
