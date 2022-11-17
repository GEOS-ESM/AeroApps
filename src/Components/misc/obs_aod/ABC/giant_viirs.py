#!/bin/env python

import os
import sys

from types     import *
from netCDF4   import Dataset
from datetime  import date, datetime, timedelta
from numpy     import savez, meshgrid, array, concatenate, zeros, ones, \
                      linspace, sqrt, load, shape, random, interp
import numpy   as np
from dateutil.parser import parse as isoparse

from pyobs.npz       import NPZ

META =  ( "Date",
          "Time",
          "Latitude",
          "Longitude",
          "SolarZenith",
          "SensorZenith",
          "ScatteringAngle",
          "GlintAngle")

ANET = (  "mean_AOD0340",
         "mean_AOD0380",
         "mean_AOD0410",
         "mean_AOD0440",
         "mean_AOD0470intrp",
         "mean_AOD0500",
         "mean_AOD0550intrp",
         "mean_AOD0660intrp",
         "mean_AOD0670",
         "mean_AOD0870",
         "mean_AOD1020",
         "mean_AOD1640",
#         "mean_AOD2100intrp",
         "mean_WaterVapor",         
         "nval_AOD0550intrp" 
         )


xDT_LAND = ("mean_AOD0480corr-l",
        "mean_AOD0550corr-l",
        "mean_AOD0670corr-l",
        "mean_AOD2250corr-l",
#        "mean_mref0412-l",
#        "mean_mref0443-l",
        "mean_mref0480-l",
#        "mean_mref0550-l",
        "mean_mref0670-l",
#        "mean_mref0745-l",
#        "mean_mref0870-l",
#        "mean_mref1200-l",
#        "mean_mref1600-l",
        "mean_mref2250-l",
        "mean_surfre0480-l",
        "mean_surfre0670-l",
        "mean_surfre2250-l",
        "mean_acfrac-l",   
        "mode_QA-l", 
#        "cval_cldpixdistavg",
#        "nval_npu0550-l",
#        "nval_AOD0550corr-l" 
        )

xDT_OCEAN = ( "mean_AOD0480ea-o",
          "mean_AOD0550ea-o",
          "mean_AOD0670ea-o",
          "mean_AOD0860ea-o",
          "mean_AOD1240ea-o",
          "mean_AOD1600ea-o",
          "mean_AOD2250ea-o",
          "mean_mref0480-o",
          "mean_mref0550-o",
          "mean_mref0670-o",
          "mean_mref0860-o",
          "mean_mref1240-o",
          "mean_mref1600-o",
          "mean_mref2250-o",
#          "mean_wspeed-o",
          "mean_acfrac-o",
          "mode_QAavg-o",
#          "cval_cldpixdistavg",
          "nval_AOD0550ea-o"
         )

xDB_LAND = ( "mean_AOD0412dpbl-l",
          "mean_AOD0488dpbl-l",
#          "mean_AOD0550bestdpbl-l",
          "mean_AOD0550dpbl-l",
          "mean_AOD0670dpbl-l",
          "mean_cfracdpbl-l",
          "mean_mref0412dpbl-l",
          "mean_mref0488dpbl-l",
          "mean_mref0550dpbl-l",
          "mean_mref0670dpbl-l",
          "mean_mref0865dpbl-l",
          "mean_mref1240dpbl-l",
          "mean_mref1640dpbl-l",
          "mean_mref2250dpbl-l",
          "mean_surfre0412dpbl-l",
          "mean_surfre0488dpbl-l",
          "mean_surfre0670dpbl-l",
          "mode_QAdpbl-l",
#          "cval_cldpixdistavg",
          "nval_AOD0550dpbl-l",
        )

xDB_OCEAN = ("mean_AOD0488dpbl-o",
#          "mean_AOD0550bestdpbl-l",
          "mean_AOD0550dpbl-o",
          "mean_AOD0670dpbl-o",
          "mean_AOD0865dpbl-o",
          "mean_AOD1240dpbl-o",
          "mean_AOD1640dpbl-o",
          "mean_AOD2250dpbl-o",
          "mean_cfracdpbl-o",
          "mean_mref0412dpbl-o",
          "mean_mref0488dpbl-o",
          "mean_mref0550dpbl-o",
          "mean_mref0670dpbl-o",
          "mean_mref0865dpbl-o",
          "mean_mref1240dpbl-o",
          "mean_mref1640dpbl-o",
          "mean_mref2250dpbl-o",
          "mode_QAdpbl-o",
#          "cval_cldpixdistavg",
        )

ALIAS = {
                "Latitude"              : 'lat',
                "Longitude"             : 'lon',
                "mean_AOD0340"          : 'aTau340',
                "mean_AOD0380"          : 'aTau380',
                "mean_AOD0410"          : 'aTau412',
                "mean_AOD0440"          : 'aTau440',
                "mean_AOD0470intrp"     : 'aTau470',
                "mean_AOD0500"          : 'aTau500',
                "mean_AOD0550intrp"     : 'aTau550',
                "mean_AOD0660intrp"     : 'aTau660',
                "mean_AOD0670"          : 'aTau670',
                "mean_AOD0870"          : 'aTau870',
                "mean_AOD1020"          : 'aTau1020',
                "mean_AOD1640"          : 'aTau1640',
#                "mean_AOD2100intrp"     : 'aTau2100',
                "mean_WaterVapor"       : 'aWaterVapor',

                "mean_AOD0480corr-l"    : 'mTau480',
                "mean_AOD0550corr-l"    : 'mTau550',
                "mean_AOD0670corr-l"    : 'mTau670',
                "mean_AOD2250corr-l"    : 'mTau2250',
#                "mean_mref0412-l"       : 'mRef412',
#                "mean_mref0443-l"       : 'mRef440',
                "mean_mref0480-l"       : 'mRef480',
#                "mean_mref0550-l"       : 'mRef550',
                "mean_mref0670-l"       : 'mRef670',
#                "mean_mref0745-l"       : 'mRef745',
#                "mean_mref0870-l"       : 'mRef870',
#                "mean_mref1200-l"       : 'mRef1200',
#                "mean_mref1600-l"       : 'mRef1600',
                "mean_mref2250-l"       : 'mRef2250',
                "mean_surfre0480-l"     : 'mSre480',
                "mean_surfre0670-l"     : 'mSre670',
                "mean_surfre2250-l"     : 'mSre2250',
                "mean_acfrac-l"         : 'cloud',    
                "mode_QA-l"             : 'qa',                

                "mean_AOD0480ea-o"      : 'mTau480',
                "mean_AOD0550ea-o"      : 'mTau550',
                "mean_AOD0670ea-o"      : 'mTau670',
                "mean_AOD0860ea-o"      : 'mTau860',
                "mean_AOD1240ea-o"      : 'mTau1240',
                "mean_AOD1600ea-o"      : 'mTau1600',
                "mean_AOD2250ea-o"      : 'mTau2250',
                "mean_mref0480-o"       : 'mRef480',
                "mean_mref0550-o"       : 'mRef550',
                "mean_mref0670-o"       : 'mRef670',
                "mean_mref0860-o"       : 'mRef860',
                "mean_mref1240-o"       : 'mRef1240',
                "mean_mref1600-o"       : 'mRef1600',
                "mean_mref2250-o"       : 'mRef2250',
                "mean_acfrac-o"         : 'cloud',
                "mode_QAavg-o"          : 'qa',

                "mean_AOD0412dpbl-l"    : 'mTau412',
                "mean_AOD0488dpbl-l"    : 'mTau488',
                "mean_AOD0550dpbl-l"    : 'mTau550',
                "mean_AOD0670dpbl-l"    : 'mTau670',
                "mean_AOD2250dpbl-l"    : 'mTau2250',
                "mean_cfracdpbl-l"      : 'cloud',    
                "mean_mref0412dpbl-l"   : 'mRef412',
                "mean_mref0488dpbl-l"   : 'mRef488',
                "mean_mref0550dpbl-l"   : 'mRef550',
                "mean_mref0670dpbl-l"   : 'mRef670',
                "mean_mref0865dpbl-l"   : 'mRef865',
                "mean_mref1240dpbl-l"   : 'mRef1240',
                "mean_mref1640dpbl-l"   : 'mRef1640',
                "mean_mref2250dpbl-l"   : 'mRef2250',
                "mean_surfre0412dpbl-l" : 'mSre412',
                "mean_surfre0488dpbl-l" : 'mSre488',
                "mean_surfre0670dpbl-l" : 'mSre670',
                "mode_QAdpbl-l"         : 'qa',

                "mean_AOD0488dpbl-o"    : 'mTau488',
                "mean_AOD0550dpbl-o"    : 'mTau550',
                "mean_AOD0670dpbl-o"    : 'mTau670',
                "mean_AOD0865dpbl-o"    : 'mTau865',
                "mean_AOD1240dpbl-o"    : 'mTau1240',
                "mean_AOD1640dpbl-o"    : 'mTau1640',
                "mean_AOD2250dpbl-o"    : 'mTau2250',
                "mean_cfracdpbl-o"      : 'cloud',
                "mean_mref0412dpbl-o"   : 'mRef412',
                "mean_mref0488dpbl-o"   : 'mRef488',
                "mean_mref0550dpbl-o"   : 'mRef550',
                "mean_mref0670dpbl-o"   : 'mRef670',
                "mean_mref0865dpbl-o"   : 'mRef865',
                "mean_mref1240dpbl-o"   : 'mRef1240',
                "mean_mref1640dpbl-o"   : 'mRef1640',
                "mean_mref2250dpbl-o"   : 'mRef2250',
                "mode_QAdpbl-o"         : 'qa',
             }

MISSING = 1.E20

# ..................................................................................
class CX_ALBEDO(object):
  """
    Container for Cox-Munk albedo 
  """
  def __init__(self,s_channels):
    for ch in s_channels:
      self.__dict__['CxAlbedo' + ch] = []


class GIANT(object):

  """
      Read Level-2 GIANT Co-located AERONET/MODIS/VIIRS aerosol files from
      GSFC MODIS group.
  """


  def __init__ (self,filename,xVars=(),only_good=True,tymemax=None):
    """
     Creates an GIANT object defining the attributes corresponding
     to the SDS of interest.
    """

    if 'SNPP' in filename:     self.sat = 'SNPP'
    else:                      self.sat = 'Unknown'    

    self.only_good = only_good
    Names = META+ANET+xVars

    # Simplify variable names
    # ------------------------
    self.ALIAS = ALIAS.copy()

    # Read in variables
    # -----------------
    print 'filename ',filename
    nc = Dataset(filename)
    Alias = self.ALIAS.keys()
    self.giantList =[]
    for name in Names:
      data = nc.variables[name][:]
      if name in Alias:
          name = self.ALIAS[name]

      # old files use -9999.0 for fill value
      # new files use masked arrays
      # convert everythong to regular array filling with -9999.0
      # make sure _fill_value is -9999.0
      self.__dict__[name] = np.array(data)
      self.giantList.append(name)
    nc.close()

    # Form python tyme
    # ----------------
    # new files have an ISO_DateTime variable
    nc = Dataset(filename)
    if 'ISO_DateTime' in nc.variables.keys():
        try:
            iso = nc.variables['ISO_DateTime'][:]
            self.tyme = array([isoparse(''.join(array(t))) for t in iso])
        except:
        # old file only have Date and Time variables
            D = self.Date[:,0:10]
            T = self.Time[:,0:5] # they didn't save the seconds
            # Bug in dataset, first field is blank
            D[0] = D[1]
            T[0] = T[1]
            self.aTau550[0] = -9999.0
            self.tyme = array([ isoparse(''.join(D[i])+'T'+''.join(t)) for i, t in enumerate(T) ])    
    nc.close()

    # Limit to the MERRA-2 time series
    #---------------------------------
    if tymemax is not None:
      tymemax = isoparse(tymemax)
      I = self.tyme < tymemax
      
      for name in Names:
        if name in Alias:
          name = self.ALIAS[name]
        self.__dict__[name] = self.__dict__[name][I]

      self.tyme = self.tyme[I]

    del self.Date, self.Time
    self.giantList.remove('Date')
    self.giantList.remove('Time')
    self.giantList.append('tyme')

    # Record number of observations
    # -----------------------------
    self.nobs = len(self.lon)

#--
  def balance(self,N):
    """
    Return indices of observations so that each species does not have more than
    N observations. This is meant to be performed with a reduced dataset.
    """
    I = zeros(self.lon.shape).astype(bool)
    random.seed(32768) # so that we get the same permutation

    for f in (self.fdu,self.fss,self.fcc,self.fsu):

      J = f>0.5                      # all obs for which species dominate
      n = len(self.lon[J])              # no. obs for this species
      P = random.permutation(n)      # randomize obs for this species
      m = min(n,N)                   # keep this many

      K = I[J]
      K[P[0:m]] = True
      I[J] = K

    return I

#---
  def reduce(self,I):
    """
    Reduce observations according to index I. 
    """
    for name in self.giantList:
      q = self.__dict__[name]
      print "{} Reducing "+name,q.shape
      self.__dict__[name] = q[I]

    self.nobs = len(self.lon)


  def getCoxMunk(self,filename='/nobackup/NNR/Misc/coxmunk_lut.npz',channel=550.):
    """
    Returns ocean albedo as a function of wind speed from look up table.
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

    self.CoxMunkLUT = albedo


  def calcCoxMunk(self,windFile,channels=[470. ,550. ,660. ,870. ,1200.,1600.,2100.],npzFile=None):
    """
    Calls VLIDORT wrapper to calculate CoxMunk Bidirectional Surface Reflectance.
    """
    import VLIDORT_BRDF_ABC_

    channeli = [200,250,300,337,400,488,515,550,633,694,860,1060,1300,1536,1800,2000,2250,2500]
    # refractive index
    mr       = [1.396,1.362,1.349,1.345,1.339,1.335,1.334,1.333,1.332,1.331,1.329,1.326,
               1.323,1.318,1.312,1.306,1.292,1.261]

    wind = np.load(windFile)
    u10m = wind['u10m'].astype('float64')
    v10m = wind['v10m'].astype('float64')
    sza = self.SolarZenith.astype('float64')
    vza = self.SensorZenith.astype('float64')

    # needs to be photon travel direction.  raa = 0 is forward scattering
    saa = self.SolarAzimuth.astype('float64')
    saa = saa + 180.0
    saa[saa>=360.0] = saa[saa>=360.0] - 360.0
    raa = np.abs(self.SensorAzimuth - saa).astype('float64')

    if type(channels) is list:
      channels = np.array(channels)

    try:
      some_object_iterator = iter(channels)
    except:
      channels = np.array([channels])

    strch = [str(int(ch)) for ch in channels]
    m  = interp(channels,channeli,mr)
    
    albedos, rc = VLIDORT_BRDF_ABC_.coxmunk(1,channels,u10m,v10m,m,sza,raa,vza,-999,1)

    self.sample = CX_ALBEDO(strch)
    for i,ch in enumerate(strch):
      self.sample.__dict__['CxAlbedo' + ch] = albedos[:,i]


    if npzFile is not None:
      savez(npzFile,**self.sample.__dict__)

#---
  def speciate(self,aer_x,FineMode=False,Verbose=False):
    """
    Use GAAS to derive fractional composition.
    """

    self.sampleFile(aer_x,onlyVars=('TOTEXTTAU',
                                    'DUEXTTAU',
                                    'SSEXTTAU',
                                    'BCEXTTAU',
                                    'OCEXTTAU',
                                    'SUEXTTAU',
                                    ),Verbose=Verbose)
    s = self.sample
    I = (s.TOTEXTTAU<=0)
    s.TOTEXTTAU[I] = 1.E30
    self.fdu  = s.DUEXTTAU / s.TOTEXTTAU
    self.fss  = s.SSEXTTAU / s.TOTEXTTAU
    self.fbc  = s.BCEXTTAU / s.TOTEXTTAU
    self.foc  = s.OCEXTTAU / s.TOTEXTTAU
    self.fcc  = self.fbc + self.foc
    self.fsu  = s.SUEXTTAU / s.TOTEXTTAU

    if FineMode:
      TOTEXTTAU = s.TOTEXTTAU[:]
      self.sampleFile(aer_x,onlyVars=('DUEXTTFM','SSEXTTFM'),Verbose=Verbose)
      self.fduf = s.DUEXTTFM / TOTEXTTAU
      self.fssf = s.SSEXTTFM / TOTEXTTAU 

    del self.sample

#---
  def sampleG5(self,gas_x=None,avk_x=None,int_x=None,slv_x=None,ext_Nc=None):
    """
    Sample key parameters from GAAS files.
    """
    from gfio import GFIOHandle

    if gas_x is not None:
      self.sampleFile(gas_x,onlyVars=('AODANA',))
      self.tau_550 = self.sample.AODANA[:]

    if avk_x is not None:
      tyme = self.tyme[:]
      self.tyme = getSyn(tyme)
      self.sampleFile(avk_x,onlyVars=('AOD',))
      self.avk = self.sample.AOD[:]
      self.tyme[:] = tyme[:]

    if int_x is not None:
      try:
        self.sampleFile(int_x,onlyVars=('TQV',)) # As in file spec
        self.tpw = self.sample.TQV[:]
      except:
        self.sampleFile(int_x,onlyVars=('TPW',)) # Larry's name
        self.tpw = self.sample.TPW[:]

    if slv_x is not None:
      self.sampleFile(slv_x,onlyVars=('U10M','V10M'))
      self.wind = sqrt(self.sample.U10M[:]**2 + self.sample.V10M[:]**2)

    if ext_Nc is not None:
      self.sampleFile(ext_Nc,onlyVars=('taod',))
      self.tau_660 = self.sample.taod[:,5] # 660


    del self.sample


  def sampleMERRA(self,slv_x='tavg1_2d_slv_Nx',aer_x='tavg1_2d_aer_Nx',
                  FineMode=False,npzFile=None,Verbose=False):
    self.sampleFile(slv_x,onlyVars=('U10M','V10M'), Verbose=Verbose)
    self.u10m = self.sample.U10M
    self.v10m = self.sample.V10M
    self.wind = sqrt(self.sample.U10M[:]**2 + self.sample.V10M[:]**2)

    del self.sample

    self.speciate(aer_x,FineMode=FineMode,Verbose=Verbose)

    if npzFile is not None:
      if FineMode:
        savez(npzFile,wind=self.wind,u10m=self.u10m,v10m=self.v10m,
                    fdu=self.fdu,fss=self.fss,fcc=self.fcc,fsu=self.fsu,
                    fduf=self.fduf,fssf=self.fssf)     
      else:
        savez(npzFile,wind=self.wind,u10m=self.u10m,v10m=self.v10m,
                    fdu=self.fdu,fss=self.fss,fcc=self.fcc,fsu=self.fsu)     

  def sampleMCD43C(self,npzFile=None,Verbose=False):
    from mcd43c import MCD43C

    brdf = MCD43C()
    brdf.sample(self,Verbose=Verbose)

    if npzFile is not None:
        savez(npzFile,**self.brdf.__dict__)     

#---
  def sampleFile(self, inFile, npzFile=None, onlyVars=None, Verbose=False, clmYear=None):
        """
        Interpolates all variables of inFile and optionally
        save them to file *npzFile*
        """
        from gfio import GFIO, GFIOctl, GFIOHandle

        # Instantiate grads and open file
        # -------------------------------
        name, ext = os.path.splitext(inFile)
        if ext in ( '.nc4', '.nc', '.hdf'):
          fh = GFIO(inFile)     # open single file
          if fh.lm == 1:
            timeInterp = False    # no time interpolation in this case
          else:
            raise ValueError, "cannot handle files with more tha 1 time, use ctl instead"
        else:
          fh = GFIOctl(inFile)  # open timeseries
          timeInterp = True     # perform time interpolation

        self.sample = GFIOHandle(inFile)
        if onlyVars is None:
            onlyVars = fh.vname

        nt = self.lon.shape
        lons = self.lon
        lats = self.lat

        if clmYear is None:
          tymes = self.tyme
        else:
          tymes = array([t + timedelta(days=365*(clmYear-t.year)) for t in self.tyme])

        print 'trange',tymes.min(),tymes.max()
        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print "<> Sampling ", v
            if timeInterp:
              var = fh.sample(v,lons,lats,tymes,Verbose=Verbose)
            else:
              var = fh.interp(v,lons,lats)
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var
            else:
                raise IndexError, 'variable <%s> has rank = %d'%(v,len(var.shape))

        if npzFile is not None:
            savez(npzFile,**self.sample.__dict__)            
       
  def sampleLoadz(self,npzFile):
        """
        Loads sample from npz file.
        """
        from grads.gahandle import GaHandle
        self.sample = GaHandle(npzFile)
        npz = load(npzFile)
        for v in npz.keys():
            self.sample.__dict__[v] = npz[v]
                
#---

class DT_LAND(GIANT):
    def __init__(self,filename,tymemax=None): 
        GIANT.__init__(self,filename,xVars=xDT_LAND,tymemax=tymemax)
        if self.sat == 'SNPP':
            self.ident = 'vsdtl'
        self.ident = self.ident + '_'+ filename.split('/')[-1].split('.')[0]
        self.surface = 'land'

        # Angstrom interpolate retrieved AOD to
        # nominal MODIS wavelengths (470, 550, 660)
        # for ease of comparison
        # VIIRS channels = 480., 550., 670., 2250.
        # ------------------------------------------------        
        self.mTau470 = aodInterpAngs(470.,self.mTau480,self.mTau550,480.,550.)
        self.mTau660 = aodInterpAngs(660.,self.mTau550,self.mTau670,550.,670.)

        self.giantList += ['mTau470','mTau660']

class DT_OCEAN(GIANT):
    def __init__(self,filename,tymemax=None): 
        GIANT.__init__(self,filename,xVars=xDT_OCEAN,tymemax=tymemax)
        if self.sat == 'SNPP':
            self.ident = 'vsdto'
        self.ident = self.ident + '_' + filename.split('/')[-1].split('.')[0]
        self.surface = 'ocean'

        # Angstrom interpolate retrieved AOD to
        # nominal MODIS wavelengths (470, 550, 660, 870, 1200, 1600, 2100)
        # for ease of comparison
        # VIIRS channels = 480., 550., 670., 860., 1240., 1600., 2250.
        # ------------------------------------------------
        self.mTau470  = aodInterpAngs(470.,self.mTau480,self.mTau550,480.,550.)
        self.mTau660  = aodInterpAngs(660.,self.mTau550,self.mTau670,550.,670.)        
        self.mTau870  = aodInterpAngs(870.,self.mTau670,self.mTau860,670.,860.)
        self.mTau2100 = aodInterpAngs(2100.,self.mTau1600,self.mTau2250,1600.,2250.)

        self.giantList +=['mTau470','mTau660','mTau870','mTau2100']

class DB_OCEAN(GIANT):
    def __init__(self,filename,tymemax=None):
        GIANT.__init__(self,filename,xVars=xDB_OCEAN,tymemax=tymemax)
        if self.sat == 'SNPP':
            self.ident = 'vsdto'
        self.ident = self.ident + '_' + filename.split('/')[-1].split('.')[0]
        self.surface = 'ocean'

        # Angstrom interpolate retrieved AOD to
        # nominal MODIS wavelengths (470, 550, 660, 870, 1200, 1600, 2100)
        # for ease of comparison
        # VIIRS channels = 488.,  550.,  670.,  865., 1240., 1640., 2250.
        # ------------------------------------------------        
        self.mTau470  = aodInterpAngs(470.,self.mTau488,self.mTau550,488.,550.)
        self.mTau660  = aodInterpAngs(660.,self.mTau550,self.mTau670,550.,670.)
        self.mTau870  = aodInterpAngs(870.,self.mTau670,self.mTau865,670.,865.)
        self.mTau2100 = aodInterpAngs(2100.,self.mTau1640,self.mTau2250,1640.,2250.)

        self.giantList += ['mTau470','mTau660','mTau870','mTau2100']

class DB_LAND(GIANT):
    def __init__(self,filename,tymemax=None): 
        GIANT.__init__(self,filename,xVars=xDB_LAND,tymemax=tymemax)
        if self.sat == 'SNPP':
            self.ident = 'vsdto'
        self.ident = self.ident + '_' + filename.split('/')[-1].split('.')[0]
        self.surface = 'land'

        # Angstrom interpolate retrieved AOD to 
        # nominal MODIS wavelengths (412, 470, 550, 660)
        # for ease of comparison
        # VIIRS channels = 412, 488, 550, 670
        # ------------------------------------------------
        self.mTau470 = aodInterpAngs(470.,self.mTau412,self.mTau488,412.,488.)
        self.mTau660 = aodInterpAngs(660.,self.mTau550,self.mTau670,550.,670.)

        self.giantList += ['mTau470','mTau660']

#---
def aodInterpAngs(lambda_,tau1,tau2,lambda1,lambda2):
    """
       Angstrom-interpolated AOD.
    """
    I = (tau1 > 0.0) & (tau2 > 0.0)

    angstrom = -np.log(tau1[I]/tau2[I])/np.log(lambda1/lambda2)
    tau = -9999. * np.ones(len(tau1))
    tau[I] = tau2[I] * (lambda2/lambda_)**angstrom

    I = (tau1 == 0.0) & (tau2 == 0.0)
    tau[I] = 0.0
    return tau
# .......................................................................................................

if __name__ == "__main__":

    # lnd = LAND('/nobackup/6/NNR/Training/giant_C6_10km_Terra_05Oct2015.nc')
    ocn = DT_OCEAN('/nobackup/6/NNR/Training/giant_C6_10km_Terra_20150921.nc')
    #xblu = DEEP('/Users/adasilva/workspace.local/Data_Analysis/C6/giant_C6_10km_9April2015.nc')

