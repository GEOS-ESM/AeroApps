#!/bin/env python

import os
import sys

from types     import *
from netCDF4   import Dataset
from datetime  import date, datetime, timedelta
from numpy     import savez, meshgrid, array, concatenate, zeros, ones, \
                      linspace, sqrt, load, shape, random, interp
from dateutil.parser import parse as isoparse

from pyobs.npz       import NPZ

META =  ( "Date",
          "Time",
          "Latitude",
          "Longitude",
          "SolarZenith",
          "SolarAzimuth",
          "SensorZenith",
          "SensorAzimuth",
          "ScatteringAngle",
          "GlintAngle")

ANET = ( # "mean_AOD0340",
#         "mean_AOD0380",
#         "mean_AOD0410",
#         "mean_AOD0440",
#         "mean_AOD0470intrp",
#         "mean_AOD0500",
#         "mean_AOD0550",
#         "mean_AOD0550img",
         "mean_AOD0550intrp",
         "mean_AOD0660intrp",
#         "mean_AOD0670",
         "mean_AOD0870",
#         "mean_AOD1020",
#         "mean_AOD1640",
         "mean_AOD2100intrp",
         )


xLAND = ("mean_AOD0470corr-l",
        "mean_AOD0550corr-l",
        "mean_AOD0660corr-l",
        "mean_AOD2100corr-l",
        "mean_QAfilteredAOD0470corr-l",
        "mean_QAfilteredAOD0550corr-l",
        "mean_QAfilteredAOD0660corr-l",
        "mean_QAfilteredAOD2100corr-l",
        "mean_mref0412-l",
        "mean_mref0443-l",
        "mean_mref0470-l",
        "mean_mref0550-l",
        "mean_mref0660-l",
        "mean_mref0745-l",
        "mean_mref0870-l",
        "mean_mref1200-l",
        "mean_mref1600-l",
        "mean_mref2100-l",
        "mean_surfre0470-l",
        "mean_surfre0660-l",
        "mean_surfre2100-l",
        "mean_acfrac-l",
        )

xOCEAN = ( "mean_AOD0470ea-o",
          "mean_AOD0550ea-o",
          "mean_AOD0550m01-o",
          "mean_AOD0550m02-o",
          "mean_AOD0550m03-o",
          "mean_AOD0550m04-o",
          "mean_AOD0550m05-o",
          "mean_AOD0550m06-o",
          "mean_AOD0550m07-o",
          "mean_AOD0550m08-o",
          "mean_AOD0550m09-o",
          "mean_AOD0660ea-o",
          "mean_AOD0870ea-o",
          "mean_AOD1200ea-o",
          "mean_AOD1600ea-o",
          "mean_AOD2100ea-o",
          "mean_QAfilteredAOD0470ea-o",
          "mean_QAfilteredAOD0550ea-o",
          "mean_QAfilteredAOD0660ea-o",
          "mean_QAfilteredAOD0870ea-o",
          "mean_QAfilteredAOD1200ea-o",
          "mean_QAfilteredAOD1600ea-o",
          "mean_QAfilteredAOD2100ea-o",
          "mean_mref0470-o",
          "mean_mref0550-o",
          "mean_mref0660-o",
          "mean_mref0870-o",
          "mean_mref1200-o",
          "mean_mref1600-o",
          "mean_mref2100-o",
          "mean_wspeed-o",
          "mean_acfrac-o",
         )

xDEEP = ( "mean_AOD0412dpbl-l",
         "mean_AOD0470dpbl-l",
#         "mean_AOD0550bestdpbl-l",
         "mean_AOD0550dpbl-l",
         "mean_AOD0660dpbl-l",
         "mean_QAfilteredAOD0412dpbl-l",
         "mean_QAfilteredAOD0470dpbl-l",
         "mean_QAfilteredAOD0550dpbl-l",
         "mean_QAfilteredAOD0660dpbl-l",
         "mean_cfracdpbl-l",
         "mean_mref0412dpbl-l",
         "mean_mref0470dpbl-l",
         "mean_mref0660dpbl-l",
         "mean_surfre0412dpbl-l",
         "mean_surfre0470dpbl-l",
         "mean_surfre0660dpbl-l",
#         "mean_uncertAOD0550dpbl-l",
         "mean_cfracdpbl-l",
        )


ALIAS = dict (
                 Latitude  = 'lat',
                Longitude  = 'lon',
                mean_AOD0470intrp = 'aTau470',
                mean_AOD0550intrp = 'aTau550',
                mean_AOD0660intrp = 'aTau660',
                mean_AOD0870      = 'aTau870',
                mean_AOD2100intrp = 'aTau2100',
             )

MISSING = 1.E20

# ..................................................................................

class GIANT(object):

  """
      Read Level-2 GIANT Co-located AERONET/MODIS/VIIRS aerosol files from
      GSFC MODIS group.
  """


  def __init__ (self,filename,xVars=(),only_good=True):
     """
       Creates an GIANT object defining the attributes corresponding
       to the SDS of interest.
     """

     self.only_good = only_good
     Names = META+ANET+xVars

     # Simpligfy variable names
     # ------------------------
     self.ALIAS = ALIAS.copy()
     for name in xVars:
         simple = name.replace('mean_AOD','mTau')\
                      .replace('mean_mref','mRef')\
                      .replace('mean_acfrac','cloud')\
                      .replace('mean_cfrac','cloud')\
                      .replace('corr','')\
                      .replace('mean_QAfilteredAOD','mQA')\
                      .replace('mean_surfre','mSref')\
                      .replace('dpbl','')\
                      .replace('-l','')\
                      .replace('-o','')
         self.ALIAS[name] = simple

     # Read in variables
     # -----------------
     nc = Dataset(filename)
     Alias = self.ALIAS.keys()
     for name in Names:
         v = nc.variables[name]
         rank = len(v.shape)
         if rank == 1:
            data = v[:]
         elif rank==2:
            data = v[:,:]
         else:
            raise ValueError, 'variable %s has invalid rank %d'%(name,rank)
         if name in Alias:
             name = self.ALIAS[name]
         self.__dict__[name] = data
     nc.close()

     # Form python tyme
     # ----------------
     D = self.Date.data[:][0:10]
     T = self.Time.data[:][0:10]
     self.tyme = array([ isoparse(''.join(d)+'T'+''.join(t)) for d, t in zip(D,T) ])
     del self.Date, self.Time

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

  def getCoxMunk(self,filename='/nobackup/NNR/Misc/c6.coxmunk_lut.npz',channel=550.):
    """
    Returns ocean albedo as a function of wid speed.
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

#---
  def reduce(self,I):
    """
    Reduce observations according to index I. 
    """
    Nicknames = ALIAS.values()
    for name in self.__dict__:
      if name in Nicknames:
        continue # alias do not get reduced
      q = self.__dict__[name]
      if type(q) is type(self.lon):
        if len(q) == self.nobs:
          # print "{} Reducing "+name
          self.__dict__[name] = q[I]

    Alias = ALIAS.keys()
    for sds in self.Names:
      if sds in Alias:
        self.__dict__[ALIAS[sds]] = self.__dict__[sds] # redefine aliases

    self.nobs = len(self.lon)

#---
  def speciate(self,aer_x,FineMode=False):
    """
    Use GAAS to derive fractional composition.
    """
    from gfio import GFIOHandle

    self.sampleFile(aer_x,onlyVars=('TOTEXTTAU',
                                    'DUEXTTAU',
                                    'SSEXTTAU',
                                    'BCEXTTAU',
                                    'OCEXTTAU',
                                    'SUEXTTAU',
                                    ))
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
      self.sampleFile(aer_x,onlyVars=('DUEXTTFM','SSEXTTFM'))
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

#---
  def sampleFile(self, inFile, npzFile=None, onlyVars=None, Verbose=False):
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
        tymes = self.tyme
        lons = self.lon
        lats = self.lat

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

class LAND(GIANT):
    def __init__(self,filename):
        GIANT.__init__(self,filename,xVars=xLAND)

class OCEAN(GIANT):
    def __init__(self,filename):
        GIANT.__init__(self,filename,xVars=xOCEAN)

class DEEP(GIANT):
    def __init__(self,filename):
        GIANT.__init__(self,filename,xVars=xDEEP)

# .......................................................................................................

if __name__ == "__main__":

    #lnd = LAND('/Users/adasilva/workspace.local/Data_Analysis/C6/giant_C6_10km_9April2015.nc')
    #ocn = OCEAN('/Users/adasilva/workspace.local/Data_Analysis/C6/giant_C6_10km_9April2015.nc')
    blu = DEEP('/Users/adasilva/workspace.local/Data_Analysis/C6/giant_C6_10km_9April2015.nc')
