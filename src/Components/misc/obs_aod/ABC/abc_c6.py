"""
   This module implements a Neural Net based MODIS Collection 6 Neural Net Retrieval.

   Arlindo da Silva, June 2015.

"""

import os, sys
sys.path.insert(0,'/home/pcastell/Enthought/Canopy_64bit/System/lib/python2.7/site-packages')
from   matplotlib.pyplot    import  cm, imshow, plot, figure
from   matplotlib.pyplot    import  xlabel, ylabel, title, grid, savefig, legend
import matplotlib.pyplot    as      plt
from   matplotlib.ticker    import  MultipleLocator
import matplotlib.patches   as      mpatches
from   numpy                import  c_ as cat
from   numpy                import  random, sort, pi, load, cos, log, std, exp
from   numpy                import  reshape, arange, ones, zeros, interp, sqrt
from   numpy                import  meshgrid, concatenate, squeeze
import numpy                as      np
from   giant                import  MISSING, LAND, OCEAN, DEEP
from   nn                   import  NN, _plotKDE
import itertools
from   sklearn.linear_model import LinearRegression
from   multiprocessing      import cpu_count
from   abc_c6_aux           import SummarizeCombinations, get_Iquartiles, get_Ispecies, get_ImRef
from   abc_c6_aux           import make_plots, make_error_pdfs, TestStats, SummaryPDFs
from   brdf                 import rtlsReflectance
from   mcd43c               import BRDF

# ------
MODVARNAMES = {'mRef470': 'MOD04 470 nm Reflectance',
               'mRef550': 'MOD04 550 nm Reflectance',
               'mRef660': 'MOD04 660 nm Reflectance',
               'mRef870': 'MOD04 870 nm Reflectance',
               'mRef1200': 'MOD04 1200 nm Reflectance',
               'mRef1600': 'MOD04 1600 nm Reflectance',
               'mRef2100': 'MOD04 2100 nm Reflectance'}

MYDVARNAMES = {'mRef470': 'MYD04 470 nm Reflectance',
               'mRef550': 'MYD04 550 nm Reflectance',
               'mRef660': 'MYD04 660 nm Reflectance',
               'mRef870': 'MYD04 870 nm Reflectance',
               'mRef1200': 'MYD04 1200 nm Reflectance',
               'mRef1600': 'MYD04 1600 nm Reflectance',
               'mRef2100': 'MYD04 2100 nm Reflectance'}

VARNAMES    = {'cloud': 'MOD04 Cloud Fraction',
               'ScatteringAngle': 'Scattering Angle',
               'GlintAngle': 'Glint Angle',
               'AMF': 'Air Mass Factor',
               'SolarZenith': 'Solar Zenith Angle',
               'CoxMunkLUT': 'Cox-Munk White Sky Albedo',
               'COxMunkBRF': 'Cox-Munk Bidirectional Surface Reflectance',
               'MOD43BClimAlbedo': 'MOD43B Albedo Climatology',
               'fdu': 'MERRA2 Fraction Dust Aerosol',
               'fcc': 'MERRA2 Fraction Carbonaceous Aerosol',
               'fsu': 'MERRA2 Fraction Sulfate Aerosol',
               'year': 'Year'}
#--------------------------------------------------------------------------------------
class SETUP(object):
  def setupNN(self,retrieval,expid,
                 nHidden=None,
                 nHLayers=1,
                 combinations=False,
                 Input_nnr  = ['mRef470','mRef550','mRef660', 'mRef870',
                                'mRef1200','mRef1600','mRef2100',
                                'ScatteringAngle', 'GlintAngle',
                                'AMF', 'SolarZenith',
                                'cloud', 'albedo','fdu','fcc','fsu' ],
                 Input_const = None,
                 Target = ['aTau550',],
                 K=None):
       
    
    self.retrieval = retrieval
    # Create outdir if it doesn't exist
    # ---------------------------------
    self.outdir = "./{}/".format(expid)
    if not os.path.exists(self.outdir):
      os.makedirs(self.outdir)

    self.plotdir = self.outdir

    # save some inputs
    # -----------------
    self.expid   = expid
    self.Target  = Target
    self.K       = K
    self.nHidden = nHidden
      
    # Balance the dataset before splitting
    # No aerosol type should make up more that 35% 
    # of the total number of obs
    # --------------------------------------
    self.iValid = self.balance(self.nobs*0.35)

    # Flatten Input_nnr into one list
    # -------------------------------
    input_nnr = flatten_list(Input_nnr)

    # Create list of combinations
    # ---------------------------
    if combinations:
      self.comblist, self.combgroups = get_combinations(Input_nnr,Input_const)
    else:
      self.comblist = [input_nnr]
          
    # Initialize arrays to hold stats
    # ------------------------------
    self.nnr  = STATS(K,self.comblist)
    self.orig = STATS(K,self.comblist)

    # Initialize K-folding
    # --------------------
    if K is None:
      self.iTest = ones([self.nobs]).astype(bool)
      self.iTrain = self.iValid
    else:
      self.kfold(K=K)

    # Create list of topologies
    # -------------------------  
    self.topology = []
    if not combinations:
      if self.nHidden is None:
        self.nHidden  = len(input_nnr)
      else:
        self.nHidden = nHidden

      self.topology.append((len(input_nnr),) + (self.nHidden,)*nHLayers + (len(Target),))

    else:
      for c,Input in enumerate(self.comblist):
        if nHidden is None:
          self.nHidden  = len(Input)
        else:
          self.nHidden = nHidden

        self.topology.append((len(Input),) + (self.nHidden,)*nHLayers + (len(Target),))

    self.combinations = combinations

#----------------------------------------------------------------------------    
class ABC(object):

    """
    Common Subroutines to all the ABC Classes
    """

    def __init__(self,fname,Albedo,coxmunk_lut=None):

        # Get Auxiliary Data
        self.setfnameRoot(fname)
        self.setWind()
        self.setAlbedo(Albedo,coxmunk_lut=coxmunk_lut)
        self.setSpec()
        self.setNDVI()

    def setfnameRoot(self,fname):
        if self.sat == 'Aqua':
            self.fnameRoot = 'myd_' + fname.split('/')[-1].split('.')[0]
        elif self.sat == 'Terra':
            self.fnameRoot = 'mod_' + fname.split('/')[-1].split('.')[0]
        

    def setWind(self):
        # Read in wind
        # ------------------------
        self.wind = load(self.fnameRoot + "_MERRA2.npz")['wind']
        self.giantList.append('wind')
        self.Wind = '' #need this for backwards compatibility

    def setAlbedo(self,Albedo,coxmunk_lut=None):
        # Define wind speed dependent ocean albedo
        # ----------------------------------------
        if Albedo is not None:
          for albedo in Albedo:
            if albedo == 'CoxMunkLUT':
              self.getCoxMunk(coxmunk_lut) 
              self.giantList.append(albedo)    
            elif albedo == 'MCD43C1':
              self.setBRDF()     
            elif 'CxAlbedo' in albedo :
              self.setCoxMunkBRF(albedo)
            else:
              self.__dict__[albedo] = squeeze(load(self.fnameRoot+'_'+albedo+'.npz')["albedo"])
              self.giantList.append(albedo)

    def setSpec(self):
        # Read in Aerosol Fractional Composition
        # --------------------------------------
        names = ('fdu','fss','fcc','fsu')
        for name in names:
            self.__dict__[name] = load(self.fnameRoot + "_MERRA2.npz")[name]
            self.giantList.append(name)

    def setCoxMunkBRF(self,albedo):
        # Read in Cox Munk Bidirectional surface reflectance
        # --------------------------------------------------
        self.__dict__[albedo] = squeeze(load(self.fnameRoot+'_CoxMunkBRF.npz')[albedo])
        self.giantList.append(albedo)

    def setBRDF(self):
        # Read in MCD43C1 BRDF
        # Calculate bidirectional surface reflectance
        # ---------------------------------------------
        brdf = BRDF(self.nobs)
        names = ('BRDFvis','BRDFnir','BRDF470','BRDF550',
                 'BRDF650','BRDF850','BRDF1200','BRDF1600',
                 'BRDF2100')
        for name in brdf.__dict__:
            brdf.__dict__[name] = load(self.fnameRoot + "_MCD43C1.npz")[name]

        for name in names:
            ch = name[4:]
            Kiso = brdf.__dict__['Riso'+ch]
            Kvol = brdf.__dict__['Rvol'+ch]
            Kgeo = brdf.__dict__['Rgeo'+ch]
            self.__dict__[name] = rtlsReflectance(Kiso,Kgeo,Kvol,
                                                  self.SolarZenith,self.SensorZenith,
                                                  self.SolarAzimuth,self.SensorAzimuth)
            self.giantList.append(name)    


    def setNDVI(self):
        # Read in NDVI
        # -------------
        names = ('NDVI','EVI','NIRref')
        for name in names:
          self.__dict__[name] = load(self.fnameRoot + "_NDVI.npz")[name]
          self.giantList.append(name)

    def outlierRemoval(self,outliers):

        # # Outlier removal based on log-transformed AOD
        # # --------------------------------------------
        if outliers > 0.:
            d = log(self.mTau550[self.iValid]+0.01) - log(self.aTau550[self.iValid]+0.01)
            if self.verbose>0:
                print "Outlier removal: %d   sig_d = %f  nGood=%d "%(-1,std(d),d.size)
            for iter in range(3):
                iValid = (abs(d)<outliers*std(d))
                self.iValid[self.iValid] = iValid
                d = log(self.mTau550[self.iValid]+0.01) - log(self.aTau550[self.iValid]+0.01)
                if self.verbose>0:
                    print "Outlier removal: %d   sig_d = %f  nGood=%d "%(iter,std(d),d.size)
              
    def angleTranform(self):            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      
        self.AMF             = (1/self.SolarZenith) + (1/self.SensorZenith)
        self.giantList.append('AMF')

    def setYear(self):
        # Year
        #-------
        self.year = np.array([t.year for t in self.tyme])
        self.giantList.append('year')

    def addFilter(self,aFilter):
        if aFilter is not None:
          filters = []
          for f in aFilter:
            filters.append(self.__dict__[f]>0)
            
          oiValid = reduce(lambda x,y: x&y,filters)
          self.iValid = self.iValid & oiValid

#---------------------------------------------------------------------------- 
class ABC_Ocean (OCEAN,NN,SETUP,ABC):

    def __init__ (self,fname, 
                  coxmunk_lut='/nobackup/NNR/Misc/coxmunk_lut.npz',
                  outliers=3., 
                  laod=True, 
                  verbose=0,
                  cloud_thresh=0.70,
                  glint_thresh=40.0,
                  Albedo=None,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Ocean algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)

        Reads in two Albedo variables
              albedo - cox munk lut that parameterizes albedo with wind speed
                                        Requires coxmunk_lut npz file.
              BRF - cox munk bidirection reflectance computed with VLIDORT
                                        and stored in npz
                                        Requires a NPZ file with the data.
              Both require a wind speed npz file.
        """

        self.verbose = verbose
        self.laod    = laod

        OCEAN.__init__(self,fname) # initialize superclass

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo,coxmunk_lut=coxmunk_lut)

        # Q/C
        # ---
        self.iValid = (self.qa>0) & \
                      (self.aTau470 > -0.01) &\
                      (self.aTau550 > -0.01) &\
                      (self.aTau660 > -0.01) &\
                      (self.aTau870 > -0.01) &\
                      (self.mTau470 > -0.01) &\
                      (self.mTau550 > -0.01) &\
                      (self.mTau660 > -0.01) &\
                      (self.mTau870 > -0.01) &\
                      (self.mRef470 > 0.0)   &\
                      (self.mRef550 > 0.0)   &\
                      (self.mRef660 > 0.0)   &\
                      (self.mRef870 > 0.0)   &\
                      (self.mRef1200 > 0.0)  &\
                      (self.mRef1600 > 0.0)  &\
                      (self.mRef2100 > 0.0)  &\
                      (self.cloud <cloud_thresh) &\
                      (self.cloud > 0)           &\
                      (self.GlintAngle != MISSING ) &\
                      (self.GlintAngle > glint_thresh) 

        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)

        # glint_thresh > 40 is a bit redundant b/c MOD04 should already give these a qa==0 or
        # does not retrieve.  However, there are a few cases (~200) where this does not happen.
        # the GlingAngle is very close to 40, greater than 38.  Not sure why these get through.

        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)
              
        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)
            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()

        # Year
        #-------
        self.setYear()

#----------------------------------------------------------------------------    

class ABC_Land (LAND,NN,SETUP,ABC):

    def __init__ (self, fname,
                  Albedo=None,
                  alb_max = 0.25,
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname)  # initialize superclass

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo)

        # Q/C: enforce QA=3 and albedo in (0,0.25), scattering angle<170
        # --------------------------------------------------------------
        self.iValid = (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.mTau2100> -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.cloud > 0)            &\
                      (self.ScatteringAngle<170.) & \
                      (self.mRef412 > 0)          & \
                      (self.mRef440 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef550 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mRef870 > 0)          & \
                      (self.mRef1200 > 0)         & \
                      (self.mRef1600 > 0)         & \
                      (self.mRef2100 > 0)         & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660 >  0.0)       & \
                      (self.mSre2100>  0.0)       

        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)

        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()    
#----------------------------------------------------------------------------    

class ABC_Deep (DEEP,NN,SETUP,ABC):

    def __init__ (self, fname,
                  useLAND=False,
                  Albedo=None,
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.verbose = verbose
        self.laod = laod

        DEEP.__init__(self,fname)  # initialize superclass

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo)

        # Q/C: enforce QA=3 and albedo in (0,0.25), scattering angle<170
        # --------------------------------------------------------------
        self.iValid = (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau412 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.cloud > 0)            &\
                      (self.ScatteringAngle<170.) & \
                      (self.mRef412 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mSre412 >  0.0)       & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660 >  0.0)      

        LANDref     = (self.mRef440 > 0)          & \
                      (self.mRef550 > 0)          & \
                      (self.mRef870 > 0)          & \
                      (self.mRef1200 > 0)         & \
                      (self.mRef1600 > 0)         & \
                      (self.mRef2100 > 0)         

        if useLAND:
          self.iValid = self.iValid & LANDref

        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()    

#----------------------------------------------------------------------------    

class ABC_DBDT (LAND,NN,SETUP,ABC):

    def __init__ (self, fname,
                  useDEEP = False,
                  Albedo=None,
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname)  # initialize superclass
        dbl = DEEP(fname)

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo)

        # Q/C: enforce QA=3 and scattering angle<170
        # Combines deep blue and dark target
        # --------------------------------------------------------------
        LANDiValid =  (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.cloud > 0)            & \
                      (self.ScatteringAngle<170.) & \
                      (self.mRef412 > 0)          & \
                      (self.mRef440 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef550 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mRef870 > 0)          & \
                      (self.mRef1200 > 0)         & \
                      (self.mRef1600 > 0)         & \
                      (self.mRef2100 > 0)         & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660>  0.0)        &\
                      (self.mSre2100 > 0.0)   

        DEEPiValid =  (dbl.qa==3)                & \
                      (dbl.aTau470 > -0.01)      & \
                      (dbl.aTau550 > -0.01)      & \
                      (dbl.aTau660 > -0.01)      & \
                      (dbl.mTau412 > -0.01)      & \
                      (dbl.mTau470 > -0.01)      & \
                      (dbl.mTau550 > -0.01)      & \
                      (dbl.mTau660 > -0.01)      & \
                      (dbl.cloud<cloud_thresh)   & \
                      (dbl.cloud > 0)           & \
                      (dbl.ScatteringAngle<170.) & \
                      (dbl.mRef412 > 0)          & \
                      (dbl.mRef440 > 0)          & \
                      (dbl.mRef470 > 0)          & \
                      (dbl.mRef550 > 0)          & \
                      (dbl.mRef660 > 0)          & \
                      (dbl.mRef870 > 0)          & \
                      (dbl.mRef1200 > 0)         & \
                      (dbl.mRef1600 > 0)         & \
                      (dbl.mRef2100 > 0)         & \
                      (dbl.mSre412 >  0.0)       & \
                      (dbl.mSre470 >  0.0)       & \
                      (dbl.mSre660>  0.0)       


        LANDref    =  (self.mRef412 > 0)          & \
                      (self.mRef440 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef550 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mRef870 > 0)          & \
                      (self.mRef1200 > 0)         & \
                      (self.mRef1600 > 0)         & \
                      (self.mRef2100 > 0)         & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660>  0.0)        &\
                      (self.mSre2100 > 0.0)  
       
        addiValid = (LANDiValid == False) & (DEEPiValid) & (LANDref)

        if useDEEP:
          replace = ['mRef412','mRef470','mRef660','mSre470','mSre660','mTau550']          
        else:
          replace = ['mTau550']
        for name in replace:
          self.__dict__[name][addiValid] = dbl.__dict__[name][addiValid]

        self.iValid = LANDiValid
        self.iValid[addiValid] = True

        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()    

#----------------------------------------------------------------------------    

class ABC_DBDT_INT (LAND,NN,SETUP,ABC):

    def __init__ (self, fname,
                  useDEEP = False,
                  testDEEP= False,
                  Albedo=['MOD43BClimAlbedo'],
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname)  # initialize superclass
        dbl = DEEP(fname)

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo)        

        # Q/C: enforce QA=3 and scattering angle<170
        # Combines deep blue and dark target
        # --------------------------------------------------------------
        LANDiValid =  (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.cloud > 0)            & \
                      (self.ScatteringAngle<170.) & \
                      (self.mRef412 > 0)          & \
                      (self.mRef440 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef550 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mRef870 > 0)          & \
                      (self.mRef1200 > 0)         & \
                      (self.mRef1600 > 0)         & \
                      (self.mRef2100 > 0)         & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660>  0.0)        &\
                      (self.mSre2100 > 0.0)  

        DEEPiValid =  (dbl.qa==3)                & \
                      (dbl.aTau470 > -0.01)      & \
                      (dbl.aTau550 > -0.01)      & \
                      (dbl.aTau660 > -0.01)      & \
                      (dbl.mTau412 > -0.01)      & \
                      (dbl.mTau470 > -0.01)      & \
                      (dbl.mTau550 > -0.01)      & \
                      (dbl.mTau660 > -0.01)      & \
                      (dbl.cloud<cloud_thresh)   & \
                      (dbl.cloud > 0)            & \
                      (dbl.ScatteringAngle<170.) & \
                      (dbl.mRef412 > 0)          & \
                      (dbl.mRef470 > 0)          & \
                      (dbl.mRef660 > 0)          & \
                      (dbl.mSre412 >  0.0)       & \
                      (dbl.mSre470 >  0.0)       & \
                      (dbl.mSre660>  0.0)       

        # Where both LAND and DEEP decide to retrieve
        intiValid = LANDiValid & DEEPiValid

        if useDEEP:
          replace = ['mRef412','mRef470','mRef660','mSre470','mSre660']
          add     = ['mSre412']
          for name in replace:
            self.__dict__[name][intiValid] = dbl.__dict__[name][intiValid]
          for name in add:
            self.__dict__[name] = dbl.__dict__[name]
            self.giantList.append(name)  

        if testDEEP:
          replace = ['mTau550']
          for name in replace:
            self.__dict__[name][intiValid] = dbl.__dict__[name][intiValid]          

        self.iValid = intiValid

        self.dbmTau550 = dbl.mTau550
        self.giantList.append('dbmTau550')  

        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()    

#----------------------------------------------------------------------------    

class ABC_LAND_COMP (LAND,NN,SETUP,ABC):

    def __init__ (self, fname,
                  useDEEP = False,
                  Albedo=['MOD43BClimAlbedo'],
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname)  # initialize superclass
        dbl = DEEP(fname)

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo)                

        # Q/C: enforce QA=3 and scattering angle<170
        # Combines deep blue and dark target
        # --------------------------------------------------------------
        LANDiValid =  (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.cloud > 0)            & \
                      (self.ScatteringAngle<170.) & \
                      (self.mRef412 > 0)          & \
                      (self.mRef440 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef550 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mRef870 > 0)          & \
                      (self.mRef1200 > 0)         & \
                      (self.mRef1600 > 0)         & \
                      (self.mRef2100 > 0)         & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660>  0.0)        &\
                      (self.mSre2100 > 0.0)  

        DEEPiValid =  (dbl.qa==3)                & \
                      (dbl.aTau470 > -0.01)      & \
                      (dbl.aTau550 > -0.01)      & \
                      (dbl.aTau660 > -0.01)      & \
                      (dbl.mTau412 > -0.01)      & \
                      (dbl.mTau470 > -0.01)      & \
                      (dbl.mTau550 > -0.01)      & \
                      (dbl.mTau660 > -0.01)      & \
                      (dbl.cloud<cloud_thresh)   & \
                      (dbl.cloud > 0)           & \
                      (dbl.ScatteringAngle<170.) & \
                      (dbl.mRef412 > 0)          & \
                      (dbl.mRef470 > 0)          & \
                      (dbl.mRef660 > 0)          & \
                      (dbl.mSre412 >  0.0)       & \
                      (dbl.mSre470 >  0.0)       & \
                      (dbl.mSre660>  0.0)       

        DEEPref    =  (dbl.mRef412 > 0)          & \
                      (dbl.mRef470 > 0)          & \
                      (dbl.mRef660 > 0)          & \
                      (dbl.mSre412 >  0.0)       & \
                      (dbl.mSre470 >  0.0)       & \
                      (dbl.mSre660>  0.0)       

        # Where LAND retrieves, and DEEP does not
        # The LAND Complement
        compiValid = LANDiValid & ~DEEPiValid

        if useDEEP:
          compiValid = compiValid & DEEPref
          replace = ['mRef412','mRef470','mRef660','mSre470','mSre660']
          add     = ['mSre412']
          for name in replace:
            self.__dict__[name][compiValid] = dbl.__dict__[name][compiValid]
          for name in add:
            self.__dict__[name] = dbl.__dict__[name]
            self.giantList.append(name)  

        self.iValid = compiValid


        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()    
#----------------------------------------------------------------------------    


class ABC_DEEP_COMP (DEEP,NN,SETUP,ABC):

    def __init__ (self, fname,
                  useLAND = False,
                  DBDT=False,
                  noLANDref=False,
                  Albedo=['MOD43BClimAlbedo'],
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  aFilter=None):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.verbose = verbose
        self.laod = laod

        DEEP.__init__(self,fname)  # initialize superclass
        lnd = LAND(fname)

        # Get Auxiliary Data
        ABC.__init__(self,fname,Albedo)           


        # Q/C: enforce QA=3 and scattering angle<170
        # Combines deep blue and dark target
        # --------------------------------------------------------------
        LANDiValid =  (lnd.qa==3)                & \
                      (lnd.aTau470 > -0.01)      & \
                      (lnd.aTau550 > -0.01)      & \
                      (lnd.aTau660 > -0.01)      & \
                      (lnd.mTau470 > -0.01)      & \
                      (lnd.mTau550 > -0.01)      & \
                      (lnd.mTau660 > -0.01)      & \
                      (lnd.cloud<cloud_thresh)   & \
                      (lnd.cloud > 0)            & \
                      (lnd.ScatteringAngle<170.) & \
                      (lnd.mRef412 > 0)          & \
                      (lnd.mRef440 > 0)          & \
                      (lnd.mRef470 > 0)          & \
                      (lnd.mRef550 > 0)          & \
                      (lnd.mRef660 > 0)          & \
                      (lnd.mRef870 > 0)          & \
                      (lnd.mRef1200 > 0)         & \
                      (lnd.mRef1600 > 0)         & \
                      (lnd.mRef2100 > 0)         & \
                      (lnd.mSre470 >  0.0)       & \
                      (lnd.mSre660>  0.0)        &\
                      (lnd.mSre2100 > 0.0)  

        DEEPiValid =  (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau412 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.cloud > 0)           & \
                      (self.ScatteringAngle<170.) & \
                      (self.mRef412 > 0)          & \
                      (self.mRef470 > 0)          & \
                      (self.mRef660 > 0)          & \
                      (self.mSre412 >  0.0)       & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660>  0.0)       

        LANDref    =  (lnd.mRef412 > 0)          & \
                      (lnd.mRef440 > 0)          & \
                      (lnd.mRef470 > 0)          & \
                      (lnd.mRef550 > 0)          & \
                      (lnd.mRef660 > 0)          & \
                      (lnd.mRef870 > 0)          & \
                      (lnd.mRef1200 > 0)         & \
                      (lnd.mRef1600 > 0)         & \
                      (lnd.mRef2100 > 0)         & \
                      (lnd.mSre470 >  0.0)       & \
                      (lnd.mSre660>  0.0)        &\
                      (lnd.mSre2100 > 0.0)  

        # Where DEEP retrieves, and LAND does not
        # The DEEP Complement
        compiValid = ~LANDiValid & DEEPiValid

        if useLAND:
          compiValid = compiValid & LANDref
          replace = ['mRef412','mRef440','mRef470','mRef550','mRef660','mRef870',
                     'mRef1200','mRef1600','mRef2100','mSre470','mSre660']
          add     = ['mSre2100']

          for name in add:
            self.__dict__[name] = lnd.__dict__[name]
            self.giantList.append(name)  
          if not DBDT:
            for name in replace:
              self.__dict__[name][compiValid] = lnd.__dict__[name][compiValid]
        elif noLANDref:
          compiValid = compiValid & ~LANDref


        self.iValid = compiValid

        # Filter by additional variables
        # ------------------------------
        self.addFilter(aFilter)
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        self.outlierRemoval(outliers)

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.angleTranform()    
#--------------------------------------------------------------------------

class STATS(object):

  def __init__ (self,K,comblist):
    c = max([len(comblist),1])
    if K is None:
      k = 1
    else:
      k = K

    self.slope     = np.ones([k,c])*-999.
    self.intercept = np.ones([k,c])*-999.
    self.R         = np.ones([k,c])*-999.
    self.rmse      = np.ones([k,c])*-999.
    self.mae       = np.ones([k,c])*-999.
    self.me        = np.ones([k,c])*-999.

#---------------------------------------------------------------------
def _train(mxd,expid,c):

  ident  = mxd.ident
  outdir = mxd.outdir

  nHidden  = mxd.nHidden
  topology = mxd.topology[c]

  Input    = mxd.comblist[c]
  Target   = mxd.Target
  
  print "-"*80
  print "--> nHidden = ", nHidden
  print "-->  Inputs = ", Input
  
  n = cpu_count()
  kwargs = {'nproc' : n}
  if mxd.K is None:
    mxd.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology,**kwargs)
    mxd.savenet(outdir+"/"+expid+'_Tau.net')    
  else:
    k = 1
    for iTrain, iTest in mxd.kf:
      I = arange(mxd.nobs)
      iValid = I[mxd.iValid]
      mxd.iTrain = iValid[iTrain]
      mxd.iTest  = iValid[iTest]

      mxd.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology,**kwargs)
      mxd.savenet(outdir+"/"+expid+'.k={}_Tau.net'.format(str(k)))
      k = k + 1
#--------------------------------------------------------------------------------------

def _test(mxd,expid,c,plotting=True):
  ident  = mxd.ident
  outdir = mxd.outdir    
  if mxd.K is None:
    if mxd.combinations:
      inputs = expid.split('.')
      found = False
      for invars in itertools.permutations(inputs):
        try: 
          netFile = outdir+"/"+".".join(invars)+'_Tau.net'
          mxd.loadnet(netFile)
          found = True
          break
        except:
          pass
        if found: break

      if not found:
        print '{} not found.  Need to train this combinatin of inputs'.format(netFile)
        raise
    else:
      netFile = outdir+"/"+expid+'_Tau.net'

    mxd.net = mxd.loadnet(netFile)
    mxd.Input = mxd.comblist[0]
    TestStats(mxd,mxd.K,c)
    if plotting: make_plots(mxd,expid,ident,I=mxd.iTest)
  else:
    k = 1
    for iTrain, iTest in mxd.kf:
      I = arange(mxd.nobs)
      iValid = I[mxd.iValid]
      mxd.iTrain = iValid[iTrain]
      mxd.iTest  = iValid[iTest]

      if mxd.combinations:
        inputs = expid.split('.')
        found = False
        for invars in itertools.permutations(inputs):
          try: 
            netFile = outdir+"/"+".".join(invars)+'.k={}_Tau.net'.format(str(k))
            mxd.loadnet(netFile)
            found = True
            print 'found file',netFile
            break
          except:
            pass

        if not found:
          print '{} not found.  Need to train this combinatin of inputs'.format(netFile)
          raise
      else:
        netFile = outdir+"/"+expid+'.k={}_Tau.net'.format(str(k))


      mxd.net = mxd.loadnet(netFile)
      mxd.Input = mxd.comblist[c]      
      TestStats(mxd,k-1,c)
      if plotting: make_plots(mxd,expid,'.k={}'.format(str(k)),I=mxd.iTest)
      k = k + 1    

#---------------------------------------------------------------------
def _trainMODIS(mxdx):

  if not mxdx.combinations:
    _train(mxdx,mxdx.expid,0)
  else:
    for c,Input in enumerate(mxdx.comblist):
      _train(mxdx,'.'.join(Input),c)

# -------------------------------------------------------------------
def _testMODIS(mxdx):

  if not mxdx.combinations:
    _test(mxdx,mxdx.expid,0,plotting=True)
  else:
    for c,Input in enumerate(mxdx.comblist):
      _test(mxdx,'.'.join(Input),c,plotting=False)


#---------------------------------------------------------------------
def get_combinations(Input_nnr,Input_const):
  comblist   = []
  combgroups = []
  for n in arange(len(Input_nnr)):
    for invars in itertools.combinations(Input_nnr,n+1):        
      b = ()
      for c in invars:
        if type(c) is list:
          b = b + tuple(c)
        else:
          b = b + (c,)
      #don't do both kinds of abledo together
      if not (('BRF' in b) and ('albedo' in b)):
        if Input_const is not None:
          comblist.append(tuple(Input_const)+b)
          combgroups.append((Input_const,)+invars)
        else:
          comblist.append(b)
          combgroups.append(invars)

  if Input_const is not None:
    comblist.insert(0,tuple(Input_const))
    combgroups.insert(0,(Input_const,))

  return comblist,combgroups
#---------------------------------------------------------------------
def flatten_list(Input_nnr):
  Input = ()
  for i in Input_nnr:
    if type(i) is list:
      Input = Input + tuple(i)
    else:
      Input = Input + (i,)

  return list(Input)
        
#------------------------------------------------------------------
  
if __name__ == "__main__":

  """
    Example Training/testing
  """
  from   abc_c6_aux           import SummarizeCombinations, SummaryPDFs

  sat          = 'Aqua'
  nHidden      = None
  nHLayers     = 1
  combinations = True
  Target       = ['aTau550',]
  Albedo       = ['CoxMunkBRF']
  K            = 2

  if sat is 'Terra':
    filename     = '/nobackup/6/NNR/Training/giant_C6_10km_Terra_20150921.nc'
    retrieval    = 'MOD_OCEAN'
  if sat is 'Aqua':
    filename     = '/nobackup/6/NNR/Training/giant_C6_10km_Aqua_20151005.nc'
    retrieval    = 'MYD_OCEAN'


  doTrain      = True
  doTest       = True
  expid        = 'CoxMunkTest_wSZA'
  Input_const  = ['SolarZenith','ScatteringAngle', 'GlintAngle','mRef470','mRef550','mRef660', 'mRef870','mRef1200','mRef1600','mRef2100']
  Input_nnr    = ['CoxMunkBRF',['fdu','fcc','fsu']]
  aFilter      = ['CoxMunkBRF']
 
  expid        = '{}_{}'.format(retrieval,expid)

  ocean = ABC_Ocean(filename,Albedo=Albedo,verbose=1,aFilter=aFilter)  
  # Initialize class for training/testing
  # ---------------------------------------------
  ocean.setupNN(retrieval, expid,
                      nHidden      = nHidden,
                      nHLayers     = nHLayers,
                      combinations = combinations,
                      Input_const  = Input_const,
                      Input_nnr    = Input_nnr,                                         
                      Target       = Target,                      
                      K            = K)

  if Input_const is not None:
    InputMaster = list((Input_const,) + tuple(Input_nnr))      
  else:
    InputMaster = Input_nnr

  # Do Training and Testing
  # ------------------------
  if doTrain:
    _trainMODIS(ocean)

  if doTest:
    _testMODIS(ocean)

    if combinations:
      SummarizeCombinations(ocean,InputMaster,yrange=None,sortname='slope')
      SummaryPDFs(ocean,varnames=['mRef870','mRef660'])  