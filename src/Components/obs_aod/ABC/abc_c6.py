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
#----------------------------------------------------------------------------    

class ABC_Ocean (OCEAN,NN):

    def __init__ (self,fname, 
                  coxmunk_lut='/nobackup/NNR/Misc/coxmunk_lut.npz',
                  outliers=3., 
                  laod=True, 
                  verbose=0,
                  cloud_thresh=0.70,
                  glint_thresh=40.0,
                  Albedo=['CoxMunkLUT'],
                  Input = ['mTAU550','mTAU470','mTAU660','mTAU870',
                           'ScatteringAngle','GlintAngle',
                           'SolarAzimuth','SolarZenith',
                           'SensorAzimuth','SensorZenith',
                           'cloud', 'wind'],
                  Target = [ 'aTau470', 'aTau550', 'aTau660', 'aTau870' ]):
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

        self.Input   = Input
        self.Target  = Target
        self.verbose = verbose
        self.laod    = laod

        OCEAN.__init__(self,fname) # initialize superclass
        if self.sat == 'Aqua':
            fnameRoot = 'myd_' + fname.split('/')[-1].split('.')[0]
        elif self.sat == 'Terra':
            fnameRoot = 'mod_' + fname.split('/')[-1].split('.')[0]

        # Read in wind
        # ------------------------
        self.wind = load(fnameRoot + "_MERRA2.npz")['wind']
        self.giantList.append('wind')
        self.Wind = '' #need this for backwards compatibility

        # Define wind speed dependent ocean albedo
        # ----------------------------------------
        for albedo in Albedo:
          if albedo == 'CoxMunkLUT':
            self.getCoxMunk(coxmunk_lut)          
          else:
            self.__dict__[albedo] = squeeze(load(fnameRoot+'_'+albedo+'.npz')["albedo"])
          self.giantList.append(albedo)
        # Read in Aerosol Fractional Composition
        # --------------------------------------
        names = ('fdu','fss','fcc','fsu')
        for name in names:
            self.__dict__[name] = load(fnameRoot + "_MERRA2.npz")[name]
            self.giantList.append(name)


        # Q/C
        # ---
        self.iValid = (self.qa>0) & \
                      (self.aTau470 > -0.01) & (self.aTau550 > -0.01) &  \
                      (self.aTau660 > -0.01) & (self.aTau870 > -0.01) &  \
                      (self.mTau470 > -0.01) & (self.mTau550 > -0.01) &  \
                      (self.mTau660 > -0.01) & (self.mTau870 > -0.01) &  \
                      (self.mRef470 > 0.0) & (self.mRef550 > 0.0) &  \
                      (self.mRef660 > 0.0) & (self.mRef870 > 0.0) &  \
                      (self.mRef1200 > 0.0) & (self.mRef1600 > 0.0) &  \
                      (self.mRef2100 > 0.0) &  \
                      (self.cloud <cloud_thresh) & (self.cloud > 0) &\
                      (self.GlintAngle != MISSING ) & (self.GlintAngle > glint_thresh) #&\
                      #(self.mNcollo >= 5)
                      #(self.aTau550 > 0.01)
                      #(self.wind>=0.0) 

        # glint_thresh > 40 is a bit redundant b/c MOD04 should already give these a qa==0 or
        # does not retrieve.  However, there are a few cases (~200) where this does not happen.
        # the GlingAngle is very close to 40, greater than 38.  Not sure why these get through.


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
              

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)
            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      
        self.AMF             = (1/self.SolarZenith) + (1/self.SensorZenith)

        # Year
        #-------
        self.year = np.array([t.year for t in self.tyme])

#----------------------------------------------------------------------------    

class ABC_Land (LAND,NN):

    def __init__ (self, fname,
                  Albedo=['MOD43BClimAlbedo'],
                  alb_max = 0.25,
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  Input = ['mTau550','mTau470','mTau660',
                           'mTAU550','mTAU470','mTAU660',
                           'ScatteringAngle',
                           'SolarAzimuth','SolarZenith',
                           'SensorAzimuth','SensorZenith',
                           'cloud' ],
                  Target = [ 'aTau470', 'aTau550', 'aTau660' ]):
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

        self.Input = Input
        self.Target = Target
        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname)  # initialize superclass

        if self.sat == 'Aqua':
            fnameRoot = 'myd_' + fname.split('/')[-1].split('.')[0]
        elif self.sat == 'Terra':
            fnameRoot = 'mod_' + fname.split('/')[-1].split('.')[0]


        # Read in desired albedo
        # ------------------------
        for albedo in Albedo:
          self.__dict__[albedo] = load(fnameRoot + "_" + albedo + ".npz")['albedo']
          self.giantList.append(albedo)

        # Read in Aerosol Fractional Composition
        # --------------------------------------
        names = ('fdu','fss','fcc','fsu')
        for name in names:
            self.__dict__[name] = load(fnameRoot + "_MERRA2.npz")[name]  
            self.giantList.append(name)      

        # Q/C: enforce QA=3 and albedo in (0,0.25), scattering angle<170
        # --------------------------------------------------------------
        self.iValid =  (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.mTau2100> -0.01)      & \
                      (self.cloud<cloud_thresh)   & \
                      (self.ScatteringAngle<170.) & \
                      (self.__dict__[Albedo[0]]>0)             & \
                      (self.__dict__[Albedo[0]]<alb_max) 

                      # (self.mSre470 >  0.0)       & \
                      # (self.mSre660 >  0.0)       & \
                      # (self.mSre2100>  0.0)       & \


        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
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

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      
        self.AMF             = (1/self.SolarZenith) + (1/self.SensorZenith)
#----------------------------------------------------------------------------    

class ABC_Deep (DEEP,NN):

    def __init__ (self, fname,
                  Albedo=['MOD43BClimAlbedo'],
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  Input = ['mTau550','mTau470','mTau660',
                           'mTAU550','mTAU470','mTAU660',
                           'ScatteringAngle',
                           'SolarAzimuth','SolarZenith',
                           'SensorAzimuth','SensorZenith',
                           'cloud' ],
                  Target = [ 'aTau470', 'aTau550', 'aTau660' ]):
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

        self.Input = Input
        self.Target = Target
        self.verbose = verbose
        self.laod = laod

        DEEP.__init__(self,fname)  # initialize superclass

        if self.sat == 'Aqua':
            fnameRoot = 'myd_' + fname.split('/')[-1].split('.')[0]
        elif self.sat == 'Terra':
            fnameRoot = 'mod_' + fname.split('/')[-1].split('.')[0]


        # Read in desired albedo
        # ------------------------
        for albedo in Albedo:
          self.__dict__[albedo] = load(fnameRoot + "_" + albedo + ".npz")['albedo']
          self.giantList.append(albedo)

        # Read in Aerosol Fractional Composition
        # --------------------------------------
        names = ('fdu','fss','fcc','fsu')
        for name in names:
            self.__dict__[name] = load(fnameRoot + "_MERRA2.npz")[name]  
            self.giantList.append(name)      

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
                      (self.ScatteringAngle<170.) & \
                      (self.__dict__[Albedo[0]]>0.00)  #&\
                      #(self.mRef660 >= 0.05) &  (self.mRef660 < 0.10)        

                      # (self.mSre470 >  0.0)       & \
                      # (self.mSre660 >  0.0)       & \
                      # (self.mSre2100>  0.0)       & \


        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
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

        # Reduce the Dataset
        # --------------------
        self.reduce(self.iValid)                    
        self.iValid = ones(self.lon.shape).astype(bool)        

            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      
        self.AMF             = (1/self.SolarZenith) + (1/self.SensorZenith)


#---------------------------------------------------------------------
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
def train_test(mxd,expid,Input,Target,K,plotting=True,c=None):

  ident  = mxd.ident
  outdir = mxd.outdir

  nHidden  = mxd.nHidden
  topology = mxd.topology
  
  print "-"*80
  print "--> nHidden = ", nHidden
  print "-->  Inputs = ", Input
  
  n = cpu_count()
  kwargs = {'nproc' : n}
  if K is None:
    # Split into training and testing sets
    # ------------------------------------
    mxd.split()

    mxd.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology,**kwargs)
    mxd.savenet(outdir+"/"+expid+"."+ident+'_Tau.net')
    TestStats(mxd,K,c)
    if plotting: make_plots(mxd,expid,ident)
    
  else:
    k = 1
    for iTrain, iTest in mxd.kf:
      I = arange(mxd.nobs)
      iValid = I[mxd.iValid]
      mxd.iTrain = iValid[iTrain]
      mxd.iTest  = iValid[iTest]

      mxd.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology,**kwargs)
      mxd.savenet(outdir+"/"+expid+"."+ident+'.k={}_Tau.net'.format(str(k)))
      TestStats(mxd,k-1,c)
      if plotting: make_plots(mxd,expid,ident+'.k={}'.format(str(k)))
      k = k + 1
  return mxd

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
    # Split into training and testing sets
    # ------------------------------------
    mxd.split()

    mxd.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology,**kwargs)
    mxd.savenet(outdir+"/"+expid+"."+ident+'_Tau.net')    
  else:
    k = 1
    for iTrain, iTest in mxd.kf:
      I = arange(mxd.nobs)
      iValid = I[mxd.iValid]
      mxd.iTrain = iValid[iTrain]
      mxd.iTest  = iValid[iTest]

      mxd.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology,**kwargs)
      mxd.savenet(outdir+"/"+expid+"."+ident+'.k={}_Tau.net'.format(str(k)))
      k = k + 1
#--------------------------------------------------------------------------------------

def _test(mxd,expid,c,plotting=True):
  ident  = mxd.ident
  outdir = mxd.outdir    
  if mxd.K is None:
    if mxdx.combinations:
      inputs = expid.split('.')
      found = False
      for invars in itertools.permutations(inputs):
        try: 
          netFile = outdir+"/"+".".join(invars)+"."+ident+'_Tau.net'
          mxdx.loadnet(netFile)
          found = True
          break
        except:
          pass
        if found: break

      if not found:
        print '{} not found.  Need to train this combinatin of inputs'.format(netFile)
        raise
    else:
      netFile = outdir+"/"+expid+"."+ident+'_Tau.net'

    mxd.net = mxd.loadnet(netFile)
    mxd.Input = mxd.comblist[0]
    TestStats(mxd,mxd.K,c)
    if plotting: make_plots(mxd,expid,ident)
  else:
    k = 1
    for iTrain, iTest in mxd.kf:
      I = arange(mxd.nobs)
      iValid = I[mxd.iValid]
      mxd.iTrain = iValid[iTrain]
      mxd.iTest  = iValid[iTest]

      if mxdx.combinations:
        inputs = expid.split('.')
        found = False
        for invars in itertools.permutations(inputs):
          try: 
            netFile = outdir+"/"+".".join(invars)+"."+ident+'.k={}_Tau.net'.format(str(k))
            mxdx.loadnet(netFile)
            found = True
            print 'found file',netFile
            break
          except:
            pass

        if not found:
          print '{} not found.  Need to train this combinatin of inputs'.format(netFile)
          raise
      else:
        netFile = outdir+"/"+expid+"."+ident+'.k={}_Tau.net'.format(str(k))

      mxd.net = mxd.loadnet(netFile)
      mxd.Input = mxd.comblist[c]      
      TestStats(mxd,k-1,c)
      if plotting: make_plots(mxd,expid,ident+'.k={}'.format(str(k)))
      k = k + 1    

#--------------------------------------------------------------------------------------

def _testMODIS(filename,retrieval,expid,
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
               Albedo=['CoxMunkLUT'],
               K=None):


  #------------------------------------------------------------
  # Read in data
  # --------------------------------------------------
  if retrieval.upper() == 'OCEAN':
    mxdx = ABC_Ocean(filename,Albedo=Albedo,verbose=1,laod=True)
  elif retrieval.upper() == 'LAND':
    mxdx = ABC_Land(filename,Albedo=Albedo,alb_max=0.25,verbose=1)
  elif retrieval.upper() == 'DEEP':
    mxdx = ABC_Deep(filename,Albedo=Albedo,verbose=1)

  mxdx.outdir = "./{}/".format(expid)
  if not os.path.exists(mxdx.outdir):
    os.makedirs(mxdx.outdir)
    
  # Balance the dataset before splitting
  # No aerosol type should make up more that 35% 
  # of the total number of obs
  # --------------------------------------
  mxdx.iValid = mxdx.balance(mxdx.nobs*0.35)


  # Create list of combinations
  # ---------------------------
  if combinations:
    mxdx.comblist, mxdx.combgroups = get_combinations(Input_nnr,Input_const)
  else:
    mxdx.comblist = []
        
  # Flatten Input_nnr into one list
  # -------------------------------
  input_nnr = flatten_list(Input_nnr)

  # Initialize arrays to hold stats
  # ------------------------------
  mxdx.nnr  = STATS(K,mxdx.comblist)
  mxdx.orig = STATS(K,mxdx.comblist)

  # Initialize K-folding
  # --------------------
  if K is not None:
    mxdx.kfold(K=K)

  if not combinations:
    if nHidden is None:
      mxdx.nHidden  = len(input_nnr)
    else:
      mxdx.nHidden = nHidden

    mxdx.topology = (len(input_nnr),) + (mxdx.nHidden,)*nHLayers + (len(Target),)

    train_test(mxdx,expid,input_nnr,Target,K)
  else:
    for c,Input in enumerate(mxdx.comblist):
      if nHidden is None:
        mxdx.nHidden  = len(Input)
      else:
        mxdx.nHidden = nHidden

      mxdx.topology = (len(Input),) + (mxdx.nHidden,)*nHLayers + (len(Target),)

      train_test(mxdx,'.'.join(Input),Input,Target,K,c=c,plotting=False)
  
  return mxdx

#--------------------------------------------------------------------------------------

def _setupMODIS(filename,retrieval,expid,
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
               Albedo=['CoxMunkLUT'],
               K=None):


  #------------------------------------------------------------
  # Read in data
  # --------------------------------------------------
  if retrieval.upper() == 'OCEAN':
    mxdx = ABC_Ocean(filename,Albedo=Albedo,verbose=1,laod=True)
  elif retrieval.upper() == 'LAND':
    mxdx = ABC_Land(filename,Albedo=Albedo,alb_max=0.25,verbose=1)
  elif retrieval.upper() == 'DEEP':
    mxdx = ABC_Deep(filename,Albedo=Albedo,verbose=1)

  # Create outdir if it doesn't exist
  # ---------------------------------
  mxdx.outdir = "./{}/".format(expid)
  if not os.path.exists(mxdx.outdir):
    os.makedirs(mxdx.outdir)

  # save some inputs
  # -----------------
  mxdx.expid   = expid
  mxdx.Target  = Target
  mxdx.K       = K
  mxdx.nHidden = nHidden
    
  # Balance the dataset before splitting
  # No aerosol type should make up more that 35% 
  # of the total number of obs
  # --------------------------------------
  mxdx.iValid = mxdx.balance(mxdx.nobs*0.35)

  # Flatten Input_nnr into one list
  # -------------------------------
  input_nnr = flatten_list(Input_nnr)

  # Create list of combinations
  # ---------------------------
  if combinations:
    mxdx.comblist, mxdx.combgroups = get_combinations(Input_nnr,Input_const)
  else:
    mxdx.comblist = [input_nnr]
        
  # Initialize arrays to hold stats
  # ------------------------------
  mxdx.nnr  = STATS(K,mxdx.comblist)
  mxdx.orig = STATS(K,mxdx.comblist)

  # Initialize K-folding
  # --------------------
  if K is not None:
    mxdx.kfold(K=K)

  # Create list of topologies
  # -------------------------  
  mxdx.topology = []
  if not combinations:
    if mxdx.nHidden is None:
      mxdx.nHidden  = len(input_nnr)
    else:
      mxdx.nHidden = nHidden

    mxdx.topology.append((len(input_nnr),) + (mxdx.nHidden,)*nHLayers + (len(Target),))

  else:
    for c,Input in enumerate(mxdx.comblist):
      if nHidden is None:
        mxdx.nHidden  = len(Input)
      else:
        mxdx.nHidden = nHidden

      mxdx.topology.append((len(Input),) + (mxdx.nHidden,)*nHLayers + (len(Target),))

  mxdx.combinations = combinations

  return mxdx

#---------------------------------------------------------------------
def _trainMODIS(mxdx):

  if not mxdx.combinations:
    _train(mxdx,mxdx.expid,0)
  else:
    for c,Input in enumerate(mxdx.comblist):
      _train(mxdx,'.'.join(Input),c)

# -------------------------------------------------------------------
def _TestMODIS(mxdx):

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
    doOcean = False
    doLand  = True
    doDeep  = False
    filename     = '/nobackup/6/NNR/Training/giant_C6_10km_Terra_20150921.nc'
    if doOcean:
      # OCEAN
      retrieval    = 'OCEAN'
      expid        = 'OCEAN_TEST'
      nHidden      = None  
      nHLayers     = 1
      combinations = True
      Input_const  = ['ScatteringAngle', 'GlintAngle', 'mRef870']
      Input_nnr    = [['mRef470','mRef550','mRef660','mRef1200','mRef1600','mRef2100']]
                    #['mRef470','mRef550','mRef660', 'mRef870',
                    #'mRef1200','mRef1600','mRef2100'],
                    #['ScatteringAngle', 'GlintAngle',
                    #'AMF', 'SolarZenith'],
                    #'cloud', 'CoxMunkLUT','CoxMunkBRF',['fdu','fcc','fsu'] ], 
      Target       = ['aTau550'] 
      Albedo       = ['CoxMunkLUT','CoxMunkBRF']      
      K            = 2           

      mxdx = _setupMODIS(filename, retrieval, expid,
                        nHidden      = nHidden,
                        nHLayers     = nHLayers,
                        combinations = combinations,
                        Input_const  = Input_const,
                        Input_nnr    = Input_nnr,                                         
                        Target       = Target,
                        Albedo       = Albedo,
                        K            = K)

    if doLand:
      # LAND
      retrieval    = 'LAND'
      expid        = 'LAND_Aux'
      nHidden      = None
      nHLayers     = 1
      combinations = True
      Input_const  = ['ScatteringAngle', 'GlintAngle', 'mRef870','mRef412','mRef440','mRef470','mRef550','mRef660','mRef1200','mRef1600','mRef2100']
      Input_nnr    = ['cloud', 'MOD43BClimAlbedo',['fdu','fcc','fsu'] ]
                    #,'mRef745'
                    #[['mRef412','mRef440','mRef470',
                    #'mRef550','mRef660', 'mRef870',
                    #'mRef1200','mRef1600','mRef2100'],
                    #['ScatteringAngle', 'GlintAngle',
                    #'AMF', 'SolarZenith']],
                    #'cloud', 'MOD43BClimAlbedo',['fdu','fcc','fsu'] ],  
      Target      = ['aTau550',]
      Albedo      = ['MOD43BClimAlbedo']
      K           = 2

      mxdx = _setupMODIS(filename, retrieval, expid,
                        nHidden      = nHidden,
                        nHLayers     = nHLayers,
                        combinations = combinations,
                        Input_const  = Input_const,
                        Input_nnr    = Input_nnr,                                         
                        Target       = Target,
                        Albedo       = Albedo,
                        K            = K)

    if doDeep:
      # DEEP BLUE
      retrieval    = 'DEEP'
      expid        = 'DEEP_AllInputsTest'
      nHidden      = None
      nHLayers     = 1
      combinations = True
      Input_const  = ['ScatteringAngle', 'GlintAngle', 'mRef660']
      Input_nnr    = [['mTau550','MOD43BClimAlbedo','fdu','fcc','fsu','mRef412','mRef470']]
                    #[['mRef412','mRef440','mRef470',
                    #'mRef550','mRef660', 'mRef870',
                    #'mRef1200','mRef1600','mRef2100'],
                    #['ScatteringAngle', 'GlintAngle',
                    #'AMF', 'SolarZenith']],
                    #'cloud', 'MOD43BClimAlbedo',['fdu','fcc','fsu'] ],  
      Target      = ['aTau550',]
      Albedo      = ['MOD43BClimAlbedo']
      K           = 2

      mxdl = _setupMODIS(filename, retrieval, expid,
                        nHidden      = nHidden,
                        nHLayers     = nHLayers,
                        combinations = combinations,
                        Input_const  = Input_const,
                        Input_nnr    = Input_nnr,                                         
                        Target       = Target,
                        Albedo       = Albedo,
                        K            = K)


    # Do Training and Testing
    # ------------------------
    _trainMODIS(mxdx)
    _TestMODIS(mxdx)

    if Input_const is not None:
      InputMaster = list((Input_const,) + tuple(Input_nnr))      
    else:
      InputMaster = Input_nnr

    if combinations:
      SummarizeCombinations(mxdx,InputMaster,yrange=None,sortname='slope')
      SummaryPDFs(mxdx)
