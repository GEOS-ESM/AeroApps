"""
   This module implements a Neural Net based MODIS Collection 6 Neural Net Retrieval.

   Arlindo da Silva, June 2015.

"""

import os, sys
sys.path.insert(0,'/home/pcastell/Enthought/Canopy_64bit/System/lib/python2.7/site-packages')
from   matplotlib.pyplot import  cm, imshow, plot, figure
from   matplotlib.pyplot import  xlabel, ylabel, title, grid, savefig, legend
import matplotlib.pyplot as      plt
from   matplotlib.ticker import  MultipleLocator
from   numpy             import  c_ as cat
from   numpy             import  random, sort, pi, load, cos, log, std, exp
from   numpy             import  reshape, arange, ones, zeros, interp, sqrt
from   numpy             import  meshgrid, concatenate, squeeze
import numpy             as      np
from   giant             import  MISSING, LAND, OCEAN
from   nn                import  NN, _plotKDE
import itertools
from   sklearn.linear_model import LinearRegression
from   multiprocessing    import cpu_count
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
               'albedo': 'Cox-Munk White Sky Albedo',
               'BRF': 'Cox-Munk Bidirectional Surface Reflectance',
               'fdu': 'MERRA2 Fraction Dust Aerosol',
               'fcc': 'MERRA2 Fraction Carbonaceous Aerosol',
               'fsu': 'MERRA2 Fraction Sulfate Aerosol',
               'year': 'Year'}
#----------------------------------------------------------------------------    

class ABC_Ocean (OCEAN,NN):

    def __init__ (self,fname, 
                  coxmunk_lut='/nobackup/NNR/Misc/coxmunk_lut.npz',
                  outliers=3., laod=True, verbose=0,
                  cloud_thresh=0.70,
                  glint_thresh=40.0,
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
        self.Wind = '' #need this for backwards compatibility

        # Define wind speed dependent ocean albedo
        # ----------------------------------------
        self.getCoxMunk(coxmunk_lut)
        self.BRF = squeeze(load(fnameRoot+'_CoxMunkAlbedo.npz')["albedo"])

        # Read in Aerosol Fractional Composition
        # --------------------------------------
        names = ('fdu','fss','fcc','fsu')
        for name in names:
            self.__dict__[name] = load(fnameRoot + "_MERRA2.npz")[name]


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
                      (self.GlintAngle != MISSING ) & (self.GlintAngle > glint_thresh) &\
                      (self.wind>=0.0) 

        # glint_thresh > 40 is a bit redundant b/c MOD04 should already give these a qa==0 or
        # does not retrieve.  However, there are a few cases (~200) where this does not happen.
        # the GlingAngle is very close to 40, greater than 38.  Not sure why these get through.


        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        if outliers > 0.:
            d = log(self.mTau550+0.01) - log(self.aTau550+0.01)
            dg = d[self.iValid]
            if self.verbose>0:
                print "Outlier removal: %d   sig_d = %f  nGood=%d "%(-1,std(dg),dg.size)
            for iter in range(3):
                self.iValid = self.iValid & (abs(d)<outliers*std(d[self.iValid]))
                if self.verbose>0:
                    dg = d[self.iValid]
                    print "Outlier removal: %d   sig_d = %f  nGood=%d "%(iter,std(dg),dg.size)

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
                  Albedo='albedo',
                  alb_min = 0.25,
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


        # Read in desired albedo
        # ------------------------
        self.albedo = load(self.ident + "_" + Albedo + ".npz")['albedo']

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
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660 >  0.0)       & \
                      (self.mSre2100>  0.0)       & \
                      (self.cloud<cloud_thresh)   & \
                      (self.ScatteringAngle<170.) & \
                      (self.albedo>0)             & \
                      (self.albedo<alb_min)
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        if outliers > 0.:
            d = log(self.mTau550+0.01) - log(self.aTau550+0.01)
            dg = d[self.iValid]
            if self.verbose>0:
                print "Outlier removal: %d   sig_d = %f  nGood=%d "%\
                      (-1,std(dg),dg.size)
            for iter in range(3):
                self.iValid = self.iValid & (abs(d)<outliers*std(d[self.iValid]))
                if self.verbose>0:
                    dg = d[self.iValid]
                    print "Outlier removal: %d   sig_d = %f  nGood=%d "\
                          %(iter,std(dg),dg.size)
            
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

def boxplot_imshow(data,blocks,masterlist,title,filename,
                   vseps=None, yrange=None,ylabel=None):
    nvars,ncomb = blocks.shape

    fig = plt.figure()
    ax  = plt.subplot(211)
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    plt.boxplot(data,showfliers=False,showbox=False,whis='range',
                whiskerprops={'linestyle':'-'})
    ax.set_xticks(np.arange(ncomb)+0.5)
    ax.set_xticklabels([])    


    if yrange is not None:
        ax.set_ylim(yrange)

    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # Minor Y Tick Marks    
    # ylim = ax.get_ylim()
    # dy   = (ylim[1] - ylim[0])/(len(yticks)-1)
    # minorLocator = MultipleLocator(0.5*dy)
    # ax.yaxis.set_minor_locator(minorLocator)

    if ylabel is not None:
      ax.set_ylabel(ylabel,fontsize=14)
    # Make vertical lines to separate bins of number of inputvars
    if vseps is not None:
        for v in vseps:
            ax.plot([v,v],np.array(ax.get_ylim()),'k-')

    plt.title(title)
    #ax.minorticks_on()
    plt.grid(True,axis='y',which='both')
    plt.tick_params(
        axis='y',          # changes apply to the y-axis
        which='major',     # major ticks are affected
        direction='out',
        right='off') 


    axblocks = plt.subplot(212)
    plt.imshow(blocks,interpolation='none',aspect='auto')
    axblocks.set_yticks(np.arange(nvars))
    axblocks.set_yticklabels(masterlist)
    axblocks.set_yticks(np.arange(nvars)+0.5,minor=True)
    axblocks.set_xticks(np.arange(ncomb)+0.5,minor=True)
    # plt.draw()  # this is needed because get_window_extent needs a renderer to work
    # yax = axblocks.get_yaxis()
    # # find the maximum width of the label on the major ticks
    # pad = max(T.label.get_window_extent().width for T in yax.majorTicks)
    # yax.set_tick_params(pad=pad)

    plt.tick_params(
        axis='both',          # changes apply to both
        which='major',     # major ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        left='off',
        right='off',
        labelbottom='off') # labels along the bottom edge are off 
    plt.grid(True,which='minor')
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.99,hspace=0.001)
    fig.set_size_inches(10.5, 5.5)   
    plt.savefig(filename, transparent='true',dpi=300)
    #plt.show()
    plt.close(fig)

# ---
def SummarizeCombinations(mxd,Input_nnr,yrange=None,sortname='slope'):
    """
    Create summary plot
    """
    outdir = mxd.outdir

    # Flatten Input_nnr into one list
    # -------------------------------
    Input = ()
    for i in Input_nnr:
      if type(i) is list:
        Input = Input + tuple(i)
      else:
        Input = Input + (i,)

    input_nnr = list(Input)

    expid      = '.'.join(input_nnr)
    ident      = mxd.ident
    nvars      = len(input_nnr)
    ngroups    = len(Input_nnr)
    ncomb      = len(mxd.comblist)
    blocks     = np.zeros([nvars,ncomb])
    masterlist = np.array(tuple(input_nnr))

    for c,comb in enumerate(mxd.comblist):
        blocks[:,c] = np.array([v == masterlist for v in comb]).sum(axis=0)

    blocks = np.insert(blocks,0,np.zeros(nvars),axis=1)  #for the original MODIS data

    newlist = []
    for var in masterlist:
        if mxd.sat == 'Terra':
          newlist.append(dict(MODVARNAMES,**VARNAMES)[var])
        else:
          newlist.append(dict(MYDVARNAMES,**VARNAMES)[var])

    masterlist = np.array(newlist)

    #nblocks = blocks.sum(axis=0)
    nblocks = [len(group) for group in mxd.combgroups]
    nblocks.insert(0,0)
    print 'nblocks',nblocks
    print "MASTERLIST", masterlist

    #--------------
    # Default sort by mean SLOPE
    #------------------
    sortvalue = mxd.nnr.__dict__[sortname]
    sortvalue  = np.insert(sortvalue,0,np.array(mxd.orig.__dict__[sortname][:,0]),axis=1)
    sortmetric = np.median(sortvalue, axis=0)
    isort = np.empty(0,dtype='int')
    vseps = np.empty(0,dtype='int')
    for i in np.sort(np.unique(nblocks)):
        print 'i',i
        istart  = np.where(nblocks == i)[0].min()
        vseps   = np.append(vseps,istart+0.5)
        isort   = np.append(isort,istart + np.argsort(sortmetric[nblocks == i]))
    blocks = blocks[:,isort]

    for i in np.arange(nvars):
        blocks[i,blocks[i,:] == 1] = i+1

    blocks = np.ma.masked_array(blocks)
    blocks.mask = [blocks == 0]


    def getplotdata(mxd,varname):
      vardata = mxd.nnr.__dict__[varname]
      vardata = np.insert(vardata,0,np.array(mxd.orig.__dict__[varname][:,0]),axis=1)
      vardata = vardata[:,isort]
      
      if vardata.shape[0] == 1:
        vardata = np.append(vardata,vardata,axis=0)

      return vardata

    boxplot_imshow(getplotdata(mxd,'slope'),
                   blocks,masterlist,
                   'Slope',
                   '{}/Slope.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=[0,1])

    boxplot_imshow(getplotdata(mxd,'R'),
                   blocks,masterlist,
                   'R',
                   '{}/R.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=[0,1])

    boxplot_imshow(getplotdata(mxd,'intercept'),
                   blocks,masterlist,
                   'Intercept',
                   '{}/Intercept.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)

    boxplot_imshow(getplotdata(mxd,'rmse'),
                   blocks,masterlist,
                   'RMSE',
                   '{}/RMSE.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)

    boxplot_imshow(getplotdata(mxd,'me'),
                   blocks,masterlist,
                   'Mean Bias',
                   '{}/ME.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)

    boxplot_imshow(getplotdata(mxd,'mae'),
                   blocks,masterlist,
                   'Mean Absolute Error',
                   '{}/MAE.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)    


#--------------------------------------------------------------------------------------

def _remove1():

    Input_all = ['mTau550','mTau470','mTau660','mTau870',
                 'mTAU550','mTAU470','mTAU660','mTAU870',
                 'ScatteringAngle','GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'ustar']

    Target=['aTau550',]

    mydo = ABC_Ocean('SUPER_ocean.Aqua.csv',ustar=False,verbose=1)
    mydo.split()

    for i in [-1,]+range(len(Input_all)):

        print "------------------------------------------------------------------------------"
        Input = Input_all[:] # make a copy of it
        if i<0:
            print "--> Excluding: (nothing)"
        else:
            print "--> Excluding: ", Input_all[i]
            del Input[i]         # delete ith item

        nHidden = len(Input)

        print "--> nHidden = ", nHidden
        print "-->  Inputs = ", Input
        
        mydo.train(Input=Input,Target=Target,nHidden=nHidden)
        out, reg = mydo.test()

#--------------------------------------------------------------------------------------

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
    mxd.kfold(K=K)
    k = 1
    for iTrain, iTest in mxd.kf:
      I = arange(mxd.nobs)
      iValid = I[mxd.iValid]
      mxd.iTrain = iValid[iTrain]
      mxd.iTest  = iValid[iTest]

      mxd.train(Input=['albedo' if x=='BRF' else x for x in Input],Target=Target,nHidden=nHidden,topology=topology,**kwargs)
      mxd.savenet(outdir+"/"+expid+"."+ident+'.k={}_Tau.net'.format(str(k)))
      TestStats(mxd,k-1,c)
      if plotting: make_plots(mxd,expid,ident+'.k={}'.format(str(k)))
      k = k + 1
  return mxd

#--------------------------------------------------------------------------------------

def _testOcean(filename,expid,
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


  #------------------------------------------------------------
  # Read in data
  # --------------------------------------------------
  mxdo = ABC_Ocean(filename,verbose=1)
  mxdo.outdir = "./{}/".format(expid)
  if not os.path.exists(mxdo.outdir):
    os.makedirs(mxdo.outdir)
    
  # Balance the dataset before splitting
  # No aerosol type should make up more that 35% 
  # of the total number of obs
  # --------------------------------------
  mxdo.iValid = mxdo.balance(mxdo.nobs*0.35)


  # Create list of combinations
  # ---------------------------
  comblist = []
  combgroups = []
  if combinations:
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

  mxdo.comblist = comblist
  mxdo.combgroups = combgroups

  # Flatten Input_nnr into one list
  # -------------------------------
  Input = ()
  for i in Input_nnr:
    if type(i) is list:
      Input = Input + tuple(i)
    else:
      Input = Input + (i,)

  input_nnr = list(Input)

  # Initialize arrays to hold stats
  # ------------------------------
  mxdo.nnr  = STATS(K,comblist)
  mxdo.orig = STATS(K,comblist)

  if not combinations:
    if nHidden is None:
      mxdo.nHidden  = len(input_nnr)
    else:
      mxdo.nHidden = nHidden

    mxdo.topology = (len(input_nnr),) + (mxdo.nHidden,)*nHLayers + (len(Target),)

    train_test(mxdo,expid,input_nnr,Target,K)
  else:
    for c,Input in enumerate(comblist):
      if nHidden is None:
        mxdo.nHidden  = len(Input)
      else:
        mxdo.nHidden = nHidden

      mxdo.topology = (len(Input),) + (mxdo.nHidden,)*nHLayers + (len(Target),)

      train_test(mxdo,'.'.join(Input),Input,Target,K,c=c,plotting=False)

    if Input_const is not None:
      SummarizeCombinations(mxdo,list((Input_const,) + tuple(Input_nnr)),yrange=None,sortname='slope')
    else:
      SummarizeCombinations(mxdo,Input_nnr,yrange=None,sortname='slope')
  
  return mxdo

#---------------------------------------------------------------------
def make_plots(mxd,expid,ident,I=None):  
  outdir = mxd.outdir
  if I is None:
    I = ones(mxd.lon.shape).astype(bool)
  # Plot KDE of corrected AOD
  # -------------------------
  # mxd.plotKDE(I=I,figfile=expid+"."+ident+"_kde-"+mxd.Target[0][1:]+"-corrected.png")
  targets  = mxd.getTargets(I).squeeze()
  results = mxd.eval(I).squeeze()
  _plotKDE(targets,results,y_label='NNR')
  title("Log("+mxd.Target[0][1:]+"+0.01)- "+ident)
  savefig(outdir+"/"+expid+"."+ident+"_kde-"+mxd.Target[0][1:]+'-corrected.png')

  # Plot KDE of uncorrected AOD
  # ---------------------------   
  original = log(mxd.mTau550[I]+0.01)
  _plotKDE(targets,original,y_label='Original MODIS')
  title("Log("+mxd.Target[0][1:]+"+0.01)- "+ident)
  savefig(outdir+"/"+expid+"."+ident+"_kde-"+mxd.Target[0][1:]+'.png')

  # Scatter diagram for testing
  # ---------------------------
  mxd.plotScat(I=I,figfile=outdir+"/"+expid+"."+ident+"_scat-"+mxd.Target[0][1:]+'.png')


#---------------------------------------------------------------------
def TestStats(mxd,K,C):
    if K is None:
      k = 0
    else:
      k = K

    if C is None:
      c = 0
    else:
      c = C

    # regression[0,2] = slope, intercept, r-value
    out, reg = mxd.test(iprint=False)

    mxd.nnr.slope[k,c]     = reg[0][0]
    mxd.nnr.intercept[k,c] = reg[0][1]
    mxd.nnr.R[k,c]         = reg[0][2]

    targets  = mxd.getTargets(mxd.iTest).squeeze()
    original = log(mxd.mTau550[mxd.iTest]+0.01)

    mxd.nnr.rmse[k,c] = rmse(out,targets)
    mxd.nnr.mae[k,c]  = mae(out,targets)
    mxd.nnr.me[k,c]   = me(out,targets)

    lm = LinearRegression()
    targets.shape = targets.shape + (1,)
    lm.fit(targets,original)
    mxd.orig.slope[k,c]     = lm.coef_[0]
    mxd.orig.intercept[k,c] = lm.intercept_
    mxd.orig.R[k,c]         = sqrt(lm.score(targets,original))

    mxd.orig.rmse[k,c] = rmse(original,targets)
    mxd.orig.mae[k,c]  = mae(original,targets)
    mxd.orig.me[k,c]   = me(original,targets)
    
# ---
def rmse(predictions, targets):
    return sqrt(((predictions - targets) ** 2).mean())
# ---
def mae(predictions, targets):
    return np.abs(predictions-targets).mean()
# ---
def me(predictions, targets):
    return (predictions-targets).mean()    
#---------------------------------------------------------------------
def _testLand(filename):

    Input_nnc = ['mTau550','mTau470','mTau660', 'mTau2100',
                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_nnr1 = ['mRef550','mRef470','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
#                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_nnr2 = ['mRef550','mRef470','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_min = ['mTau550',
                 'ScatteringAngle', 'cloud', 'albedo' ]

#    Target=['aTau_c',]
    Target=['aTau550',]
#    Target = [ 'aTau550', 'aTau470', 'aTau660', ]


    # Read and split dataset in training/testing subsets
    # --------------------------------------------------
    mxdl = ABC_Land(filename,alb_min=0.25,verbose=1)
    mxdl.split()

    ident = mxdl.ident
    expid = 'nnr_002'
    for Input in (Input_nnr2,):

        nHidden = len(Input)
        
        print "-"*80
        print "--> nHidden = ", nHidden
        print "-->  Inputs = ", Input
        
        mxdl.train(Input=Input,Target=Target,nHidden=nHidden)
        out, reg = mxdl.test()

        mxdl.savenet(expid+"."+ident+'_Tau.net')

        # Plot KDE of corrected AOD
        # -------------------------
        mxdl.plotKDE(figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-corrected.png")

        # Plot KDE of uncorrected AOD
        # ---------------------------
        targets = mxdl.getTargets(mxdl.iValid).squeeze()
        original = log(mxdl.mTau550[mxdl.iValid]+0.01)
        _plotKDE(targets,original,y_label='Original MODIS')
        title("Log("+Target[0][1:]+"+0.01)- "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'.png')

        # Scatter diagram for testing
        # ---------------------------
        mxdl.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

    return mxdl


def doAlbedo():
    from anet import LAND
    from grads import GrADS
    ga = GrADS(Echo=False,Window=False)
    ga('open albedo_clim.ctl')

    modl = LAND('SUPER_land.Terra.csv')
    modl.addVar(ga,'modl_albedo.npz',expr='albedo',clmYear=2000)
        
    mydl = LAND('SUPER_land.Aqua.csv')
    mydl.addVar(ga,'mydl_albedo.npz',expr='albedo',clmYear=2000)
        
def doWind():
    from anet import OCEAN
    from grads import GrADS
    ga = GrADS(Echo=False,Window=False)

    ga('xdfopen merra_slv-hourly.ddf')

    modo = OCEAN('SUPER2_combo.Terra.csv',)
    modo.addVar(ga,'modo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')

#    mydo = OCEAN('SUPER2_combo.Aqua.csv')
#    mydo.addVar(ga,'mydo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')

        
#------------------------------------------------------------------
  
if __name__ == "__main__":

    #modl = _svrLand('SUPER_land.Terra.csv')

    # modo = _testOcean('SUPER2_combo.Terra.csv')
    modo = _testOcean('/nobackup/6/NNR/Training/giant_C6_10km_Terra_20150921.nc','SA_GA_870_2100',
                      nHidden=None,
                      nHLayers=1,
                      combinations=True,
                      Input_nnr  =  ['mRef2100'],
                                    #['mRef470','mRef550','mRef660', 'mRef870',
                                    #'mRef1200','mRef1600','mRef2100'],
                                    #['ScatteringAngle', 'GlintAngle',
                                    #'AMF', 'SolarZenith'],
                                    #'cloud', 'albedo','BRF',['fdu','fcc','fsu'] ],  
                      Input_const = ['ScatteringAngle', 'GlintAngle', 'mRef870'],                  
                      Target = ['aTau550',],
                      K=3)
    # mydl = _testLand('/nobackup/6/NNR/Training/giant_C6_10km_9April2015.nc')
    # modl = _testLand('SUPER_land.Terra.csv')

def hold():

    doWind()
    doAlbedo()

    mxdx = _testOcean('SUPER2_combo.Terra.csv')

    mydo = _testOcean('SUPER2_combo.Aqua.csv')
        
    mxdl = _testLand('SUPER_land.Aqua.csv')

#    mxdx = _testOcean('SUPER2_combo.Aqua.csv')
    mxdx = _testOcean('SUPER2_combo.Terra.csv')

    mxdx = _testLand('SUPER_land.Aqua.csv')


