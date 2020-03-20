"""
   This module implements a Neural Net based MODIS AOD bias correction.  

   Important: The new SUPER2_combo datasets cannot be used for *land* because
              it lacks QA flads; the QAdark_l and QAdpbl_l columns are blanl.

   Arlindo da Silva, October 2010.
"""
import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
#import ffnet as nn
import pyobs.sknet as nn

from warnings import  warn
from pylab    import  cm, imshow, plot,figure
from pylab    import  xlabel, ylabel, title, grid, savefig, legend
from numpy    import  c_ as cat
from numpy    import  random, sort, pi, load, cos, log, std, exp
from numpy    import  reshape, arange, ones, zeros, interp
from numpy    import  meshgrid, concatenate
from matplotlib import ticker
from scipy    import  stats, mgrid
from anet     import  MISSING, LAND, OCEAN

#..............................................................
class aodFormat(ticker.Formatter):
    def __call__(self,x,pos=None):
        y = exp(x)-0.01
        return '%4.2f'%y

class NN(object):

    def train (self,Input=None,Target=None,nHidden=200,maxfun=2550,biases=True,
               topology=None, **kwargs):
        """
        Train the Neural Net, using a maximum of *maxfun* iterations.
        On input,
            Input   ---  string list with the name of the predictors;
                         if specified dataset specific default is
                         redefined.
            Target  ---  string list with the name of the targets;
                         if specified dataset specific default is
                         redefined.
           nHidden  ---  number of hidden nodes
           maxfun   ---  max number of iterations
           biases   ---  whether to include bias nodes
         topology   ---  Network topology; default is (nInput,nHidden,nTarget)
         
         Returns:
            Nothing.
        """
            
        # Possibly redefine Input/Targets
        # -------------------------------
        if Input != None:
            self.Input = Input
        if Target != None:
            self.Target = Target

        # Instantiate Neural Net
        # ----------------------
        if topology==None:
            topology = (len(self.Input), nHidden,len(self.Target))
        #self.net = nn.ffnet(nn.mlgraph(topology,biases=biases))
        self.net = nn.SKNET(nn.mlgraph(topology,biases=biases))

        # Add these attributes to net so that later on
        # we now how to apply it to regular MODIS data
        # --------------------------------------------
        self.net.InputNames = self.Input
        self.net.TargetNames = self.Target
        self.net.laod = self.laod
        if self.surface == 'ocean':
            self.net.Wind = self.Wind

        # Indices for training set
        # ------------------------
        try:
            iTrain = self.iTrain
        except AttributeError:
            iTrain = self.iValid # good QC marks
            
        # Prepare inputs and targets
        # --------------------------
        inputs  = self.getInputs(iTrain)
        targets = self.getTargets(iTrain) 

        # Train
        # -----
        if self.verbose>0:
            print "Starting training with %s inputs and %s targets"\
                  %(str(inputs.shape),str(targets.shape))
        self.net.train_tnc(inputs,targets, maxfun=maxfun, **kwargs)
#        self.net.train_bfgs(inputs,targets, maxfun=maxfun)


    def test(self,iprint=1,fname=None):

        # Indices for training set
        # ------------------------
        try:
            iTest = self.iTest
        except AttributeError:
            iTest = self.iValid
            
        # Prepare inputs and targets
        # --------------------------
        inputs  = self.getInputs(iTest)
        targets = self.getTargets(iTest) 

        return self.net.test(inputs,targets,iprint=iprint,filename=fname)
        
    def eval(self,I=None):
        if I == None: I = self.iValid
        return self.net(self.getInputs(I))

    __call__ = eval

    def derivative(self,I=None):
        if I == None: I = self.iValid
        return self.net.derivative(self.getInputs(I))
    
    def savenet(self,fname):
        nn.savenet(self.net,fname)

    def exportnet(self,fname):
        nn.exportnet(self.net,fname)
        
    def split (self,fTrain=0.9):
        """
        Splits the input dataset in training and testing subsets. No data is
        actually moved only attributes with indices iTrain/iTest are created;
        only data with an iValid Q/C flag is considered. On input, *fTrain* is
        the fraction of the dataset to be used for training.
        Returns: (nothing)
        """
        n = self.lon.size
        nTrain = int(fTrain * n)
        random.seed(32768) # so that we get the same permutation
        i = random.permutation(n)
        iValid = self.iValid[i]
        self.iTrain = i[0:nTrain][iValid[0:nTrain]] # Keep only good obs
        self.iTest  = i[nTrain:][iValid[nTrain:]]   # Keep only good obs

    def getInputs(self,I,Input=None):
        """
        Given a set of indices *I*, returns the corresponding
        inputs for a neural net evaluation.
        Returns: inputs
        """
        if self.verbose:
            print " "
            print "       Feature          Min      Max"
            print "  ------------------  -------  -------"
        if Input==None:
            Input = self.Input
        inputs = self.__dict__[Input[0]][I]
        if self.verbose:
            print "%20s %8.4f %8.4f"%(Input[0],inputs.min(),inputs.max())
        for var in Input[1:]:
            q = self.__dict__[var][I]
            inputs = cat[inputs,q]
            if self.verbose:
                print "%20s %8.4f %8.4f"%(var,q.min(),q.max())
        if self.verbose:
            print "  ------------------  -------  -------"
            print ""
        return inputs
    
    def getTargets(self,I):
        """
        Given a set of indices *I*, return the corresponding
        targets for a neural net evaluation:
        Returns: tagets
        """
        targets = self.__dict__[self.Target[0]][I]
        for var in self.Target[1:]:
            targets = cat[targets,self.__dict__[var][I]]
        if self.laod:
            targets = log(targets + 0.01)
        return targets
 
    def plotKDE(self,bins=None,I=None,figfile=None,
                x_label='AERONET'):
        """
        Plot Target vs Model using a 2D Kernel Density Estime.
        """
        if I==None: I = self.iValid # All data by default
        results = self.eval(I)
        targets = self.getTargets(I)
        if self.laod:
            formatter = aodFormat()
        else:
            formatter = None
        if bins == None:
            if self.laod:
                bins = arange(-5., 1., 0.1 )
            else:
                bins = arange(0., 0.6, 0.01 )
        x_bins = bins
        y_bins = bins
        if len(targets.shape) == 1:
            x_values = targets
            y_values = results.squeeze()
        else:
            x_values = targets[:,0]            
            y_values = results[:,0]
        _plotKDE(x_values,y_values,x_bins,y_bins,y_label='NNR',
                 formatter=formatter,x_label=x_label)        
        title("Log("+self.Target[0][1:]+"+0.01) - "+self.ident)
        if figfile != None:
            savefig(figfile)
            
    def plotScat(self,bins=None,I=None,figfile=None):
        """
        Plot Target vs Model using a 2D Kernel Density Estime.
        """
        if I==None: I = self.iTest # Testing data by default
        results = self.eval(I)
        targets = self.getTargets(I)
        original = log(self.__dict__['m'+self.Target[0][1:]][I] + 0.01)
        if bins == None:
            bins = arange(-5., 1., 0.1 )

        figure()
        plot(targets,original,'bo',label='Original')
        plot(targets,results,'ro',label='Corrected')
        legend(loc='upper left')
        plot(bins,bins,'k')
        grid()
        xlabel('AERONET')
        ylabel('MODIS')
        title("Log("+self.Target[0][1:]+"+0.01) - "+self.ident)
        if figfile != None:
            savefig(figfile)
            
#----------------------------------------------------------------------------    

class ABC_Ocean (OCEAN,NN):

    def __init__ (self,fname, Wind=None,
                  coxmunk_lut='/nobackup/NNR/Misc/coxmunk_lut.npz',
                  outliers=3., laod=True, verbose=0,
                  cloud_thresh=0.70, csvVersion = 1,
                  Input = ['mTau550','mTau470','mTau660','mTau870',
                           'mTAU550','mTAU470','mTAU660','mTAU870',
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
        Wind    ---  Type of wind related parameter to be read from a NPZ file. Typical
                     values are:
                                  merra_ustar
                                  merra_wind
                     Requires a NPZ file with the data.
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.Input = Input
        self.Target = Target
        self.verbose = verbose
        self.laod = laod
        self.Wind = Wind

        OCEAN.__init__(self,fname,csvVersion=csvVersion) # initialize superclass

        # Read in wind if desired
        # ------------------------
        if Wind != None:
            if csvVersion==1:
                self.wind = load(self.ident + "_" + Wind + ".npz")['var']
            else:
                self.wind = load(self.ident + "2_" + Wind + ".npz")['var']
        else:
            self.wind = zeros(self.N)
            print "WARNING: No longer using *ncep_windspd* because of too many undefs" 
            print "WARNING: wind attribute being set to zero"

        # Define wind speed dependent ocean albedo
        # ----------------------------------------
        self.getCoxMunk(coxmunk_lut)

        # Q/C
        # ---
        self.iValid =  (self.qa>0) & \
                      (self.aTau470 > 0.0) & (self.aTau550 > 0.0) &  \
                      (self.aTau660 > 0.0) & (self.aTau870 > 0.0) &  \
                      (self.mTau470 > 0.0) & (self.mTau550 > 0.0) &  \
                      (self.mTau660 > 0.0) & (self.mTau870 > 0.0) &  \
                      (self.mtau470 > 0.0) & (self.mtau550 > 0.0) &  \
                      (self.mtau660 > 0.0) & (self.mtau870 > 0.0) &  \
                      (self.cloud <cloud_thresh) & (self.wind>=0.0) &  \
                      (self.GlintAngle != MISSING )

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
            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      

#----------------------------------------------------------------------------    

class ABC_Land (LAND,NN):

    def __init__ (self, fname,
                  Albedo='albedo',
                  alb_min = 0.25,
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  csvVersion = 1,
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

        if csvVersion != 1:
            raise ValueError, 'must use CVS Version 1 for land because of QA flags'

        self.Input = Input
        self.Target = Target
        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname,csvVersion=csvVersion) # initialize superclass

        # Read in wind if desired
        # ------------------------
        self.albedo = load(self.ident + "_" + Albedo + ".npz")['var']

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

        print self.qa[self.iValid].shape
        
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

    def getAlbedo(self,npz_file):
        from grads import GrADS
        ga = GrADS(Echo=False,Window=False)
        ga('open albedo_clim.ctl')
        self.addVar(ga,npz_file,expr='albedo',clmYear=2000)

#---------------------------------------------------------------------------------
def _cat2 (X, Y):
    """
    Given 2 arrays of same shape, returns array of shape (2,N),
    where N = X.size = Y.size
    """
    xy = concatenate((X.ravel(),Y.ravel())) # shape is (N+N)
    return reshape(xy,(2,X.size))         # shape is (2,N)

def _plotKDE(x_values,y_values,x_bins=None,y_bins=None,
             x_label='AERONET', y_label='MODIS',formatter=None):
        """
        Plot Target vs Model using a 2D Kernel Density Estimate.
        """

        if x_bins == None: x_bins = arange(-5., 1., 0.1 )
        if y_bins == None: y_bins = x_bins

        Nx = len(x_bins)
        Ny = len(y_bins)

        print "Evaluating 2D kernel on grid with (Nx,Ny)=(%d,%d) ..."%(Nx,Ny)
        kernel = stats.kde.gaussian_kde(_cat2(x_values,y_values))
        X, Y = meshgrid(x_bins,y_bins)   # each has shape (Ny,Nx)
        Z = kernel(_cat2(X,Y))           # shape is (Ny*Nx)
        Z = reshape(Z,X.shape)

        fig = figure()
        # ax = fig.add_axes([0.1,0.1,0.75,0.75])
        ax = fig.add_axes([0.1,0.1,0.75,0.75])
        if formatter != None:
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
        imshow(Z, cmap=cm.gist_earth_r, origin='lower', 
               extent=(x_bins[0],x_bins[-1],y_bins[0],y_bins[-1]) )
        plot([x_bins[0],x_bins[-1]],[y_bins[0],y_bins[-1]],'k')
        xlabel(x_label)
        ylabel(y_label)
        grid()

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

def _testOcean(filename):

    Input_nnc = ['mTau550','mTau470','mTau660','mTau870',
                 'mTAU550','mTAU470','mTAU660','mTAU870',
                 'ScatteringAngle','GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'wind' ]

    Input_nnr1 = ['mRef470','mRef550','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'ScatteringAngle', 'GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'wind' ]

    Input_nnr2 = ['mRef470','mRef550','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'ScatteringAngle', 'GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_min = ['mTau550','mTAU550',
                 'GlintAngle', 'cloud', 'wind' ]

#    Target=['aTau_c',]
    Target=['aTau550',]
#    Target = [ 'aTau550','aTau470','aTau660', 'aTau870' ]

    # Read and split dataset in training/testing subsets
    # --------------------------------------------------
    mxdo = ABC_Ocean(filename,Wind='merra_wind',
                     verbose=1,csvVersion=2)
    mxdo.split()
    ident = mxdo.ident

    expid = 'nnr_002'
    for Input in (Input_nnr2,):

        nHidden = len(Input)
        topology = (len(Input), nHidden, len(Target))
        
        print "-"*80
        print "--> nHidden = ", nHidden
        print "-->  Inputs = ", Input
        
        mxdo.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology)
        out, reg = mxdo.test()

        mxdo.savenet(expid+"."+ident+'_Tau.net')

        # Plot KDE of corrected AOD
        # -------------------------
        mxdo.plotKDE(figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-corrected.png")

        # Plot KDE of uncorrected AOD
        # ---------------------------
        targets = mxdo.getTargets(mxdo.iValid).squeeze()
        original = log(mxdo.mTau550[mxdo.iValid]+0.01)
        _plotKDE(targets,original,y_label='Original MODIS')
        title("Log("+Target[0][1:]+"+0.01)- "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'.png')

        # Scatter diagram for testing
        # ---------------------------
        mxdo.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

    return mxdo

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
    mxdl = ABC_Land(filename,alb_min=0.25,verbose=1,csvVersion=1)
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

#---------------------------------------------------------------------
def _svrLand(filename):
    from sklearn.svm import SVR
    
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
    mxdl = ABC_Land(filename,alb_min=0.25,verbose=1,csvVersion=1)
    mxdl.split()

    ident = mxdl.ident
    expid = 'svr_002'
    for Input in (Input_nnr2,):

        X = mxdl.getInputs(I=mxdl.iValid,Input=Input)
        y = log(mxdl.aTau550[mxdl.iValid]+0.01)

        print X.shape, y.shape

        # SVR fitting
        # -----------
        svr = SVR(kernel='rbf', C=1e3, gamma=0.1,verbose=True)
        svr = svr.fit(X, y)

        # Plot KDE of corrected AOD
        # -------------------------
        targets = y
        results = svr.predict(X)
        _plotKDE(targets,results,y_label='SVR Fit')
        title("SVR Log("+Target[0][1:]+"+0.01)- "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'.png')

        # Scatter diagram for testing
        # ---------------------------
        #mxdl.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

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

    modl = _svrLand('SUPER_land.Terra.csv')

    # modo = _testOcean('SUPER2_combo.Terra.csv')
    # mydo = _testOcean('SUPER2_combo.Aqua.csv')
    # mydl = _testLand('SUPER_land.Aqua.csv')
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

    
